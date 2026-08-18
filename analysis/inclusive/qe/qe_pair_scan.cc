//******************************************************************
//*  qe_pair_scan.cc
//*
//*  Purpose:
//*    Many-run QE tensor pair scanner for RGC ND3 studies.
//*
//*    This is a diagnostic/ranking tool, not the final Azz extractor.
//*    It scans QE kin0 skim ROOT files, builds a run-level table, and
//*    ranks all possible run pairs by tensor lever arm and run-condition
//*    compatibility.
//*
//*  Intended input:
//*    rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1kin0tf1/spring23/*.root
//*
//*  Required ROOT trees:
//*    inclusive
//*    scaler
//*
//*  Optional ROOT tree:
//*    run_summary
//*
//*  Main outputs:
//*    outputs/qe_pair_scan_spring23_kin0.root
//*    outputs/qe_pair_scan_spring23_kin0_run_summary.csv
//*    outputs/qe_pair_scan_spring23_kin0_pair_summary.csv
//******************************************************************

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <set>
#include <limits>
#include <regex>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>

using namespace std;
namespace fs = std::filesystem;

const double PI = acos(-1.0);
const double RAD_TO_DEG = 180.0 / PI;
const string QE_PAIR_SCAN_VERSION = "qe_pair_scan_v2_cut_profile_2026_08_13";

//============================================================
// Utility helpers
//============================================================

bool isFinite(double x){
    return std::isfinite(x);
}

double radToDeg(double x){
    return x * RAD_TO_DEG;
}

bool strToBool01(const string &s){
    if(s == "1" || s == "true" || s == "TRUE" || s == "yes" || s == "YES") return true;
    if(s == "0" || s == "false" || s == "FALSE" || s == "no" || s == "NO") return false;

    cerr << "Error: expected 0/1 or true/false but got: " << s << endl;
    exit(1);
}

double tensorFromVector(double P){
    const double inside = 4.0 - 3.0 * P * P;
    if(inside < 0.0) return -9999.0;
    return 2.0 - sqrt(inside);
}

double tensorErrFromVector(double P, double Perr){
    const double inside = 4.0 - 3.0 * P * P;
    if(inside <= 0.0) return 0.0;
    return fabs(3.0 * P / sqrt(inside)) * fabs(Perr);
}

double safeRatio(double num, double den, double bad=-9999.0){
    if(!isFinite(num) || !isFinite(den) || fabs(den) <= 1.0e-30) return bad;
    return num / den;
}

double safeLogRatio(double num, double den, double bad=999.0){
    if(!isFinite(num) || !isFinite(den) || num <= 0.0 || den <= 0.0) return bad;
    return log(num / den);
}

int targetVectorSign(double P){
    if(P > 0.0) return 1;
    if(P < 0.0) return -1;
    return 0;
}

int targetPairClass(double P1, double P2){
    const int s1 = targetVectorSign(P1);
    const int s2 = targetVectorSign(P2);
    if(s1 == 0 || s2 == 0) return 0;
    if(s1 > 0 && s2 > 0) return 1;    // positive-positive target vector polarization
    if(s1 < 0 && s2 < 0) return -1;   // negative-negative target vector polarization
    return 2;                         // opposite-sign target vector polarization
}

string targetPairClassName(int cls){
    if(cls == 1) return "pos_pos";
    if(cls == -1) return "neg_neg";
    if(cls == 2) return "pos_neg";
    return "unknown";
}

void addTriggerBitsToMap(map<int, Long64_t> &m, Long64_t trigger){
    const unsigned long long trig = static_cast<unsigned long long>(trigger);
    for(int bit = 0; bit <= 62; ++bit){
        const unsigned long long mask = (1ULL << bit);
        if((trig & mask) != 0ULL){
            m[bit]++;
        }
    }
}

int dominantKey(const map<int, Long64_t> &m, int bad=-9999){
    if(m.empty()) return bad;
    int bestKey = bad;
    Long64_t bestVal = -1;
    for(const auto &kv : m){
        if(kv.second > bestVal){
            bestKey = kv.first;
            bestVal = kv.second;
        }
    }
    return bestKey;
}

bool mapHasKey(const map<int, Long64_t> &m, int key){
    return m.find(key) != m.end();
}

int extractRunNumberFromPath(const string &path){
    const string base = fs::path(path).filename().string();

    // Prefer the common RGC form *_017492.root or sidisdvcs_017492.hipo.
    regex sixDigitRun("([0-9]{6})");
    sregex_iterator it(base.begin(), base.end(), sixDigitRun);
    sregex_iterator end;

    int run = -1;
    for(; it != end; ++it){
        run = atoi((*it)[1].str().c_str());
    }

    return run;
}

vector<int> parseRunList(const string &s){
    vector<int> runs;
    if(s.empty() || s == "all" || s == "ALL") return runs;

    string token;
    stringstream ss(s);
    while(getline(ss, token, ',')){
        if(token.empty()) continue;
        runs.push_back(atoi(token.c_str()));
    }

    sort(runs.begin(), runs.end());
    runs.erase(unique(runs.begin(), runs.end()), runs.end());
    return runs;
}

bool runAllowed(int run, const vector<int> &allowed){
    if(allowed.empty()) return true;
    return binary_search(allowed.begin(), allowed.end(), run);
}

//============================================================
// Command-line configuration
//============================================================

struct Options {
    string indir = "";
    string contains = ".root";
    string runs = "all";

    string output_root = "outputs/qe_pair_scan_spring23_kin0.root";
    string run_csv = "outputs/qe_pair_scan_spring23_kin0_run_summary.csv";
    string pair_csv = "outputs/qe_pair_scan_spring23_kin0_pair_summary.csv";

    // Cut-profile convention:
    //   kin0 = inclusive QE detector/PID cuts only.
    //   kin1 = inclusive QE detector/PID cuts plus tight QE kinematic cuts.
    string cut_profile = "kin1";
    bool apply_kinematic_cuts = true;

    string charge_branch = "fcupgated_delta";
    double charge_rel_err = 0.002;

    double theta_min_deg = 7.8;
    double theta_max_deg = 8.2;
    double vz_min = -5.758;
    double vz_max = 1.5165;
    double chi2_abs_max = 3.0;
    double detector_pe_min = 0.0;
    double pe_min = 2.0;
    double lv_min = 14.0;
    double lw_min = 14.0;
    double pcal_min = 0.07;
    double sf_max = 0.28;
    bool use_sf_cut = true;
    bool use_ecin_pcal_cut = true;
    double ecin_pcal_slope = -0.625;
    double ecin_pcal_intercept = 0.15;
    double q2_min = 1.9433;
    double q2_max = 2.0574;
    double w_min = 0.0;
    double w_max = 1.073;
    double xmin = 0.8;
    double xmax = 1.6;

    bool require_target_match = true;
    bool require_is_qe = true;
    bool require_fd_status = false;

    double min_abs_delta_q = 0.0;
    double min_charge = 1.0;
    Long64_t min_final = 1;
    int max_pairs_to_write = -1; // -1 means all pairs
};

void printUsage(const char *prog){
    cout << endl;
    cout << "Usage:" << endl;
    cout << "  " << prog << " --indir <skim_root_dir> [options]" << endl;
    cout << endl;
    cout << "Required:" << endl;
    cout << "  --indir <dir>                  Directory containing per-run QE skim ROOT files" << endl;
    cout << endl;
    cout << "Common options:" << endl;
    cout << "  --contains <text>              Only process ROOT filenames containing this text. Default .root" << endl;
    cout << "  --runs <all|r1,r2,...>         Optional run filter. Default all" << endl;
    cout << "  --output-root <root>           Default outputs/qe_pair_scan_spring23_kin0.root" << endl;
    cout << "  --run-csv <csv>                Default outputs/qe_pair_scan_spring23_kin0_run_summary.csv" << endl;
    cout << "  --pair-csv <csv>               Default outputs/qe_pair_scan_spring23_kin0_pair_summary.csv" << endl;
    cout << "  --cut-profile <kin0|kin1>      kin0=detector/PID only, kin1=detector/PID + tight kinematics. Default kin1" << endl;
    cout << "  --charge-branch <name>         Default fcupgated_delta" << endl;
    cout << "  --require-is-qe 0|1            Default 1 for qe kin0 skims" << endl;
    cout << "  --require-target-match 0|1     Default 1" << endl;
    cout << "  --require-fd-status 0|1        Default 0" << endl;
    cout << "  --min-abs-delta-q <value>      Default 0" << endl;
    cout << "  --min-final <N>                Default 1" << endl;
    cout << "  --max-pairs <N>                Default -1, write all pairs" << endl;
    cout << endl;
    cout << "Cut convention:" << endl;
    cout << "  Detector/PID cuts are always active for kin0 and kin1:" << endl;
    cout << "    |chi2pid|<3, lv/lw>14, pcal>0.07, SF_QE<0.28, ECIN/PCAL cut." << endl;
    cout << "  Tight kinematic cuts are active only for --cut-profile kin1:" << endl;
    cout << "    theta 7.8-8.2 deg, vz -5.758 to 1.5165, pe>2," << endl;
    cout << "    Q2 1.9433-2.0574, W 0-1.073. xB is always the counting range." << endl;
    cout << endl;
    cout << "Example:" << endl;
    cout << "  " << prog << " \\" << endl;
    cout << "    --indir ../../rgc-workflow/skim/rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1kin0tf1/spring23 \\" << endl;
    cout << "    --output-root outputs/qe_pair_scan_spring23_kin0.root \\" << endl;
    cout << "    --run-csv outputs/qe_pair_scan_spring23_kin0_run_summary.csv \\" << endl;
    cout << "    --pair-csv outputs/qe_pair_scan_spring23_kin0_pair_summary.csv" << endl;
    cout << endl;
}

Options parseOptions(int argc, char **argv){
    Options opt;

    if(argc < 2){
        printUsage(argv[0]);
        exit(1);
    }

    for(int i = 1; i < argc; ++i){
        string a = argv[i];
        auto needValue = [&](const string &name){
            if(i + 1 >= argc){
                cerr << "Error: " << name << " requires a value." << endl;
                exit(1);
            }
        };

        if(a == "--indir"){
            needValue(a); opt.indir = argv[++i];
        }else if(a == "--contains"){
            needValue(a); opt.contains = argv[++i];
        }else if(a == "--runs"){
            needValue(a); opt.runs = argv[++i];
        }else if(a == "--output-root"){
            needValue(a); opt.output_root = argv[++i];
        }else if(a == "--run-csv"){
            needValue(a); opt.run_csv = argv[++i];
        }else if(a == "--pair-csv"){
            needValue(a); opt.pair_csv = argv[++i];
        }else if(a == "--cut-profile"){
            needValue(a); opt.cut_profile = argv[++i];
        }else if(a == "--charge-branch"){
            needValue(a); opt.charge_branch = argv[++i];
        }else if(a == "--charge-rel-err"){
            needValue(a); opt.charge_rel_err = atof(argv[++i]);
        }else if(a == "--theta-min"){
            needValue(a); opt.theta_min_deg = atof(argv[++i]);
        }else if(a == "--theta-max"){
            needValue(a); opt.theta_max_deg = atof(argv[++i]);
        }else if(a == "--vz-min"){
            needValue(a); opt.vz_min = atof(argv[++i]);
        }else if(a == "--vz-max"){
            needValue(a); opt.vz_max = atof(argv[++i]);
        }else if(a == "--chi2-max"){
            needValue(a); opt.chi2_abs_max = atof(argv[++i]);
        }else if(a == "--detector-pe-min"){
            needValue(a); opt.detector_pe_min = atof(argv[++i]);
        }else if(a == "--pe-min"){
            needValue(a); opt.pe_min = atof(argv[++i]);
        }else if(a == "--lv-min"){
            needValue(a); opt.lv_min = atof(argv[++i]);
        }else if(a == "--lw-min"){
            needValue(a); opt.lw_min = atof(argv[++i]);
        }else if(a == "--pcal-min"){
            needValue(a); opt.pcal_min = atof(argv[++i]);
        }else if(a == "--sf-max"){
            needValue(a); opt.sf_max = atof(argv[++i]);
        }else if(a == "--use-sf-cut"){
            needValue(a); opt.use_sf_cut = strToBool01(argv[++i]);
        }else if(a == "--use-ecin-pcal-cut"){
            needValue(a); opt.use_ecin_pcal_cut = strToBool01(argv[++i]);
        }else if(a == "--ecin-slope"){
            needValue(a); opt.ecin_pcal_slope = atof(argv[++i]);
        }else if(a == "--ecin-intercept"){
            needValue(a); opt.ecin_pcal_intercept = atof(argv[++i]);
        }else if(a == "--q2-min"){
            needValue(a); opt.q2_min = atof(argv[++i]);
        }else if(a == "--q2-max"){
            needValue(a); opt.q2_max = atof(argv[++i]);
        }else if(a == "--w-min"){
            needValue(a); opt.w_min = atof(argv[++i]);
        }else if(a == "--w-max"){
            needValue(a); opt.w_max = atof(argv[++i]);
        }else if(a == "--xmin"){
            needValue(a); opt.xmin = atof(argv[++i]);
        }else if(a == "--xmax"){
            needValue(a); opt.xmax = atof(argv[++i]);
        }else if(a == "--require-target-match"){
            needValue(a); opt.require_target_match = strToBool01(argv[++i]);
        }else if(a == "--require-is-qe"){
            needValue(a); opt.require_is_qe = strToBool01(argv[++i]);
        }else if(a == "--require-fd-status"){
            needValue(a); opt.require_fd_status = strToBool01(argv[++i]);
        }else if(a == "--min-abs-delta-q"){
            needValue(a); opt.min_abs_delta_q = atof(argv[++i]);
        }else if(a == "--min-charge"){
            needValue(a); opt.min_charge = atof(argv[++i]);
        }else if(a == "--min-final"){
            needValue(a); opt.min_final = atoll(argv[++i]);
        }else if(a == "--max-pairs"){
            needValue(a); opt.max_pairs_to_write = atoi(argv[++i]);
        }else if(a == "--help" || a == "-h"){
            printUsage(argv[0]);
            exit(0);
        }else{
            cerr << "Error: unknown option: " << a << endl;
            printUsage(argv[0]);
            exit(1);
        }
    }

    if(opt.indir.empty()){
        cerr << "Error: --indir is required." << endl;
        exit(1);
    }

    if(!fs::exists(opt.indir) || !fs::is_directory(opt.indir)){
        cerr << "Error: --indir is not a directory: " << opt.indir << endl;
        exit(1);
    }

    if(opt.cut_profile == "kin0"){
        opt.apply_kinematic_cuts = false;
    }else if(opt.cut_profile == "kin1"){
        opt.apply_kinematic_cuts = true;
    }else{
        cerr << "Error: unsupported --cut-profile: " << opt.cut_profile << endl;
        cerr << "       Allowed values are kin0 and kin1." << endl;
        exit(1);
    }

    return opt;
}

//============================================================
// Input branch containers
//============================================================

struct EventBranches {
    Int_t run = 0;
    Int_t event = 0;
    Int_t event_file_index = -9999;
    Int_t input_file_index = -9999;
    Int_t input_file_number = -9999;

    Long64_t trigger = -9999;

    Int_t helicity = -9999;
    Int_t helicity_valid = 0;

    Int_t target_match = 0;
    Float_t target_pol_offline = -9999.0;
    Float_t target_pol_offline_err = -9999.0;
    Float_t vector_pol = -9999.0;
    Float_t vector_pol_err = -9999.0;
    Float_t tensor_pol = -9999.0;
    Float_t tensor_pol_err = -9999.0;

    Int_t pid_e = 0;
    Int_t e_status = 0;
    Int_t e_det_status = 0;

    Float_t pe = -9999.0;
    Float_t thetae = -9999.0;
    Float_t vze = -9999.0;

    Float_t e_chi2pid = -9999.0;
    Float_t e_pcal_e = -9999.0;
    Float_t e_ecin_e = -9999.0;
    Float_t e_cal_e = -9999.0;
    Float_t e_sampling_fraction = -9999.0;
    Float_t e_lv = -9999.0;
    Float_t e_lw = -9999.0;

    Float_t Q2 = -9999.0;
    Float_t W = -9999.0;
    Float_t xB = -9999.0;

    Int_t is_qe = -9999;
};

struct BranchAvailability {
    bool event_file_index = false;
    bool input_file_index = false;
    bool input_file_number = false;
    bool trigger = false;
    bool is_qe = false;
    bool e_sampling_fraction = false;
    bool e_cal_e = false;
    bool vector_pol_err = false;
    bool target_pol_offline_err = false;
};

bool setBranchIfExists(TTree *tree, const char *name, void *address, bool required){
    if(!tree){
        cerr << "Error: null TTree in setBranchIfExists." << endl;
        exit(1);
    }

    TBranch *br = tree->GetBranch(name);
    if(br){
        tree->SetBranchStatus(name, 1);
        int status = tree->SetBranchAddress(name, address);
        if(status < 0){
            cerr << "Error: SetBranchAddress failed for branch " << name
                 << " status = " << status << endl;
            exit(1);
        }
        return true;
    }

    if(required){
        cerr << "Error: required branch not found: " << name << endl;
        cerr << "Available branches:" << endl;
        tree->GetListOfBranches()->Print();
        exit(1);
    }

    return false;
}

BranchAvailability setupInclusiveBranches(TTree *tree, EventBranches &e){
    BranchAvailability avail;

    tree->SetBranchStatus("*", 0);

    const bool req = true;
    const bool opt = false;

    setBranchIfExists(tree, "run", &e.run, req);
    setBranchIfExists(tree, "event", &e.event, req);
    avail.event_file_index = setBranchIfExists(tree, "event_file_index", &e.event_file_index, opt);
    avail.input_file_index = setBranchIfExists(tree, "input_file_index", &e.input_file_index, opt);
    avail.input_file_number = setBranchIfExists(tree, "input_file_number", &e.input_file_number, opt);
    avail.trigger = setBranchIfExists(tree, "trigger", &e.trigger, opt);

    setBranchIfExists(tree, "helicity", &e.helicity, opt);
    setBranchIfExists(tree, "helicity_valid", &e.helicity_valid, opt);

    setBranchIfExists(tree, "target_match", &e.target_match, opt);
    setBranchIfExists(tree, "target_pol_offline", &e.target_pol_offline, opt);
    avail.target_pol_offline_err = setBranchIfExists(tree, "target_pol_offline_err", &e.target_pol_offline_err, opt);
    setBranchIfExists(tree, "vector_pol", &e.vector_pol, opt);
    avail.vector_pol_err = setBranchIfExists(tree, "vector_pol_err", &e.vector_pol_err, opt);
    setBranchIfExists(tree, "tensor_pol", &e.tensor_pol, opt);
    setBranchIfExists(tree, "tensor_pol_err", &e.tensor_pol_err, opt);

    setBranchIfExists(tree, "pid_e", &e.pid_e, req);
    setBranchIfExists(tree, "e_status", &e.e_status, opt);
    setBranchIfExists(tree, "e_det_status", &e.e_det_status, opt);

    setBranchIfExists(tree, "pe", &e.pe, req);
    setBranchIfExists(tree, "thetae", &e.thetae, req);
    setBranchIfExists(tree, "vze", &e.vze, req);

    setBranchIfExists(tree, "e_chi2pid", &e.e_chi2pid, req);
    setBranchIfExists(tree, "e_pcal_e", &e.e_pcal_e, req);
    setBranchIfExists(tree, "e_ecin_e", &e.e_ecin_e, req);
    avail.e_cal_e = setBranchIfExists(tree, "e_cal_e", &e.e_cal_e, opt);
    avail.e_sampling_fraction = setBranchIfExists(tree, "e_sampling_fraction", &e.e_sampling_fraction, opt);
    setBranchIfExists(tree, "e_lv", &e.e_lv, req);
    setBranchIfExists(tree, "e_lw", &e.e_lw, req);

    setBranchIfExists(tree, "Q2", &e.Q2, req);
    setBranchIfExists(tree, "W", &e.W, req);
    setBranchIfExists(tree, "xB", &e.xB, req);

    avail.is_qe = setBranchIfExists(tree, "is_qe", &e.is_qe, opt);

    return avail;
}

//============================================================
// Counting, charge, and cut-flow containers
//============================================================

struct ChargeSummary {
    Long64_t n_scaler_rows = 0;
    double q_plus = 0.0;
    double q_minus = 0.0;
    double q_zero = 0.0;
    double q_total = 0.0;
};

struct CutFlow {
    Long64_t n_all = 0;
    Long64_t n_target_match = 0;
    Long64_t n_is_qe = 0;
    Long64_t n_pid = 0;
    Long64_t n_status = 0;
    Long64_t n_theta = 0;
    Long64_t n_vz = 0;
    Long64_t n_chi2 = 0;
    Long64_t n_lv_lw = 0;
    Long64_t n_pe = 0;
    Long64_t n_pcal = 0;
    Long64_t n_sf = 0;
    Long64_t n_ecin_pcal = 0;
    Long64_t n_q2 = 0;
    Long64_t n_w = 0;
    Long64_t n_xb_range = 0;
};

struct RunSummary {
    string file = "";
    int run = -1;

    Long64_t n_inclusive_entries = 0;
    Long64_t n_final = 0;
    Long64_t n_plus = 0;
    Long64_t n_minus = 0;
    Long64_t n_zero = 0;

    ChargeSummary charge;
    CutFlow cutflow;

    double mean_vector_pol = -9999.0;
    double mean_vector_pol_err = 0.0;
    double mean_target_pol_offline = -9999.0;
    double mean_target_pol_offline_err = 0.0;
    double mean_tensor_pol_branch = -9999.0;
    double P = -9999.0;
    double Perr = 0.0;
    double Q = -9999.0;
    double Qerr = 0.0;
    int P_sign = 0;

    Long64_t skim_n_events_all = -1;
    Long64_t skim_n_inclusive_written = -1;
    Long64_t skim_n_sidis_rows_written = -1;
    Long64_t skim_n_dihadron_rows_written = -1;

    map<int, Long64_t> all_by_trigger_bit;
    map<int, Long64_t> selected_by_trigger_bit;
    int dominant_trigger_bit_all = -9999;
    int dominant_trigger_bit_selected = -9999;
    int has_bit29_all = 0;
    int has_bit29_selected = 0;

    double raw_yield = -9999.0;
    double inclusive_yield = -9999.0;
    double final_yield = -9999.0;
    double yield_plus = -9999.0;
    double yield_minus = -9999.0;
    double theta_survival = -9999.0;
    double vz_survival = -9999.0;
    double chi2_survival = -9999.0;
    double lv_lw_survival = -9999.0;
    double ecin_pcal_survival = -9999.0;
    double q2_survival = -9999.0;
    double w_survival = -9999.0;
    double xb_survival = -9999.0;

    int usable = 0;
};

struct PairSummary {
    int run1 = 0;
    int run2 = 0;
    double P1 = -9999.0;
    double P2 = -9999.0;
    double Q1 = -9999.0;
    double Q2 = -9999.0;
    double deltaQ = -9999.0;
    double abs_deltaQ = -9999.0;
    int target_pair_class = 0;
    string target_pair_class_name = "unknown";

    double C1 = -9999.0;
    double C2 = -9999.0;
    double N1 = -9999.0;
    double N2 = -9999.0;
    double Y1 = -9999.0;
    double Y2 = -9999.0;
    double yield_ratio_12 = -9999.0;
    double rel_yield_diff_12 = -9999.0;
    double R2Q_screen = -9999.0;
    double R2P_screen = -9999.0;

    double raw_yield_ratio_12 = -9999.0;
    double inclusive_yield_ratio_12 = -9999.0;
    double theta_survival_ratio_12 = -9999.0;
    double lv_lw_survival_ratio_12 = -9999.0;
    double ecin_pcal_survival_ratio_12 = -9999.0;

    int dominant_trigger_bit_all_1 = -9999;
    int dominant_trigger_bit_all_2 = -9999;
    int dominant_trigger_bit_selected_1 = -9999;
    int dominant_trigger_bit_selected_2 = -9999;
    int same_dominant_trigger_bit_all = 0;
    int same_dominant_trigger_bit_selected = 0;
    int bit29_pair_flag = 0;

    double pair_quality_score = 0.0;
};

bool passFDStatus(Int_t status){
    return (status > -4000 && status <= -2000);
}

bool getQESamplingFraction(const EventBranches &e, double &sf){
    // Inclusive-QE detector/PID definition used by the skim:
    //   SF_QE = (E_PCAL + E_ECIN) / pe
    // Do not use the generic total-calorimeter sampling-fraction branch here.
    if(e.pe <= 0.0) {
        sf = -9999.0;
        return false;
    }
    if(e.e_pcal_e <= -900.0 || e.e_ecin_e <= -900.0) {
        sf = -9999.0;
        return false;
    }
    if(!isFinite(e.e_pcal_e) || !isFinite(e.e_ecin_e) || !isFinite(e.pe)){
        sf = -9999.0;
        return false;
    }
    sf = (e.e_pcal_e + e.e_ecin_e) / e.pe;
    return isFinite(sf);
}

bool passEventCuts(
    const EventBranches &e,
    const BranchAvailability &avail,
    const Options &opt,
    CutFlow &cf
){
    cf.n_all++;

    if(opt.require_target_match && e.target_match != 1) return false;
    cf.n_target_match++;

    if(opt.require_is_qe && avail.is_qe && e.is_qe != 1) return false;
    cf.n_is_qe++;

    if(e.pid_e != 11) return false;
    cf.n_pid++;

    if(opt.require_fd_status && !passFDStatus(e.e_status)) return false;
    cf.n_status++;

    // Detector/PID cut block: active for both kin0 and kin1.
    if(fabs(e.e_chi2pid) >= opt.chi2_abs_max) return false;
    cf.n_chi2++;

    if(e.e_lv <= opt.lv_min || e.e_lw <= opt.lw_min) return false;
    cf.n_lv_lw++;

    if(e.pe <= opt.detector_pe_min) return false;

    if(e.e_pcal_e <= opt.pcal_min) return false;
    cf.n_pcal++;

    if(opt.use_sf_cut){
        double sf_qe = -9999.0;
        if(!getQESamplingFraction(e, sf_qe)) return false;
        if(sf_qe >= opt.sf_max) return false;
    }
    cf.n_sf++;

    if(opt.use_ecin_pcal_cut){
        if(e.pe <= 0.0) return false;
        const double pcal_over_p = e.e_pcal_e / e.pe;
        const double ecin_over_p = e.e_ecin_e / e.pe;
        const double threshold = opt.ecin_pcal_slope * pcal_over_p + opt.ecin_pcal_intercept;
        if(ecin_over_p < threshold) return false;
    }
    cf.n_ecin_pcal++;

    // Tight QE kinematic cut block: active only for kin1.
    if(opt.apply_kinematic_cuts){
        const double theta_deg = radToDeg(e.thetae);
        if(theta_deg <= opt.theta_min_deg || theta_deg >= opt.theta_max_deg) return false;
        cf.n_theta++;

        if(e.vze <= opt.vz_min || e.vze >= opt.vz_max) return false;
        cf.n_vz++;

        if(e.pe <= opt.pe_min) return false;
        cf.n_pe++;

        if(e.Q2 <= opt.q2_min || e.Q2 >= opt.q2_max) return false;
        cf.n_q2++;

        if(e.W <= opt.w_min || e.W >= opt.w_max) return false;
        cf.n_w++;
    }else{
        cf.n_theta++;
        cf.n_vz++;
        cf.n_pe++;
        cf.n_q2++;
        cf.n_w++;
    }

    // xB is always the plotting/counting interval.
    if(e.xB < opt.xmin || e.xB >= opt.xmax) return false;
    cf.n_xb_range++;

    return true;
}

ChargeSummary readChargeSummary(const string &fileName, int run, const string &chargeBranch){
    ChargeSummary charge;

    TFile fin(fileName.c_str(), "READ");
    if(fin.IsZombie()){
        cerr << "Error: could not open ROOT file for scaler read: " << fileName << endl;
        exit(1);
    }

    TTree *scaler = dynamic_cast<TTree*>(fin.Get("scaler"));
    if(!scaler){
        cerr << "Error: scaler tree not found in file: " << fileName << endl;
        exit(1);
    }

    Int_t br_run = -9999;
    Int_t br_helicity = -9999;
    Float_t br_charge = 0.0;

    if(!scaler->GetBranch("run") || !scaler->GetBranch("helicity") || !scaler->GetBranch(chargeBranch.c_str())){
        cerr << "Error: scaler tree is missing one of these branches: run, helicity, "
             << chargeBranch << " in file " << fileName << endl;
        scaler->GetListOfBranches()->Print();
        exit(1);
    }

    scaler->SetBranchAddress("run", &br_run);
    scaler->SetBranchAddress("helicity", &br_helicity);
    scaler->SetBranchAddress(chargeBranch.c_str(), &br_charge);

    const Long64_t n = scaler->GetEntries();

    for(Long64_t i = 0; i < n; ++i){
        scaler->GetEntry(i);
        if(br_run != run) continue;
        if(!isFinite(br_charge)) continue;

        charge.n_scaler_rows++;
        charge.q_total += br_charge;

        if(br_helicity == 1){
            charge.q_plus += br_charge;
        }else if(br_helicity == -1){
            charge.q_minus += br_charge;
        }else{
            charge.q_zero += br_charge;
        }
    }

    return charge;
}

void readSkimRunSummary(TFile &fin, RunSummary &r){
    TTree *t = dynamic_cast<TTree*>(fin.Get("run_summary"));
    if(!t || t->GetEntries() <= 0) return;

    Long64_t n_events_all = -1;
    Long64_t n_inclusive_written = -1;
    Long64_t n_sidis_rows_written = -1;
    Long64_t n_dihadron_rows_written = -1;

    if(t->GetBranch("n_events_all")){
        t->SetBranchAddress("n_events_all", &n_events_all);
    }
    if(t->GetBranch("n_inclusive_written")){
        t->SetBranchAddress("n_inclusive_written", &n_inclusive_written);
    }
    if(t->GetBranch("n_sidis_rows_written")){
        t->SetBranchAddress("n_sidis_rows_written", &n_sidis_rows_written);
    }
    if(t->GetBranch("n_dihadron_rows_written")){
        t->SetBranchAddress("n_dihadron_rows_written", &n_dihadron_rows_written);
    }

    t->GetEntry(0);

    r.skim_n_events_all = n_events_all;
    r.skim_n_inclusive_written = n_inclusive_written;
    r.skim_n_sidis_rows_written = n_sidis_rows_written;
    r.skim_n_dihadron_rows_written = n_dihadron_rows_written;
}

void finalizeRunDerivedQuantities(RunSummary &r){
    if(r.charge.q_total > 0.0){
        r.raw_yield = (r.skim_n_events_all > 0) ? static_cast<double>(r.skim_n_events_all) / r.charge.q_total : -9999.0;
        r.inclusive_yield = static_cast<double>(r.n_inclusive_entries) / r.charge.q_total;
        r.final_yield = static_cast<double>(r.n_final) / r.charge.q_total;
    }
    r.yield_plus = safeRatio(static_cast<double>(r.n_plus), r.charge.q_plus);
    r.yield_minus = safeRatio(static_cast<double>(r.n_minus), r.charge.q_minus);

    r.theta_survival = safeRatio(static_cast<double>(r.cutflow.n_theta), static_cast<double>(r.cutflow.n_status));
    r.vz_survival = safeRatio(static_cast<double>(r.cutflow.n_vz), static_cast<double>(r.cutflow.n_theta));
    r.chi2_survival = safeRatio(static_cast<double>(r.cutflow.n_chi2), static_cast<double>(r.cutflow.n_vz));
    r.lv_lw_survival = safeRatio(static_cast<double>(r.cutflow.n_lv_lw), static_cast<double>(r.cutflow.n_chi2));
    r.ecin_pcal_survival = safeRatio(static_cast<double>(r.cutflow.n_ecin_pcal), static_cast<double>(r.cutflow.n_sf));
    r.q2_survival = safeRatio(static_cast<double>(r.cutflow.n_q2), static_cast<double>(r.cutflow.n_ecin_pcal));
    r.w_survival = safeRatio(static_cast<double>(r.cutflow.n_w), static_cast<double>(r.cutflow.n_q2));
    r.xb_survival = safeRatio(static_cast<double>(r.cutflow.n_xb_range), static_cast<double>(r.cutflow.n_w));

    r.P = r.mean_vector_pol;
    r.Perr = r.mean_vector_pol_err;
    r.Q = tensorFromVector(r.P);
    r.Qerr = tensorErrFromVector(r.P, r.Perr);
    r.P_sign = targetVectorSign(r.P);

    r.dominant_trigger_bit_all = dominantKey(r.all_by_trigger_bit);
    r.dominant_trigger_bit_selected = dominantKey(r.selected_by_trigger_bit);
    r.has_bit29_all = mapHasKey(r.all_by_trigger_bit, 29) ? 1 : 0;
    r.has_bit29_selected = mapHasKey(r.selected_by_trigger_bit, 29) ? 1 : 0;

    if(r.run > 0 && r.charge.q_total > 0.0 && r.n_final >= 0 && r.Q > -900.0){
        r.usable = 1;
    }else{
        r.usable = 0;
    }
}

RunSummary processRunFile(const string &fileName, const Options &opt){
    RunSummary r;
    r.file = fileName;
    r.run = extractRunNumberFromPath(fileName);

    if(r.run <= 0){
        cerr << "Warning: could not parse run number from filename, skipping: " << fileName << endl;
        r.usable = 0;
        return r;
    }

    TFile fin(fileName.c_str(), "READ");
    if(fin.IsZombie()){
        cerr << "Warning: could not open ROOT file, skipping: " << fileName << endl;
        r.usable = 0;
        return r;
    }

    readSkimRunSummary(fin, r);

    TTree *tree = dynamic_cast<TTree*>(fin.Get("inclusive"));
    if(!tree){
        cerr << "Warning: inclusive tree not found, skipping: " << fileName << endl;
        r.usable = 0;
        return r;
    }

    r.charge = readChargeSummary(fileName, r.run, opt.charge_branch);

    EventBranches e;
    BranchAvailability avail = setupInclusiveBranches(tree, e);

    const Long64_t nEntries = tree->GetEntries();
    r.n_inclusive_entries = nEntries;

    double sum_vector_pol = 0.0;
    double sum_vector_pol_err = 0.0;
    double sum_target_pol = 0.0;
    double sum_target_pol_err = 0.0;
    double sum_tensor_pol = 0.0;
    Long64_t n_vector_pol = 0;
    Long64_t n_vector_pol_err = 0;
    Long64_t n_target_pol = 0;
    Long64_t n_target_pol_err = 0;
    Long64_t n_tensor_pol = 0;

    for(Long64_t i = 0; i < nEntries; ++i){
        tree->GetEntry(i);

        if(e.run != r.run) continue;

        if(avail.trigger){
            addTriggerBitsToMap(r.all_by_trigger_bit, e.trigger);
        }

        if(e.vector_pol > -900.0 && isFinite(e.vector_pol)){
            sum_vector_pol += e.vector_pol;
            n_vector_pol++;
        }
        if(e.vector_pol_err > -900.0 && isFinite(e.vector_pol_err)){
            sum_vector_pol_err += fabs(e.vector_pol_err);
            n_vector_pol_err++;
        }
        if(e.target_pol_offline > -900.0 && isFinite(e.target_pol_offline)){
            sum_target_pol += e.target_pol_offline;
            n_target_pol++;
        }
        if(e.target_pol_offline_err > -900.0 && isFinite(e.target_pol_offline_err)){
            sum_target_pol_err += fabs(e.target_pol_offline_err);
            n_target_pol_err++;
        }
        if(e.tensor_pol > -900.0 && isFinite(e.tensor_pol)){
            sum_tensor_pol += e.tensor_pol;
            n_tensor_pol++;
        }

        if(!passEventCuts(e, avail, opt, r.cutflow)) continue;

        if(avail.trigger){
            addTriggerBitsToMap(r.selected_by_trigger_bit, e.trigger);
        }

        if(e.helicity == 1){
            r.n_plus++;
        }else if(e.helicity == -1){
            r.n_minus++;
        }else{
            r.n_zero++;
        }
        r.n_final++;
    }

    if(n_vector_pol > 0) r.mean_vector_pol = sum_vector_pol / static_cast<double>(n_vector_pol);
    if(n_vector_pol_err > 0) r.mean_vector_pol_err = sum_vector_pol_err / static_cast<double>(n_vector_pol_err);
    if(n_target_pol > 0) r.mean_target_pol_offline = sum_target_pol / static_cast<double>(n_target_pol);
    if(n_target_pol_err > 0) r.mean_target_pol_offline_err = sum_target_pol_err / static_cast<double>(n_target_pol_err);
    if(n_tensor_pol > 0) r.mean_tensor_pol_branch = sum_tensor_pol / static_cast<double>(n_tensor_pol);

    finalizeRunDerivedQuantities(r);

    return r;
}

vector<string> findInputFiles(const Options &opt){
    vector<string> files;
    const vector<int> allowedRuns = parseRunList(opt.runs);

    for(const auto &entry : fs::directory_iterator(opt.indir)){
        if(!entry.is_regular_file()) continue;

        const string path = entry.path().string();
        const string base = entry.path().filename().string();

        if(base.size() < 5 || base.substr(base.size() - 5) != ".root") continue;
        if(!opt.contains.empty() && base.find(opt.contains) == string::npos) continue;

        const int run = extractRunNumberFromPath(path);
        if(run <= 0) continue;
        if(!runAllowed(run, allowedRuns)) continue;

        files.push_back(path);
    }

    sort(files.begin(), files.end(), [](const string &a, const string &b){
        return extractRunNumberFromPath(a) < extractRunNumberFromPath(b);
    });

    return files;
}

PairSummary makePairSummary(const RunSummary &a, const RunSummary &b){
    PairSummary p;

    p.run1 = a.run;
    p.run2 = b.run;
    p.P1 = a.P;
    p.P2 = b.P;
    p.Q1 = a.Q;
    p.Q2 = b.Q;
    p.deltaQ = a.Q - b.Q;
    p.abs_deltaQ = fabs(p.deltaQ);
    p.target_pair_class = targetPairClass(a.P, b.P);
    p.target_pair_class_name = targetPairClassName(p.target_pair_class);

    p.C1 = a.charge.q_total;
    p.C2 = b.charge.q_total;
    p.N1 = static_cast<double>(a.n_final);
    p.N2 = static_cast<double>(b.n_final);
    p.Y1 = safeRatio(p.N1, p.C1);
    p.Y2 = safeRatio(p.N2, p.C2);
    p.yield_ratio_12 = safeRatio(p.Y1, p.Y2);
    p.rel_yield_diff_12 = (p.yield_ratio_12 > -9000.0) ? p.yield_ratio_12 - 1.0 : -9999.0;

    const double denomQ = p.Q1 * p.Y2 - p.Q2 * p.Y1;
    const double denomP = p.P1 * p.Y2 - p.P2 * p.Y1;
    p.R2Q_screen = safeRatio(p.Y1 - p.Y2, denomQ);
    p.R2P_screen = safeRatio(p.Y1 - p.Y2, denomP);

    p.raw_yield_ratio_12 = safeRatio(a.raw_yield, b.raw_yield);
    p.inclusive_yield_ratio_12 = safeRatio(a.inclusive_yield, b.inclusive_yield);
    p.theta_survival_ratio_12 = safeRatio(a.theta_survival, b.theta_survival);
    p.lv_lw_survival_ratio_12 = safeRatio(a.lv_lw_survival, b.lv_lw_survival);
    p.ecin_pcal_survival_ratio_12 = safeRatio(a.ecin_pcal_survival, b.ecin_pcal_survival);

    p.dominant_trigger_bit_all_1 = a.dominant_trigger_bit_all;
    p.dominant_trigger_bit_all_2 = b.dominant_trigger_bit_all;
    p.dominant_trigger_bit_selected_1 = a.dominant_trigger_bit_selected;
    p.dominant_trigger_bit_selected_2 = b.dominant_trigger_bit_selected;
    p.same_dominant_trigger_bit_all = (a.dominant_trigger_bit_all == b.dominant_trigger_bit_all) ? 1 : 0;
    p.same_dominant_trigger_bit_selected = (a.dominant_trigger_bit_selected == b.dominant_trigger_bit_selected) ? 1 : 0;
    p.bit29_pair_flag = (a.has_bit29_all || b.has_bit29_all) ? 1 : 0;

    // Quality score: larger is better. It rewards tensor lever arm and penalizes
    // yield/cut/trigger incompatibility. It is only a ranking guide; use the
    // individual fields for final judgment.
    const double lever = min(p.abs_deltaQ / 0.15, 2.0);
    const double yPenalty = exp(-fabs(safeLogRatio(p.Y1, p.Y2, 20.0)));
    const double rawPenalty = exp(-0.5 * fabs(safeLogRatio(a.raw_yield, b.raw_yield, 20.0)));
    const double thetaPenalty = exp(-0.5 * fabs(safeLogRatio(a.theta_survival, b.theta_survival, 20.0)));
    const double lvPenalty = exp(-0.5 * fabs(safeLogRatio(a.lv_lw_survival, b.lv_lw_survival, 20.0)));
    const double ecinPenalty = exp(-0.5 * fabs(safeLogRatio(a.ecin_pcal_survival, b.ecin_pcal_survival, 20.0)));
    const double trigPenalty = p.same_dominant_trigger_bit_all ? 1.0 : 0.5;

    p.pair_quality_score = 100.0 * lever * yPenalty * rawPenalty * thetaPenalty * lvPenalty * ecinPenalty * trigPenalty;

    return p;
}

vector<PairSummary> buildPairs(const vector<RunSummary> &runs, const Options &opt){
    vector<PairSummary> pairs;

    for(size_t i = 0; i < runs.size(); ++i){
        const RunSummary &a = runs[i];
        if(!a.usable) continue;
        if(a.charge.q_total < opt.min_charge) continue;
        if(a.n_final < opt.min_final) continue;

        for(size_t j = i + 1; j < runs.size(); ++j){
            const RunSummary &b = runs[j];
            if(!b.usable) continue;
            if(b.charge.q_total < opt.min_charge) continue;
            if(b.n_final < opt.min_final) continue;

            PairSummary p = makePairSummary(a, b);
            if(p.abs_deltaQ < opt.min_abs_delta_q) continue;
            if(p.target_pair_class == 0) continue;
            pairs.push_back(p);
        }
    }

    sort(pairs.begin(), pairs.end(), [](const PairSummary &a, const PairSummary &b){
        if(a.pair_quality_score != b.pair_quality_score) return a.pair_quality_score > b.pair_quality_score;
        if(a.abs_deltaQ != b.abs_deltaQ) return a.abs_deltaQ > b.abs_deltaQ;
        return fabs(a.rel_yield_diff_12) < fabs(b.rel_yield_diff_12);
    });

    if(opt.max_pairs_to_write > 0 && static_cast<int>(pairs.size()) > opt.max_pairs_to_write){
        pairs.resize(opt.max_pairs_to_write);
    }

    return pairs;
}

void ensureOutputDir(const string &path){
    fs::path p(path);
    if(p.has_parent_path()){
        fs::create_directories(p.parent_path());
    }
}

void writeRunCSV(const string &csv, const vector<RunSummary> &runs, const Options &opt){
    ensureOutputDir(csv);
    ofstream f(csv.c_str());

    f << "cut_profile,apply_kinematic_cuts,run,file,P,Perr,P_sign,Q,Qerr,charge_total,charge_plus,charge_minus,n_scaler_rows,skim_n_events_all,n_inclusive_entries,n_final,n_plus,n_minus,n_zero,raw_yield,inclusive_yield,final_yield,yield_plus,yield_minus,theta_survival,vz_survival,chi2_survival,lv_lw_survival,ecin_pcal_survival,q2_survival,w_survival,xb_survival,dominant_trigger_bit_all,dominant_trigger_bit_selected,has_bit29_all,has_bit29_selected,usable\n";

    for(const RunSummary &r : runs){
        f << opt.cut_profile << "," << (opt.apply_kinematic_cuts ? 1 : 0) << ","
          << r.run << ","
          << r.file << ","
          << r.P << "," << r.Perr << "," << r.P_sign << ","
          << r.Q << "," << r.Qerr << ","
          << r.charge.q_total << "," << r.charge.q_plus << "," << r.charge.q_minus << "," << r.charge.n_scaler_rows << ","
          << r.skim_n_events_all << "," << r.n_inclusive_entries << "," << r.n_final << ","
          << r.n_plus << "," << r.n_minus << "," << r.n_zero << ","
          << r.raw_yield << "," << r.inclusive_yield << "," << r.final_yield << ","
          << r.yield_plus << "," << r.yield_minus << ","
          << r.theta_survival << "," << r.vz_survival << "," << r.chi2_survival << ","
          << r.lv_lw_survival << "," << r.ecin_pcal_survival << ","
          << r.q2_survival << "," << r.w_survival << "," << r.xb_survival << ","
          << r.dominant_trigger_bit_all << "," << r.dominant_trigger_bit_selected << ","
          << r.has_bit29_all << "," << r.has_bit29_selected << "," << r.usable << "\n";
    }
}

void writePairCSV(const string &csv, const vector<PairSummary> &pairs, const Options &opt){
    ensureOutputDir(csv);
    ofstream f(csv.c_str());

    f << "rank,cut_profile,apply_kinematic_cuts,run1,run2,P1,P2,Q1,Q2,deltaQ,abs_deltaQ,target_pair_class,target_pair_class_name,C1,C2,N1,N2,Y1,Y2,yield_ratio_12,rel_yield_diff_12,R2Q_screen,R2P_screen,raw_yield_ratio_12,inclusive_yield_ratio_12,theta_survival_ratio_12,lv_lw_survival_ratio_12,ecin_pcal_survival_ratio_12,dominant_trigger_bit_all_1,dominant_trigger_bit_all_2,dominant_trigger_bit_selected_1,dominant_trigger_bit_selected_2,same_dominant_trigger_bit_all,same_dominant_trigger_bit_selected,bit29_pair_flag,pair_quality_score\n";

    for(size_t i = 0; i < pairs.size(); ++i){
        const PairSummary &p = pairs[i];
        f << (i + 1) << "," << opt.cut_profile << "," << (opt.apply_kinematic_cuts ? 1 : 0) << ","
          << p.run1 << "," << p.run2 << ","
          << p.P1 << "," << p.P2 << ","
          << p.Q1 << "," << p.Q2 << ","
          << p.deltaQ << "," << p.abs_deltaQ << ","
          << p.target_pair_class << "," << p.target_pair_class_name << ","
          << p.C1 << "," << p.C2 << ","
          << p.N1 << "," << p.N2 << ","
          << p.Y1 << "," << p.Y2 << ","
          << p.yield_ratio_12 << "," << p.rel_yield_diff_12 << ","
          << p.R2Q_screen << "," << p.R2P_screen << ","
          << p.raw_yield_ratio_12 << "," << p.inclusive_yield_ratio_12 << ","
          << p.theta_survival_ratio_12 << "," << p.lv_lw_survival_ratio_12 << "," << p.ecin_pcal_survival_ratio_12 << ","
          << p.dominant_trigger_bit_all_1 << "," << p.dominant_trigger_bit_all_2 << ","
          << p.dominant_trigger_bit_selected_1 << "," << p.dominant_trigger_bit_selected_2 << ","
          << p.same_dominant_trigger_bit_all << "," << p.same_dominant_trigger_bit_selected << ","
          << p.bit29_pair_flag << "," << p.pair_quality_score << "\n";
    }
}

void writeRootOutput(const string &rootFile, const vector<RunSummary> &runs, const vector<PairSummary> &pairs){
    ensureOutputDir(rootFile);

    TFile fout(rootFile.c_str(), "RECREATE");
    if(fout.IsZombie()){
        cerr << "Error: could not create ROOT output: " << rootFile << endl;
        exit(1);
    }

    // Run tree
    TTree tr("qe_pair_scan_run_summary", "Run-level QE pair-scan diagnostics");

    Int_t run = 0;
    string file = "";
    Double_t P = 0.0, Perr = 0.0, Q = 0.0, Qerr = 0.0;
    Int_t P_sign = 0;
    Double_t charge_total = 0.0, charge_plus = 0.0, charge_minus = 0.0;
    Long64_t n_scaler_rows = 0, skim_n_events_all = 0, n_inclusive_entries = 0, n_final = 0;
    Long64_t n_plus = 0, n_minus = 0, n_zero = 0;
    Double_t raw_yield = 0.0, inclusive_yield = 0.0, final_yield = 0.0;
    Double_t yield_plus = 0.0, yield_minus = 0.0;
    Double_t theta_survival = 0.0, vz_survival = 0.0, chi2_survival = 0.0;
    Double_t lv_lw_survival = 0.0, ecin_pcal_survival = 0.0, q2_survival = 0.0, w_survival = 0.0, xb_survival = 0.0;
    Int_t dominant_trigger_bit_all = 0, dominant_trigger_bit_selected = 0, has_bit29_all = 0, has_bit29_selected = 0, usable = 0;

    tr.Branch("run", &run, "run/I");
    tr.Branch("file", &file);
    tr.Branch("P", &P, "P/D");
    tr.Branch("Perr", &Perr, "Perr/D");
    tr.Branch("P_sign", &P_sign, "P_sign/I");
    tr.Branch("Q", &Q, "Q/D");
    tr.Branch("Qerr", &Qerr, "Qerr/D");
    tr.Branch("charge_total", &charge_total, "charge_total/D");
    tr.Branch("charge_plus", &charge_plus, "charge_plus/D");
    tr.Branch("charge_minus", &charge_minus, "charge_minus/D");
    tr.Branch("n_scaler_rows", &n_scaler_rows, "n_scaler_rows/L");
    tr.Branch("skim_n_events_all", &skim_n_events_all, "skim_n_events_all/L");
    tr.Branch("n_inclusive_entries", &n_inclusive_entries, "n_inclusive_entries/L");
    tr.Branch("n_final", &n_final, "n_final/L");
    tr.Branch("n_plus", &n_plus, "n_plus/L");
    tr.Branch("n_minus", &n_minus, "n_minus/L");
    tr.Branch("n_zero", &n_zero, "n_zero/L");
    tr.Branch("raw_yield", &raw_yield, "raw_yield/D");
    tr.Branch("inclusive_yield", &inclusive_yield, "inclusive_yield/D");
    tr.Branch("final_yield", &final_yield, "final_yield/D");
    tr.Branch("yield_plus", &yield_plus, "yield_plus/D");
    tr.Branch("yield_minus", &yield_minus, "yield_minus/D");
    tr.Branch("theta_survival", &theta_survival, "theta_survival/D");
    tr.Branch("vz_survival", &vz_survival, "vz_survival/D");
    tr.Branch("chi2_survival", &chi2_survival, "chi2_survival/D");
    tr.Branch("lv_lw_survival", &lv_lw_survival, "lv_lw_survival/D");
    tr.Branch("ecin_pcal_survival", &ecin_pcal_survival, "ecin_pcal_survival/D");
    tr.Branch("q2_survival", &q2_survival, "q2_survival/D");
    tr.Branch("w_survival", &w_survival, "w_survival/D");
    tr.Branch("xb_survival", &xb_survival, "xb_survival/D");
    tr.Branch("dominant_trigger_bit_all", &dominant_trigger_bit_all, "dominant_trigger_bit_all/I");
    tr.Branch("dominant_trigger_bit_selected", &dominant_trigger_bit_selected, "dominant_trigger_bit_selected/I");
    tr.Branch("has_bit29_all", &has_bit29_all, "has_bit29_all/I");
    tr.Branch("has_bit29_selected", &has_bit29_selected, "has_bit29_selected/I");
    tr.Branch("usable", &usable, "usable/I");

    for(const RunSummary &r : runs){
        run = r.run;
        file = r.file;
        P = r.P; Perr = r.Perr; P_sign = r.P_sign; Q = r.Q; Qerr = r.Qerr;
        charge_total = r.charge.q_total; charge_plus = r.charge.q_plus; charge_minus = r.charge.q_minus;
        n_scaler_rows = r.charge.n_scaler_rows;
        skim_n_events_all = r.skim_n_events_all;
        n_inclusive_entries = r.n_inclusive_entries;
        n_final = r.n_final; n_plus = r.n_plus; n_minus = r.n_minus; n_zero = r.n_zero;
        raw_yield = r.raw_yield; inclusive_yield = r.inclusive_yield; final_yield = r.final_yield;
        yield_plus = r.yield_plus; yield_minus = r.yield_minus;
        theta_survival = r.theta_survival; vz_survival = r.vz_survival; chi2_survival = r.chi2_survival;
        lv_lw_survival = r.lv_lw_survival; ecin_pcal_survival = r.ecin_pcal_survival;
        q2_survival = r.q2_survival; w_survival = r.w_survival; xb_survival = r.xb_survival;
        dominant_trigger_bit_all = r.dominant_trigger_bit_all;
        dominant_trigger_bit_selected = r.dominant_trigger_bit_selected;
        has_bit29_all = r.has_bit29_all;
        has_bit29_selected = r.has_bit29_selected;
        usable = r.usable;
        tr.Fill();
    }

    tr.Write();

    // Pair tree
    TTree tp("qe_pair_scan_pair_summary", "Pair-level QE tensor pair ranking diagnostics");

    Int_t rank = 0, run1 = 0, run2 = 0, target_pair_class = 0;
    string target_pair_class_name = "";
    Double_t P1 = 0.0, P2 = 0.0, Q1 = 0.0, Q2 = 0.0, deltaQ = 0.0, abs_deltaQ = 0.0;
    Double_t C1 = 0.0, C2 = 0.0, N1 = 0.0, N2 = 0.0, Y1 = 0.0, Y2 = 0.0;
    Double_t yield_ratio_12 = 0.0, rel_yield_diff_12 = 0.0, R2Q_screen = 0.0, R2P_screen = 0.0;
    Double_t raw_yield_ratio_12 = 0.0, inclusive_yield_ratio_12 = 0.0, theta_survival_ratio_12 = 0.0;
    Double_t lv_lw_survival_ratio_12 = 0.0, ecin_pcal_survival_ratio_12 = 0.0;
    Int_t dominant_trigger_bit_all_1 = 0, dominant_trigger_bit_all_2 = 0;
    Int_t dominant_trigger_bit_selected_1 = 0, dominant_trigger_bit_selected_2 = 0;
    Int_t same_dominant_trigger_bit_all = 0, same_dominant_trigger_bit_selected = 0, bit29_pair_flag = 0;
    Double_t pair_quality_score = 0.0;

    tp.Branch("rank", &rank, "rank/I");
    tp.Branch("run1", &run1, "run1/I");
    tp.Branch("run2", &run2, "run2/I");
    tp.Branch("P1", &P1, "P1/D");
    tp.Branch("P2", &P2, "P2/D");
    tp.Branch("Q1", &Q1, "Q1/D");
    tp.Branch("Q2", &Q2, "Q2/D");
    tp.Branch("deltaQ", &deltaQ, "deltaQ/D");
    tp.Branch("abs_deltaQ", &abs_deltaQ, "abs_deltaQ/D");
    tp.Branch("target_pair_class", &target_pair_class, "target_pair_class/I");
    tp.Branch("target_pair_class_name", &target_pair_class_name);
    tp.Branch("C1", &C1, "C1/D");
    tp.Branch("C2", &C2, "C2/D");
    tp.Branch("N1", &N1, "N1/D");
    tp.Branch("N2", &N2, "N2/D");
    tp.Branch("Y1", &Y1, "Y1/D");
    tp.Branch("Y2", &Y2, "Y2/D");
    tp.Branch("yield_ratio_12", &yield_ratio_12, "yield_ratio_12/D");
    tp.Branch("rel_yield_diff_12", &rel_yield_diff_12, "rel_yield_diff_12/D");
    tp.Branch("R2Q_screen", &R2Q_screen, "R2Q_screen/D");
    tp.Branch("R2P_screen", &R2P_screen, "R2P_screen/D");
    tp.Branch("raw_yield_ratio_12", &raw_yield_ratio_12, "raw_yield_ratio_12/D");
    tp.Branch("inclusive_yield_ratio_12", &inclusive_yield_ratio_12, "inclusive_yield_ratio_12/D");
    tp.Branch("theta_survival_ratio_12", &theta_survival_ratio_12, "theta_survival_ratio_12/D");
    tp.Branch("lv_lw_survival_ratio_12", &lv_lw_survival_ratio_12, "lv_lw_survival_ratio_12/D");
    tp.Branch("ecin_pcal_survival_ratio_12", &ecin_pcal_survival_ratio_12, "ecin_pcal_survival_ratio_12/D");
    tp.Branch("dominant_trigger_bit_all_1", &dominant_trigger_bit_all_1, "dominant_trigger_bit_all_1/I");
    tp.Branch("dominant_trigger_bit_all_2", &dominant_trigger_bit_all_2, "dominant_trigger_bit_all_2/I");
    tp.Branch("dominant_trigger_bit_selected_1", &dominant_trigger_bit_selected_1, "dominant_trigger_bit_selected_1/I");
    tp.Branch("dominant_trigger_bit_selected_2", &dominant_trigger_bit_selected_2, "dominant_trigger_bit_selected_2/I");
    tp.Branch("same_dominant_trigger_bit_all", &same_dominant_trigger_bit_all, "same_dominant_trigger_bit_all/I");
    tp.Branch("same_dominant_trigger_bit_selected", &same_dominant_trigger_bit_selected, "same_dominant_trigger_bit_selected/I");
    tp.Branch("bit29_pair_flag", &bit29_pair_flag, "bit29_pair_flag/I");
    tp.Branch("pair_quality_score", &pair_quality_score, "pair_quality_score/D");

    for(size_t i = 0; i < pairs.size(); ++i){
        const PairSummary &p = pairs[i];
        rank = static_cast<Int_t>(i + 1);
        run1 = p.run1; run2 = p.run2;
        P1 = p.P1; P2 = p.P2; Q1 = p.Q1; Q2 = p.Q2; deltaQ = p.deltaQ; abs_deltaQ = p.abs_deltaQ;
        target_pair_class = p.target_pair_class; target_pair_class_name = p.target_pair_class_name;
        C1 = p.C1; C2 = p.C2; N1 = p.N1; N2 = p.N2; Y1 = p.Y1; Y2 = p.Y2;
        yield_ratio_12 = p.yield_ratio_12; rel_yield_diff_12 = p.rel_yield_diff_12;
        R2Q_screen = p.R2Q_screen; R2P_screen = p.R2P_screen;
        raw_yield_ratio_12 = p.raw_yield_ratio_12; inclusive_yield_ratio_12 = p.inclusive_yield_ratio_12;
        theta_survival_ratio_12 = p.theta_survival_ratio_12;
        lv_lw_survival_ratio_12 = p.lv_lw_survival_ratio_12;
        ecin_pcal_survival_ratio_12 = p.ecin_pcal_survival_ratio_12;
        dominant_trigger_bit_all_1 = p.dominant_trigger_bit_all_1;
        dominant_trigger_bit_all_2 = p.dominant_trigger_bit_all_2;
        dominant_trigger_bit_selected_1 = p.dominant_trigger_bit_selected_1;
        dominant_trigger_bit_selected_2 = p.dominant_trigger_bit_selected_2;
        same_dominant_trigger_bit_all = p.same_dominant_trigger_bit_all;
        same_dominant_trigger_bit_selected = p.same_dominant_trigger_bit_selected;
        bit29_pair_flag = p.bit29_pair_flag;
        pair_quality_score = p.pair_quality_score;
        tp.Fill();
    }

    tp.Write();
    fout.Close();
}

void printTopPairs(const vector<PairSummary> &pairs, const string &className, int maxRows=10){
    cout << endl;
    cout << "Top " << maxRows << " pairs for class " << className << ":" << endl;
    cout << "----------------------------------------------------------------------------------------------------------------" << endl;
    cout << setw(6) << "rank"
         << setw(8) << "run1"
         << setw(8) << "run2"
         << setw(10) << "class"
         << setw(12) << "abs_dQ"
         << setw(12) << "Y1/Y2"
         << setw(12) << "R2Q"
         << setw(12) << "rawRat"
         << setw(12) << "thetaRat"
         << setw(8) << "trig"
         << setw(12) << "score" << endl;
    cout << "----------------------------------------------------------------------------------------------------------------" << endl;

    int shown = 0;
    for(size_t i = 0; i < pairs.size(); ++i){
        const PairSummary &p = pairs[i];
        if(className != "all" && p.target_pair_class_name != className) continue;
        shown++;
        cout << setw(6) << (i + 1)
             << setw(8) << p.run1
             << setw(8) << p.run2
             << setw(10) << p.target_pair_class_name
             << setw(12) << p.abs_deltaQ
             << setw(12) << p.yield_ratio_12
             << setw(12) << p.R2Q_screen
             << setw(12) << p.raw_yield_ratio_12
             << setw(12) << p.theta_survival_ratio_12
             << setw(8) << p.same_dominant_trigger_bit_all
             << setw(12) << p.pair_quality_score << endl;
        if(shown >= maxRows) break;
    }

    if(shown == 0){
        cout << "  No pairs found." << endl;
    }
}

int main(int argc, char **argv){
    Options opt = parseOptions(argc, argv);

    cout << endl;
    cout << "============================================================" << endl;
    cout << "  QE many-run pair scanner" << endl;
    cout << "  Version: " << QE_PAIR_SCAN_VERSION << endl;
    cout << "============================================================" << endl;
    cout << "Input dir:           " << opt.indir << endl;
    cout << "Filename contains:   " << opt.contains << endl;
    cout << "Run filter:          " << opt.runs << endl;
    cout << "Output ROOT:         " << opt.output_root << endl;
    cout << "Run CSV:             " << opt.run_csv << endl;
    cout << "Pair CSV:            " << opt.pair_csv << endl;
    cout << "Cut profile:         " << opt.cut_profile << endl;
    cout << "Apply kin cuts:      " << (opt.apply_kinematic_cuts ? 1 : 0) << endl;
    cout << "Charge branch:       " << opt.charge_branch << endl;
    cout << "Require is_qe:       " << (opt.require_is_qe ? 1 : 0) << endl;
    cout << "theta range [deg]:   " << opt.theta_min_deg << " to " << opt.theta_max_deg << endl;
    cout << "Q2 range [GeV2]:     " << opt.q2_min << " to " << opt.q2_max << endl;
    cout << "W range [GeV]:       " << opt.w_min << " to " << opt.w_max << endl;
    cout << "xB range:            " << opt.xmin << " to " << opt.xmax << endl;
    cout << "============================================================" << endl;

    vector<string> files = findInputFiles(opt);
    cout << endl;
    cout << "Found " << files.size() << " input ROOT files." << endl;

    if(files.empty()){
        cerr << "Error: no input files found." << endl;
        return 1;
    }

    vector<RunSummary> runs;
    runs.reserve(files.size());

    for(size_t i = 0; i < files.size(); ++i){
        cout << "[" << (i + 1) << " / " << files.size() << "] " << files[i] << endl;
        RunSummary r = processRunFile(files[i], opt);
        runs.push_back(r);

        cout << "  run " << r.run
             << " P=" << r.P
             << " Q=" << r.Q
             << " charge=" << r.charge.q_total
             << " inclusive=" << r.n_inclusive_entries
             << " selected=" << r.n_final
             << " Y=" << r.final_yield
             << " rawY=" << r.raw_yield
             << " trigAll=" << r.dominant_trigger_bit_all
             << " usable=" << r.usable << endl;
    }

    vector<PairSummary> pairs = buildPairs(runs, opt);

    cout << endl;
    cout << "Built " << pairs.size() << " candidate pairs." << endl;

    writeRunCSV(opt.run_csv, runs, opt);
    writePairCSV(opt.pair_csv, pairs, opt);
    writeRootOutput(opt.output_root, runs, pairs);

    cout << endl;
    cout << "Wrote outputs:" << endl;
    cout << "  " << opt.output_root << endl;
    cout << "  " << opt.run_csv << endl;
    cout << "  " << opt.pair_csv << endl;

    printTopPairs(pairs, "all", 10);
    printTopPairs(pairs, "pos_pos", 10);
    printTopPairs(pairs, "neg_neg", 10);
    printTopPairs(pairs, "pos_neg", 10);

    cout << endl;
    cout << "Done." << endl;
    return 0;
}
