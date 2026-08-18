//******************************************************************
//*  qe_analysis.cc
//*
//*  Purpose:
//*    Two-run inclusive quasi-elastic tensor-asymmetry extraction
//*    for RGC ND3 data using QE skim ROOT files.
//*
//*  Method:
//*    For two selected runs with different tensor polarizations Q1 and Q2,
//*    count inclusive QE electrons in xB bins, normalize each run by its
//*    integrated gated FCup charge, and calculate
//*
//*        R_2^Q(xB) = (Y1 - Y2) / (Q1*Y2 - Q2*Y1)
//*
//*    with
//*
//*        Yj = Nj / Cj,
//*
//*    where Nj is the event count in a given xB bin and Cj is the
//*    run-integrated FCup gated charge.
//*
//*  Tensor polarization:
//*    If Q1/Q2 are not provided directly, they are calculated from the
//*    vector polarization P using the spin-temperature relation
//*
//*        Q = 2 - sqrt(4 - 3*P^2)
//*
//*    with P supplied from --p1/--p2 or read from the skim vector_pol branch.
//*
//*  Error treatment:
//*    Event-counting uncertainty:       sigma_N = sqrt(N)
//*    FCup relative uncertainty:        sigma_C/C = --charge-rel-err
//*    Yield uncertainty:
//*        sigma_Y_stat   = sqrt(N)/C
//*        sigma_Y_charge = Y * charge_rel_err
//*        sigma_Y_total  = sqrt(sigma_Y_stat^2 + sigma_Y_charge^2)
//*
//*    R error is propagated analytically from Y1, Y2, Q1, and Q2.
//*
//*  Required input ROOT trees:
//*    inclusive
//*    scaler
//*
//*  Main output:
//*    qe_tensor_bin_summary
//*    qe_tensor_run_summary
//*    qe_tensor_cutflow
//*    h_R2Q_vs_xB_total
//*    h_Azz_approx_vs_xB_total, where Azz_approx = R_2^Q
//*    h_Azz_times2_diagnostic_vs_xB_total, where Azz_times2_diagnostic = 2*R_2^Q
//*    g_R2Q_vs_xB_total
//*    g_Azz_approx_vs_xB_total
//*    g_Azz_times2_diagnostic_vs_xB_total
//*    qe_tensor_trigger_bit_summary
//*    qe_tensor_input_file_summary
//******************************************************************

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <set>
#include <limits>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TSystem.h>
#include <TString.h>

using namespace std;

const double PI = acos(-1.0);
const double RAD_TO_DEG = 180.0 / PI;
const string QE_ANALYSIS_VERSION = "qe_analysis_diagnostics_v5_cut_profile_2026_08_13";

// Hardcoded plot y-axis range. Edit these constants directly when you want
// a different presentation range. No command-line option is used.
const bool USE_HARDCODED_Y_RANGE = true;
const double R2Q_YMIN = -5.0;
const double R2Q_YMAX =  5.0;
const bool SCALE_TIMES2_DIAGNOSTIC_Y_RANGE = true;

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

string removeRootSuffix(const string &s){
    const string suffix = ".root";
    if(s.size() >= suffix.size() && s.substr(s.size() - suffix.size()) == suffix){
        return s.substr(0, s.size() - suffix.size());
    }
    return s;
}

string csvPrefixFromRootOutput(const string &s){
    return removeRootSuffix(s);
}

double safeRatio(double num, double den, double bad=-9999.0){
    if(!isFinite(num) || !isFinite(den) || fabs(den) <= 1.0e-30) return bad;
    return num / den;
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

void addTriggerBitsToMap(map<int, Long64_t> &m, int trigger){
    // Decode the trigger word into individual set bits.
    // The skim branch is Int_t; we use bits 0..30, which includes high
    // positive bits such as bit 29 seen in the 17792 trigger words.
    const unsigned int trig = static_cast<unsigned int>(trigger);
    for(int bit = 0; bit <= 30; ++bit){
        const unsigned int mask = (1u << bit);
        if((trig & mask) != 0u){
            m[bit]++;
        }
    }
}

string targetPairClassName(int cls){
    if(cls == 1) return "pos_pos";
    if(cls == -1) return "neg_neg";
    if(cls == 2) return "pos_neg";
    return "unknown";
}

void updateMinMaxInt(bool &seen, Int_t v, Int_t &vmin, Int_t &vmax){
    if(!seen){
        vmin = v;
        vmax = v;
        seen = true;
        return;
    }
    if(v < vmin) vmin = v;
    if(v > vmax) vmax = v;
}

void updateMinMaxLong(bool &seen, Long64_t v, Long64_t &vmin, Long64_t &vmax){
    if(!seen){
        vmin = v;
        vmax = v;
        seen = true;
        return;
    }
    if(v < vmin) vmin = v;
    if(v > vmax) vmax = v;
}

//============================================================
// Command-line configuration
//============================================================

struct Options {
    string file1 = "";
    string file2 = "";
    string output = "qe_tensor_R.root";

    // Cut-profile convention used in this analysis workflow:
    //   kin0 = inclusive QE detector/PID cuts only.
    //   kin1 = inclusive QE detector/PID cuts plus tight QE kinematic cuts.
    // The xB range is always the plotting/counting range, not part of the
    // skim kincut definition.
    string cut_profile = "kin1";
    bool apply_kinematic_cuts = true;

    int run1 = 0;
    int run2 = 0;

    // Direct tensor-polarization input. If not provided, Q is derived from P.
    bool q1_user = false;
    bool q2_user = false;
    double q1 = -9999.0;
    double q2 = -9999.0;

    bool qerr1_user = false;
    bool qerr2_user = false;
    double qerr1 = 0.0;
    double qerr2 = 0.0;

    // Optional vector-polarization input. Used to derive Q if Q not supplied.
    bool p1_user = false;
    bool p2_user = false;
    double p1 = -9999.0;
    double p2 = -9999.0;

    bool perr1_user = false;
    bool perr2_user = false;
    double perr1 = 0.0;
    double perr2 = 0.0;

    int xbins = 120;
    double xmin = 0.8;
    double xmax = 1.6;

    // Default cuts reflect the narrow QE check used in the slides.
    // They can be loosened from the command line.
    double theta_min_deg = 7.8;
    double theta_max_deg = 8.2;

    double vz_min = -5.758;
    double vz_max = 1.5165;

    double chi2_abs_max = 3.0;

    // Detector/PID-level momentum requirement used only to make calorimeter
    // ratios well-defined. The tight pe > 2 GeV cut is pe_min below and is
    // applied only for cut_profile=kin1.
    double detector_pe_min = 0.0;
    double pe_min = 2.0;

    double lv_min = 14.0;
    double lw_min = 14.0;

    double pcal_min = 0.07;
    double sf_max = 0.28;

    bool use_sf_cut = true;
    bool use_ecin_pcal_cut = true;

    // ECIN/PCAL diagonal cut:
    //   ECIN/p >= slope*(PCAL/p) + intercept
    double ecin_pcal_slope = -0.625;
    double ecin_pcal_intercept = 0.15;

    double q2_min = 1.9433;
    double q2_max = 2.0574;

    double w_min = 0.0;
    double w_max = 1.073;

    double charge_rel_err = 0.002;
    string charge_branch = "fcupgated_delta";

    bool require_target_match = true;
    bool require_is_qe = true;
    bool require_fd_status = false;

    bool save_plots = true;

    // Diagnostic output is on by default. It writes extra ROOT trees and CSV files.
    bool write_diagnostics = true;
};

void printUsage(const char *prog){
    cout << endl;
    cout << "Usage:" << endl;
    cout << "  " << prog << " --file1 run_highQ.root --file2 run_lowQ.root --run1 RUN1 --run2 RUN2 [options]" << endl;
    cout << endl;
    cout << "Required:" << endl;
    cout << "  --file1 <root>             QE skim ROOT file for run 1" << endl;
    cout << "  --file2 <root>             QE skim ROOT file for run 2" << endl;
    cout << "  --run1 <run>               Run number for file1" << endl;
    cout << "  --run2 <run>               Run number for file2" << endl;
    cout << "  --cut-profile <kin0|kin1>  kin0=detector/PID only, kin1=detector/PID + tight kinematics. Default kin1" << endl;
    cout << endl;
    cout << "Tensor/vector polarization options:" << endl;
    cout << "  --q1 <Q1>                  Tensor polarization for run1, as fraction, e.g. 0.2074" << endl;
    cout << "  --q2 <Q2>                  Tensor polarization for run2, as fraction, e.g. 0.0525" << endl;
    cout << "  --qerr1 <dQ1>              Absolute tensor-polarization uncertainty" << endl;
    cout << "  --qerr2 <dQ2>              Absolute tensor-polarization uncertainty" << endl;
    cout << "  --p1 <P1>                  Vector polarization for run1. Q is calculated from P." << endl;
    cout << "  --p2 <P2>                  Vector polarization for run2. Q is calculated from P." << endl;
    cout << "  --perr1 <dP1>              Vector-polarization uncertainty for run1" << endl;
    cout << "  --perr2 <dP2>              Vector-polarization uncertainty for run2" << endl;
    cout << endl;
    cout << "Binning:" << endl;
    cout << "  --xbins <N>                Default 120" << endl;
    cout << "  --xmin <x>                 Default 0.8" << endl;
    cout << "  --xmax <x>                 Default 1.6" << endl;
    cout << endl;
    cout << "Analysis cuts:" << endl;
    cout << "  --theta-min <deg>          Default 7.8" << endl;
    cout << "  --theta-max <deg>          Default 8.2" << endl;
    cout << "  --vz-min <cm>              Default -5.758" << endl;
    cout << "  --vz-max <cm>              Default 1.5165" << endl;
    cout << "  --chi2-max <value>         Default 3.0, applied as detector/PID abs(e_chi2pid) < chi2-max" << endl;
    cout << "  --detector-pe-min <GeV>    Default 0.0, detector/PID pe threshold for calorimeter ratios" << endl;
    cout << "  --pe-min <GeV>             Default 2.0, tight kinematic cut applied only for cut-profile kin1" << endl;
    cout << "  --lv-min <value>           Default 14.0" << endl;
    cout << "  --lw-min <value>           Default 14.0" << endl;
    cout << "  --pcal-min <GeV>           Default 0.07" << endl;
    cout << "  --sf-max <value>           Default 0.28" << endl;
    cout << "  --use-sf-cut 0|1           Default 1" << endl;
    cout << "  --use-ecin-pcal-cut 0|1    Default 1" << endl;
    cout << "  --ecin-slope <value>       Default -0.625" << endl;
    cout << "  --ecin-intercept <value>   Default 0.15" << endl;
    cout << "  --q2-min <GeV2>            Default 1.9433" << endl;
    cout << "  --q2-max <GeV2>            Default 2.0574" << endl;
    cout << "  --w-min <GeV>              Default 0.0" << endl;
    cout << "  --w-max <GeV>              Default 1.073" << endl;
    cout << endl;
    cout << "Normalization and other options:" << endl;
    cout << "  --charge-branch <name>     Default fcupgated_delta" << endl;
    cout << "  --charge-rel-err <value>   Default 0.002" << endl;
    cout << "  --require-target-match 0|1 Default 1" << endl;
    cout << "  --require-is-qe 0|1        Default 1 if branch exists" << endl;
    cout << "  --require-fd-status 0|1    Default 0. If 1, require -4000 < e_status <= -2000" << endl;
    cout << "  --save-plots 0|1           Default 1" << endl;
    cout << "  --write-diagnostics 0|1    Default 1. Write pair/run/file/trigger diagnostic trees and CSVs" << endl;
    cout << "  --output <root>            Default qe_tensor_R.root" << endl;
    cout << endl;
    cout << "Example:" << endl;
    cout << "  " << prog << " \\" << endl;
    cout << "    --file1 rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1kin1tf1/spring23/sidisdvcs_inclusive_qe_eT0pid1kin1tf1_017492.root \\" << endl;
    cout << "    --file2 rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1kin1tf1/spring23/sidisdvcs_inclusive_qe_eT0pid1kin1tf1_017792.root \\" << endl;
    cout << "    --run1 17492 --run2 17792 \\" << endl;
    cout << "    --p1 0.5120999 --p2 -0.262899 \\" << endl;
    cout << "    --output qe_tensor_R_17492_17792.root" << endl;
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

        if(a == "--file1"){
            needValue(a); opt.file1 = argv[++i];
        }else if(a == "--file2"){
            needValue(a); opt.file2 = argv[++i];
        }else if(a == "--run1"){
            needValue(a); opt.run1 = atoi(argv[++i]);
        }else if(a == "--run2"){
            needValue(a); opt.run2 = atoi(argv[++i]);
        }else if(a == "--cut-profile"){
            needValue(a); opt.cut_profile = argv[++i];
        }else if(a == "--output" || a == "-o"){
            needValue(a); opt.output = argv[++i];
        }else if(a == "--q1"){
            needValue(a); opt.q1 = atof(argv[++i]); opt.q1_user = true;
        }else if(a == "--q2"){
            needValue(a); opt.q2 = atof(argv[++i]); opt.q2_user = true;
        }else if(a == "--qerr1"){
            needValue(a); opt.qerr1 = atof(argv[++i]); opt.qerr1_user = true;
        }else if(a == "--qerr2"){
            needValue(a); opt.qerr2 = atof(argv[++i]); opt.qerr2_user = true;
        }else if(a == "--p1"){
            needValue(a); opt.p1 = atof(argv[++i]); opt.p1_user = true;
        }else if(a == "--p2"){
            needValue(a); opt.p2 = atof(argv[++i]); opt.p2_user = true;
        }else if(a == "--perr1"){
            needValue(a); opt.perr1 = atof(argv[++i]); opt.perr1_user = true;
        }else if(a == "--perr2"){
            needValue(a); opt.perr2 = atof(argv[++i]); opt.perr2_user = true;
        }else if(a == "--xbins"){
            needValue(a); opt.xbins = atoi(argv[++i]);
        }else if(a == "--xmin"){
            needValue(a); opt.xmin = atof(argv[++i]);
        }else if(a == "--xmax"){
            needValue(a); opt.xmax = atof(argv[++i]);
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
        }else if(a == "--charge-rel-err"){
            needValue(a); opt.charge_rel_err = atof(argv[++i]);
        }else if(a == "--charge-branch"){
            needValue(a); opt.charge_branch = argv[++i];
        }else if(a == "--require-target-match"){
            needValue(a); opt.require_target_match = strToBool01(argv[++i]);
        }else if(a == "--require-is-qe"){
            needValue(a); opt.require_is_qe = strToBool01(argv[++i]);
        }else if(a == "--require-fd-status"){
            needValue(a); opt.require_fd_status = strToBool01(argv[++i]);
        }else if(a == "--save-plots"){
            needValue(a); opt.save_plots = strToBool01(argv[++i]);
        }else if(a == "--write-diagnostics"){
            needValue(a); opt.write_diagnostics = strToBool01(argv[++i]);
        }else if(a == "--help" || a == "-h"){
            printUsage(argv[0]);
            exit(0);
        }else{
            cerr << "Error: unknown option: " << a << endl;
            printUsage(argv[0]);
            exit(1);
        }
    }

    if(opt.file1.empty() || opt.file2.empty()){
        cerr << "Error: --file1 and --file2 are required." << endl;
        exit(1);
    }

    if(opt.run1 <= 0 || opt.run2 <= 0){
        cerr << "Error: --run1 and --run2 must be positive run numbers." << endl;
        exit(1);
    }

    if(opt.xbins <= 0 || opt.xmax <= opt.xmin){
        cerr << "Error: invalid xB binning." << endl;
        exit(1);
    }

    if(opt.charge_rel_err < 0.0){
        cerr << "Error: --charge-rel-err must be non-negative." << endl;
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

    Int_t run_config = -9999;
    Int_t event_config = -9999;
    Long64_t trigger = -9999;
    Long64_t timestamp = -9999;
    Long64_t unixtime = -9999;

    Int_t helicity = -9999;
    Int_t helicity_valid = 0;

    Int_t target_match = 0;
    Float_t target_pol = -9999.0;
    Float_t target_pol_err = -9999.0;
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
    Float_t phie = -9999.0;

    Float_t vze = -9999.0;

    Float_t e_chi2pid = -9999.0;
    Float_t e_pcal_e = -9999.0;
    Float_t e_ecin_e = -9999.0;
    Float_t e_ecout_e = -9999.0;
    Float_t e_cal_e = -9999.0;
    Float_t e_sampling_fraction = -9999.0;
    Float_t e_lu = -9999.0;
    Float_t e_lv = -9999.0;
    Float_t e_lw = -9999.0;

    Float_t Q2 = -9999.0;
    Float_t nu = -9999.0;
    Float_t W = -9999.0;
    Float_t xB = -9999.0;
    Float_t y = -9999.0;

    Int_t is_qe = -9999;
};

struct BranchAvailability {
    bool is_qe = false;
    bool e_sampling_fraction = false;
    bool e_cal_e = false;
    bool vector_pol_err = false;
    bool target_pol_offline_err = false;

    bool event_file_index = false;
    bool input_file_index = false;
    bool input_file_number = false;
    bool run_config = false;
    bool event_config = false;
    bool trigger = false;
    bool timestamp = false;
    bool unixtime = false;
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

    cerr << "Warning: optional branch not found: " << name
         << ". Default value will be used." << endl;
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

    avail.run_config = setBranchIfExists(tree, "run_config", &e.run_config, opt);
    avail.event_config = setBranchIfExists(tree, "event_config", &e.event_config, opt);
    avail.trigger = setBranchIfExists(tree, "trigger", &e.trigger, opt);
    avail.timestamp = setBranchIfExists(tree, "timestamp", &e.timestamp, opt);
    avail.unixtime = setBranchIfExists(tree, "unixtime", &e.unixtime, opt);

    setBranchIfExists(tree, "helicity", &e.helicity, opt);
    setBranchIfExists(tree, "helicity_valid", &e.helicity_valid, opt);

    setBranchIfExists(tree, "target_match", &e.target_match, opt);
    setBranchIfExists(tree, "target_pol", &e.target_pol, opt);
    setBranchIfExists(tree, "target_pol_err", &e.target_pol_err, opt);
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
    setBranchIfExists(tree, "phie", &e.phie, opt);
    setBranchIfExists(tree, "vze", &e.vze, req);

    setBranchIfExists(tree, "e_chi2pid", &e.e_chi2pid, req);
    setBranchIfExists(tree, "e_pcal_e", &e.e_pcal_e, req);
    setBranchIfExists(tree, "e_ecin_e", &e.e_ecin_e, req);
    setBranchIfExists(tree, "e_ecout_e", &e.e_ecout_e, opt);
    avail.e_cal_e = setBranchIfExists(tree, "e_cal_e", &e.e_cal_e, opt);
    avail.e_sampling_fraction = setBranchIfExists(tree, "e_sampling_fraction", &e.e_sampling_fraction, opt);
    setBranchIfExists(tree, "e_lu", &e.e_lu, opt);
    setBranchIfExists(tree, "e_lv", &e.e_lv, req);
    setBranchIfExists(tree, "e_lw", &e.e_lw, req);

    setBranchIfExists(tree, "Q2", &e.Q2, req);
    setBranchIfExists(tree, "nu", &e.nu, opt);
    setBranchIfExists(tree, "W", &e.W, req);
    setBranchIfExists(tree, "xB", &e.xB, req);
    setBranchIfExists(tree, "y", &e.y, opt);

    avail.is_qe = setBranchIfExists(tree, "is_qe", &e.is_qe, opt);

    return avail;
}

int chooseInputFileKey(const EventBranches &e, const BranchAvailability &avail){
    // Prefer a true HIPO file number if it exists, then the command-level
    // input file index, and only fall back to the old local event counter.
    if(avail.input_file_number && e.input_file_number > -9000){
        return e.input_file_number;
    }

    if(avail.input_file_index && e.input_file_index > -9000){
        return e.input_file_index;
    }

    if(avail.event_file_index && e.event_file_index > -9000){
        return e.event_file_index;
    }

    return -9999;
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

struct BinCounts {
    Long64_t total = 0;
    Long64_t plus = 0;
    Long64_t minus = 0;
    Long64_t zero = 0;
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

struct RunResult {
    string file = "";
    int run = 0;
    Long64_t n_inclusive_entries = 0;

    vector<BinCounts> bins;
    CutFlow cutflow;
    ChargeSummary charge;

    Long64_t final_plus = 0;
    Long64_t final_minus = 0;
    Long64_t final_zero = 0;
    Long64_t final_total = 0;

    double mean_vector_pol = -9999.0;
    double mean_vector_pol_err = -9999.0;
    double mean_target_pol_offline = -9999.0;
    double mean_target_pol_offline_err = -9999.0;
    double mean_tensor_pol_branch = -9999.0;

    bool has_is_qe = false;

    // Diagnostic bookkeeping. These are filled before and after analysis cuts.
    map<int, Long64_t> all_by_file_index;
    map<int, Long64_t> selected_by_file_index;
    map<int, Long64_t> all_by_input_file;
    map<int, Long64_t> selected_by_input_file;
    map<int, Long64_t> all_by_trigger;
    map<int, Long64_t> selected_by_trigger;
    map<int, Long64_t> all_by_trigger_bit;
    map<int, Long64_t> selected_by_trigger_bit;
    map<int, Long64_t> all_by_run_config;
    map<int, Long64_t> selected_by_run_config;
    map<int, Long64_t> all_by_event_config;
    map<int, Long64_t> selected_by_event_config;

    // Optional run_summary tree information from the skim.
    bool has_skim_run_summary = false;
    Int_t skim_run_summary_run = -9999;
    Int_t skim_input_file_index = -9999;
    Int_t skim_input_file_number = -9999;
    Long64_t skim_n_events_all = 0;
    Long64_t skim_n_inclusive_written = 0;
    Long64_t skim_n_sidis_rows_written = 0;
    Long64_t skim_n_dihadron_rows_written = 0;
    Long64_t skim_n_scaler_rows = 0;
    Int_t skim_first_event = -9999;
    Int_t skim_last_event = -9999;
    Long64_t skim_first_timestamp = -9999;
    Long64_t skim_last_timestamp = -9999;

    Int_t event_min_all = 0;
    Int_t event_max_all = 0;
    Int_t event_min_selected = 0;
    Int_t event_max_selected = 0;
    bool seen_event_all = false;
    bool seen_event_selected = false;

    Long64_t timestamp_min_all = 0;
    Long64_t timestamp_max_all = 0;
    Long64_t timestamp_min_selected = 0;
    Long64_t timestamp_max_selected = 0;
    bool seen_timestamp_all = false;
    bool seen_timestamp_selected = false;

    Long64_t unixtime_min_all = 0;
    Long64_t unixtime_max_all = 0;
    Long64_t unixtime_min_selected = 0;
    Long64_t unixtime_max_selected = 0;
    bool seen_unixtime_all = false;
    bool seen_unixtime_selected = false;
};

bool passFDStatus(Int_t status){
    return (status > -4000 && status <= -2000);
}

bool getQESamplingFraction(const EventBranches &e, double &sf){
    // Inclusive-QE detector/PID definition used by the skim:
    //   SF_QE = (E_PCAL + E_ECIN) / pe
    // This intentionally does NOT use the generic e_sampling_fraction branch
    // based on total calorimeter energy.
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

    // Detector/PID cut block. This block is active for both kin0 and kin1
    // profiles. It mirrors the inclusive QE detpidcut=1 definition used in
    // the skim as closely as the available skim branches allow.
    if(fabs(e.e_chi2pid) >= opt.chi2_abs_max) return false;
    cf.n_chi2++;

    if(e.e_lv <= opt.lv_min || e.e_lw <= opt.lw_min) return false;
    cf.n_lv_lw++;

    // Detector-level pe requirement used only so PCAL/ECIN ratios are defined.
    // The tight pe > 2 GeV kinematic cut is applied below only for kin1.
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

    // Tight QE kinematic block. This is active only for cut_profile=kin1.
    // For cut_profile=kin0, these counters are advanced as pass-through
    // stages so the cumulative cut-flow remains readable.
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

    // xB is treated as the plotting/counting interval, not as part of the
    // skim kincut definition.
    if(e.xB < opt.xmin || e.xB >= opt.xmax) return false;
    cf.n_xb_range++;

    return true;
}

void readSkimRunSummary(TFile &fin, RunResult &result){
    TTree *rs = dynamic_cast<TTree*>(fin.Get("run_summary"));
    if(!rs){
        return;
    }

    Int_t br_run = -9999;
    Int_t br_input_file_index = -9999;
    Int_t br_input_file_number = -9999;
    Long64_t br_n_events_all = 0;
    Long64_t br_n_inclusive_written = 0;
    Long64_t br_n_sidis_rows_written = 0;
    Long64_t br_n_dihadron_rows_written = 0;
    Long64_t br_n_scaler_rows = 0;
    Int_t br_first_event = -9999;
    Int_t br_last_event = -9999;
    Long64_t br_first_timestamp = -9999;
    Long64_t br_last_timestamp = -9999;

    if(rs->GetBranch("run")) rs->SetBranchAddress("run", &br_run);
    if(rs->GetBranch("input_file_index")) rs->SetBranchAddress("input_file_index", &br_input_file_index);
    if(rs->GetBranch("input_file_number")) rs->SetBranchAddress("input_file_number", &br_input_file_number);
    if(rs->GetBranch("n_events_all")) rs->SetBranchAddress("n_events_all", &br_n_events_all);
    if(rs->GetBranch("n_inclusive_written")) rs->SetBranchAddress("n_inclusive_written", &br_n_inclusive_written);
    if(rs->GetBranch("n_sidis_rows_written")) rs->SetBranchAddress("n_sidis_rows_written", &br_n_sidis_rows_written);
    if(rs->GetBranch("n_dihadron_rows_written")) rs->SetBranchAddress("n_dihadron_rows_written", &br_n_dihadron_rows_written);
    if(rs->GetBranch("n_scaler_rows")) rs->SetBranchAddress("n_scaler_rows", &br_n_scaler_rows);
    if(rs->GetBranch("first_event")) rs->SetBranchAddress("first_event", &br_first_event);
    if(rs->GetBranch("last_event")) rs->SetBranchAddress("last_event", &br_last_event);
    if(rs->GetBranch("first_timestamp")) rs->SetBranchAddress("first_timestamp", &br_first_timestamp);
    if(rs->GetBranch("last_timestamp")) rs->SetBranchAddress("last_timestamp", &br_last_timestamp);

    const Long64_t n = rs->GetEntries();
    Long64_t n_used = 0;

    for(Long64_t i = 0; i < n; ++i){
        // Reset per-entry values so missing branches do not retain stale data.
        br_run = -9999;
        br_input_file_index = -9999;
        br_input_file_number = -9999;
        br_n_events_all = 0;
        br_n_inclusive_written = 0;
        br_n_sidis_rows_written = 0;
        br_n_dihadron_rows_written = 0;
        br_n_scaler_rows = 0;
        br_first_event = -9999;
        br_last_event = -9999;
        br_first_timestamp = -9999;
        br_last_timestamp = -9999;

        rs->GetEntry(i);

        // Some current skim files have run_summary.run = 0 because the last
        // raw event in the file had a missing RUN::config. If this ROOT file is
        // a single-run file, still use the row.
        const bool use_row = (br_run == result.run) || (br_run == 0) || (br_run < 0) || (n == 1);
        if(!use_row) continue;

        if(n_used == 0){
            result.skim_run_summary_run = br_run;
            result.skim_input_file_index = br_input_file_index;
            result.skim_input_file_number = br_input_file_number;
            result.skim_first_event = br_first_event;
            result.skim_last_event = br_last_event;
            result.skim_first_timestamp = br_first_timestamp;
            result.skim_last_timestamp = br_last_timestamp;
        }else{
            if(br_input_file_index > -9000 && result.skim_input_file_index <= -9000){
                result.skim_input_file_index = br_input_file_index;
            }
            if(br_input_file_number > -9000 && result.skim_input_file_number <= -9000){
                result.skim_input_file_number = br_input_file_number;
            }
            if(br_first_event >= 0 && (result.skim_first_event < 0 || br_first_event < result.skim_first_event)){
                result.skim_first_event = br_first_event;
            }
            if(br_last_event >= 0 && (result.skim_last_event < 0 || br_last_event > result.skim_last_event)){
                result.skim_last_event = br_last_event;
            }
            if(br_first_timestamp > 0 && (result.skim_first_timestamp < 0 || br_first_timestamp < result.skim_first_timestamp)){
                result.skim_first_timestamp = br_first_timestamp;
            }
            if(br_last_timestamp > 0 && (result.skim_last_timestamp < 0 || br_last_timestamp > result.skim_last_timestamp)){
                result.skim_last_timestamp = br_last_timestamp;
            }
        }

        result.skim_n_events_all += br_n_events_all;
        result.skim_n_inclusive_written += br_n_inclusive_written;
        result.skim_n_sidis_rows_written += br_n_sidis_rows_written;
        result.skim_n_dihadron_rows_written += br_n_dihadron_rows_written;
        result.skim_n_scaler_rows += br_n_scaler_rows;
        n_used++;
    }

    result.has_skim_run_summary = (n_used > 0);
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
             << chargeBranch << endl;
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

    if(charge.n_scaler_rows <= 0 || charge.q_total <= 0.0){
        cerr << "Error: no usable scaler charge found for run " << run
             << " in file " << fileName << endl;
        exit(1);
    }

    return charge;
}

RunResult processRun(const string &fileName, int run, const Options &opt){
    RunResult result;
    result.file = fileName;
    result.run = run;
    result.bins.resize(opt.xbins);

    result.charge = readChargeSummary(fileName, run, opt.charge_branch);

    TFile fin(fileName.c_str(), "READ");
    if(fin.IsZombie()){
        cerr << "Error: could not open ROOT file: " << fileName << endl;
        exit(1);
    }

    readSkimRunSummary(fin, result);

    TTree *tree = dynamic_cast<TTree*>(fin.Get("inclusive"));
    if(!tree){
        cerr << "Error: inclusive tree not found in file: " << fileName << endl;
        exit(1);
    }

    EventBranches e;
    BranchAvailability avail = setupInclusiveBranches(tree, e);
    result.has_is_qe = avail.is_qe;

    const Long64_t nEntries = tree->GetEntries();
    result.n_inclusive_entries = nEntries;

    double sum_vector_pol = 0.0;
    double sum_vector_pol_err = 0.0;
    double sum_target_pol_offline = 0.0;
    double sum_target_pol_offline_err = 0.0;
    double sum_tensor_pol = 0.0;

    Long64_t n_vector_pol = 0;
    Long64_t n_vector_pol_err = 0;
    Long64_t n_target_pol = 0;
    Long64_t n_target_pol_err = 0;
    Long64_t n_tensor_pol = 0;

    const double dx = (opt.xmax - opt.xmin) / static_cast<double>(opt.xbins);

    for(Long64_t i = 0; i < nEntries; ++i){
        tree->GetEntry(i);

        if(e.run != run) continue;

        updateMinMaxInt(result.seen_event_all, e.event, result.event_min_all, result.event_max_all);
        if(avail.event_file_index && e.event_file_index > -9000){
            result.all_by_file_index[e.event_file_index]++;
        }
        {
            const int input_key = chooseInputFileKey(e, avail);
            if(input_key > -9000){
                result.all_by_input_file[input_key]++;
            }
        }
        if(avail.trigger){
            result.all_by_trigger[e.trigger]++;
            addTriggerBitsToMap(result.all_by_trigger_bit, e.trigger);
        }
        if(avail.run_config){
            result.all_by_run_config[e.run_config]++;
        }
        if(avail.event_config){
            result.all_by_event_config[e.event_config]++;
        }
        if(avail.timestamp){
            updateMinMaxLong(result.seen_timestamp_all, e.timestamp, result.timestamp_min_all, result.timestamp_max_all);
        }
        if(avail.unixtime){
            updateMinMaxLong(result.seen_unixtime_all, e.unixtime, result.unixtime_min_all, result.unixtime_max_all);
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
            sum_target_pol_offline += e.target_pol_offline;
            n_target_pol++;
        }

        if(e.target_pol_offline_err > -900.0 && isFinite(e.target_pol_offline_err)){
            sum_target_pol_offline_err += fabs(e.target_pol_offline_err);
            n_target_pol_err++;
        }

        if(e.tensor_pol > -900.0 && isFinite(e.tensor_pol)){
            sum_tensor_pol += e.tensor_pol;
            n_tensor_pol++;
        }

        if(!passEventCuts(e, avail, opt, result.cutflow)) continue;

        updateMinMaxInt(result.seen_event_selected, e.event, result.event_min_selected, result.event_max_selected);
        if(avail.event_file_index && e.event_file_index > -9000){
            result.selected_by_file_index[e.event_file_index]++;
        }
        {
            const int input_key = chooseInputFileKey(e, avail);
            if(input_key > -9000){
                result.selected_by_input_file[input_key]++;
            }
        }
        if(avail.trigger){
            result.selected_by_trigger[e.trigger]++;
            addTriggerBitsToMap(result.selected_by_trigger_bit, e.trigger);
        }
        if(avail.run_config){
            result.selected_by_run_config[e.run_config]++;
        }
        if(avail.event_config){
            result.selected_by_event_config[e.event_config]++;
        }
        if(avail.timestamp){
            updateMinMaxLong(result.seen_timestamp_selected, e.timestamp, result.timestamp_min_selected, result.timestamp_max_selected);
        }
        if(avail.unixtime){
            updateMinMaxLong(result.seen_unixtime_selected, e.unixtime, result.unixtime_min_selected, result.unixtime_max_selected);
        }

        int ibin = static_cast<int>(floor((e.xB - opt.xmin) / dx));
        if(ibin < 0 || ibin >= opt.xbins) continue;

        BinCounts &bc = result.bins[ibin];
        bc.total++;

        if(e.helicity == 1){
            bc.plus++;
            result.final_plus++;
        }else if(e.helicity == -1){
            bc.minus++;
            result.final_minus++;
        }else{
            bc.zero++;
            result.final_zero++;
        }

        result.final_total++;
    }

    if(n_vector_pol > 0) result.mean_vector_pol = sum_vector_pol / static_cast<double>(n_vector_pol);
    if(n_vector_pol_err > 0) result.mean_vector_pol_err = sum_vector_pol_err / static_cast<double>(n_vector_pol_err);
    if(n_target_pol > 0) result.mean_target_pol_offline = sum_target_pol_offline / static_cast<double>(n_target_pol);
    if(n_target_pol_err > 0) result.mean_target_pol_offline_err = sum_target_pol_offline_err / static_cast<double>(n_target_pol_err);
    if(n_tensor_pol > 0) result.mean_tensor_pol_branch = sum_tensor_pol / static_cast<double>(n_tensor_pol);

    return result;
}

//============================================================
// Printing helpers
//============================================================

void printCutFlow(const string &label, const CutFlow &cf){
    cout << endl;
    cout << "Cut-flow for " << label << ":" << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << left << setw(28) << "Step" << right << setw(18) << "Remaining" << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << left << setw(28) << "all inclusive entries" << right << setw(18) << cf.n_all << endl;
    cout << left << setw(28) << "+ target_match" << right << setw(18) << cf.n_target_match << endl;
    cout << left << setw(28) << "+ is_qe" << right << setw(18) << cf.n_is_qe << endl;
    cout << left << setw(28) << "+ pid_e == 11" << right << setw(18) << cf.n_pid << endl;
    cout << left << setw(28) << "+ status" << right << setw(18) << cf.n_status << endl;
    cout << left << setw(28) << "+ det chi2pid" << right << setw(18) << cf.n_chi2 << endl;
    cout << left << setw(28) << "+ det lv/lw" << right << setw(18) << cf.n_lv_lw << endl;
    cout << left << setw(28) << "+ det pcal" << right << setw(18) << cf.n_pcal << endl;
    cout << left << setw(28) << "+ det SF_QE" << right << setw(18) << cf.n_sf << endl;
    cout << left << setw(28) << "+ det ecin/pcal" << right << setw(18) << cf.n_ecin_pcal << endl;
    cout << left << setw(28) << "+ kin theta" << right << setw(18) << cf.n_theta << endl;
    cout << left << setw(28) << "+ kin vz" << right << setw(18) << cf.n_vz << endl;
    cout << left << setw(28) << "+ kin pe" << right << setw(18) << cf.n_pe << endl;
    cout << left << setw(28) << "+ kin Q2" << right << setw(18) << cf.n_q2 << endl;
    cout << left << setw(28) << "+ kin W" << right << setw(18) << cf.n_w << endl;
    cout << left << setw(28) << "+ xB range" << right << setw(18) << cf.n_xb_range << endl;
    cout << "------------------------------------------------------------" << endl;
}

void printRunSummary(const RunResult &r){
    cout << endl;
    cout << "Run summary for run " << r.run << endl;
    cout << "  file:                         " << r.file << endl;
    cout << "  inclusive entries:            " << r.n_inclusive_entries << endl;
    cout << "  final selected entries:        " << r.final_total << endl;
    cout << "    helicity +1:                 " << r.final_plus << endl;
    cout << "    helicity -1:                 " << r.final_minus << endl;
    cout << "    helicity  0/other:           " << r.final_zero << endl;
    cout << "  scaler rows:                  " << r.charge.n_scaler_rows << endl;
    cout << "  FCup gated charge total:       " << r.charge.q_total << endl;
    cout << "    charge +1:                   " << r.charge.q_plus << endl;
    cout << "    charge -1:                   " << r.charge.q_minus << endl;
    cout << "    charge  0/other:             " << r.charge.q_zero << endl;
    cout << "  mean vector_pol:               " << r.mean_vector_pol << endl;
    cout << "  mean vector_pol_err:           " << r.mean_vector_pol_err << endl;
    cout << "  mean target_pol_offline:       " << r.mean_target_pol_offline << endl;
    cout << "  mean target_pol_offline_err:   " << r.mean_target_pol_offline_err << endl;
    cout << "  mean tensor_pol branch:        " << r.mean_tensor_pol_branch << endl;
    cout << "  event_file_index all count:    " << r.all_by_file_index.size() << endl;
    cout << "  event_file_index selected cnt: " << r.selected_by_file_index.size() << endl;
    cout << "  input-file key all count:      " << r.all_by_input_file.size() << endl;
    cout << "  input-file key selected cnt:   " << r.selected_by_input_file.size() << endl;
    if(r.has_skim_run_summary){
        cout << "  skim run_summary run:          " << r.skim_run_summary_run << endl;
        cout << "  skim n_events_all:             " << r.skim_n_events_all << endl;
        cout << "  skim n_inclusive_written:      " << r.skim_n_inclusive_written << endl;
    }
    if(r.seen_event_all){
        cout << "  event range all:               " << r.event_min_all << " to " << r.event_max_all << endl;
    }
    if(r.seen_event_selected){
        cout << "  event range selected:          " << r.event_min_selected << " to " << r.event_max_selected << endl;
    }
}

//============================================================
// Resolve tensor-polarization values
//============================================================

struct TensorInputs {
    double P1 = -9999.0;
    double P2 = -9999.0;
    double Perr1 = 0.0;
    double Perr2 = 0.0;
    double Q1 = -9999.0;
    double Q2 = -9999.0;
    double Qerr1 = 0.0;
    double Qerr2 = 0.0;
};

TensorInputs resolveTensorInputs(const Options &opt, const RunResult &r1, const RunResult &r2){
    TensorInputs t;

    // Resolve P values first.
    t.P1 = opt.p1_user ? opt.p1 : r1.mean_vector_pol;
    t.P2 = opt.p2_user ? opt.p2 : r2.mean_vector_pol;

    if(opt.perr1_user){
        t.Perr1 = opt.perr1;
    }else if(r1.mean_vector_pol_err > -900.0){
        t.Perr1 = r1.mean_vector_pol_err;
    }

    if(opt.perr2_user){
        t.Perr2 = opt.perr2;
    }else if(r2.mean_vector_pol_err > -900.0){
        t.Perr2 = r2.mean_vector_pol_err;
    }

    // Resolve Q values.
    if(opt.q1_user){
        t.Q1 = opt.q1;
    }else{
        t.Q1 = tensorFromVector(t.P1);
    }

    if(opt.q2_user){
        t.Q2 = opt.q2;
    }else{
        t.Q2 = tensorFromVector(t.P2);
    }

    // Resolve Q errors.
    if(opt.qerr1_user){
        t.Qerr1 = opt.qerr1;
    }else if(!opt.q1_user){
        t.Qerr1 = tensorErrFromVector(t.P1, t.Perr1);
    }

    if(opt.qerr2_user){
        t.Qerr2 = opt.qerr2;
    }else if(!opt.q2_user){
        t.Qerr2 = tensorErrFromVector(t.P2, t.Perr2);
    }

    if(t.Q1 <= -900.0 || t.Q2 <= -900.0){
        cerr << "Error: could not resolve tensor polarizations Q1 and Q2." << endl;
        cerr << "       Provide --q1/--q2 or --p1/--p2." << endl;
        exit(1);
    }

    if(fabs(t.Q1 - t.Q2) <= 1.0e-12){
        cerr << "Error: Q1 and Q2 are effectively equal; two-run extraction is singular." << endl;
        exit(1);
    }

    return t;
}

//============================================================
// ROOT output writer
//============================================================

void writeCutFlowTree(TFile *fout, const RunResult &r1, const RunResult &r2){
    fout->cd();

    TTree *t = new TTree("qe_tensor_cutflow", "QE tensor analysis cumulative cut-flow by run");

    Int_t run = 0;
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

    t->Branch("run", &run, "run/I");
    t->Branch("n_all", &n_all, "n_all/L");
    t->Branch("n_target_match", &n_target_match, "n_target_match/L");
    t->Branch("n_is_qe", &n_is_qe, "n_is_qe/L");
    t->Branch("n_pid", &n_pid, "n_pid/L");
    t->Branch("n_status", &n_status, "n_status/L");
    t->Branch("n_theta", &n_theta, "n_theta/L");
    t->Branch("n_vz", &n_vz, "n_vz/L");
    t->Branch("n_chi2", &n_chi2, "n_chi2/L");
    t->Branch("n_lv_lw", &n_lv_lw, "n_lv_lw/L");
    t->Branch("n_pe", &n_pe, "n_pe/L");
    t->Branch("n_pcal", &n_pcal, "n_pcal/L");
    t->Branch("n_sf", &n_sf, "n_sf/L");
    t->Branch("n_ecin_pcal", &n_ecin_pcal, "n_ecin_pcal/L");
    t->Branch("n_q2", &n_q2, "n_q2/L");
    t->Branch("n_w", &n_w, "n_w/L");
    t->Branch("n_xb_range", &n_xb_range, "n_xb_range/L");

    auto fillOne = [&](const RunResult &r){
        run = r.run;
        n_all = r.cutflow.n_all;
        n_target_match = r.cutflow.n_target_match;
        n_is_qe = r.cutflow.n_is_qe;
        n_pid = r.cutflow.n_pid;
        n_status = r.cutflow.n_status;
        n_theta = r.cutflow.n_theta;
        n_vz = r.cutflow.n_vz;
        n_chi2 = r.cutflow.n_chi2;
        n_lv_lw = r.cutflow.n_lv_lw;
        n_pe = r.cutflow.n_pe;
        n_pcal = r.cutflow.n_pcal;
        n_sf = r.cutflow.n_sf;
        n_ecin_pcal = r.cutflow.n_ecin_pcal;
        n_q2 = r.cutflow.n_q2;
        n_w = r.cutflow.n_w;
        n_xb_range = r.cutflow.n_xb_range;
        t->Fill();
    };

    fillOne(r1);
    fillOne(r2);

    t->Write();
}

void writeRunSummaryTree(TFile *fout, const RunResult &r1, const RunResult &r2, const TensorInputs &tin){
    fout->cd();

    TTree *t = new TTree("qe_tensor_run_summary", "Run-level inputs and FCup-normalized totals");

    Int_t run = 0;
    Long64_t n_inclusive_entries = 0;
    Long64_t n_final = 0;
    Long64_t n_plus = 0;
    Long64_t n_minus = 0;
    Long64_t n_zero = 0;
    Long64_t n_scaler_rows = 0;

    Double_t charge_total = 0.0;
    Double_t charge_plus = 0.0;
    Double_t charge_minus = 0.0;
    Double_t charge_zero = 0.0;
    Double_t yield_total = -9999.0;
    Double_t yield_plus = -9999.0;
    Double_t yield_minus = -9999.0;
    Double_t yield_zero = -9999.0;
    Int_t target_vector_sign = 0;
    Long64_t n_event_files_all = 0;
    Long64_t n_event_files_selected = 0;
    Long64_t n_event_file_index_all = 0;
    Long64_t n_event_file_index_selected = 0;
    // Optional run_summary tree information from the skim.
    Int_t has_skim_run_summary = 0;
    Int_t skim_run_summary_run = -9999;
    Int_t skim_input_file_index = -9999;
    Int_t skim_input_file_number = -9999;
    Long64_t skim_n_events_all = 0;
    Long64_t skim_n_inclusive_written = 0;
    Long64_t skim_n_sidis_rows_written = 0;
    Long64_t skim_n_dihadron_rows_written = 0;
    Long64_t skim_n_scaler_rows = 0;
    Int_t skim_first_event = -9999;
    Int_t skim_last_event = -9999;
    Long64_t skim_first_timestamp = -9999;
    Long64_t skim_last_timestamp = -9999;

    Int_t event_min_all = 0;
    Int_t event_max_all = 0;
    Int_t event_min_selected = 0;
    Int_t event_max_selected = 0;
    Long64_t timestamp_min_all = 0;
    Long64_t timestamp_max_all = 0;
    Long64_t timestamp_min_selected = 0;
    Long64_t timestamp_max_selected = 0;
    Long64_t unixtime_min_all = 0;
    Long64_t unixtime_max_all = 0;
    Long64_t unixtime_min_selected = 0;
    Long64_t unixtime_max_selected = 0;

    Double_t mean_vector_pol = -9999.0;
    Double_t mean_vector_pol_err = -9999.0;
    Double_t mean_target_pol_offline = -9999.0;
    Double_t mean_target_pol_offline_err = -9999.0;
    Double_t mean_tensor_pol_branch = -9999.0;

    Double_t vector_pol_used = -9999.0;
    Double_t vector_pol_err_used = 0.0;
    Double_t tensor_pol_used = -9999.0;
    Double_t tensor_pol_err_used = 0.0;

    t->Branch("run", &run, "run/I");
    t->Branch("n_inclusive_entries", &n_inclusive_entries, "n_inclusive_entries/L");
    t->Branch("n_final", &n_final, "n_final/L");
    t->Branch("n_plus", &n_plus, "n_plus/L");
    t->Branch("n_minus", &n_minus, "n_minus/L");
    t->Branch("n_zero", &n_zero, "n_zero/L");
    t->Branch("n_scaler_rows", &n_scaler_rows, "n_scaler_rows/L");
    t->Branch("charge_total", &charge_total, "charge_total/D");
    t->Branch("charge_plus", &charge_plus, "charge_plus/D");
    t->Branch("charge_minus", &charge_minus, "charge_minus/D");
    t->Branch("charge_zero", &charge_zero, "charge_zero/D");
    t->Branch("yield_total", &yield_total, "yield_total/D");
    t->Branch("yield_plus", &yield_plus, "yield_plus/D");
    t->Branch("yield_minus", &yield_minus, "yield_minus/D");
    t->Branch("yield_zero", &yield_zero, "yield_zero/D");
    t->Branch("target_vector_sign", &target_vector_sign, "target_vector_sign/I");
    t->Branch("n_event_files_all", &n_event_files_all, "n_event_files_all/L");
    t->Branch("n_event_files_selected", &n_event_files_selected, "n_event_files_selected/L");
    t->Branch("n_event_file_index_all", &n_event_file_index_all, "n_event_file_index_all/L");
    t->Branch("n_event_file_index_selected", &n_event_file_index_selected, "n_event_file_index_selected/L");
    t->Branch("has_skim_run_summary", &has_skim_run_summary, "has_skim_run_summary/I");
    t->Branch("skim_run_summary_run", &skim_run_summary_run, "skim_run_summary_run/I");
    t->Branch("skim_input_file_index", &skim_input_file_index, "skim_input_file_index/I");
    t->Branch("skim_input_file_number", &skim_input_file_number, "skim_input_file_number/I");
    t->Branch("skim_n_events_all", &skim_n_events_all, "skim_n_events_all/L");
    t->Branch("skim_n_inclusive_written", &skim_n_inclusive_written, "skim_n_inclusive_written/L");
    t->Branch("skim_n_sidis_rows_written", &skim_n_sidis_rows_written, "skim_n_sidis_rows_written/L");
    t->Branch("skim_n_dihadron_rows_written", &skim_n_dihadron_rows_written, "skim_n_dihadron_rows_written/L");
    t->Branch("skim_n_scaler_rows", &skim_n_scaler_rows, "skim_n_scaler_rows/L");
    t->Branch("skim_first_event", &skim_first_event, "skim_first_event/I");
    t->Branch("skim_last_event", &skim_last_event, "skim_last_event/I");
    t->Branch("skim_first_timestamp", &skim_first_timestamp, "skim_first_timestamp/L");
    t->Branch("skim_last_timestamp", &skim_last_timestamp, "skim_last_timestamp/L");
    t->Branch("event_min_all", &event_min_all, "event_min_all/I");
    t->Branch("event_max_all", &event_max_all, "event_max_all/I");
    t->Branch("event_min_selected", &event_min_selected, "event_min_selected/I");
    t->Branch("event_max_selected", &event_max_selected, "event_max_selected/I");
    t->Branch("timestamp_min_all", &timestamp_min_all, "timestamp_min_all/L");
    t->Branch("timestamp_max_all", &timestamp_max_all, "timestamp_max_all/L");
    t->Branch("timestamp_min_selected", &timestamp_min_selected, "timestamp_min_selected/L");
    t->Branch("timestamp_max_selected", &timestamp_max_selected, "timestamp_max_selected/L");
    t->Branch("unixtime_min_all", &unixtime_min_all, "unixtime_min_all/L");
    t->Branch("unixtime_max_all", &unixtime_max_all, "unixtime_max_all/L");
    t->Branch("unixtime_min_selected", &unixtime_min_selected, "unixtime_min_selected/L");
    t->Branch("unixtime_max_selected", &unixtime_max_selected, "unixtime_max_selected/L");
    t->Branch("mean_vector_pol", &mean_vector_pol, "mean_vector_pol/D");
    t->Branch("mean_vector_pol_err", &mean_vector_pol_err, "mean_vector_pol_err/D");
    t->Branch("mean_target_pol_offline", &mean_target_pol_offline, "mean_target_pol_offline/D");
    t->Branch("mean_target_pol_offline_err", &mean_target_pol_offline_err, "mean_target_pol_offline_err/D");
    t->Branch("mean_tensor_pol_branch", &mean_tensor_pol_branch, "mean_tensor_pol_branch/D");
    t->Branch("vector_pol_used", &vector_pol_used, "vector_pol_used/D");
    t->Branch("vector_pol_err_used", &vector_pol_err_used, "vector_pol_err_used/D");
    t->Branch("tensor_pol_used", &tensor_pol_used, "tensor_pol_used/D");
    t->Branch("tensor_pol_err_used", &tensor_pol_err_used, "tensor_pol_err_used/D");

    auto fillOne = [&](const RunResult &r, int idx){
        run = r.run;
        n_inclusive_entries = r.n_inclusive_entries;
        n_final = r.final_total;
        n_plus = r.final_plus;
        n_minus = r.final_minus;
        n_zero = r.final_zero;
        n_scaler_rows = r.charge.n_scaler_rows;
        charge_total = r.charge.q_total;
        charge_plus = r.charge.q_plus;
        charge_minus = r.charge.q_minus;
        charge_zero = r.charge.q_zero;
        yield_total = (charge_total > 0.0) ? static_cast<double>(n_final) / charge_total : -9999.0;
        yield_plus = (charge_plus > 0.0) ? static_cast<double>(n_plus) / charge_plus : -9999.0;
        yield_minus = (charge_minus > 0.0) ? static_cast<double>(n_minus) / charge_minus : -9999.0;
        yield_zero = (charge_zero > 0.0) ? static_cast<double>(n_zero) / charge_zero : -9999.0;
        target_vector_sign = targetVectorSign((idx == 1) ? tin.P1 : tin.P2);
        n_event_files_all = static_cast<Long64_t>(r.all_by_input_file.size());
        n_event_files_selected = static_cast<Long64_t>(r.selected_by_input_file.size());
        n_event_file_index_all = static_cast<Long64_t>(r.all_by_file_index.size());
        n_event_file_index_selected = static_cast<Long64_t>(r.selected_by_file_index.size());
        has_skim_run_summary = r.has_skim_run_summary ? 1 : 0;
        skim_run_summary_run = r.skim_run_summary_run;
        skim_input_file_index = r.skim_input_file_index;
        skim_input_file_number = r.skim_input_file_number;
        skim_n_events_all = r.skim_n_events_all;
        skim_n_inclusive_written = r.skim_n_inclusive_written;
        skim_n_sidis_rows_written = r.skim_n_sidis_rows_written;
        skim_n_dihadron_rows_written = r.skim_n_dihadron_rows_written;
        skim_n_scaler_rows = r.skim_n_scaler_rows;
        skim_first_event = r.skim_first_event;
        skim_last_event = r.skim_last_event;
        skim_first_timestamp = r.skim_first_timestamp;
        skim_last_timestamp = r.skim_last_timestamp;
        event_min_all = r.seen_event_all ? r.event_min_all : 0;
        event_max_all = r.seen_event_all ? r.event_max_all : 0;
        event_min_selected = r.seen_event_selected ? r.event_min_selected : 0;
        event_max_selected = r.seen_event_selected ? r.event_max_selected : 0;
        timestamp_min_all = r.seen_timestamp_all ? r.timestamp_min_all : 0;
        timestamp_max_all = r.seen_timestamp_all ? r.timestamp_max_all : 0;
        timestamp_min_selected = r.seen_timestamp_selected ? r.timestamp_min_selected : 0;
        timestamp_max_selected = r.seen_timestamp_selected ? r.timestamp_max_selected : 0;
        unixtime_min_all = r.seen_unixtime_all ? r.unixtime_min_all : 0;
        unixtime_max_all = r.seen_unixtime_all ? r.unixtime_max_all : 0;
        unixtime_min_selected = r.seen_unixtime_selected ? r.unixtime_min_selected : 0;
        unixtime_max_selected = r.seen_unixtime_selected ? r.unixtime_max_selected : 0;
        mean_vector_pol = r.mean_vector_pol;
        mean_vector_pol_err = r.mean_vector_pol_err;
        mean_target_pol_offline = r.mean_target_pol_offline;
        mean_target_pol_offline_err = r.mean_target_pol_offline_err;
        mean_tensor_pol_branch = r.mean_tensor_pol_branch;

        if(idx == 1){
            vector_pol_used = tin.P1;
            vector_pol_err_used = tin.Perr1;
            tensor_pol_used = tin.Q1;
            tensor_pol_err_used = tin.Qerr1;
        }else{
            vector_pol_used = tin.P2;
            vector_pol_err_used = tin.Perr2;
            tensor_pol_used = tin.Q2;
            tensor_pol_err_used = tin.Qerr2;
        }

        t->Fill();
    };

    fillOne(r1, 1);
    fillOne(r2, 2);

    t->Write();
}


void writePairSummaryTree(TFile *fout, const RunResult &r1, const RunResult &r2, const TensorInputs &tin){
    fout->cd();

    TTree *t = new TTree("qe_tensor_pair_summary", "Integrated two-run pair screening diagnostics");

    Int_t run1 = r1.run;
    Int_t run2 = r2.run;
    Double_t P1 = tin.P1;
    Double_t P2 = tin.P2;
    Double_t Q1 = tin.Q1;
    Double_t Q2 = tin.Q2;
    Double_t deltaQ = Q1 - Q2;
    Int_t target_sign1 = targetVectorSign(P1);
    Int_t target_sign2 = targetVectorSign(P2);
    Int_t target_pair_class = targetPairClass(P1, P2);

    Double_t C1 = r1.charge.q_total;
    Double_t C2 = r2.charge.q_total;
    Double_t N1 = static_cast<Double_t>(r1.final_total);
    Double_t N2 = static_cast<Double_t>(r2.final_total);
    Double_t Y1 = safeRatio(N1, C1);
    Double_t Y2 = safeRatio(N2, C2);
    Double_t yield_ratio_12 = safeRatio(Y1, Y2);
    Double_t rel_yield_diff_12 = (yield_ratio_12 > -9000.0) ? (yield_ratio_12 - 1.0) : -9999.0;

    Double_t C1_plus = r1.charge.q_plus;
    Double_t C1_minus = r1.charge.q_minus;
    Double_t C2_plus = r2.charge.q_plus;
    Double_t C2_minus = r2.charge.q_minus;
    Double_t Y1_plus = safeRatio(static_cast<Double_t>(r1.final_plus), C1_plus);
    Double_t Y1_minus = safeRatio(static_cast<Double_t>(r1.final_minus), C1_minus);
    Double_t Y2_plus = safeRatio(static_cast<Double_t>(r2.final_plus), C2_plus);
    Double_t Y2_minus = safeRatio(static_cast<Double_t>(r2.final_minus), C2_minus);
    Double_t yield_ratio_plus_12 = safeRatio(Y1_plus, Y2_plus);
    Double_t yield_ratio_minus_12 = safeRatio(Y1_minus, Y2_minus);

    Double_t denomQ = Q1 * Y2 - Q2 * Y1;
    Double_t denomP = P1 * Y2 - P2 * Y1;
    Double_t R2Q_screen = safeRatio(Y1 - Y2, denomQ);
    Double_t R2P_screen = safeRatio(Y1 - Y2, denomP);

    Long64_t n_event_files_all_1 = static_cast<Long64_t>(r1.all_by_input_file.size());
    Long64_t n_event_files_all_2 = static_cast<Long64_t>(r2.all_by_input_file.size());
    Long64_t n_event_files_selected_1 = static_cast<Long64_t>(r1.selected_by_input_file.size());
    Long64_t n_event_files_selected_2 = static_cast<Long64_t>(r2.selected_by_input_file.size());
    Long64_t n_event_file_index_all_1 = static_cast<Long64_t>(r1.all_by_file_index.size());
    Long64_t n_event_file_index_all_2 = static_cast<Long64_t>(r2.all_by_file_index.size());
    Long64_t skim_n_events_all_1 = r1.skim_n_events_all;
    Long64_t skim_n_events_all_2 = r2.skim_n_events_all;
    Double_t skim_raw_yield_1 = safeRatio(static_cast<Double_t>(r1.skim_n_events_all), C1);
    Double_t skim_raw_yield_2 = safeRatio(static_cast<Double_t>(r2.skim_n_events_all), C2);
    Double_t skim_raw_yield_ratio_12 = safeRatio(skim_raw_yield_1, skim_raw_yield_2);
    Double_t skim_selected_fraction_1 = safeRatio(static_cast<Double_t>(r1.n_inclusive_entries), static_cast<Double_t>(r1.skim_n_events_all));
    Double_t skim_selected_fraction_2 = safeRatio(static_cast<Double_t>(r2.n_inclusive_entries), static_cast<Double_t>(r2.skim_n_events_all));

    t->Branch("run1", &run1, "run1/I");
    t->Branch("run2", &run2, "run2/I");
    t->Branch("P1", &P1, "P1/D");
    t->Branch("P2", &P2, "P2/D");
    t->Branch("Q1", &Q1, "Q1/D");
    t->Branch("Q2", &Q2, "Q2/D");
    t->Branch("deltaQ", &deltaQ, "deltaQ/D");
    t->Branch("target_sign1", &target_sign1, "target_sign1/I");
    t->Branch("target_sign2", &target_sign2, "target_sign2/I");
    t->Branch("target_pair_class", &target_pair_class, "target_pair_class/I");
    t->Branch("N1", &N1, "N1/D");
    t->Branch("N2", &N2, "N2/D");
    t->Branch("C1", &C1, "C1/D");
    t->Branch("C2", &C2, "C2/D");
    t->Branch("Y1", &Y1, "Y1/D");
    t->Branch("Y2", &Y2, "Y2/D");
    t->Branch("yield_ratio_12", &yield_ratio_12, "yield_ratio_12/D");
    t->Branch("rel_yield_diff_12", &rel_yield_diff_12, "rel_yield_diff_12/D");
    t->Branch("Y1_plus", &Y1_plus, "Y1_plus/D");
    t->Branch("Y1_minus", &Y1_minus, "Y1_minus/D");
    t->Branch("Y2_plus", &Y2_plus, "Y2_plus/D");
    t->Branch("Y2_minus", &Y2_minus, "Y2_minus/D");
    t->Branch("yield_ratio_plus_12", &yield_ratio_plus_12, "yield_ratio_plus_12/D");
    t->Branch("yield_ratio_minus_12", &yield_ratio_minus_12, "yield_ratio_minus_12/D");
    t->Branch("R2Q_screen", &R2Q_screen, "R2Q_screen/D");
    t->Branch("R2P_screen", &R2P_screen, "R2P_screen/D");
    t->Branch("n_event_files_all_1", &n_event_files_all_1, "n_event_files_all_1/L");
    t->Branch("n_event_files_all_2", &n_event_files_all_2, "n_event_files_all_2/L");
    t->Branch("n_event_files_selected_1", &n_event_files_selected_1, "n_event_files_selected_1/L");
    t->Branch("n_event_files_selected_2", &n_event_files_selected_2, "n_event_files_selected_2/L");
    t->Branch("n_event_file_index_all_1", &n_event_file_index_all_1, "n_event_file_index_all_1/L");
    t->Branch("n_event_file_index_all_2", &n_event_file_index_all_2, "n_event_file_index_all_2/L");
    t->Branch("skim_n_events_all_1", &skim_n_events_all_1, "skim_n_events_all_1/L");
    t->Branch("skim_n_events_all_2", &skim_n_events_all_2, "skim_n_events_all_2/L");
    t->Branch("skim_raw_yield_1", &skim_raw_yield_1, "skim_raw_yield_1/D");
    t->Branch("skim_raw_yield_2", &skim_raw_yield_2, "skim_raw_yield_2/D");
    t->Branch("skim_raw_yield_ratio_12", &skim_raw_yield_ratio_12, "skim_raw_yield_ratio_12/D");
    t->Branch("skim_selected_fraction_1", &skim_selected_fraction_1, "skim_selected_fraction_1/D");
    t->Branch("skim_selected_fraction_2", &skim_selected_fraction_2, "skim_selected_fraction_2/D");

    t->Fill();
    t->Write();

    cout << endl;
    cout << "Integrated pair diagnostics:" << endl;
    cout << "  target pair class:        " << targetPairClassName(target_pair_class) << endl;
    cout << "  Y1 = N1/C1:               " << Y1 << endl;
    cout << "  Y2 = N2/C2:               " << Y2 << endl;
    cout << "  Y1/Y2:                    " << yield_ratio_12 << endl;
    cout << "  relative yield diff:      " << rel_yield_diff_12 << endl;
    cout << "  Delta Q:                  " << deltaQ << endl;
    cout << "  R2Q screen:               " << R2Q_screen << endl;
    cout << "  R2P screen:               " << R2P_screen << endl;
    if(r1.has_skim_run_summary || r2.has_skim_run_summary){
        cout << "  skim raw yield 1:         " << skim_raw_yield_1 << endl;
        cout << "  skim raw yield 2:         " << skim_raw_yield_2 << endl;
        cout << "  skim raw yield ratio:     " << skim_raw_yield_ratio_12 << endl;
        cout << "  skim selected fraction 1: " << skim_selected_fraction_1 << endl;
        cout << "  skim selected fraction 2: " << skim_selected_fraction_2 << endl;
    }
}

void writeMapSummaryTree(
    TFile *fout,
    const string &treeName,
    const string &keyName,
    const RunResult &r1,
    const RunResult &r2,
    const map<int, Long64_t> RunResult::*allMap,
    const map<int, Long64_t> RunResult::*selMap
){
    fout->cd();

    TTree *t = new TTree(treeName.c_str(), treeName.c_str());

    Int_t run = 0;
    Int_t key = 0;
    Long64_t n_all = 0;
    Long64_t n_selected = 0;
    Double_t selected_fraction = 0.0;

    t->Branch("run", &run, "run/I");
    t->Branch(keyName.c_str(), &key, Form("%s/I", keyName.c_str()));
    t->Branch("n_all", &n_all, "n_all/L");
    t->Branch("n_selected", &n_selected, "n_selected/L");
    t->Branch("selected_fraction", &selected_fraction, "selected_fraction/D");

    auto fillRun = [&](const RunResult &r){
        set<int> keys;
        for(auto const &kv : r.*allMap) keys.insert(kv.first);
        for(auto const &kv : r.*selMap) keys.insert(kv.first);

        for(int k : keys){
            run = r.run;
            key = k;
            auto itAll = (r.*allMap).find(k);
            auto itSel = (r.*selMap).find(k);
            n_all = (itAll == (r.*allMap).end()) ? 0 : itAll->second;
            n_selected = (itSel == (r.*selMap).end()) ? 0 : itSel->second;
            selected_fraction = (n_all > 0) ? static_cast<Double_t>(n_selected) / static_cast<Double_t>(n_all) : 0.0;
            t->Fill();
        }
    };

    fillRun(r1);
    fillRun(r2);
    t->Write();
}


void writeTriggerBitSummaryTree(TFile *fout, const RunResult &r1, const RunResult &r2){
    fout->cd();

    TTree *t = new TTree("qe_tensor_trigger_bit_summary", "Decoded trigger-bit summary by run");

    Int_t run = 0;
    Int_t trigger_bit = 0;
    UInt_t trigger_mask = 0;
    Long64_t n_all_with_bit = 0;
    Long64_t n_selected_with_bit = 0;
    Double_t selected_fraction = 0.0;
    Double_t fraction_all = 0.0;
    Double_t fraction_selected = 0.0;

    t->Branch("run", &run, "run/I");
    t->Branch("trigger_bit", &trigger_bit, "trigger_bit/I");
    t->Branch("trigger_mask", &trigger_mask, "trigger_mask/i");
    t->Branch("n_all_with_bit", &n_all_with_bit, "n_all_with_bit/L");
    t->Branch("n_selected_with_bit", &n_selected_with_bit, "n_selected_with_bit/L");
    t->Branch("selected_fraction", &selected_fraction, "selected_fraction/D");
    t->Branch("fraction_all", &fraction_all, "fraction_all/D");
    t->Branch("fraction_selected", &fraction_selected, "fraction_selected/D");

    auto fillRun = [&](const RunResult &r){
        set<int> bits;
        for(auto const &kv : r.all_by_trigger_bit) bits.insert(kv.first);
        for(auto const &kv : r.selected_by_trigger_bit) bits.insert(kv.first);

        for(int bit : bits){
            run = r.run;
            trigger_bit = bit;
            trigger_mask = (1u << bit);
            auto itAll = r.all_by_trigger_bit.find(bit);
            auto itSel = r.selected_by_trigger_bit.find(bit);
            n_all_with_bit = (itAll == r.all_by_trigger_bit.end()) ? 0 : itAll->second;
            n_selected_with_bit = (itSel == r.selected_by_trigger_bit.end()) ? 0 : itSel->second;
            selected_fraction = (n_all_with_bit > 0) ? static_cast<Double_t>(n_selected_with_bit) / static_cast<Double_t>(n_all_with_bit) : 0.0;
            fraction_all = (r.n_inclusive_entries > 0) ? static_cast<Double_t>(n_all_with_bit) / static_cast<Double_t>(r.n_inclusive_entries) : 0.0;
            fraction_selected = (r.final_total > 0) ? static_cast<Double_t>(n_selected_with_bit) / static_cast<Double_t>(r.final_total) : 0.0;
            t->Fill();
        }
    };

    fillRun(r1);
    fillRun(r2);
    t->Write();
}

void writeDiagnosticCSVs(const Options &opt, const RunResult &r1, const RunResult &r2, const TensorInputs &tin){
    const string prefix = csvPrefixFromRootOutput(opt.output);
    const string runCsv = prefix + "_run_diagnostics.csv";
    const string pairCsv = prefix + "_pair_diagnostics.csv";

    ofstream fr(runCsv.c_str());
    fr << "run,P,Q,N_total,C_total,Y_total,N_plus,C_plus,Y_plus,N_minus,C_minus,Y_minus,target_vector_sign,"
       << "n_event_files_all,n_event_files_selected,n_event_file_index_all,n_event_file_index_selected,"
       << "skim_n_events_all,skim_n_inclusive_written,skim_raw_yield,skim_selected_fraction,"
       << "event_min_all,event_max_all,event_min_selected,event_max_selected\n";

    auto writeRun = [&](const RunResult &r, double P, double Q){
        double Ytotal = safeRatio(static_cast<double>(r.final_total), r.charge.q_total);
        double Yplus = safeRatio(static_cast<double>(r.final_plus), r.charge.q_plus);
        double Yminus = safeRatio(static_cast<double>(r.final_minus), r.charge.q_minus);
        fr << r.run << ","
           << P << "," << Q << ","
           << r.final_total << "," << r.charge.q_total << "," << Ytotal << ","
           << r.final_plus << "," << r.charge.q_plus << "," << Yplus << ","
           << r.final_minus << "," << r.charge.q_minus << "," << Yminus << ","
           << targetVectorSign(P) << ","
           << r.all_by_input_file.size() << "," << r.selected_by_input_file.size() << ","
           << r.all_by_file_index.size() << "," << r.selected_by_file_index.size() << ","
           << r.skim_n_events_all << "," << r.skim_n_inclusive_written << ","
           << safeRatio(static_cast<double>(r.skim_n_events_all), r.charge.q_total) << ","
           << safeRatio(static_cast<double>(r.n_inclusive_entries), static_cast<double>(r.skim_n_events_all)) << ","
           << (r.seen_event_all ? r.event_min_all : 0) << "," << (r.seen_event_all ? r.event_max_all : 0) << ","
           << (r.seen_event_selected ? r.event_min_selected : 0) << "," << (r.seen_event_selected ? r.event_max_selected : 0)
           << "\n";
    };

    writeRun(r1, tin.P1, tin.Q1);
    writeRun(r2, tin.P2, tin.Q2);
    fr.close();

    const double Y1 = safeRatio(static_cast<double>(r1.final_total), r1.charge.q_total);
    const double Y2 = safeRatio(static_cast<double>(r2.final_total), r2.charge.q_total);
    const double ratio = safeRatio(Y1, Y2);
    const double rel = (ratio > -9000.0) ? ratio - 1.0 : -9999.0;
    const double denomQ = tin.Q1 * Y2 - tin.Q2 * Y1;
    const double denomP = tin.P1 * Y2 - tin.P2 * Y1;
    const double R2Qscreen = safeRatio(Y1 - Y2, denomQ);
    const double R2Pscreen = safeRatio(Y1 - Y2, denomP);
    const double skimRawYield1 = safeRatio(static_cast<double>(r1.skim_n_events_all), r1.charge.q_total);
    const double skimRawYield2 = safeRatio(static_cast<double>(r2.skim_n_events_all), r2.charge.q_total);
    const double skimRawRatio = safeRatio(skimRawYield1, skimRawYield2);
    const double skimSelectedFrac1 = safeRatio(static_cast<double>(r1.n_inclusive_entries), static_cast<double>(r1.skim_n_events_all));
    const double skimSelectedFrac2 = safeRatio(static_cast<double>(r2.n_inclusive_entries), static_cast<double>(r2.skim_n_events_all));

    ofstream fp(pairCsv.c_str());
    fp << "run1,run2,P1,P2,Q1,Q2,deltaQ,target_pair_class,Y1,Y2,yield_ratio_12,rel_yield_diff_12,R2Q_screen,R2P_screen\n";
    fp << r1.run << "," << r2.run << ","
       << tin.P1 << "," << tin.P2 << ","
       << tin.Q1 << "," << tin.Q2 << "," << (tin.Q1 - tin.Q2) << ","
       << targetPairClass(tin.P1, tin.P2) << ","
       << Y1 << "," << Y2 << "," << ratio << "," << rel << ","
       << R2Qscreen << "," << R2Pscreen << "\n";
    fp.close();

    cout << endl;
    cout << "Wrote diagnostic CSV files:" << endl;
    cout << "  " << runCsv << endl;
    cout << "  " << pairCsv << endl;
}

void writeBinSummaryAndPlots(
    TFile *fout,
    const RunResult &r1,
    const RunResult &r2,
    const TensorInputs &tin,
    const Options &opt
){
    fout->cd();

    TTree *t = new TTree("qe_tensor_bin_summary", "Binned two-run QE tensor extraction summary");

    Int_t bin = 0;
    Double_t xB_low = 0.0;
    Double_t xB_high = 0.0;
    Double_t xB_center = 0.0;
    Double_t xB_width = 0.0;

    Long64_t N1 = 0;
    Long64_t N2 = 0;
    Long64_t N1_plus = 0;
    Long64_t N1_minus = 0;
    Long64_t N1_zero = 0;
    Long64_t N2_plus = 0;
    Long64_t N2_minus = 0;
    Long64_t N2_zero = 0;

    Double_t C1_total = r1.charge.q_total;
    Double_t C2_total = r2.charge.q_total;
    Double_t C1_plus = r1.charge.q_plus;
    Double_t C1_minus = r1.charge.q_minus;
    Double_t C1_zero = r1.charge.q_zero;
    Double_t C2_plus = r2.charge.q_plus;
    Double_t C2_minus = r2.charge.q_minus;
    Double_t C2_zero = r2.charge.q_zero;

    Double_t Y1 = -9999.0;
    Double_t Y2 = -9999.0;
    Double_t Y1_err_stat = 0.0;
    Double_t Y2_err_stat = 0.0;
    Double_t Y1_err_charge = 0.0;
    Double_t Y2_err_charge = 0.0;
    Double_t Y1_err_total = 0.0;
    Double_t Y2_err_total = 0.0;

    Double_t P1_vector = tin.P1;
    Double_t P2_vector = tin.P2;
    Double_t P1_vector_err = tin.Perr1;
    Double_t P2_vector_err = tin.Perr2;
    Double_t Q1_tensor = tin.Q1;
    Double_t Q2_tensor = tin.Q2;
    Double_t Q1_tensor_err = tin.Qerr1;
    Double_t Q2_tensor_err = tin.Qerr2;
    Double_t delta_Q = tin.Q1 - tin.Q2;

    Double_t denom = -9999.0;
    Double_t R2Q = -9999.0;
    Double_t R2Q_err_stat = 0.0;
    Double_t R2Q_err_y_total = 0.0;
    Double_t R2Q_err_total = 0.0;

    // Main presentation convention:
    //   Azz_approx is identical to R2Q.
    // The doubled value is kept only as a diagnostic branch/histogram.
    Double_t Azz_approx = -9999.0;
    Double_t Azz_approx_err_stat = 0.0;
    Double_t Azz_approx_err_y_total = 0.0;
    Double_t Azz_approx_err_total = 0.0;

    Double_t Azz_times2_diagnostic = -9999.0;
    Double_t Azz_times2_diagnostic_err_stat = 0.0;
    Double_t Azz_times2_diagnostic_err_y_total = 0.0;
    Double_t Azz_times2_diagnostic_err_total = 0.0;

    Int_t valid = 0;

    t->Branch("bin", &bin, "bin/I");
    t->Branch("xB_low", &xB_low, "xB_low/D");
    t->Branch("xB_high", &xB_high, "xB_high/D");
    t->Branch("xB_center", &xB_center, "xB_center/D");
    t->Branch("xB_width", &xB_width, "xB_width/D");

    t->Branch("N1", &N1, "N1/L");
    t->Branch("N2", &N2, "N2/L");
    t->Branch("N1_plus", &N1_plus, "N1_plus/L");
    t->Branch("N1_minus", &N1_minus, "N1_minus/L");
    t->Branch("N1_zero", &N1_zero, "N1_zero/L");
    t->Branch("N2_plus", &N2_plus, "N2_plus/L");
    t->Branch("N2_minus", &N2_minus, "N2_minus/L");
    t->Branch("N2_zero", &N2_zero, "N2_zero/L");

    t->Branch("C1_total", &C1_total, "C1_total/D");
    t->Branch("C2_total", &C2_total, "C2_total/D");
    t->Branch("C1_plus", &C1_plus, "C1_plus/D");
    t->Branch("C1_minus", &C1_minus, "C1_minus/D");
    t->Branch("C1_zero", &C1_zero, "C1_zero/D");
    t->Branch("C2_plus", &C2_plus, "C2_plus/D");
    t->Branch("C2_minus", &C2_minus, "C2_minus/D");
    t->Branch("C2_zero", &C2_zero, "C2_zero/D");

    t->Branch("Y1", &Y1, "Y1/D");
    t->Branch("Y2", &Y2, "Y2/D");
    t->Branch("Y1_err_stat", &Y1_err_stat, "Y1_err_stat/D");
    t->Branch("Y2_err_stat", &Y2_err_stat, "Y2_err_stat/D");
    t->Branch("Y1_err_charge", &Y1_err_charge, "Y1_err_charge/D");
    t->Branch("Y2_err_charge", &Y2_err_charge, "Y2_err_charge/D");
    t->Branch("Y1_err_total", &Y1_err_total, "Y1_err_total/D");
    t->Branch("Y2_err_total", &Y2_err_total, "Y2_err_total/D");

    t->Branch("P1_vector", &P1_vector, "P1_vector/D");
    t->Branch("P2_vector", &P2_vector, "P2_vector/D");
    t->Branch("P1_vector_err", &P1_vector_err, "P1_vector_err/D");
    t->Branch("P2_vector_err", &P2_vector_err, "P2_vector_err/D");
    t->Branch("Q1_tensor", &Q1_tensor, "Q1_tensor/D");
    t->Branch("Q2_tensor", &Q2_tensor, "Q2_tensor/D");
    t->Branch("Q1_tensor_err", &Q1_tensor_err, "Q1_tensor_err/D");
    t->Branch("Q2_tensor_err", &Q2_tensor_err, "Q2_tensor_err/D");
    t->Branch("delta_Q", &delta_Q, "delta_Q/D");

    t->Branch("denom", &denom, "denom/D");
    t->Branch("R2Q", &R2Q, "R2Q/D");
    t->Branch("R2Q_err_stat", &R2Q_err_stat, "R2Q_err_stat/D");
    t->Branch("R2Q_err_y_total", &R2Q_err_y_total, "R2Q_err_y_total/D");
    t->Branch("R2Q_err_total", &R2Q_err_total, "R2Q_err_total/D");
    t->Branch("Azz_approx", &Azz_approx, "Azz_approx/D");
    t->Branch("Azz_approx_err_stat", &Azz_approx_err_stat, "Azz_approx_err_stat/D");
    t->Branch("Azz_approx_err_y_total", &Azz_approx_err_y_total, "Azz_approx_err_y_total/D");
    t->Branch("Azz_approx_err_total", &Azz_approx_err_total, "Azz_approx_err_total/D");

    t->Branch("Azz_times2_diagnostic", &Azz_times2_diagnostic, "Azz_times2_diagnostic/D");
    t->Branch("Azz_times2_diagnostic_err_stat", &Azz_times2_diagnostic_err_stat, "Azz_times2_diagnostic_err_stat/D");
    t->Branch("Azz_times2_diagnostic_err_y_total", &Azz_times2_diagnostic_err_y_total, "Azz_times2_diagnostic_err_y_total/D");
    t->Branch("Azz_times2_diagnostic_err_total", &Azz_times2_diagnostic_err_total, "Azz_times2_diagnostic_err_total/D");
    t->Branch("valid", &valid, "valid/I");

    TH1D *h_N1 = new TH1D("h_N1_vs_xB", Form("Run %d counts;x_{B};counts", r1.run), opt.xbins, opt.xmin, opt.xmax);
    TH1D *h_N2 = new TH1D("h_N2_vs_xB", Form("Run %d counts;x_{B};counts", r2.run), opt.xbins, opt.xmin, opt.xmax);
    TH1D *h_Y1 = new TH1D("h_Y1_vs_xB", Form("Run %d FCup-normalized yield;x_{B};N/C", r1.run), opt.xbins, opt.xmin, opt.xmax);
    TH1D *h_Y2 = new TH1D("h_Y2_vs_xB", Form("Run %d FCup-normalized yield;x_{B};N/C", r2.run), opt.xbins, opt.xmin, opt.xmax);

    TH1D *h_R2Q_stat = new TH1D("h_R2Q_vs_xB_stat", "R_{2}^{Q} vs x_{B} with statistical errors;x_{B};R_{2}^{Q}", opt.xbins, opt.xmin, opt.xmax);
    TH1D *h_R2Q_total = new TH1D("h_R2Q_vs_xB_total", "R_{2}^{Q} vs x_{B} with total errors;x_{B};R_{2}^{Q}", opt.xbins, opt.xmin, opt.xmax);
    TH1D *h_Azz_approx_total = new TH1D("h_Azz_approx_vs_xB_total", "A_{zz}^{approx}=R_{2}^{Q} vs x_{B} with total errors;x_{B};A_{zz}^{approx}", opt.xbins, opt.xmin, opt.xmax);
    TH1D *h_Azz_times2_diagnostic_total = new TH1D("h_Azz_times2_diagnostic_vs_xB_total", "Diagnostic 2R_{2}^{Q} vs x_{B} with total errors;x_{B};2R_{2}^{Q}", opt.xbins, opt.xmin, opt.xmax);

    if(USE_HARDCODED_Y_RANGE){
        h_R2Q_stat->SetMinimum(R2Q_YMIN);
        h_R2Q_stat->SetMaximum(R2Q_YMAX);

        h_R2Q_total->SetMinimum(R2Q_YMIN);
        h_R2Q_total->SetMaximum(R2Q_YMAX);

        h_Azz_approx_total->SetMinimum(R2Q_YMIN);
        h_Azz_approx_total->SetMaximum(R2Q_YMAX);

        const double y2min = SCALE_TIMES2_DIAGNOSTIC_Y_RANGE ? 2.0 * R2Q_YMIN : R2Q_YMIN;
        const double y2max = SCALE_TIMES2_DIAGNOSTIC_Y_RANGE ? 2.0 * R2Q_YMAX : R2Q_YMAX;
        h_Azz_times2_diagnostic_total->SetMinimum(y2min);
        h_Azz_times2_diagnostic_total->SetMaximum(y2max);
    }

    vector<double> gx;
    vector<double> gex;
    vector<double> gR;
    vector<double> gRerrStat;
    vector<double> gRerrTotal;
    vector<double> gAzzApprox;
    vector<double> gAzzApproxErrTotal;
    vector<double> gAzzTimes2Diagnostic;
    vector<double> gAzzTimes2DiagnosticErrTotal;

    const double dx = (opt.xmax - opt.xmin) / static_cast<double>(opt.xbins);

    for(int ib = 0; ib < opt.xbins; ++ib){
        bin = ib;
        xB_low = opt.xmin + ib * dx;
        xB_high = xB_low + dx;
        xB_center = 0.5 * (xB_low + xB_high);
        xB_width = dx;

        const BinCounts &b1 = r1.bins[ib];
        const BinCounts &b2 = r2.bins[ib];

        N1 = b1.total;
        N2 = b2.total;
        N1_plus = b1.plus;
        N1_minus = b1.minus;
        N1_zero = b1.zero;
        N2_plus = b2.plus;
        N2_minus = b2.minus;
        N2_zero = b2.zero;

        Y1 = -9999.0;
        Y2 = -9999.0;
        Y1_err_stat = 0.0;
        Y2_err_stat = 0.0;
        Y1_err_charge = 0.0;
        Y2_err_charge = 0.0;
        Y1_err_total = 0.0;
        Y2_err_total = 0.0;
        denom = -9999.0;
        R2Q = -9999.0;
        R2Q_err_stat = 0.0;
        R2Q_err_y_total = 0.0;
        R2Q_err_total = 0.0;
        Azz_approx = -9999.0;
        Azz_approx_err_stat = 0.0;
        Azz_approx_err_y_total = 0.0;
        Azz_approx_err_total = 0.0;
        Azz_times2_diagnostic = -9999.0;
        Azz_times2_diagnostic_err_stat = 0.0;
        Azz_times2_diagnostic_err_y_total = 0.0;
        Azz_times2_diagnostic_err_total = 0.0;
        valid = 0;

        h_N1->SetBinContent(ib + 1, N1);
        h_N2->SetBinContent(ib + 1, N2);
        h_N1->SetBinError(ib + 1, (N1 > 0) ? sqrt(static_cast<double>(N1)) : 0.0);
        h_N2->SetBinError(ib + 1, (N2 > 0) ? sqrt(static_cast<double>(N2)) : 0.0);

        if(N1 > 0 && N2 > 0 && C1_total > 0.0 && C2_total > 0.0){
            Y1 = static_cast<double>(N1) / C1_total;
            Y2 = static_cast<double>(N2) / C2_total;

            Y1_err_stat = sqrt(static_cast<double>(N1)) / C1_total;
            Y2_err_stat = sqrt(static_cast<double>(N2)) / C2_total;

            Y1_err_charge = fabs(Y1) * opt.charge_rel_err;
            Y2_err_charge = fabs(Y2) * opt.charge_rel_err;

            Y1_err_total = sqrt(Y1_err_stat * Y1_err_stat + Y1_err_charge * Y1_err_charge);
            Y2_err_total = sqrt(Y2_err_stat * Y2_err_stat + Y2_err_charge * Y2_err_charge);

            h_Y1->SetBinContent(ib + 1, Y1);
            h_Y2->SetBinContent(ib + 1, Y2);
            h_Y1->SetBinError(ib + 1, Y1_err_total);
            h_Y2->SetBinError(ib + 1, Y2_err_total);

            denom = Q1_tensor * Y2 - Q2_tensor * Y1;

            if(fabs(denom) > 1.0e-30){
                const double numerator = Y1 - Y2;
                R2Q = numerator / denom;

                const double dR_dY1 = Y2 * (Q1_tensor - Q2_tensor) / (denom * denom);
                const double dR_dY2 = -Y1 * (Q1_tensor - Q2_tensor) / (denom * denom);
                const double dR_dQ1 = -numerator * Y2 / (denom * denom);
                const double dR_dQ2 = numerator * Y1 / (denom * denom);

                R2Q_err_stat = sqrt(
                    dR_dY1 * dR_dY1 * Y1_err_stat * Y1_err_stat +
                    dR_dY2 * dR_dY2 * Y2_err_stat * Y2_err_stat
                );

                R2Q_err_y_total = sqrt(
                    dR_dY1 * dR_dY1 * Y1_err_total * Y1_err_total +
                    dR_dY2 * dR_dY2 * Y2_err_total * Y2_err_total
                );

                R2Q_err_total = sqrt(
                    R2Q_err_y_total * R2Q_err_y_total +
                    dR_dQ1 * dR_dQ1 * Q1_tensor_err * Q1_tensor_err +
                    dR_dQ2 * dR_dQ2 * Q2_tensor_err * Q2_tensor_err
                );

                Azz_approx = R2Q;
                Azz_approx_err_stat = R2Q_err_stat;
                Azz_approx_err_y_total = R2Q_err_y_total;
                Azz_approx_err_total = R2Q_err_total;

                Azz_times2_diagnostic = 2.0 * R2Q;
                Azz_times2_diagnostic_err_stat = 2.0 * R2Q_err_stat;
                Azz_times2_diagnostic_err_y_total = 2.0 * R2Q_err_y_total;
                Azz_times2_diagnostic_err_total = 2.0 * R2Q_err_total;

                valid = 1;

                h_R2Q_stat->SetBinContent(ib + 1, R2Q);
                h_R2Q_stat->SetBinError(ib + 1, R2Q_err_stat);
                h_R2Q_total->SetBinContent(ib + 1, R2Q);
                h_R2Q_total->SetBinError(ib + 1, R2Q_err_total);
                h_Azz_approx_total->SetBinContent(ib + 1, Azz_approx);
                h_Azz_approx_total->SetBinError(ib + 1, Azz_approx_err_total);
                h_Azz_times2_diagnostic_total->SetBinContent(ib + 1, Azz_times2_diagnostic);
                h_Azz_times2_diagnostic_total->SetBinError(ib + 1, Azz_times2_diagnostic_err_total);

                gx.push_back(xB_center);
                gex.push_back(0.5 * dx);
                gR.push_back(R2Q);
                gRerrStat.push_back(R2Q_err_stat);
                gRerrTotal.push_back(R2Q_err_total);
                gAzzApprox.push_back(Azz_approx);
                gAzzApproxErrTotal.push_back(Azz_approx_err_total);
                gAzzTimes2Diagnostic.push_back(Azz_times2_diagnostic);
                gAzzTimes2DiagnosticErrTotal.push_back(Azz_times2_diagnostic_err_total);
            }
        }

        t->Fill();
    }

    t->Write();

    TDirectory *dHists = fout->mkdir("hists_qe_tensor");
    dHists->cd();
    h_N1->Write();
    h_N2->Write();
    h_Y1->Write();
    h_Y2->Write();
    h_R2Q_stat->Write();
    h_R2Q_total->Write();
    h_Azz_approx_total->Write();
    h_Azz_times2_diagnostic_total->Write();
    fout->cd();

    const int nGraph = static_cast<int>(gx.size());

    if(nGraph > 0){
        TGraphErrors *gR_stat = new TGraphErrors(nGraph, gx.data(), gR.data(), gex.data(), gRerrStat.data());
        gR_stat->SetName("g_R2Q_vs_xB_stat");
        gR_stat->SetTitle("R_{2}^{Q} vs x_{B};x_{B};R_{2}^{Q}");
        gR_stat->SetMarkerStyle(20);

        TGraphErrors *gR_total = new TGraphErrors(nGraph, gx.data(), gR.data(), gex.data(), gRerrTotal.data());
        gR_total->SetName("g_R2Q_vs_xB_total");
        gR_total->SetTitle("R_{2}^{Q} vs x_{B};x_{B};R_{2}^{Q}");
        gR_total->SetMarkerStyle(20);

        TGraphErrors *gAzz_approx_total = new TGraphErrors(nGraph, gx.data(), gAzzApprox.data(), gex.data(), gAzzApproxErrTotal.data());
        gAzz_approx_total->SetName("g_Azz_approx_vs_xB_total");
        gAzz_approx_total->SetTitle("A_{zz}^{approx}=R_{2}^{Q} vs x_{B};x_{B};A_{zz}^{approx}");
        gAzz_approx_total->SetMarkerStyle(20);

        TGraphErrors *gAzz_times2_diagnostic_total = new TGraphErrors(nGraph, gx.data(), gAzzTimes2Diagnostic.data(), gex.data(), gAzzTimes2DiagnosticErrTotal.data());
        gAzz_times2_diagnostic_total->SetName("g_Azz_times2_diagnostic_vs_xB_total");
        gAzz_times2_diagnostic_total->SetTitle("Diagnostic 2R_{2}^{Q} vs x_{B};x_{B};2R_{2}^{Q}");
        gAzz_times2_diagnostic_total->SetMarkerStyle(20);

        if(USE_HARDCODED_Y_RANGE){
            gR_stat->SetMinimum(R2Q_YMIN);
            gR_stat->SetMaximum(R2Q_YMAX);

            gR_total->SetMinimum(R2Q_YMIN);
            gR_total->SetMaximum(R2Q_YMAX);

            gAzz_approx_total->SetMinimum(R2Q_YMIN);
            gAzz_approx_total->SetMaximum(R2Q_YMAX);

            const double y2min = SCALE_TIMES2_DIAGNOSTIC_Y_RANGE ? 2.0 * R2Q_YMIN : R2Q_YMIN;
            const double y2max = SCALE_TIMES2_DIAGNOSTIC_Y_RANGE ? 2.0 * R2Q_YMAX : R2Q_YMAX;
            gAzz_times2_diagnostic_total->SetMinimum(y2min);
            gAzz_times2_diagnostic_total->SetMaximum(y2max);
        }

        gR_stat->Write();
        gR_total->Write();
        gAzz_approx_total->Write();
        gAzz_times2_diagnostic_total->Write();

        TCanvas *cR = new TCanvas("c_R2Q_vs_xB", "R2Q vs xB", 900, 700);
        gR_total->Draw("AP");
        cR->Write();

        TCanvas *cA = new TCanvas("c_Azz_approx_vs_xB", "Azz approx vs xB", 900, 700);
        gAzz_approx_total->Draw("AP");
        cA->Write();

        TCanvas *cA2 = new TCanvas("c_Azz_times2_diagnostic_vs_xB", "2R2Q diagnostic vs xB", 900, 700);
        gAzz_times2_diagnostic_total->Draw("AP");
        cA2->Write();

        if(opt.save_plots){
            const string prefix = removeRootSuffix(opt.output);
            cR->SaveAs((prefix + "_R2Q_vs_xB.png").c_str());
            cR->SaveAs((prefix + "_R2Q_vs_xB.pdf").c_str());
            cA->SaveAs((prefix + "_Azz_approx_vs_xB.png").c_str());
            cA->SaveAs((prefix + "_Azz_approx_vs_xB.pdf").c_str());
            cA2->SaveAs((prefix + "_Azz_times2_diagnostic_vs_xB.png").c_str());
            cA2->SaveAs((prefix + "_Azz_times2_diagnostic_vs_xB.pdf").c_str());
        }
    }else{
        cerr << "Warning: no valid xB bins for R2Q graph. Check cuts and binning." << endl;
    }
}

//============================================================
// Main
//============================================================

int main(int argc, char **argv){
    Options opt = parseOptions(argc, argv);

    cout << endl;
    cout << "============================================================" << endl;
    cout << "  QE tensor two-run analysis" << endl;
    cout << "  Version: " << QE_ANALYSIS_VERSION << endl;
    cout << "============================================================" << endl;
    cout << "File 1:              " << opt.file1 << endl;
    cout << "File 2:              " << opt.file2 << endl;
    cout << "Run 1:               " << opt.run1 << endl;
    cout << "Run 2:               " << opt.run2 << endl;
    cout << "Cut profile:         " << opt.cut_profile << endl;
    cout << "Apply kin cuts:      " << (opt.apply_kinematic_cuts ? 1 : 0) << endl;
    cout << "Output:              " << opt.output << endl;
    cout << "Charge branch:       " << opt.charge_branch << endl;
    cout << "Charge rel. error:   " << opt.charge_rel_err << endl;
    cout << "xB bins:             " << opt.xbins << endl;
    cout << "xB range:            " << opt.xmin << " to " << opt.xmax << endl;
    if(USE_HARDCODED_Y_RANGE){
        cout << "Hardcoded y range:   " << R2Q_YMIN << " to " << R2Q_YMAX << endl;
        if(SCALE_TIMES2_DIAGNOSTIC_Y_RANGE){
            cout << "2R2Q y range:        " << 2.0 * R2Q_YMIN << " to " << 2.0 * R2Q_YMAX << endl;
        }else{
            cout << "2R2Q y range:        " << R2Q_YMIN << " to " << R2Q_YMAX << endl;
        }
    }
    cout << "theta range [deg]:   " << opt.theta_min_deg << " to " << opt.theta_max_deg << endl;
    cout << "vz range [cm]:       " << opt.vz_min << " to " << opt.vz_max << endl;
    cout << "Q2 range [GeV2]:     " << opt.q2_min << " to " << opt.q2_max << endl;
    cout << "W range [GeV]:       " << opt.w_min << " to " << opt.w_max << endl;
    cout << "Write diagnostics:   " << (opt.write_diagnostics ? 1 : 0) << endl;
    cout << "============================================================" << endl;

    RunResult r1 = processRun(opt.file1, opt.run1, opt);
    RunResult r2 = processRun(opt.file2, opt.run2, opt);

    printRunSummary(r1);
    printRunSummary(r2);

    printCutFlow(Form("run %d", r1.run), r1.cutflow);
    printCutFlow(Form("run %d", r2.run), r2.cutflow);

    TensorInputs tin = resolveTensorInputs(opt, r1, r2);

    cout << endl;
    cout << "Resolved tensor inputs:" << endl;
    cout << "  P1 = " << tin.P1 << " +/- " << tin.Perr1 << endl;
    cout << "  P2 = " << tin.P2 << " +/- " << tin.Perr2 << endl;
    cout << "  Q1 = " << tin.Q1 << " +/- " << tin.Qerr1 << endl;
    cout << "  Q2 = " << tin.Q2 << " +/- " << tin.Qerr2 << endl;
    cout << "  Delta Q = " << (tin.Q1 - tin.Q2) << endl;
    cout << endl;

    TFile *fout = new TFile(opt.output.c_str(), "RECREATE");
    if(!fout || fout->IsZombie()){
        cerr << "Error: could not create output file: " << opt.output << endl;
        return 1;
    }

    writeRunSummaryTree(fout, r1, r2, tin);
    writeCutFlowTree(fout, r1, r2);
    if(opt.write_diagnostics){
        writePairSummaryTree(fout, r1, r2, tin);
        writeMapSummaryTree(fout, "qe_tensor_file_index_summary", "event_file_index", r1, r2,
                            &RunResult::all_by_file_index, &RunResult::selected_by_file_index);
        writeMapSummaryTree(fout, "qe_tensor_input_file_summary", "input_file_key", r1, r2,
                            &RunResult::all_by_input_file, &RunResult::selected_by_input_file);
        writeMapSummaryTree(fout, "qe_tensor_trigger_summary", "trigger", r1, r2,
                            &RunResult::all_by_trigger, &RunResult::selected_by_trigger);
        writeTriggerBitSummaryTree(fout, r1, r2);
        writeMapSummaryTree(fout, "qe_tensor_run_config_summary", "run_config", r1, r2,
                            &RunResult::all_by_run_config, &RunResult::selected_by_run_config);
        writeMapSummaryTree(fout, "qe_tensor_event_config_summary", "event_config", r1, r2,
                            &RunResult::all_by_event_config, &RunResult::selected_by_event_config);
        writeDiagnosticCSVs(opt, r1, r2, tin);
    }
    writeBinSummaryAndPlots(fout, r1, r2, tin, opt);

    fout->Close();

    cout << endl;
    cout << "Wrote output ROOT file:" << endl;
    cout << "  " << opt.output << endl;
    if(opt.write_diagnostics){
        cout << "Diagnostic ROOT trees added:" << endl;
        cout << "  qe_tensor_pair_summary" << endl;
        cout << "  qe_tensor_file_index_summary" << endl;
        cout << "  qe_tensor_input_file_summary" << endl;
        cout << "  qe_tensor_trigger_summary" << endl;
        cout << "  qe_tensor_trigger_bit_summary" << endl;
        cout << "  qe_tensor_run_config_summary" << endl;
        cout << "  qe_tensor_event_config_summary" << endl;
    }

    if(opt.save_plots){
        const string prefix = removeRootSuffix(opt.output);
        cout << "Saved plots:" << endl;
        cout << "  " << prefix << "_R2Q_vs_xB.png" << endl;
        cout << "  " << prefix << "_R2Q_vs_xB.pdf" << endl;
        cout << "  " << prefix << "_Azz_approx_vs_xB.png" << endl;
        cout << "  " << prefix << "_Azz_approx_vs_xB.pdf" << endl;
        cout << "  " << prefix << "_Azz_times2_diagnostic_vs_xB.png" << endl;
        cout << "  " << prefix << "_Azz_times2_diagnostic_vs_xB.pdf" << endl;
    }

    cout << endl;
    cout << "Done." << endl;

    return 0;
}