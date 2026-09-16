//******************************************************************
//*  ██╗  ██╗██╗██████╗  ██████╗     ██╗  ██╗    ██████╗
//*  ██║  ██║██║██╔══██╗██╔═══██╗    ██║  ██║   ██╔═████╗
//*  ███████║██║██████╔╝██║   ██║    ███████║   ██║██╔██║
//*  ██╔══██║██║██╔═══╝ ██║   ██║    ╚════██║   ████╔╝██║
//*  ██║  ██║██║██║     ╚██████╔╝         ██║██╗╚██████╔╝
//*  ╚═╝  ╚═╝╚═╝╚═╝      ╚═════╝          ╚═╝╚═╝ ╚═════╝
//************************ Jefferson National Lab (2017) ***********
//******************************************************************

#include <cstdlib>
#include <iostream>
#include "reader.h"
#include "writer.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <math.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <iomanip>
#include <time.h>
#include <TVector3.h>
#include <TMath.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TLegend.h"
#include <unordered_map>
#include <float.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCutG.h>
#include <TCut.h>
#include <TPad.h>
#include <TGLabel.h>
#include <TGaxis.h>
#include "TTree.h"
#include "TBranch.h"
#include <algorithm>
#include <cctype>

using namespace std;

// Original value used in earlier skim production.
// const double beam_enrg = 10.6;       // beam energy in GeV

// Use the same beam energy as FCup_reading_summer2026_V1.cc.
// This affects W, Q2, nu, xB, and y.
const double beam_enrg = 10.6041;    // beam energy in GeV

// Inclusive comparison switch.
// true  = inclusive mode uses REC::Particle row 0, matching FCup_reading_summer2026_V1.cc
// false = inclusive mode goes back to your original highest momentum electron choice
const bool USE_ROW0_ELECTRON_FOR_INCLUSIVE = true;

const double M_d = 1.875612928;      // deuteron mass in GeV
const double M_p = 0.938272088;      // proton mass in GeV
const double M_pion = 0.13957039;    // charged pion mass in GeV
const double M_electron = 0.000510999; // electron mass in GeV

// Dedicated inclusive DIS electron cuts.
// These are used only for:
//   --mode inclusive --region dis --detpidcut 1
// Detector cuts are selected with detpidcut.
// Kinematic cuts are selected with kincut.

const double DEDICATED_DIS_Q2_MIN = 1.0;
const double DEDICATED_DIS_W_MIN = 2.0;
const double DEDICATED_DIS_EP_MIN = 2.6;

const double DEDICATED_DIS_THETA_MIN_DEG = 5.0;
const double DEDICATED_DIS_THETA_MAX_DEG = 40.0;

const double DEDICATED_DIS_VZ_MIN = -5.758;
const double DEDICATED_DIS_VZ_MAX = 1.515;

const double DEDICATED_DIS_NPHE_MIN = 2.0;
const double DEDICATED_DIS_PCAL_E_MIN = 0.06;

const double DEDICATED_DIS_PCAL_LVW_S1 = 22.5;
const double DEDICATED_DIS_PCAL_LVW_S26 = 13.5;

const int DEDICATED_DIS_DC_DETECTOR = 6;
const int DEDICATED_DIS_DC_LAYER_R1 = 6;
const int DEDICATED_DIS_DC_LAYER_R2 = 18;
const int DEDICATED_DIS_DC_LAYER_R3 = 36;

const double DEDICATED_DIS_DC_EDGE_R1 = 4.0;
const double DEDICATED_DIS_DC_EDGE_R2 = 5.0;
const double DEDICATED_DIS_DC_EDGE_R3 = 8.0;

const double DEDICATED_DIS_SF_LOW[6][3] = {
    {0.188243, +0.004928, -0.000498},
    {0.183701, +0.008171, -0.000765},
    {0.185788, +0.007974, -0.000773},
    {0.177398, +0.010710, -0.000954},
    {0.181192, +0.007219, -0.000601},
    {0.187191, +0.005478, -0.000496}
};

const double DEDICATED_DIS_SF_HIGH[6][3] = {
    {0.308340, -0.007150, +0.000487},
    {0.310391, -0.003058, -0.000184},
    {0.310277, -0.004496, -0.000111},
    {0.314889, -0.003913, -0.000101},
    {0.312211, -0.007023, +0.000140},
    {0.309825, -0.004853, -0.000123}
};

// Dedicated inclusive QE electron cuts.
// These are used only for inclusive QE when detpidcut or kincut is selected.
// They are kept separate from the broad QE flag.

const double DEDICATED_QE_THETA_MIN_DEG = 7.80;
const double DEDICATED_QE_THETA_MAX_DEG = 8.20;

const double DEDICATED_QE_VZ_MIN = -5.758;
const double DEDICATED_QE_VZ_MAX = 1.5165;

const double DEDICATED_QE_ABS_CHI2PID_MAX = 3.0;

const double DEDICATED_QE_LV_MIN = 14.0;
const double DEDICATED_QE_LW_MIN = 14.0;

const double DEDICATED_QE_SF_MAX = 0.28;

const double DEDICATED_QE_ECIN_PCAL_SLOPE = -0.625;
const double DEDICATED_QE_ECIN_PCAL_INTERCEPT = 0.15;

const double DEDICATED_QE_P_MIN = 2.0;
const double DEDICATED_QE_PCAL_E_MIN = 0.07;

const double DEDICATED_QE_Q2_MIN = 1.9433;
const double DEDICATED_QE_Q2_MAX = 2.0574;

const double DEDICATED_QE_W_MIN = 0.0;
const double DEDICATED_QE_W_MAX = 1.073;

typedef unordered_map<int,vector<int>> vectorMap;

// Container for detector/PID information connected to one REC:Particle row.
// This is used only internally.
// The ROOT output uses flat scalar branches.
struct DetectorInfo {
    int charge = -9999;
    float beta = -9999;
    float chi2pid = -9999;
    int has_htcc = 0;
    float htcc_nphe = -9999;
    int has_cal = 0;
    int cal_sector = -9999;
    float pcal_e = 0.0;
    float ecin_e = 0.0;
    float ecout_e = 0.0;
    float cal_e = 0.0;
    float sampling_fraction = -9999;
    float pcal_x = -9999;
    float pcal_y = -9999;
    float lu = -9999;
    float lv = -9999;
    float lw = -9999;
};

// Container for DC fiducial edge information from REC::Traj.
// This is needed for the FCup-style inclusive DIS electron cut.
struct DCEdgeInfo {
    float edge_r1 = -9999;
    float edge_r2 = -9999;
    float edge_r3 = -9999;

    int has_r1 = 0;
    int has_r2 = 0;
    int has_r3 = 0;
};

// Container for event level information.
// This is shared by all future skim modes:
//   SIDIS
//   inclusive
//   exclusive, later
// Phase 1 reads the most important quantities from:
//   REC::Event
//   RUN::config
// Target/tensor polarization quantities are placeholders for now.
// They will be filled later from a run dependent target table.

struct EventInfo {
    // True run/event metadata from RUN::config when available.
    int run = -9999;
    int event = -9999;

    // Local event counter inside this input HIPO file.
    // This preserves the old counter style event indexing if needed.
    int event_file_index = -9999;

    // True input-HIPO-file provenance.
    // input_file_index is the index in the current command's hipoFiles vector.
    // input_file_number is parsed from the HIPO filename when possible.
    int input_file_index = -9999;
    int input_file_number = -9999;

    // Beam helicity from REC::Event.helicity.
    // Expected valid values are usually -1 and +1.
    // Anything else is treated as invalid.
    int helicity = -9999;
    int helicity_valid = 0;

    // RUN::config metadata.
    Long64_t trigger = -9999;
    Long64_t timestamp = -9999;
    int unixtime = -9999;

    // Placeholders for tensor target / b1 / Azz analysis.
    // These are not filled from HIPO yet.
    // Later we will load them from a run table or target log.
    int target_state = -9999;

    float target_pol = -9999;
    float target_pol_err = -9999;

    float tensor_pol = -9999;
    float tensor_pol_err = -9999;

    float vector_pol = -9999;
    float vector_pol_err = -9999;

    // Raw target map information.
    int target_match = 0;

    string target_species = "";
    int target_species_id = -9999;

    int target_cell = -9999;

    float target_pol_online = -9999;
    float target_pol_offline = -9999;
    float target_pol_offline_err = -9999;

    float target_run_dose = -9999;

    string target_start_time = "";
    string target_stop_time = "";
};

// Container for one RUN::scaler row written to the skim ROOT file.
// This is needed for FCup normalized helicity asymmetries.
// The tree is independent of SIDIS/inclusive selections and is written for every skim mode.
struct ScalerBranchVars {
    int run = -9999;
    int event = -9999;
    int event_file_index = -9999;

    int input_file_index = -9999;
    int input_file_number = -9999;

    int helicity = -9999;
    int helicity_valid = 0;

    // Raw cumulative scaler values from RUN::scaler.
    float fcup = 0.0;
    float fcupgated = 0.0;

    // Difference between consecutive scaler reads.
    // These are the values that should be integrated for Q+ and Q-.
    float fcup_delta = 0.0;
    float fcupgated_delta = 0.0;
};

// Per input file summary. This keeps the all trigger helicity counts that are lost if we only look at selected skim rows later.
struct RunSummaryVars {
    int run = -9999;

    int input_file_index = -9999;
    int input_file_number = -9999;
    string input_file_name = "";

    Long64_t n_events_all = 0;
    Long64_t n_inclusive_written = 0;
    Long64_t n_sidis_rows_written = 0;
    Long64_t n_dihadron_rows_written = 0;

    int first_event = -9999;
    int last_event = -9999;
    Long64_t first_timestamp = -9999;
    Long64_t last_timestamp = -9999;
    Long64_t n_hel_plus_all = 0;
    Long64_t n_hel_minus_all = 0;
    Long64_t n_hel_zero_all = 0;

    Long64_t n_scaler_rows = 0;

    double fcup_plus = 0.0;
    double fcup_minus = 0.0;
    double fcup_zero = 0.0;

    double fcupgated_plus = 0.0;
    double fcupgated_minus = 0.0;
    double fcupgated_zero = 0.0;
};

// Container for run level target/polarization information loaded from CSV.
// CSV columns expected:
//   run
//   start_time
//   stop_time
//   species
//   cell
//   charge_avg_online
//   charge_avg_offline
//   charge_avg_offline_err
//   run_dose(Pe/cm2)
// Important:
//   charge_avg_offline is treated as the default analysis polarization.
//   tensor polarization is NOT derived here.

struct TargetInfo {

    int target_match = 0;

    int run = -9999;

    string start_time = "";
    string stop_time = "";

    string species = "";
    int species_id = -9999;

    int cell = -9999;

    float pol_online = -9999;
    float pol_offline = -9999;
    float pol_offline_err = -9999;

    float run_dose = -9999;
};

// Flat ROOT branch variables for target metadata.
// We use this helper so the same branch block can be attached to:
//   sidis
//   electrons
//   inclusive
struct TargetBranchVars {

    Int_t target_match = 0;

    string target_species = "";
    Int_t target_species_id = -9999;

    Int_t target_cell = -9999;

    Float_t target_pol_online = -9999;
    Float_t target_pol_offline = -9999;
    Float_t target_pol_offline_err = -9999;

    Float_t target_run_dose = -9999;

    string target_start_time = "";
    string target_stop_time = "";
};

// Container for inclusive electron scattering kinematics.
// Used by inclusive mode only.
// One output row = one selected electron event.
struct InclusiveKin {

    float Q2 = -9999;
    float nu = -9999;

    // nucleon level invariant mass.
    // This is the W used for DIS/resonance separation.
    float W2 = -9999;
    float W = -9999;

    // Deuteron system invariant mass.
    // Useful to keep for tensor deuteron studies, but not used for W < 2 resonance cuts.
    float Wd2 = -9999;
    float Wd = -9999;

    // Nucleon Bjorken x:
    //   xB = Q2 / (2 M_N nu)
    float xB = -9999;

    // Deuteron x:
    //   xD = Q2 / (2 M_D nu)
    float xD = -9999;

    // In the target rest frame, y = nu / beam_energy.
    float y = -9999;

    int is_dis = 0;
    int is_resonance = 0;

    // Broad inclusive quasi elastic candidate flag.
    // This is intentionally loose at the skim level.
    // Final QE physics cuts will be applied later in rgc_analysis.cc.
    int is_qe = 0;
};

// Forward declaration.
// The full ParticleCandidate struct is defined later before main().
struct ParticleCandidate;

// Functions used
vectorMap loadMapIndex(hipo::bank bankName, int detectorID);

EventInfo getEventInfo(
    hipo::bank bankEvent,
    hipo::bank bankRunConfig,
    int event_file_index,
    int fallback_run,
    int input_file_index,
    int input_file_number
);

bool passElectronDetPidCuts(const DetectorInfo &eDet);

InclusiveKin calculateInclusiveKinematics(
    const ParticleCandidate &electron
);

bool passInclusiveRegion(
    const InclusiveKin &kin,
    const string &inclusiveRegion
);

DCEdgeInfo getDCEdgeInfo(
    int pindex,
    hipo::bank bankTraj
);

bool passInclusiveDISSamplingFraction(
    int sector,
    double p,
    double sampling_fraction
);

bool passInclusiveDISPCALFiducial(
    int sector,
    double lv,
    double lw
);

bool passInclusiveDISDCFiducial(
    const DCEdgeInfo &dcInfo
);

bool passInclusiveDISDetectorCuts(
    const ParticleCandidate &electron,
    const DetectorInfo &eDet,
    const DCEdgeInfo &dcInfo
);

bool passInclusiveDISKinematicCuts(
    const ParticleCandidate &electron,
    const InclusiveKin &kin
);

bool passInclusiveQEDetectorCuts(
    const ParticleCandidate &electron,
    const DetectorInfo &eDet
);

bool passInclusiveQEKinematicCuts(
    const ParticleCandidate &electron,
    const InclusiveKin &kin
);

vector<string> splitCSVLine(
    const string &line
);

string trimString(
    const string &input
);

int getTargetSpeciesID(
    const string &species
);

bool loadTargetMap(
    const string &targetMapFile,
    unordered_map<int, TargetInfo> &targetMap
);

TargetInfo findTargetInfo(
    int run,
    const unordered_map<int, TargetInfo> &targetMap
);

void applyTargetInfoToEventInfo(
    EventInfo &evInfo,
    const TargetInfo &targetInfo,
    const string &polSource
);

void addTargetBranches(
    TTree *tree,
    TargetBranchVars &targetVars
);

void copyTargetInfoToBranches(
    const EventInfo &evInfo,
    TargetBranchVars &targetVars
);

DetectorInfo getDetectorInfo(
    int pindex,
    hipo::bank bankParticle,
    hipo::bank bankCalorimeter,
    hipo::bank bankCherenkov,
    const vectorMap &Calomap,
    const vectorMap &Htccmap,
    float momentum
);

bool passDetPidCuts(
    const DetectorInfo &eDet,
    const DetectorInfo &hDet,
    int selected_hadron_pid
);

// Generic electron hadron SIDIS event test.
// This allows the same skim code to be used for e pi+, e pi-, e K+, e K-, e p, etc.
void eh_sidis_event_test(hipo::bank bankParticle, int selected_hadron_pid, bool &eh_sidis);

void epippim_dihadron_event_test(
    hipo::bank bankParticle,
    bool &dihadron_event
);

void fillScalerTree(
    hipo::bank bankHelScaler,
    const EventInfo &evInfo,
    TTree *scalerTree,
    ScalerBranchVars &scalerVars,
    RunSummaryVars &runSummary
);

// Small helper used only to make readable output ROOT file names.
string getHadronTag(int pid);

int extractInputFileNumber(const string &baseName);
bool stringIsAllDigits(const string &s);

// Container for one selected particle candidate.
// This is not written directly as a vector or object to the ROOT file.
// It is only used internally to remember the best trigger electron in the event.
// The actual ROOT tree branches remain simple scalar branches like pxe, pye, pze, etc.
struct ParticleCandidate {
  bool found = false;

  int pindex = -1;
  int pid = 0;
  int status = 0;
  int det_status = 0;

  float px = -9999;
  float py = -9999;
  float pz = -9999;

  float vx = -9999;
  float vy = -9999;
  float vz = -9999;
  float vt = -9999;

  float p = -9999;
  float theta = -9999;
  float phi = -9999;
};
 

//#################  main program   ###################
//#####################################################

int main(int argc, char** argv) {

   std::cout << " reading HIPO file  "  << __cplusplus << std::endl;

   char inputFile[256];
      
    // Counters used for the final skim summary.
    int total_events=0;
    int n_file=0;
    int sidis_events = 0;

    // Inclusive counters. Inclusive_events_seen counts accepted best electron events before region filtering.
    // inclusive_events counts rows actually written to the inclusive tree.
    int inclusive_events_seen = 0;
    int inclusive_events = 0;

    int inclusive_dis_events = 0;
    int inclusive_res_events = 0;
    int inclusive_qe_events = 0;
    int inclusive_region_rejected = 0;

    int dihadron_events = 0;
    long long dihadron_pairs_written = 0;
    long long hipo_events_written = 0;

    int target_map_entries = 0;

    long long target_matched_input_events = 0;
    long long target_unmatched_input_events = 0;

    long long target_species_p_input_events = 0;
    long long target_species_d_input_events = 0;

    
    if(argc < 2){
        cout << "Usage:" << endl;
        cout << endl;

        cout << "Old SIDIS style:" << endl;
        cout << "  " << argv[0] << " <hadron_pid> <hipo_file1> <hipo_file2> ... [options]" << endl;
        cout << endl;

        cout << "New mode-based style:" << endl;
        cout << "  " << argv[0] << " --mode sidis --pid <hadron_pid> <hipo_file1> <hipo_file2> ... [options]" << endl;
        cout << "  " << argv[0] << " --mode inclusive <hipo_file1> <hipo_file2> ... [options]" << endl;
        cout << endl;

        cout << "Modes:" << endl;
        cout << "  sidis       electron + selected hadron" << endl;
        cout << "  inclusive   electron only" << endl;
        cout << "  dihadron    electron + pi+ + pi-" << endl;
        cout << endl;

        cout << "SIDIS hadron_pid:" << endl;
        cout << "   211   = e pi+" << endl;
        cout << "  -211   = e pi-" << endl;
        cout << "   321   = e K+" << endl;
        cout << "  -321   = e K-" << endl;
        cout << "   2212  = e proton" << endl;
        cout << endl;

        cout << "Options:" << endl;
        cout << "  --electrontree 0 or 1    Write the extra electrons tree in SIDIS mode. Default = 0" << endl;
        cout << "  --detpidcut    0 or 1    Apply detector/PID cuts. Default = 0" << endl;
        cout << "  --kincut       0 or 1    Apply dedicated region kinematic cuts. Default = 0" << endl;
        cout << "  --outformat    root, hipo, or both. Default = root" << endl;
        cout << endl;

        cout << "Examples:" << endl;
        cout << "  " << argv[0] << " --mode sidis --pid 211 file.hipo --electrontree 0 --detpidcut 1" << endl;
        cout << "  " << argv[0] << " --mode sidis --pid 211 file.hipo --detpidcut 1" << endl;
        cout << "  " << argv[0] << " --mode sidis --pid 211 file.hipo --electrontree 1 --detpidcut 1" << endl;
        cout << "  " << argv[0] << " --mode inclusive file.hipo --detpidcut 1 --kincut 1" << endl;

        exit(0);
    }else{
        int n_argc = argc;
        int Run_num = 0;

        // ------------------------------------------------------------
        // Mode parser
        // Backward compatibility:
        //   ./rgcsidis 211 file.hipo ...
        // New style:
        //   ./rgcsidis --mode sidis --pid 211 file.hipo ...
        //   ./rgcsidis --mode inclusive file.hipo ...
        // ------------------------------------------------------------

        string skimMode = "sidis";
        int selected_hadron_pid = 0;

        // inclusive region selector.
        // Used only when skimMode == "inclusive".
        // Allowed values:
        //   all
        //   dis
        //   res
        //   resonance
        //   qe
        //   quasi
        //   quasielastic
        //   quasi elastic
        // Internally, aliases are normalized:
        //   resonance     -> res
        //   quasi*        -> qe
        // Internally, "resonance" is converted to "res".
        string inclusiveRegion = "all";

        // Optional run level target/polarization map.
        // If empty, target map branches stay at default placeholder values.
        string targetMapFile = "";

        // Polarization source used to fill the generic target_pol/vector_pol branches.
        // Allowed values:
        //   offline
        //   online
        // Both raw online and offline values are still stored when a target map is provided.
        string polSource = "offline";

        // Output format.
        // Allowed values:
        //   root
        //   hipo
        //   both
        string outputFormat = "root";
        
        // Output path
        // Optional argument to set custom output paht.
        string outputPath = "./"; //default value

        // Optional diagnostic electron candidate tree.
        // Default = 0:
        //   do not write the extra electrons tree.
        // Use:
        //   --electrontree 1
        // to write all good electron candidates for QA/debugging.
        int writeElectronTree = 0;

        int applyDetPidCut = 0;
        int applyKinCut = 0;

        vector<string> hipoFiles;

        // If first argument is not an option, assume old SIDIS style.
        int iarg_start = 1;

        if(n_argc >= 2){
            string firstArg = argv[1];

            if(firstArg.rfind("--", 0) != 0){
                skimMode = "sidis";
                selected_hadron_pid = atoi(argv[1]);
                iarg_start = 2;
            }
        }

        for(int iarg = iarg_start; iarg < n_argc; ++iarg){

            string arg = argv[iarg];

            if(arg == "--mode"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: --mode requires sidis or inclusive." << endl;
                    exit(1);
                }

                skimMode = argv[iarg + 1];
                iarg++;

            }else if(arg == "--pid" || arg == "--hadronpid"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: " << arg << " requires a hadron PID." << endl;
                    exit(1);
                }

                selected_hadron_pid = atoi(argv[iarg + 1]);
                iarg++;

            }else if(arg == "--electrontree" || arg == "--electronTree"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: " << arg << " requires 0 or 1." << endl;
                    exit(1);
                }

                writeElectronTree = atoi(argv[iarg + 1]);
                iarg++;

            }else if(arg == "--detpidcut"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: --detpidcut requires 0 or 1." << endl;
                    exit(1);
                }

                applyDetPidCut = atoi(argv[iarg + 1]);
                iarg++;

            }else if(arg == "--kincut"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: --kincut requires 0 or 1." << endl;
                    exit(1);
                }

                applyKinCut = atoi(argv[iarg + 1]);
                iarg++;

            }else if(arg == "--region"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: --region requires all, dis, res, or qe." << endl;
                    exit(1);
                }

                inclusiveRegion = argv[iarg + 1];
                iarg++;
            
            }else if(arg == "--targetmap"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: --targetmap requires a CSV file path." << endl;
                    exit(1);
                }

                targetMapFile = argv[iarg + 1];
                iarg++;

            }else if(arg == "--polsource"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: --polsource requires offline or online." << endl;
                    exit(1);
                }

                polSource = argv[iarg + 1];
                iarg++;

            }else if(arg == "--outformat" || arg == "--outputformat"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: " << arg << " requires root, hipo, or both." << endl;
                    exit(1);
                }

                outputFormat = argv[iarg + 1];
                iarg++;

            }else if(arg == "--outpath"){

                if(iarg + 1 >= n_argc){
                    cout << "Error: " << arg << " requires output path." << endl;
                    exit(1);
                }

                outputPath = argv[iarg + 1];
                iarg++;

            }
            else{

                hipoFiles.push_back(arg);
            }
        }

        if(skimMode != "sidis" && skimMode != "inclusive" && skimMode != "dihadron"){
            cout << "Error: unknown mode: " << skimMode << endl;
            cout << "Allowed modes: sidis, inclusive, dihadron" << endl;
            exit(1);
        }

        // Normalize inclusive region option.
        if(inclusiveRegion == "resonance"){
            inclusiveRegion = "res";
        }

        if(inclusiveRegion == "quasi" || inclusiveRegion == "quasielastic" || inclusiveRegion == "quasi-elastic"){
            inclusiveRegion = "qe";
        }

        if(inclusiveRegion != "all" &&
            inclusiveRegion != "dis" &&
            inclusiveRegion != "res" &&
            inclusiveRegion != "qe"){
            cout << "Error: unknown inclusive region: " << inclusiveRegion << endl;
            cout << "Allowed regions: all, dis, res, resonance, qe, quasi, quasielastic, quasi-elastic" << endl;
            exit(1);
        }

        if(skimMode == "sidis" && inclusiveRegion != "all"){
            cout << "Error: --region is only meaningful for --mode inclusive." << endl;
            cout << "For SIDIS mode, do not pass --region." << endl;
            exit(1);
        }

        if(skimMode == "sidis" && selected_hadron_pid == 0){
            cout << "Error: SIDIS mode requires a selected hadron PID." << endl;
            cout << "Example:" << endl;
            cout << "  " << argv[0] << " --mode sidis --pid 211 file.hipo" << endl;
            exit(1);
        }

        if(writeElectronTree != 0 && writeElectronTree != 1){
            cout << "Error: --electrontree must be 0 or 1." << endl;
            exit(1);
        }

        if(applyDetPidCut != 0 && applyDetPidCut != 1){
            cout << "Error: --detpidcut must be 0 or 1." << endl;
            exit(1);
        }

        if(applyKinCut != 0 && applyKinCut != 1){
            cout << "Error: --kincut must be 0 or 1." << endl;
            exit(1);
        }

        if(polSource != "offline" && polSource != "online"){
            cout << "Error: --polsource must be offline or online." << endl;
            exit(1);
        }

        if(outputFormat != "root" && outputFormat != "hipo" && outputFormat != "both"){
            cout << "Error: --outformat must be root, hipo, or both." << endl;
            exit(1);
        }

        bool writeRootOutput = (outputFormat == "root" || outputFormat == "both");
        bool writeHipoOutput = (outputFormat == "hipo" || outputFormat == "both");

        if(hipoFiles.size() == 0){
            cout << "Error: no input HIPO files were provided." << endl;
            exit(1);
        }

        unordered_map<int, TargetInfo> targetMap;

        if(targetMapFile != ""){

            bool loadedTargetMap = loadTargetMap(targetMapFile, targetMap);

            if(!loadedTargetMap){
                cout << "Error: failed to load target map:" << endl;
                cout << "  " << targetMapFile << endl;
                exit(1);
            }

            target_map_entries = targetMap.size();

            cout << "Loaded target map:" << endl;
            cout << "  file    = " << targetMapFile << endl;
            cout << "  entries = " << target_map_entries << endl;
            cout << "  polsource = " << polSource << endl;
        }

        cout << "Options:" << endl;
        cout << "  --electrontree 0 or 1    Write the extra electrons tree in SIDIS mode. Default = 0" << endl;
        cout << "  --detpidcut    0 or 1    Apply detector/PID cuts. Default = 0" << endl;
        cout << "  --kincut       0 or 1    Apply dedicated region kinematic cuts. Default = 0" << endl;
        cout << "  --region all/dis/res/qe  Inclusive mode only. Default = all" << endl;    
        cout << "  --targetmap <csv>       Optional run-level target/polarization CSV" << endl;
        cout << "  --polsource offline|online  Polarization source for target_pol. Default = offline" << endl;
        cout << "  --outformat root|hipo|both  Output format. Default = root" << endl;
        cout << endl;

        cout << "  " << argv[0] << " --mode inclusive file.hipo --region all --detpidcut 1 --kincut 0" << endl;
        cout << "  " << argv[0] << " --mode inclusive file.hipo --region dis --detpidcut 1 --kincut 1" << endl;
        cout << "  " << argv[0] << " --mode inclusive file.hipo --region res --detpidcut 1 --kincut 0" << endl;
        cout << "  " << argv[0] << " --mode inclusive file.hipo --region qe --detpidcut 1 --kincut 1" << endl;

        string hadronTag = "none";

        if(skimMode == "sidis"){
            hadronTag = getHadronTag(selected_hadron_pid);
        }

        cout << "Selected skim mode: " << skimMode << endl;

        if(skimMode == "sidis"){
            cout << "Selected SIDIS channel: e + " << hadronTag
                << "  with hadron PID = " << selected_hadron_pid << endl;
        }else{
            cout << "Selected inclusive electron skim." << endl;
            cout << "Inclusive region: " << inclusiveRegion << endl;
        }

        cout << "Option --electrontree = " << writeElectronTree << endl;
        cout << "Option --detpidcut    = " << applyDetPidCut << endl;
        cout << "Option --kincut       = " << applyKinCut << endl;
        cout << "Option --targetmap    = " << targetMapFile << endl;
        cout << "Option --polsource    = " << polSource << endl;
        cout << "Option --outformat    = " << outputFormat << endl;

        //Loop to include all the input files to analyze
        for(size_t n_argv = 0; n_argv < hipoFiles.size(); ++n_argv){

            sprintf(inputFile,"%s",hipoFiles[n_argv].c_str());
            sscanf(inputFile,"%*87c%d%*c",&Run_num);

            // Old hard cutoff removed. RUN::config.run is the reliable run number and is read later by getEventInfo().
            // if(Run_num > 16432) break;

            cout << "Reading run number: " << Run_num << endl;

            // ---------- Create one ROOT output file per input HIPO file ----------
            // Output file name example:
            //   input HIPO:  clas_016000.evio.00123.hipo
            //   channel:     e pi+
            //   output:      rootfiles/skim_epip_clas_016000.evio.00123.root
            // This keeps the skim files channel specific and file specific.

            std::string inName(inputFile);
            std::string baseName = inName.substr(inName.find_last_of("/\\") + 1);

            const int input_file_index_this = static_cast<int>(n_argv);
            const int input_file_number_this = extractInputFileNumber(baseName);

            cout << "Input HIPO provenance:" << endl;
            cout << "  input_file_index  = " << input_file_index_this << endl;
            cout << "  input_file_number = " << input_file_number_this << endl;
            cout << "  input_file_name   = " << baseName << endl;

            std::string outRootName = baseName;
            size_t hipoPos = outRootName.rfind(".hipo");
            if(hipoPos != std::string::npos){
                outRootName.replace(hipoPos, 5, ".root");
            }else{
                outRootName += ".root";
            }

            std::string outHipoName = baseName;
            hipoPos = outHipoName.rfind(".hipo");
            if(hipoPos == std::string::npos){
                outHipoName += ".hipo";
            }

            if(skimMode == "sidis"){
                outRootName = outputPath + "rootfiles/skim_e" + hadronTag + "_" + outRootName;
                outHipoName = outputPath + "hipofiles/skim_e" + hadronTag + "_" + outHipoName;
            }else if(skimMode == "inclusive"){
                outRootName = outputPath + "rootfiles/skim_inclusive_" + outRootName;
                outHipoName = "hipofiles/skim_inclusive_" + outHipoName;
            }else if(skimMode == "dihadron"){
                outRootName = outputPath + "rootfiles/skim_epippim_" + outRootName;
                outHipoName = outputPath + "hipofiles/skim_epippim_" + outHipoName;
            }

            TFile *fout = nullptr;
            if(writeRootOutput){
                fout = new TFile(outRootName.c_str(), "RECREATE");
            }

            // ------------------------------------------------------------
            // FCup/scaler tree and run summary tree
            // ------------------------------------------------------------
            // These are written for every skim mode. They are independent of electron/hadron cuts, so the analysis stage can reproduce the FCup normalized raw asymmetry used in the comparison macro.

            ScalerBranchVars scalerVars;
            RunSummaryVars runSummary;

            runSummary.input_file_index = input_file_index_this;
            runSummary.input_file_number = input_file_number_this;
            runSummary.input_file_name = baseName;

            // State used to convert cumulative FCup scaler values into per read deltas.
            // These reset for each output skim file.
            // These were used for RUN::scaler cumulative difference logic.
            // The comparison macro does not use that method for Q+ and Q-.
            // double lastScalerFcup = 0.0;
            // double lastScalerFcupGated = 0.0;
            // bool haveLastScalerRead = false;

            TTree *scalerTree = new TTree("scaler", "RUN::scaler FCup information");

            scalerTree->Branch("run",              &scalerVars.run,              "run/I");
            scalerTree->Branch("event",            &scalerVars.event,            "event/I");
            scalerTree->Branch("event_file_index", &scalerVars.event_file_index, "event_file_index/I");

            scalerTree->Branch("input_file_index",  &scalerVars.input_file_index,  "input_file_index/I");
            scalerTree->Branch("input_file_number", &scalerVars.input_file_number, "input_file_number/I");

            scalerTree->Branch("helicity",         &scalerVars.helicity,         "helicity/I");
            scalerTree->Branch("helicity_valid",   &scalerVars.helicity_valid,   "helicity_valid/I");

            scalerTree->Branch("fcup",             &scalerVars.fcup,             "fcup/F");
            scalerTree->Branch("fcupgated",        &scalerVars.fcupgated,        "fcupgated/F");

            scalerTree->Branch("fcup_delta",       &scalerVars.fcup_delta,       "fcup_delta/F");
            scalerTree->Branch("fcupgated_delta",  &scalerVars.fcupgated_delta,  "fcupgated_delta/F");

            TTree *runSummaryTree = new TTree("run_summary", "All-trigger run/file summary");

            runSummaryTree->Branch("run",               &runSummary.run,               "run/I");

            runSummaryTree->Branch("input_file_index",  &runSummary.input_file_index,  "input_file_index/I");
            runSummaryTree->Branch("input_file_number", &runSummary.input_file_number, "input_file_number/I");
            runSummaryTree->Branch("input_file_name",   &runSummary.input_file_name);

            runSummaryTree->Branch("n_events_all",          &runSummary.n_events_all,          "n_events_all/L");
            runSummaryTree->Branch("n_inclusive_written",   &runSummary.n_inclusive_written,   "n_inclusive_written/L");
            runSummaryTree->Branch("n_sidis_rows_written",  &runSummary.n_sidis_rows_written,  "n_sidis_rows_written/L");
            runSummaryTree->Branch("n_dihadron_rows_written",&runSummary.n_dihadron_rows_written,"n_dihadron_rows_written/L");

            runSummaryTree->Branch("first_event",      &runSummary.first_event,      "first_event/I");
            runSummaryTree->Branch("last_event",       &runSummary.last_event,       "last_event/I");
            runSummaryTree->Branch("first_timestamp",  &runSummary.first_timestamp,  "first_timestamp/L");
            runSummaryTree->Branch("last_timestamp",   &runSummary.last_timestamp,   "last_timestamp/L");
            runSummaryTree->Branch("n_hel_plus_all",   &runSummary.n_hel_plus_all,   "n_hel_plus_all/L");
            runSummaryTree->Branch("n_hel_minus_all",  &runSummary.n_hel_minus_all,  "n_hel_minus_all/L");
            runSummaryTree->Branch("n_hel_zero_all",   &runSummary.n_hel_zero_all,   "n_hel_zero_all/L");

            runSummaryTree->Branch("n_scaler_rows",    &runSummary.n_scaler_rows,    "n_scaler_rows/L");

            runSummaryTree->Branch("fcup_plus",        &runSummary.fcup_plus,        "fcup_plus/D");
            runSummaryTree->Branch("fcup_minus",       &runSummary.fcup_minus,       "fcup_minus/D");
            runSummaryTree->Branch("fcup_zero",        &runSummary.fcup_zero,        "fcup_zero/D");

            runSummaryTree->Branch("fcupgated_plus",   &runSummary.fcupgated_plus,   "fcupgated_plus/D");
            runSummaryTree->Branch("fcupgated_minus",  &runSummary.fcupgated_minus,  "fcupgated_minus/D");
            runSummaryTree->Branch("fcupgated_zero",   &runSummary.fcupgated_zero,   "fcupgated_zero/D");

            // One row in this tree corresponds to one electron hadron SIDIS pair.
            TTree *skimTree = new TTree("sidis", "Flat electron hadron SIDIS skim");

            // ---------- Flat scalar branch variables ----------
            // No vectors are used.
            // Each TTree entry stores exactly one electron hadron pair.
            // For example, in an e pi+ skim:
            //   pxe, pye, pze = electron momentum
            //   pxh, pyh, pzh = pi+ momentum
            // In an e K+ skim:
            //   pxe, pye, pze = electron momentum
            //   pxh, pyh, pzh = K+ momentum

            // ---------- Shared event level variables ----------
            // These branches should eventually exist in every skim mode:
            //   sidis
            //   inclusive
            //   exclusive
            // run/event are filled from RUN::config if available.
            // event_file_index is the local counter inside this input file.

            Int_t out_run = 0;
            Int_t out_event = 0;
            Int_t out_event_file_index = 0;

            Int_t out_input_file_index = -9999;
            Int_t out_input_file_number = -9999;

            Int_t out_run_config = -9999;
            Int_t out_event_config = -9999;

            Long64_t out_trigger = -9999;
            Long64_t out_timestamp = -9999;
            Int_t out_unixtime = -9999;

            Int_t out_helicity = -9999;
            Int_t out_helicity_valid = 0;

            // Target/tensor placeholders.
            // These will be filled later from an RGC target state map.
            Int_t out_target_state = -9999;

            Float_t out_target_pol = -9999;
            Float_t out_target_pol_err = -9999;

            Float_t out_tensor_pol = -9999;
            Float_t out_tensor_pol_err = -9999;

            Float_t out_vector_pol = -9999;
            Float_t out_vector_pol_err = -9999;

            TargetBranchVars out_target_extra;

            Int_t out_pair_id = 0;

            // Event level multiplicities for this skim channel.
            // n_good_e = number of good trigger electron candidates in the event.
            // n_good_h = number of selected hadrons in the event.
            // Example:
            //   If selected_hadron_pid = 211, then n_good_h counts good pi+ candidates.
            //   If selected_hadron_pid = -211, then n_good_h counts good pi- candidates.
            //   If selected_hadron_pid = 321, then n_good_h counts good K+ candidates.
            Int_t out_n_good_e = 0;
            Int_t out_n_good_h = 0;

            Int_t out_e_pindex = -1;
            Int_t out_h_pindex = -1;

            Int_t out_pid_e = 11;
            Int_t out_pid_h = selected_hadron_pid;

            Int_t out_e_status = 0;
            Int_t out_h_status = 0;

            Int_t out_e_det_status = 0;
            Int_t out_h_det_status = 0;

            // Electron variables
            Float_t pxe = -9999;
            Float_t pye = -9999;
            Float_t pze = -9999;

            Float_t vxe = -9999;
            Float_t vye = -9999;
            Float_t vze = -9999;
            Float_t vte = -9999;

            Float_t pe = -9999;
            Float_t thetae = -9999;
            Float_t phie = -9999;

            // Hadron variables
            Float_t pxh = -9999;
            Float_t pyh = -9999;
            Float_t pzh = -9999;

            Float_t vxh = -9999;
            Float_t vyh = -9999;
            Float_t vzh = -9999;
            Float_t vth = -9999;

            Float_t ph = -9999;
            Float_t thetah = -9999;
            Float_t phih = -9999;

            // ---------- Electron detector/PID variables ----------
            // These are stored so later analysis can apply cuts without reopening HIPO files.

            Int_t e_charge = -9999;

            Float_t e_beta = -9999;
            Float_t e_chi2pid = -9999;

            Int_t e_has_htcc = 0;
            Float_t e_htcc_nphe = -9999;

            Int_t e_has_cal = 0;
            Int_t e_cal_sector = -9999;

            Float_t e_pcal_e = 0.0;
            Float_t e_ecin_e = 0.0;
            Float_t e_ecout_e = 0.0;
            Float_t e_cal_e = 0.0;
            Float_t e_sampling_fraction = -9999;

            Float_t e_pcal_x = -9999;
            Float_t e_pcal_y = -9999;

            Float_t e_lu = -9999;
            Float_t e_lv = -9999;
            Float_t e_lw = -9999;

            // ---------- Hadron detector/PID variables ----------
            // Here h means the selected hadron channel.
            // For e pi+ skim, h = pi+.
            // For e pi- skim, h = pi-.
            // For e K+ skim, h = K+.

            Int_t h_charge = -9999;

            Float_t h_beta = -9999;
            Float_t h_chi2pid = -9999;

            Int_t h_has_htcc = 0;
            Float_t h_htcc_nphe = -9999;

            Int_t h_has_cal = 0;
            Int_t h_cal_sector = -9999;

            Float_t h_pcal_e = 0.0;
            Float_t h_ecin_e = 0.0;
            Float_t h_ecout_e = 0.0;
            Float_t h_cal_e = 0.0;
            Float_t h_sampling_fraction = -9999;

            Float_t h_pcal_x = -9999;
            Float_t h_pcal_y = -9999;

            Float_t h_lu = -9999;
            Float_t h_lv = -9999;
            Float_t h_lw = -9999;

            // ---------- Create branches ----------

            // ---------- Shared event level branches ----------
            skimTree->Branch("run",              &out_run,              "run/I");
            skimTree->Branch("event",            &out_event,            "event/I");
            skimTree->Branch("event_file_index", &out_event_file_index, "event_file_index/I");

            skimTree->Branch("input_file_index",  &out_input_file_index,  "input_file_index/I");
            skimTree->Branch("input_file_number", &out_input_file_number, "input_file_number/I");

            skimTree->Branch("run_config",   &out_run_config,   "run_config/I");
            skimTree->Branch("event_config", &out_event_config, "event_config/I");

            skimTree->Branch("trigger",   &out_trigger,   "trigger/L");
            skimTree->Branch("timestamp", &out_timestamp, "timestamp/L");
            skimTree->Branch("unixtime",  &out_unixtime,  "unixtime/I");

            skimTree->Branch("helicity",       &out_helicity,       "helicity/I");
            skimTree->Branch("helicity_valid", &out_helicity_valid, "helicity_valid/I");

            // Target/tensor placeholders for b1/Azz.
            // These are intentionally present now so the tree structure is future compatible.
            skimTree->Branch("target_state", &out_target_state, "target_state/I");

            skimTree->Branch("target_pol",     &out_target_pol,     "target_pol/F");
            skimTree->Branch("target_pol_err", &out_target_pol_err, "target_pol_err/F");

            skimTree->Branch("tensor_pol",     &out_tensor_pol,     "tensor_pol/F");
            skimTree->Branch("tensor_pol_err", &out_tensor_pol_err, "tensor_pol_err/F");

            skimTree->Branch("vector_pol",     &out_vector_pol,     "vector_pol/F");
            skimTree->Branch("vector_pol_err", &out_vector_pol_err, "vector_pol_err/F");

            addTargetBranches(skimTree, out_target_extra);

            // SIDIS pair index inside this accepted event.
            skimTree->Branch("pair_id", &out_pair_id, "pair_id/I");

            skimTree->Branch("n_good_e", &out_n_good_e, "n_good_e/I");
            skimTree->Branch("n_good_h", &out_n_good_h, "n_good_h/I");

            skimTree->Branch("e_pindex", &out_e_pindex, "e_pindex/I");
            skimTree->Branch("h_pindex", &out_h_pindex, "h_pindex/I");

            skimTree->Branch("pid_e", &out_pid_e, "pid_e/I");
            skimTree->Branch("pid_h", &out_pid_h, "pid_h/I");

            skimTree->Branch("e_status", &out_e_status, "e_status/I");
            skimTree->Branch("h_status", &out_h_status, "h_status/I");

            skimTree->Branch("e_det_status", &out_e_det_status, "e_det_status/I");
            skimTree->Branch("h_det_status", &out_h_det_status, "h_det_status/I");

            skimTree->Branch("pxe", &pxe, "pxe/F");
            skimTree->Branch("pye", &pye, "pye/F");
            skimTree->Branch("pze", &pze, "pze/F");

            skimTree->Branch("vxe", &vxe, "vxe/F");
            skimTree->Branch("vye", &vye, "vye/F");
            skimTree->Branch("vze", &vze, "vze/F");
            skimTree->Branch("vte", &vte, "vte/F");

            skimTree->Branch("pe",     &pe,     "pe/F");
            skimTree->Branch("thetae", &thetae, "thetae/F");
            skimTree->Branch("phie",   &phie,   "phie/F");

            skimTree->Branch("pxh", &pxh, "pxh/F");
            skimTree->Branch("pyh", &pyh, "pyh/F");
            skimTree->Branch("pzh", &pzh, "pzh/F");

            skimTree->Branch("vxh", &vxh, "vxh/F");
            skimTree->Branch("vyh", &vyh, "vyh/F");
            skimTree->Branch("vzh", &vzh, "vzh/F");
            skimTree->Branch("vth", &vth, "vth/F");

            skimTree->Branch("ph",     &ph,     "ph/F");
            skimTree->Branch("thetah", &thetah, "thetah/F");
            skimTree->Branch("phih",   &phih,   "phih/F");

            // Electron detector/PID branches
            skimTree->Branch("e_charge", &e_charge, "e_charge/I");

            skimTree->Branch("e_beta",    &e_beta,    "e_beta/F");
            skimTree->Branch("e_chi2pid", &e_chi2pid, "e_chi2pid/F");

            skimTree->Branch("e_has_htcc",  &e_has_htcc,  "e_has_htcc/I");
            skimTree->Branch("e_htcc_nphe", &e_htcc_nphe, "e_htcc_nphe/F");

            skimTree->Branch("e_has_cal",   &e_has_cal,   "e_has_cal/I");
            skimTree->Branch("e_cal_sector",&e_cal_sector,"e_cal_sector/I");

            skimTree->Branch("e_pcal_e", &e_pcal_e, "e_pcal_e/F");
            skimTree->Branch("e_ecin_e", &e_ecin_e, "e_ecin_e/F");
            skimTree->Branch("e_ecout_e",&e_ecout_e,"e_ecout_e/F");
            skimTree->Branch("e_cal_e",  &e_cal_e,  "e_cal_e/F");

            skimTree->Branch("e_sampling_fraction", &e_sampling_fraction, "e_sampling_fraction/F");

            skimTree->Branch("e_pcal_x", &e_pcal_x, "e_pcal_x/F");
            skimTree->Branch("e_pcal_y", &e_pcal_y, "e_pcal_y/F");

            skimTree->Branch("e_lu", &e_lu, "e_lu/F");
            skimTree->Branch("e_lv", &e_lv, "e_lv/F");
            skimTree->Branch("e_lw", &e_lw, "e_lw/F");

            // Hadron detector/PID branches
            skimTree->Branch("h_charge", &h_charge, "h_charge/I");

            skimTree->Branch("h_beta",    &h_beta,    "h_beta/F");
            skimTree->Branch("h_chi2pid", &h_chi2pid, "h_chi2pid/F");

            skimTree->Branch("h_has_htcc",  &h_has_htcc,  "h_has_htcc/I");
            skimTree->Branch("h_htcc_nphe", &h_htcc_nphe, "h_htcc_nphe/F");

            skimTree->Branch("h_has_cal",   &h_has_cal,   "h_has_cal/I");
            skimTree->Branch("h_cal_sector",&h_cal_sector,"h_cal_sector/I");

            skimTree->Branch("h_pcal_e", &h_pcal_e, "h_pcal_e/F");
            skimTree->Branch("h_ecin_e", &h_ecin_e, "h_ecin_e/F");
            skimTree->Branch("h_ecout_e",&h_ecout_e,"h_ecout_e/F");
            skimTree->Branch("h_cal_e",  &h_cal_e,  "h_cal_e/F");

            skimTree->Branch("h_sampling_fraction", &h_sampling_fraction, "h_sampling_fraction/F");

            skimTree->Branch("h_pcal_x", &h_pcal_x, "h_pcal_x/F");
            skimTree->Branch("h_pcal_y", &h_pcal_y, "h_pcal_y/F");

            skimTree->Branch("h_lu", &h_lu, "h_lu/F");
            skimTree->Branch("h_lv", &h_lv, "h_lv/F");
            skimTree->Branch("h_lw", &h_lw, "h_lw/F");


            // ---------- Optional electron candidate diagnostic tree ----------
            // The main sidis tree keeps only the best electron paired with each selected hadron.
            // This optional tree stores all good trigger electron candidates from accepted
            // SIDIS events. It is useful for QA checks, especially for events with more than one electron candidate.
            // Default behavior:
            //   --electrontree 0
            //  The tree is not filled or written
            // To enable:
            //   --electrontree 1

            TTree *electronTree = new TTree("electrons", "All good electron candidates in accepted SIDIS events");

            // ---------- Shared event level variables for the electron candidate tree ----------
            Int_t elec_run = 0;
            Int_t elec_event = 0;
            Int_t elec_event_file_index = 0;

            Int_t elec_input_file_index = -9999;
            Int_t elec_input_file_number = -9999;

            Int_t elec_run_config = -9999;
            Int_t elec_event_config = -9999;

            Long64_t elec_trigger = -9999;
            Long64_t elec_timestamp = -9999;
            Int_t elec_unixtime = -9999;

            Int_t elec_helicity = -9999;
            Int_t elec_helicity_valid = 0;

            // Target/tensor placeholders.
            Int_t elec_target_state = -9999;

            Float_t elec_target_pol = -9999;
            Float_t elec_target_pol_err = -9999;

            Float_t elec_tensor_pol = -9999;
            Float_t elec_tensor_pol_err = -9999;

            Float_t elec_vector_pol = -9999;
            Float_t elec_vector_pol_err = -9999;

            TargetBranchVars elec_target_extra;

            Int_t elec_candidate_id = 0;
            Int_t elec_is_best = 0;

            Int_t elec_n_good_e = 0;
            Int_t elec_n_good_h = 0;

            Int_t elec_pindex = -1;
            Int_t elec_pid = 11;
            Int_t elec_status = 0;
            Int_t elec_det_status = 0;

            Float_t elec_px = -9999;
            Float_t elec_py = -9999;
            Float_t elec_pz = -9999;

            Float_t elec_vx = -9999;
            Float_t elec_vy = -9999;
            Float_t elec_vz = -9999;
            Float_t elec_vt = -9999;

            Float_t elec_p = -9999;
            Float_t elec_theta = -9999;
            Float_t elec_phi = -9999;

            Int_t elec_charge = -9999;

            Float_t elec_beta = -9999;
            Float_t elec_chi2pid = -9999;

            Int_t elec_has_htcc = 0;
            Float_t elec_htcc_nphe = -9999;

            Int_t elec_has_cal = 0;
            Int_t elec_cal_sector = -9999;

            Float_t elec_pcal_e = 0.0;
            Float_t elec_ecin_e = 0.0;
            Float_t elec_ecout_e = 0.0;
            Float_t elec_cal_e = 0.0;
            Float_t elec_sampling_fraction = -9999;

            Float_t elec_pcal_x = -9999;
            Float_t elec_pcal_y = -9999;

            Float_t elec_lu = -9999;
            Float_t elec_lv = -9999;
            Float_t elec_lw = -9999;

            // ---------- Shared event level branches for the electron candidate tree ----------
            electronTree->Branch("run",              &elec_run,              "run/I");
            electronTree->Branch("event",            &elec_event,            "event/I");
            electronTree->Branch("event_file_index", &elec_event_file_index, "event_file_index/I");

            electronTree->Branch("input_file_index",  &elec_input_file_index,  "input_file_index/I");
            electronTree->Branch("input_file_number", &elec_input_file_number, "input_file_number/I");

            electronTree->Branch("run_config",   &elec_run_config,   "run_config/I");
            electronTree->Branch("event_config", &elec_event_config, "event_config/I");

            electronTree->Branch("trigger",   &elec_trigger,   "trigger/L");
            electronTree->Branch("timestamp", &elec_timestamp, "timestamp/L");
            electronTree->Branch("unixtime",  &elec_unixtime,  "unixtime/I");

            electronTree->Branch("helicity",       &elec_helicity,       "helicity/I");
            electronTree->Branch("helicity_valid", &elec_helicity_valid, "helicity_valid/I");

            electronTree->Branch("target_state", &elec_target_state, "target_state/I");

            electronTree->Branch("target_pol",     &elec_target_pol,     "target_pol/F");
            electronTree->Branch("target_pol_err", &elec_target_pol_err, "target_pol_err/F");

            electronTree->Branch("tensor_pol",     &elec_tensor_pol,     "tensor_pol/F");
            electronTree->Branch("tensor_pol_err", &elec_tensor_pol_err, "tensor_pol_err/F");

            electronTree->Branch("vector_pol",     &elec_vector_pol,     "vector_pol/F");
            electronTree->Branch("vector_pol_err", &elec_vector_pol_err, "vector_pol_err/F");

            addTargetBranches(electronTree, elec_target_extra);

            electronTree->Branch("e_candidate_id", &elec_candidate_id, "e_candidate_id/I");
            electronTree->Branch("is_best_e",      &elec_is_best,      "is_best_e/I");

            electronTree->Branch("n_good_e", &elec_n_good_e, "n_good_e/I");
            electronTree->Branch("n_good_h", &elec_n_good_h, "n_good_h/I");

            electronTree->Branch("pindex",     &elec_pindex,     "pindex/I");
            electronTree->Branch("pid",        &elec_pid,        "pid/I");
            electronTree->Branch("status",     &elec_status,     "status/I");
            electronTree->Branch("det_status", &elec_det_status, "det_status/I");

            electronTree->Branch("px", &elec_px, "px/F");
            electronTree->Branch("py", &elec_py, "py/F");
            electronTree->Branch("pz", &elec_pz, "pz/F");

            electronTree->Branch("vx", &elec_vx, "vx/F");
            electronTree->Branch("vy", &elec_vy, "vy/F");
            electronTree->Branch("vz", &elec_vz, "vz/F");
            electronTree->Branch("vt", &elec_vt, "vt/F");

            electronTree->Branch("p",     &elec_p,     "p/F");
            electronTree->Branch("theta", &elec_theta, "theta/F");
            electronTree->Branch("phi",   &elec_phi,   "phi/F");

            electronTree->Branch("charge", &elec_charge, "charge/I");

            electronTree->Branch("beta",    &elec_beta,    "beta/F");
            electronTree->Branch("chi2pid", &elec_chi2pid, "chi2pid/F");

            electronTree->Branch("has_htcc",  &elec_has_htcc,  "has_htcc/I");
            electronTree->Branch("htcc_nphe", &elec_htcc_nphe, "htcc_nphe/F");

            electronTree->Branch("has_cal",   &elec_has_cal,   "has_cal/I");
            electronTree->Branch("cal_sector",&elec_cal_sector,"cal_sector/I");

            electronTree->Branch("pcal_e", &elec_pcal_e, "pcal_e/F");
            electronTree->Branch("ecin_e", &elec_ecin_e, "ecin_e/F");
            electronTree->Branch("ecout_e",&elec_ecout_e,"ecout_e/F");
            electronTree->Branch("cal_e",  &elec_cal_e,  "cal_e/F");

            electronTree->Branch("sampling_fraction", &elec_sampling_fraction, "sampling_fraction/F");

            electronTree->Branch("pcal_x", &elec_pcal_x, "pcal_x/F");
            electronTree->Branch("pcal_y", &elec_pcal_y, "pcal_y/F");

            electronTree->Branch("lu", &elec_lu, "lu/F");
            electronTree->Branch("lv", &elec_lv, "lv/F");
            electronTree->Branch("lw", &elec_lw, "lw/F");

            // ============================================================
            // Dihadron electron pi+pi- tree
            // Used only when:
            //   --mode dihadron
            // One row = one accepted electron + one pi+ pi- pair.
            // Important:
            //   This tree is separate from the existing "sidis" tree.
            //   The sidis tree is for e + one selected hadron.
            //   The dihadron tree is for e + pi+ + pi-.
            // First goal:
            //   Store enough information to count statistics and make
            //   basic control plots versus Q2, W, xB, z, Mh, and target state.
            // ============================================================

            TTree *dihadronTree = new TTree("dihadron", "Flat electron pi+ pi- SIDIS skim");

            // ---------- Shared event level variables ----------
            // These mirror the sidis/inclusive event metadata branches.

            Int_t di_run = 0;
            Int_t di_event = 0;
            Int_t di_event_file_index = 0;

            Int_t di_input_file_index = -9999;
            Int_t di_input_file_number = -9999;

            Int_t di_run_config = -9999;
            Int_t di_event_config = -9999;

            Long64_t di_trigger = -9999;
            Long64_t di_timestamp = -9999;
            Int_t di_unixtime = -9999;

            Int_t di_helicity = -9999;
            Int_t di_helicity_valid = 0;

            // ---------- Target / tensor polarization variables ----------
            // These are needed later for tensor asymmetry studies.

            Int_t di_target_state = -9999;

            Float_t di_target_pol = -9999;
            Float_t di_target_pol_err = -9999;

            Float_t di_tensor_pol = -9999;
            Float_t di_tensor_pol_err = -9999;

            Float_t di_vector_pol = -9999;
            Float_t di_vector_pol_err = -9999;

            TargetBranchVars di_target_extra;

            // ---------- Pair and multiplicity information ----------
            // pair_id is the pi+pi- pair index inside this accepted event.
            // Example:
            //   If one event has 2 pi+ and 1 pi-, then it gives two pairs:
            //      pair_id = 0
            //      pair_id = 1

            Int_t di_pair_id = 0;

            Int_t di_n_good_e = 0;
            Int_t di_n_good_pip = 0;
            Int_t di_n_good_pim = 0;
            Int_t di_n_dihadron_pairs = 0;

            // ---------- Particle indices ----------
            // These are REC::Particle row indices.

            Int_t di_e_pindex = -1;
            Int_t di_pip_pindex = -1;
            Int_t di_pim_pindex = -1;

            Int_t di_pid_e = 11;
            Int_t di_pid_pip = 211;
            Int_t di_pid_pim = -211;

            Int_t di_e_status = 0;
            Int_t di_pip_status = 0;
            Int_t di_pim_status = 0;

            Int_t di_e_det_status = 0;
            Int_t di_pip_det_status = 0;
            Int_t di_pim_det_status = 0;

            // ============================================================
            // Electron kinematics
            // ============================================================

            Float_t di_pxe = -9999;
            Float_t di_pye = -9999;
            Float_t di_pze = -9999;

            Float_t di_vxe = -9999;
            Float_t di_vye = -9999;
            Float_t di_vze = -9999;
            Float_t di_vte = -9999;

            Float_t di_pe = -9999;
            Float_t di_thetae = -9999;
            Float_t di_phie = -9999;

            // ---------- Electron detector/PID information ----------

            Int_t di_e_charge = -9999;

            Float_t di_e_beta = -9999;
            Float_t di_e_chi2pid = -9999;

            Int_t di_e_has_htcc = 0;
            Float_t di_e_htcc_nphe = -9999;

            Int_t di_e_has_cal = 0;
            Int_t di_e_cal_sector = -9999;

            Float_t di_e_pcal_e = 0.0;
            Float_t di_e_ecin_e = 0.0;
            Float_t di_e_ecout_e = 0.0;
            Float_t di_e_cal_e = 0.0;
            Float_t di_e_sampling_fraction = -9999;

            Float_t di_e_pcal_x = -9999;
            Float_t di_e_pcal_y = -9999;

            Float_t di_e_lu = -9999;
            Float_t di_e_lv = -9999;
            Float_t di_e_lw = -9999;

            // ============================================================
            // pi+ kinematics
            // ============================================================

            Float_t di_px_pip = -9999;
            Float_t di_py_pip = -9999;
            Float_t di_pz_pip = -9999;

            Float_t di_vx_pip = -9999;
            Float_t di_vy_pip = -9999;
            Float_t di_vz_pip = -9999;
            Float_t di_vt_pip = -9999;

            Float_t di_p_pip = -9999;
            Float_t di_theta_pip = -9999;
            Float_t di_phi_pip = -9999;

            // ---------- pi+ detector/PID information ----------

            Int_t di_pip_charge = -9999;

            Float_t di_pip_beta = -9999;
            Float_t di_pip_chi2pid = -9999;

            Int_t di_pip_has_htcc = 0;
            Float_t di_pip_htcc_nphe = -9999;

            Int_t di_pip_has_cal = 0;
            Int_t di_pip_cal_sector = -9999;

            Float_t di_pip_pcal_e = 0.0;
            Float_t di_pip_ecin_e = 0.0;
            Float_t di_pip_ecout_e = 0.0;
            Float_t di_pip_cal_e = 0.0;
            Float_t di_pip_sampling_fraction = -9999;

            Float_t di_pip_pcal_x = -9999;
            Float_t di_pip_pcal_y = -9999;

            Float_t di_pip_lu = -9999;
            Float_t di_pip_lv = -9999;
            Float_t di_pip_lw = -9999;

            // ============================================================
            // pi- kinematics
            // ============================================================

            Float_t di_px_pim = -9999;
            Float_t di_py_pim = -9999;
            Float_t di_pz_pim = -9999;

            Float_t di_vx_pim = -9999;
            Float_t di_vy_pim = -9999;
            Float_t di_vz_pim = -9999;
            Float_t di_vt_pim = -9999;

            Float_t di_p_pim = -9999;
            Float_t di_theta_pim = -9999;
            Float_t di_phi_pim = -9999;

            // ---------- pi- detector/PID information ----------

            Int_t di_pim_charge = -9999;

            Float_t di_pim_beta = -9999;
            Float_t di_pim_chi2pid = -9999;

            Int_t di_pim_has_htcc = 0;
            Float_t di_pim_htcc_nphe = -9999;

            Int_t di_pim_has_cal = 0;
            Int_t di_pim_cal_sector = -9999;

            Float_t di_pim_pcal_e = 0.0;
            Float_t di_pim_ecin_e = 0.0;
            Float_t di_pim_ecout_e = 0.0;
            Float_t di_pim_cal_e = 0.0;
            Float_t di_pim_sampling_fraction = -9999;

            Float_t di_pim_pcal_x = -9999;
            Float_t di_pim_pcal_y = -9999;

            Float_t di_pim_lu = -9999;
            Float_t di_pim_lv = -9999;
            Float_t di_pim_lw = -9999;

            // ============================================================
            // Inclusive DIS kinematics from the electron
            // ============================================================

            Float_t di_Q2 = -9999;
            Float_t di_nu = -9999;

            // nucleon level W.
            Float_t di_W2 = -9999;
            Float_t di_W = -9999;

            // Deuteron level W.
            // Keep this separately for tensor deuteron studies.
            Float_t di_Wd2 = -9999;
            Float_t di_Wd = -9999;

            Float_t di_xB = -9999;
            Float_t di_xD = -9999;
            Float_t di_y = -9999;

            Int_t di_is_dis = 0;
            Int_t di_is_resonance = 0;

            // ============================================================
            // Dihadron pair variables
            // ============================================================

            Float_t di_pair_px = -9999;
            Float_t di_pair_py = -9999;
            Float_t di_pair_pz = -9999;
            Float_t di_pair_E = -9999;
            Float_t di_pair_p = -9999;

            Float_t di_pair_theta = -9999;
            Float_t di_pair_phi = -9999;

            // Invariant mass of pi+pi-.
            Float_t di_Mh = -9999;

            // SIDIS energy fractions.
            Float_t di_z_pip = -9999;
            Float_t di_z_pim = -9999;
            Float_t di_z_pair = -9999;

            // Transverse momentum of pair relative to virtual photon.
            // We can fill this in the next step after defining q vector logic.
            Float_t di_PhT = -9999;

            // Dihadron azimuthal variables.
            // For now these are placeholders. We will fill them carefully later.
            Float_t di_phi_h = -9999;
            Float_t di_phi_R = -9999;

            // ============================================================
            // Create dihadron branches
            // ============================================================

            // ---------- Shared event level branches ----------
            dihadronTree->Branch("run",              &di_run,              "run/I");
            dihadronTree->Branch("event",            &di_event,            "event/I");
            dihadronTree->Branch("event_file_index", &di_event_file_index, "event_file_index/I");

            dihadronTree->Branch("input_file_index",  &di_input_file_index,  "input_file_index/I");
            dihadronTree->Branch("input_file_number", &di_input_file_number, "input_file_number/I");

            dihadronTree->Branch("run_config",   &di_run_config,   "run_config/I");
            dihadronTree->Branch("event_config", &di_event_config, "event_config/I");

            dihadronTree->Branch("trigger",   &di_trigger,   "trigger/L");
            dihadronTree->Branch("timestamp", &di_timestamp, "timestamp/L");
            dihadronTree->Branch("unixtime",  &di_unixtime,  "unixtime/I");

            dihadronTree->Branch("helicity",       &di_helicity,       "helicity/I");
            dihadronTree->Branch("helicity_valid", &di_helicity_valid, "helicity_valid/I");

            // ---------- Target/tensor branches ----------
            dihadronTree->Branch("target_state", &di_target_state, "target_state/I");

            dihadronTree->Branch("target_pol",     &di_target_pol,     "target_pol/F");
            dihadronTree->Branch("target_pol_err", &di_target_pol_err, "target_pol_err/F");

            dihadronTree->Branch("tensor_pol",     &di_tensor_pol,     "tensor_pol/F");
            dihadronTree->Branch("tensor_pol_err", &di_tensor_pol_err, "tensor_pol_err/F");

            dihadronTree->Branch("vector_pol",     &di_vector_pol,     "vector_pol/F");
            dihadronTree->Branch("vector_pol_err", &di_vector_pol_err, "vector_pol_err/F");

            addTargetBranches(dihadronTree, di_target_extra);

            // ---------- Pair/multiplicity branches ----------
            dihadronTree->Branch("pair_id", &di_pair_id, "pair_id/I");

            dihadronTree->Branch("n_good_e",          &di_n_good_e,          "n_good_e/I");
            dihadronTree->Branch("n_good_pip",        &di_n_good_pip,        "n_good_pip/I");
            dihadronTree->Branch("n_good_pim",        &di_n_good_pim,        "n_good_pim/I");
            dihadronTree->Branch("n_dihadron_pairs",  &di_n_dihadron_pairs,  "n_dihadron_pairs/I");

            // ---------- Particle index and status branches ----------
            dihadronTree->Branch("e_pindex",   &di_e_pindex,   "e_pindex/I");
            dihadronTree->Branch("pip_pindex", &di_pip_pindex, "pip_pindex/I");
            dihadronTree->Branch("pim_pindex", &di_pim_pindex, "pim_pindex/I");

            dihadronTree->Branch("pid_e",   &di_pid_e,   "pid_e/I");
            dihadronTree->Branch("pid_pip", &di_pid_pip, "pid_pip/I");
            dihadronTree->Branch("pid_pim", &di_pid_pim, "pid_pim/I");

            dihadronTree->Branch("e_status",   &di_e_status,   "e_status/I");
            dihadronTree->Branch("pip_status", &di_pip_status, "pip_status/I");
            dihadronTree->Branch("pim_status", &di_pim_status, "pim_status/I");

            dihadronTree->Branch("e_det_status",   &di_e_det_status,   "e_det_status/I");
            dihadronTree->Branch("pip_det_status", &di_pip_det_status, "pip_det_status/I");
            dihadronTree->Branch("pim_det_status", &di_pim_det_status, "pim_det_status/I");

            // ---------- Electron branches ----------
            dihadronTree->Branch("pxe", &di_pxe, "pxe/F");
            dihadronTree->Branch("pye", &di_pye, "pye/F");
            dihadronTree->Branch("pze", &di_pze, "pze/F");

            dihadronTree->Branch("vxe", &di_vxe, "vxe/F");
            dihadronTree->Branch("vye", &di_vye, "vye/F");
            dihadronTree->Branch("vze", &di_vze, "vze/F");
            dihadronTree->Branch("vte", &di_vte, "vte/F");

            dihadronTree->Branch("pe",     &di_pe,     "pe/F");
            dihadronTree->Branch("thetae", &di_thetae, "thetae/F");
            dihadronTree->Branch("phie",   &di_phie,   "phie/F");

            dihadronTree->Branch("e_charge", &di_e_charge, "e_charge/I");

            dihadronTree->Branch("e_beta",    &di_e_beta,    "e_beta/F");
            dihadronTree->Branch("e_chi2pid", &di_e_chi2pid, "e_chi2pid/F");

            dihadronTree->Branch("e_has_htcc",  &di_e_has_htcc,  "e_has_htcc/I");
            dihadronTree->Branch("e_htcc_nphe", &di_e_htcc_nphe, "e_htcc_nphe/F");

            dihadronTree->Branch("e_has_cal",    &di_e_has_cal,    "e_has_cal/I");
            dihadronTree->Branch("e_cal_sector", &di_e_cal_sector, "e_cal_sector/I");

            dihadronTree->Branch("e_pcal_e",  &di_e_pcal_e,  "e_pcal_e/F");
            dihadronTree->Branch("e_ecin_e",  &di_e_ecin_e,  "e_ecin_e/F");
            dihadronTree->Branch("e_ecout_e", &di_e_ecout_e, "e_ecout_e/F");
            dihadronTree->Branch("e_cal_e",   &di_e_cal_e,   "e_cal_e/F");

            dihadronTree->Branch("e_sampling_fraction", &di_e_sampling_fraction, "e_sampling_fraction/F");

            dihadronTree->Branch("e_pcal_x", &di_e_pcal_x, "e_pcal_x/F");
            dihadronTree->Branch("e_pcal_y", &di_e_pcal_y, "e_pcal_y/F");

            dihadronTree->Branch("e_lu", &di_e_lu, "e_lu/F");
            dihadronTree->Branch("e_lv", &di_e_lv, "e_lv/F");
            dihadronTree->Branch("e_lw", &di_e_lw, "e_lw/F");

            // ---------- pi+ branches ----------
            dihadronTree->Branch("px_pip", &di_px_pip, "px_pip/F");
            dihadronTree->Branch("py_pip", &di_py_pip, "py_pip/F");
            dihadronTree->Branch("pz_pip", &di_pz_pip, "pz_pip/F");

            dihadronTree->Branch("vx_pip", &di_vx_pip, "vx_pip/F");
            dihadronTree->Branch("vy_pip", &di_vy_pip, "vy_pip/F");
            dihadronTree->Branch("vz_pip", &di_vz_pip, "vz_pip/F");
            dihadronTree->Branch("vt_pip", &di_vt_pip, "vt_pip/F");

            dihadronTree->Branch("p_pip",     &di_p_pip,     "p_pip/F");
            dihadronTree->Branch("theta_pip", &di_theta_pip, "theta_pip/F");
            dihadronTree->Branch("phi_pip",   &di_phi_pip,   "phi_pip/F");

            dihadronTree->Branch("pip_charge", &di_pip_charge, "pip_charge/I");

            dihadronTree->Branch("pip_beta",    &di_pip_beta,    "pip_beta/F");
            dihadronTree->Branch("pip_chi2pid", &di_pip_chi2pid, "pip_chi2pid/F");

            dihadronTree->Branch("pip_has_htcc",  &di_pip_has_htcc,  "pip_has_htcc/I");
            dihadronTree->Branch("pip_htcc_nphe", &di_pip_htcc_nphe, "pip_htcc_nphe/F");

            dihadronTree->Branch("pip_has_cal",    &di_pip_has_cal,    "pip_has_cal/I");
            dihadronTree->Branch("pip_cal_sector", &di_pip_cal_sector, "pip_cal_sector/I");

            dihadronTree->Branch("pip_pcal_e",  &di_pip_pcal_e,  "pip_pcal_e/F");
            dihadronTree->Branch("pip_ecin_e",  &di_pip_ecin_e,  "pip_ecin_e/F");
            dihadronTree->Branch("pip_ecout_e", &di_pip_ecout_e, "pip_ecout_e/F");
            dihadronTree->Branch("pip_cal_e",   &di_pip_cal_e,   "pip_cal_e/F");

            dihadronTree->Branch("pip_sampling_fraction", &di_pip_sampling_fraction, "pip_sampling_fraction/F");

            dihadronTree->Branch("pip_pcal_x", &di_pip_pcal_x, "pip_pcal_x/F");
            dihadronTree->Branch("pip_pcal_y", &di_pip_pcal_y, "pip_pcal_y/F");

            dihadronTree->Branch("pip_lu", &di_pip_lu, "pip_lu/F");
            dihadronTree->Branch("pip_lv", &di_pip_lv, "pip_lv/F");
            dihadronTree->Branch("pip_lw", &di_pip_lw, "pip_lw/F");

            // ---------- pi- branches ----------
            dihadronTree->Branch("px_pim", &di_px_pim, "px_pim/F");
            dihadronTree->Branch("py_pim", &di_py_pim, "py_pim/F");
            dihadronTree->Branch("pz_pim", &di_pz_pim, "pz_pim/F");

            dihadronTree->Branch("vx_pim", &di_vx_pim, "vx_pim/F");
            dihadronTree->Branch("vy_pim", &di_vy_pim, "vy_pim/F");
            dihadronTree->Branch("vz_pim", &di_vz_pim, "vz_pim/F");
            dihadronTree->Branch("vt_pim", &di_vt_pim, "vt_pim/F");

            dihadronTree->Branch("p_pim",     &di_p_pim,     "p_pim/F");
            dihadronTree->Branch("theta_pim", &di_theta_pim, "theta_pim/F");
            dihadronTree->Branch("phi_pim",   &di_phi_pim,   "phi_pim/F");

            dihadronTree->Branch("pim_charge", &di_pim_charge, "pim_charge/I");

            dihadronTree->Branch("pim_beta",    &di_pim_beta,    "pim_beta/F");
            dihadronTree->Branch("pim_chi2pid", &di_pim_chi2pid, "pim_chi2pid/F");

            dihadronTree->Branch("pim_has_htcc",  &di_pim_has_htcc,  "pim_has_htcc/I");
            dihadronTree->Branch("pim_htcc_nphe", &di_pim_htcc_nphe, "pim_htcc_nphe/F");

            dihadronTree->Branch("pim_has_cal",    &di_pim_has_cal,    "pim_has_cal/I");
            dihadronTree->Branch("pim_cal_sector", &di_pim_cal_sector, "pim_cal_sector/I");

            dihadronTree->Branch("pim_pcal_e",  &di_pim_pcal_e,  "pim_pcal_e/F");
            dihadronTree->Branch("pim_ecin_e",  &di_pim_ecin_e,  "pim_ecin_e/F");
            dihadronTree->Branch("pim_ecout_e", &di_pim_ecout_e, "pim_ecout_e/F");
            dihadronTree->Branch("pim_cal_e",   &di_pim_cal_e,   "pim_cal_e/F");

            dihadronTree->Branch("pim_sampling_fraction", &di_pim_sampling_fraction, "pim_sampling_fraction/F");

            dihadronTree->Branch("pim_pcal_x", &di_pim_pcal_x, "pim_pcal_x/F");
            dihadronTree->Branch("pim_pcal_y", &di_pim_pcal_y, "pim_pcal_y/F");

            dihadronTree->Branch("pim_lu", &di_pim_lu, "pim_lu/F");
            dihadronTree->Branch("pim_lv", &di_pim_lv, "pim_lv/F");
            dihadronTree->Branch("pim_lw", &di_pim_lw, "pim_lw/F");

            // ---------- DIS kinematic branches ----------
            dihadronTree->Branch("Q2", &di_Q2, "Q2/F");
            dihadronTree->Branch("nu", &di_nu, "nu/F");

            dihadronTree->Branch("W2", &di_W2, "W2/F");
            dihadronTree->Branch("W",  &di_W,  "W/F");

            dihadronTree->Branch("Wd2", &di_Wd2, "Wd2/F");
            dihadronTree->Branch("Wd",  &di_Wd,  "Wd/F");

            dihadronTree->Branch("xB", &di_xB, "xB/F");
            dihadronTree->Branch("xD", &di_xD, "xD/F");
            dihadronTree->Branch("y",  &di_y,  "y/F");

            dihadronTree->Branch("is_dis",       &di_is_dis,       "is_dis/I");
            dihadronTree->Branch("is_resonance", &di_is_resonance, "is_resonance/I");

            // ---------- Dihadron pair branches ----------
            dihadronTree->Branch("pair_px", &di_pair_px, "pair_px/F");
            dihadronTree->Branch("pair_py", &di_pair_py, "pair_py/F");
            dihadronTree->Branch("pair_pz", &di_pair_pz, "pair_pz/F");
            dihadronTree->Branch("pair_E",  &di_pair_E,  "pair_E/F");
            dihadronTree->Branch("pair_p",  &di_pair_p,  "pair_p/F");

            dihadronTree->Branch("pair_theta", &di_pair_theta, "pair_theta/F");
            dihadronTree->Branch("pair_phi",   &di_pair_phi,   "pair_phi/F");

            dihadronTree->Branch("Mh", &di_Mh, "Mh/F");

            dihadronTree->Branch("z_pip",  &di_z_pip,  "z_pip/F");
            dihadronTree->Branch("z_pim",  &di_z_pim,  "z_pim/F");
            dihadronTree->Branch("z_pair", &di_z_pair, "z_pair/F");

            dihadronTree->Branch("PhT", &di_PhT, "PhT/F");

            dihadronTree->Branch("phi_h", &di_phi_h, "phi_h/F");
            dihadronTree->Branch("phi_R", &di_phi_R, "phi_R/F");

            // ============================================================
            // Inclusive electron only tree
            // Used only when:
            //   --mode inclusive
            // One row = one selected best electron event.
            // No hadron is required.
            // DIS/resonance region filtering is applied in Phase 3 through --region.
            // Regions are based on nucleon level W:
            //   all = write all accepted inclusive electron events
            //   dis = write W > 2
            //   res = write 0 < W < 2
            // Deuteron system W is stored separately as Wd/Wd2.
            // ============================================================

            TTree *inclusiveTree = new TTree("inclusive", "Flat inclusive electron skim");

            // Shared event level variables
            Int_t inc_run = 0;
            Int_t inc_event = 0;
            Int_t inc_event_file_index = 0;

            Int_t inc_input_file_index = -9999;
            Int_t inc_input_file_number = -9999;

            Int_t inc_run_config = -9999;
            Int_t inc_event_config = -9999;

            Long64_t inc_trigger = -9999;
            Long64_t inc_timestamp = -9999;
            Int_t inc_unixtime = -9999;

            Int_t inc_helicity = -9999;
            Int_t inc_helicity_valid = 0;

            Int_t inc_target_state = -9999;

            Float_t inc_target_pol = -9999;
            Float_t inc_target_pol_err = -9999;

            Float_t inc_tensor_pol = -9999;
            Float_t inc_tensor_pol_err = -9999;

            Float_t inc_vector_pol = -9999;
            Float_t inc_vector_pol_err = -9999;

            TargetBranchVars inc_target_extra;

            // Electron candidate information
            Int_t inc_n_good_e = 0;

            Int_t inc_e_pindex = -1;
            Int_t inc_pid_e = 11;
            Int_t inc_e_status = 0;
            Int_t inc_e_det_status = 0;

            Float_t inc_pxe = -9999;
            Float_t inc_pye = -9999;
            Float_t inc_pze = -9999;

            Float_t inc_vxe = -9999;
            Float_t inc_vye = -9999;
            Float_t inc_vze = -9999;
            Float_t inc_vte = -9999;

            Float_t inc_pe = -9999;
            Float_t inc_thetae = -9999;
            Float_t inc_phie = -9999;

            // Electron detector/PID information
            Int_t inc_e_charge = -9999;

            Float_t inc_e_beta = -9999;
            Float_t inc_e_chi2pid = -9999;

            Int_t inc_e_has_htcc = 0;
            Float_t inc_e_htcc_nphe = -9999;

            Int_t inc_e_has_cal = 0;
            Int_t inc_e_cal_sector = -9999;

            Float_t inc_e_pcal_e = 0.0;
            Float_t inc_e_ecin_e = 0.0;
            Float_t inc_e_ecout_e = 0.0;
            Float_t inc_e_cal_e = 0.0;
            Float_t inc_e_sampling_fraction = -9999;

            Float_t inc_e_pcal_x = -9999;
            Float_t inc_e_pcal_y = -9999;

            Float_t inc_e_lu = -9999;
            Float_t inc_e_lv = -9999;
            Float_t inc_e_lw = -9999;

            // Inclusive kinematics
            Float_t inc_Q2 = -9999;
            Float_t inc_nu = -9999;

            // nucleon level W, used for DIS/resonance separation.
            Float_t inc_W2 = -9999;
            Float_t inc_W = -9999;

            // deuteron system W, kept for tensor deuteron studies.
            Float_t inc_Wd2 = -9999;
            Float_t inc_Wd = -9999;

            Float_t inc_xB = -9999;
            Float_t inc_xD = -9999;
            Float_t inc_y = -9999;

            Int_t inc_is_dis = 0;
            Int_t inc_is_resonance = 0;
            Int_t inc_is_qe = 0;

            // Branches
            inclusiveTree->Branch("run",              &inc_run,              "run/I");
            inclusiveTree->Branch("event",            &inc_event,            "event/I");
            inclusiveTree->Branch("event_file_index", &inc_event_file_index, "event_file_index/I");

            inclusiveTree->Branch("input_file_index",  &inc_input_file_index,  "input_file_index/I");
            inclusiveTree->Branch("input_file_number", &inc_input_file_number, "input_file_number/I");

            inclusiveTree->Branch("run_config",   &inc_run_config,   "run_config/I");
            inclusiveTree->Branch("event_config", &inc_event_config, "event_config/I");

            inclusiveTree->Branch("trigger",   &inc_trigger,   "trigger/L");
            inclusiveTree->Branch("timestamp", &inc_timestamp, "timestamp/L");
            inclusiveTree->Branch("unixtime",  &inc_unixtime,  "unixtime/I");

            inclusiveTree->Branch("helicity",       &inc_helicity,       "helicity/I");
            inclusiveTree->Branch("helicity_valid", &inc_helicity_valid, "helicity_valid/I");

            inclusiveTree->Branch("target_state", &inc_target_state, "target_state/I");

            inclusiveTree->Branch("target_pol",     &inc_target_pol,     "target_pol/F");
            inclusiveTree->Branch("target_pol_err", &inc_target_pol_err, "target_pol_err/F");

            inclusiveTree->Branch("tensor_pol",     &inc_tensor_pol,     "tensor_pol/F");
            inclusiveTree->Branch("tensor_pol_err", &inc_tensor_pol_err, "tensor_pol_err/F");

            inclusiveTree->Branch("vector_pol",     &inc_vector_pol,     "vector_pol/F");
            inclusiveTree->Branch("vector_pol_err", &inc_vector_pol_err, "vector_pol_err/F");

            addTargetBranches(inclusiveTree, inc_target_extra);

            inclusiveTree->Branch("n_good_e", &inc_n_good_e, "n_good_e/I");

            inclusiveTree->Branch("e_pindex",     &inc_e_pindex,     "e_pindex/I");
            inclusiveTree->Branch("pid_e",        &inc_pid_e,        "pid_e/I");
            inclusiveTree->Branch("e_status",     &inc_e_status,     "e_status/I");
            inclusiveTree->Branch("e_det_status", &inc_e_det_status, "e_det_status/I");

            inclusiveTree->Branch("pxe", &inc_pxe, "pxe/F");
            inclusiveTree->Branch("pye", &inc_pye, "pye/F");
            inclusiveTree->Branch("pze", &inc_pze, "pze/F");

            inclusiveTree->Branch("vxe", &inc_vxe, "vxe/F");
            inclusiveTree->Branch("vye", &inc_vye, "vye/F");
            inclusiveTree->Branch("vze", &inc_vze, "vze/F");
            inclusiveTree->Branch("vte", &inc_vte, "vte/F");

            inclusiveTree->Branch("pe",     &inc_pe,     "pe/F");
            inclusiveTree->Branch("thetae", &inc_thetae, "thetae/F");
            inclusiveTree->Branch("phie",   &inc_phie,   "phie/F");

            inclusiveTree->Branch("e_charge", &inc_e_charge, "e_charge/I");

            inclusiveTree->Branch("e_beta",    &inc_e_beta,    "e_beta/F");
            inclusiveTree->Branch("e_chi2pid", &inc_e_chi2pid, "e_chi2pid/F");

            inclusiveTree->Branch("e_has_htcc",  &inc_e_has_htcc,  "e_has_htcc/I");
            inclusiveTree->Branch("e_htcc_nphe", &inc_e_htcc_nphe, "e_htcc_nphe/F");

            inclusiveTree->Branch("e_has_cal",    &inc_e_has_cal,    "e_has_cal/I");
            inclusiveTree->Branch("e_cal_sector", &inc_e_cal_sector, "e_cal_sector/I");

            inclusiveTree->Branch("e_pcal_e",  &inc_e_pcal_e,  "e_pcal_e/F");
            inclusiveTree->Branch("e_ecin_e",  &inc_e_ecin_e,  "e_ecin_e/F");
            inclusiveTree->Branch("e_ecout_e", &inc_e_ecout_e, "e_ecout_e/F");
            inclusiveTree->Branch("e_cal_e",   &inc_e_cal_e,   "e_cal_e/F");

            inclusiveTree->Branch("e_sampling_fraction", &inc_e_sampling_fraction, "e_sampling_fraction/F");

            inclusiveTree->Branch("e_pcal_x", &inc_e_pcal_x, "e_pcal_x/F");
            inclusiveTree->Branch("e_pcal_y", &inc_e_pcal_y, "e_pcal_y/F");

            inclusiveTree->Branch("e_lu", &inc_e_lu, "e_lu/F");
            inclusiveTree->Branch("e_lv", &inc_e_lv, "e_lv/F");
            inclusiveTree->Branch("e_lw", &inc_e_lw, "e_lw/F");

            inclusiveTree->Branch("Q2", &inc_Q2, "Q2/F");
            inclusiveTree->Branch("nu", &inc_nu, "nu/F");

            // nucleon level W.
            // This is the W used for DIS/resonance classification.
            inclusiveTree->Branch("W2", &inc_W2, "W2/F");
            inclusiveTree->Branch("W",  &inc_W,  "W/F");

            // deuteron system invariant mass.
            // Kept separately so we do not confuse deuteron W with nucleon W.
            inclusiveTree->Branch("Wd2", &inc_Wd2, "Wd2/F");
            inclusiveTree->Branch("Wd",  &inc_Wd,  "Wd/F");

            inclusiveTree->Branch("xB", &inc_xB, "xB/F");
            inclusiveTree->Branch("xD", &inc_xD, "xD/F");
            inclusiveTree->Branch("y",  &inc_y,  "y/F");

            inclusiveTree->Branch("is_dis",       &inc_is_dis,       "is_dis/I");
            inclusiveTree->Branch("is_resonance", &inc_is_resonance, "is_resonance/I");

            // Broad inclusive quasi elastic flag.
            // Final tensor QE cuts are applied later in analysis.
            inclusiveTree->Branch("is_qe", &inc_is_qe, "is_qe/I");

            //#### necessary lines with hipo:: to read hipo files from CLAS12 ######    
            hipo::reader  reader;
            reader.open(inputFile);
            hipo::dictionary  factory;
            reader.readDictionary(factory);

            hipo::writer hipoWriter;
            long long hipo_events_written_this_file = 0;

            if(writeHipoOutput){
                hipoWriter.getDictionary().addSchema(factory.getSchema("REC::Particle"));
                hipoWriter.getDictionary().addSchema(factory.getSchema("REC::Calorimeter"));
                hipoWriter.getDictionary().addSchema(factory.getSchema("REC::Cherenkov"));
                hipoWriter.getDictionary().addSchema(factory.getSchema("REC::Traj"));
                hipoWriter.getDictionary().addSchema(factory.getSchema("REC::Event"));
                hipoWriter.getDictionary().addSchema(factory.getSchema("RUN::config"));
                hipoWriter.getDictionary().addSchema(factory.getSchema("HEL::scaler"));
                hipoWriter.open(outHipoName.c_str());
            }

            hipo::bank  particles(factory.getSchema("REC::Particle"));
            hipo::bank  calorimeter(factory.getSchema("REC::Calorimeter"));
            hipo::bank  cherenkov(factory.getSchema("REC::Cherenkov"));
            hipo::bank  trajectory(factory.getSchema("REC::Traj"));

            // event level banks.
            // REC::Event provides helicity and event level reconstruction quantities.
            // RUN::config provides run/event/trigger/timestamp metadata.
            hipo::bank  recEvent(factory.getSchema("REC::Event"));
            hipo::bank  runConfig(factory.getSchema("RUN::config"));

            // HEL::scaler is used by FCup_reading_summer2026_V1.cc
            // for helicity dependent FCup charge normalization.
            hipo::bank  helScaler(factory.getSchema("HEL::scaler"));

            hipo::event      event;
    
            ++n_file;
            cout << "***** Reading file: " << inputFile << "  ***** " << "file processing: " << n_file << " out of " << hipoFiles.size() << endl;

            int counter = 0;
    
            //####################################################
            //##### Begining loop running through events #########

            while(reader.next()==true){
                reader.read(event);

                event.getStructure(particles);
                event.getStructure(calorimeter);
                event.getStructure(cherenkov);
                event.getStructure(trajectory);

                // event level banks.
                // These are needed for helicity and run/event metadata.
                event.getStructure(recEvent);
                event.getStructure(runConfig);
                event.getStructure(helScaler);

                bool hipo_event_written_this_event = false;

                EventInfo evInfo = getEventInfo(
                    recEvent,
                    runConfig,
                    counter,
                    Run_num,
                    input_file_index_this,
                    input_file_number_this
                );

                // All trigger helicity counters for comparison with the FCup macro.
                // These are counted before any electron, SIDIS, detector, or region selection.
                if(evInfo.run > 0 && runSummary.run <= 0){
                    runSummary.run = evInfo.run;
                }
                runSummary.n_events_all++;

                if(runSummary.first_event < 0 || evInfo.event < runSummary.first_event){
                    runSummary.first_event = evInfo.event;
                }

                if(runSummary.last_event < 0 || evInfo.event > runSummary.last_event){
                    runSummary.last_event = evInfo.event;
                }

                if(evInfo.timestamp > 0){
                    if(runSummary.first_timestamp < 0 || evInfo.timestamp < runSummary.first_timestamp){
                        runSummary.first_timestamp = evInfo.timestamp;
                    }

                    if(runSummary.last_timestamp < 0 || evInfo.timestamp > runSummary.last_timestamp){
                        runSummary.last_timestamp = evInfo.timestamp;
                    }
                }

                if(evInfo.helicity == 1){
                    runSummary.n_hel_plus_all++;
                }else if(evInfo.helicity == -1){
                    runSummary.n_hel_minus_all++;
                }else{
                    runSummary.n_hel_zero_all++;
                }

                // Store all scaler rows from this HIPO event.
                // This is independent of whether the event later passes the skim selection.
                if(writeRootOutput){
                    fillScalerTree(
                        helScaler,
                        evInfo,
                        scalerTree,
                        scalerVars,
                        runSummary
                    );
                }

                // Keep scaler events in HIPO output so FCup/helicity normalization survives HIPO readback.
                if(writeHipoOutput && helScaler.getRows() > 0){
                    hipoWriter.addEvent(event);
                    hipo_events_written++;
                    hipo_events_written_this_file++;
                    hipo_event_written_this_event = true;
                }

                if(targetMapFile != ""){

                    TargetInfo targetInfo = findTargetInfo(evInfo.run, targetMap);

                    if(targetInfo.target_match == 1){

                        target_matched_input_events++;

                        if(targetInfo.species == "p"){
                            target_species_p_input_events++;
                        }else if(targetInfo.species == "d"){
                            target_species_d_input_events++;
                        }

                    }else{

                        target_unmatched_input_events++;
                    }

                    applyTargetInfoToEventInfo(
                        evInfo,
                        targetInfo,
                        polSource
                    );
                }
      
	            if(total_events%10000==0) cout<<"Total events processed: "<<total_events<<"\r"<<flush;
                //cout<<"new event:"<<counter<<"\n";
                //if(total_events==100) break; //debugging
	
                // ------------------------------------------------------------
                // Mode dependent event preselection
                // SIDIS mode:
                //   require one trigger electron and at least one selected hadron.
                // Inclusive mode:
                //   do not require a hadron.
                //   We only need at least one good electron, which is checked below
                //   after the best electron search.
                // ------------------------------------------------------------

                if(skimMode == "sidis"){
                    bool eh_sidis_event = false;
                    eh_sidis_event_test(particles, selected_hadron_pid, eh_sidis_event);

                    if(!eh_sidis_event) continue;

                    ++sidis_events;
                }

                if(skimMode == "dihadron"){
                    bool dihadron_event = false;
                    epippim_dihadron_event_test(particles, dihadron_event);

                    if(!dihadron_event) continue;

                    ++dihadron_events;
                }
     
                //particles.show();
                //calorimeter.show();
                vectorMap Calomap = loadMapIndex(calorimeter,7);
                vectorMap Htccmap = loadMapIndex(cherenkov,15);
       

                int nrows = particles.getRows();

                // ---------- First pass over REC::Particle ----------
                // SIDIS mode:
                //   Original behavior is preserved.
                //   We choose the highest momentum electron with pid = 11 and det_status = -2.
                // Inclusive mode:
                //   For FCup comparison, use REC::Particle row 0 only.
                //   This matches FCup_reading_summer2026_V1.cc, which reads:
                //      PARTICLES.getInt("pid", 0)
                //      PARTICLES.getInt("status", 0)
                //   and uses pindex = 0 for detector matching.
                // To switch inclusive mode back to the old behavior, set:
                //   USE_ROW0_ELECTRON_FOR_INCLUSIVE = false near the top of this file.

                ParticleCandidate bestElectron;
                DetectorInfo bestElectronDet;

                // Count good electrons and selected hadrons in this event.
                // SIDIS mode:
                //   n_good_h counts the selected single hadron PID.
                // Dihadron mode:
                //   n_good_pip counts pi+
                //   n_good_pim counts pi-
                // These are status level counts before optional detector/PID cuts.
                int n_good_e = 0;
                int n_good_h = 0;

                int n_good_pip = 0;
                int n_good_pim = 0;

                for(int erow = 0; erow < nrows; erow++){

                    int pid_now = particles.getInt("pid", erow);
                    short raw_status_now = particles.getShort("status", erow);
                    int det_status_now = raw_status_now / 1000;

                    // Count selected hadrons for this skim channel.
                    // Example:
                    //   selected_hadron_pid = 211  counts pi+
                    //   selected_hadron_pid = -211 counts pi-
                    //   selected_hadron_pid = 321  counts K+
                    if(skimMode == "sidis"){
                        if(pid_now == selected_hadron_pid && det_status_now == 2){
                            n_good_h++;
                        }
                    }

                    // Count pion candidates for dihadron mode.
                    // We keep this separate from n_good_h because dihadron mode
                    // always needs both pi+ and pi- in the same event.
                    if(skimMode == "dihadron"){
                        if(pid_now == 211 && det_status_now == 2){
                            n_good_pip++;
                        }
                        if(pid_now == -211 && det_status_now == 2){
                            n_good_pim++;
                        }
                    }

                    // From here onward, only electron candidates are considered.
                    if(pid_now != 11) continue;

                    // Original status cut:
                    // if(det_status_now != -2) continue;
                    // For SIDIS, keep your original detector status convention.
                    // For inclusive FCup comparison, use the exact raw status window used there:
                    //      -4000 < status <= -2000
                    bool pass_electron_status = false;

                    if(skimMode == "inclusive" && USE_ROW0_ELECTRON_FOR_INCLUSIVE){
                        pass_electron_status = (raw_status_now > -4000 && raw_status_now <= -2000);
                    }else{
                        pass_electron_status = (det_status_now == -2);
                    }

                    if(!pass_electron_status) continue;

                    n_good_e++;

                    float px_e = particles.getFloat("px", erow);
                    float py_e = particles.getFloat("py", erow);
                    float pz_e = particles.getFloat("pz", erow);

                    float p_e = sqrt(px_e*px_e + py_e*py_e + pz_e*pz_e);

                    bool choose_this_electron = false;

                    if(skimMode == "inclusive" && USE_ROW0_ELECTRON_FOR_INCLUSIVE){

                        // FCup comparison convention:
                        // use only REC::Particle row 0.
                        // Do not replace row 0 with a later higher momentum electron.
                        choose_this_electron = (erow == 0);

                    }else{

                        // Original behavior:
                        // choose the highest momentum trigger electron.
                        // Original line was:
                        // if(!bestElectron.found || p_e > bestElectron.p){
                        choose_this_electron = (!bestElectron.found || p_e > bestElectron.p);
                    }

                    if(choose_this_electron){

                        TVector3 eMom(px_e, py_e, pz_e);

                        bestElectron.found = true;

                        bestElectron.pindex = erow;
                        bestElectron.pid = pid_now;
                        bestElectron.status = raw_status_now;
                        bestElectron.det_status = det_status_now;

                        bestElectron.px = px_e;
                        bestElectron.py = py_e;
                        bestElectron.pz = pz_e;

                        bestElectron.vx = particles.getFloat("vx", erow);
                        bestElectron.vy = particles.getFloat("vy", erow);
                        bestElectron.vz = particles.getFloat("vz", erow);
                        bestElectron.vt = particles.getFloat("vt", erow);

                        bestElectron.p = p_e;
                        bestElectron.theta = eMom.Theta();
                        bestElectron.phi = eMom.Phi();

                        bestElectronDet = getDetectorInfo(
                            erow,
                            particles,
                            calorimeter,
                            cherenkov,
                            Calomap,
                            Htccmap,
                            p_e
                        );
                    }
                }

                // This should normally not happen because eh_sidis_event_test already checked it.
                // Keep it anyway as protection.
                if(!bestElectron.found) continue;

                // ============================================================
                // Inclusive mode fill
                // One row = one selected best electron event.
                // We do this before the SIDIS electron candidate tree and hadron loop.
                // Then we continue to the next HIPO event.
                // ============================================================

                if(skimMode == "inclusive"){

                    InclusiveKin incKin = calculateInclusiveKinematics(bestElectron);

                    // Optional detector/PID cuts.
                    // For inclusive DIS and QE, use the dedicated detector cuts.
                    // For all and res, keep the original loose electron detector/PID cuts.
                    if(applyDetPidCut == 1){
                        bool good_e = true;

                        if(inclusiveRegion == "dis"){

                            DCEdgeInfo bestElectronDC = getDCEdgeInfo(
                                bestElectron.pindex,
                                trajectory
                            );

                            good_e = passInclusiveDISDetectorCuts(
                                bestElectron,
                                bestElectronDet,
                                bestElectronDC
                            );

                        }else if(inclusiveRegion == "qe"){

                            good_e = passInclusiveQEDetectorCuts(
                                bestElectron,
                                bestElectronDet
                            );

                        }else{

                            good_e = passElectronDetPidCuts(bestElectronDet);
                        }

                        if(!good_e){
                            continue;
                        }
                    }

                    // Optional dedicated kinematic cuts.
                    // These are applied only for regions that have a dedicated kinematic cut set.
                    if(applyKinCut == 1){
                        bool good_kin = true;

                        if(inclusiveRegion == "dis"){

                            good_kin = passInclusiveDISKinematicCuts(
                                bestElectron,
                                incKin
                            );

                        }else if(inclusiveRegion == "qe"){

                            good_kin = passInclusiveQEKinematicCuts(
                                bestElectron,
                                incKin
                            );
                        }

                        if(!good_kin){
                            continue;
                        }
                    }

                    inclusive_events_seen++;

                    if(incKin.is_dis == 1){
                        inclusive_dis_events++;
                    }

                    if(incKin.is_resonance == 1){
                        inclusive_res_events++;
                    }

                    if(incKin.is_qe == 1){
                        inclusive_qe_events++;
                    }

                    // Apply inclusive region filter.
                    //   all = write all accepted inclusive electron events
                    //   dis = write only W > 2
                    //   res = write only 0 < W < 2
                    if(!passInclusiveRegion(incKin, inclusiveRegion)){
                        inclusive_region_rejected++;
                        counter++;
                        total_events++;
                        continue;
                    }

                    // Shared event level information
                    inc_run = evInfo.run;
                    inc_event = evInfo.event;
                    inc_event_file_index = evInfo.event_file_index;

                    inc_input_file_index = evInfo.input_file_index;
                    inc_input_file_number = evInfo.input_file_number;

                    inc_run_config = evInfo.run;
                    inc_event_config = evInfo.event;

                    inc_trigger = evInfo.trigger;
                    inc_timestamp = evInfo.timestamp;
                    inc_unixtime = evInfo.unixtime;

                    inc_helicity = evInfo.helicity;
                    inc_helicity_valid = evInfo.helicity_valid;

                    inc_target_state = evInfo.target_state;

                    inc_target_pol = evInfo.target_pol;
                    inc_target_pol_err = evInfo.target_pol_err;

                    inc_tensor_pol = evInfo.tensor_pol;
                    inc_tensor_pol_err = evInfo.tensor_pol_err;

                    inc_vector_pol = evInfo.vector_pol;
                    inc_vector_pol_err = evInfo.vector_pol_err;

                    copyTargetInfoToBranches(evInfo, inc_target_extra);

                    // Electron information
                    inc_n_good_e = n_good_e;

                    inc_e_pindex = bestElectron.pindex;
                    inc_pid_e = bestElectron.pid;
                    inc_e_status = bestElectron.status;
                    inc_e_det_status = bestElectron.det_status;

                    inc_pxe = bestElectron.px;
                    inc_pye = bestElectron.py;
                    inc_pze = bestElectron.pz;

                    inc_vxe = bestElectron.vx;
                    inc_vye = bestElectron.vy;
                    inc_vze = bestElectron.vz;
                    inc_vte = bestElectron.vt;

                    inc_pe = bestElectron.p;
                    inc_thetae = bestElectron.theta;
                    inc_phie = bestElectron.phi;

                    // Electron detector/PID information
                    inc_e_charge = bestElectronDet.charge;

                    inc_e_beta = bestElectronDet.beta;
                    inc_e_chi2pid = bestElectronDet.chi2pid;

                    inc_e_has_htcc = bestElectronDet.has_htcc;
                    inc_e_htcc_nphe = bestElectronDet.htcc_nphe;

                    inc_e_has_cal = bestElectronDet.has_cal;
                    inc_e_cal_sector = bestElectronDet.cal_sector;

                    inc_e_pcal_e = bestElectronDet.pcal_e;
                    inc_e_ecin_e = bestElectronDet.ecin_e;
                    inc_e_ecout_e = bestElectronDet.ecout_e;
                    inc_e_cal_e = bestElectronDet.cal_e;
                    inc_e_sampling_fraction = bestElectronDet.sampling_fraction;

                    inc_e_pcal_x = bestElectronDet.pcal_x;
                    inc_e_pcal_y = bestElectronDet.pcal_y;

                    inc_e_lu = bestElectronDet.lu;
                    inc_e_lv = bestElectronDet.lv;
                    inc_e_lw = bestElectronDet.lw;

                    // Inclusive kinematics
                    inc_Q2 = incKin.Q2;
                    inc_nu = incKin.nu;

                    // nucleon level W
                    inc_W2 = incKin.W2;
                    inc_W = incKin.W;

                    // deuteron system W
                    inc_Wd2 = incKin.Wd2;
                    inc_Wd = incKin.Wd;

                    inc_xB = incKin.xB;
                    inc_xD = incKin.xD;
                    inc_y = incKin.y;

                    inc_is_dis = incKin.is_dis;
                    inc_is_resonance = incKin.is_resonance;
                    inc_is_qe = incKin.is_qe;

                    if(writeRootOutput){
                        inclusiveTree->Fill();
                        runSummary.n_inclusive_written++;
                    }

                    if(writeHipoOutput && !hipo_event_written_this_event){
                        hipoWriter.addEvent(event);
                        hipo_events_written++;
                        hipo_events_written_this_file++;
                        hipo_event_written_this_event = true;
                    }

                    inclusive_events++;

                    counter++;
                    total_events++;

                    continue;
                }

                // ============================================================
                // Dihadron mode fill
                // ============================================================
                // Used only when:
                //   --mode dihadron
                // One output row = best electron + one pi+ pi- pair.
                // Important:
                //   This block is placed before the old SIDIS electronTree
                //   and single hadron loop. Therefore, after filling the
                //   dihadron tree, we continue to the next event.
                // ============================================================

                if(skimMode == "dihadron"){

                    vector<ParticleCandidate> pipCandidates;
                    vector<ParticleCandidate> pimCandidates;

                    vector<DetectorInfo> pipDets;
                    vector<DetectorInfo> pimDets;

                    // ------------------------------------------------------------
                    // Optional electron detector/PID check.
                    // If --detpidcut 1 is requested and the best electron fails
                    // the loose electron QA cuts, then this event cannot produce
                    // a clean dihadron row.
                    // ------------------------------------------------------------

                    if(applyDetPidCut == 1){
                        bool good_e_det_pid = passElectronDetPidCuts(bestElectronDet);

                        if(!good_e_det_pid){
                            counter++;
                            total_events++;
                            continue;
                        }
                    }

                    // ------------------------------------------------------------
                    // Collect pi+ and pi- candidates.
                    // We do not cut on Q2, W, z, Mh, or missing mass here.
                    // Those cuts should be studied later in the analysis stage.
                    // ------------------------------------------------------------

                    for(int row = 0; row < nrows; row++){

                        int pid = particles.getInt("pid", row);
                        short raw_status = particles.getShort("status", row);
                        int det_status = raw_status / 1000;

                        // Dihadron channel currently means pi+ pi- only.
                        if(pid != 211 && pid != -211) continue;

                        // For now use forward detector positive status hadrons, matching the existing SIDIS convention.
                        if(det_status != 2) continue;

                        float px = particles.getFloat("px", row);
                        float py = particles.getFloat("py", row);
                        float pz = particles.getFloat("pz", row);

                        float p_mag = sqrt(px*px + py*py + pz*pz);
                        TVector3 mom(px, py, pz);

                        ParticleCandidate pion;

                        pion.found = true;
                        pion.pindex = row;
                        pion.pid = pid;
                        pion.status = raw_status;
                        pion.det_status = det_status;

                        pion.px = px;
                        pion.py = py;
                        pion.pz = pz;

                        pion.vx = particles.getFloat("vx", row);
                        pion.vy = particles.getFloat("vy", row);
                        pion.vz = particles.getFloat("vz", row);
                        pion.vt = particles.getFloat("vt", row);

                        pion.p = p_mag;
                        pion.theta = mom.Theta();
                        pion.phi = mom.Phi();

                        DetectorInfo pionDet = getDetectorInfo(
                            row,
                            particles,
                            calorimeter,
                            cherenkov,
                            Calomap,
                            Htccmap,
                            p_mag
                        );

                        // Optional loose detector/PID sanity cuts.
                        // Reuse the existing single hadron PID helper.
                        if(applyDetPidCut == 1){
                            bool good_pion_det_pid = passDetPidCuts(
                                bestElectronDet,
                                pionDet,
                                pid
                            );

                            if(!good_pion_det_pid){
                                continue;
                            }
                        }

                        if(pid == 211){
                            pipCandidates.push_back(pion);
                            pipDets.push_back(pionDet);
                        }

                        if(pid == -211){
                            pimCandidates.push_back(pion);
                            pimDets.push_back(pionDet);
                        }
                    }

                    // Number of possible pi+pi- pairs after optional detector/PID cuts.
                    int n_dihadron_pairs_this_event =
                        static_cast<int>(pipCandidates.size() * pimCandidates.size());

                    // If no pair survives, move to the next event.
                    if(n_dihadron_pairs_this_event <= 0){
                        counter++;
                        total_events++;
                        continue;
                    }

                    // ------------------------------------------------------------
                    // Electron DIS kinematics.
                    // Use the same helper already used by inclusive mode so that
                    // Q2, W, Wd, xB, xD, and y are consistent across skim modes.
                    // ------------------------------------------------------------

                    InclusiveKin diKin = calculateInclusiveKinematics(bestElectron);

                    TLorentzVector e_in4vect_di(0.0, 0.0, beam_enrg, beam_enrg);

                    TLorentzVector e_out4vect_di(
                        bestElectron.px,
                        bestElectron.py,
                        bestElectron.pz,
                        sqrt(bestElectron.p * bestElectron.p + M_electron * M_electron)
                    );

                    TLorentzVector q_4vect_di = e_in4vect_di - e_out4vect_di;
                    TVector3 q_vec_di = q_4vect_di.Vect();

                    // ------------------------------------------------------------
                    // Build all pi+ pi- combinations.
                    // If an event has:
                    //   2 pi+ and 1 pi-
                    // it produces:
                    //   2 rows in the dihadron tree.
                    // ------------------------------------------------------------

                    int pair_id_di = 0;

                    for(size_t ipip = 0; ipip < pipCandidates.size(); ipip++){

                        for(size_t ipim = 0; ipim < pimCandidates.size(); ipim++){

                            const ParticleCandidate &pip = pipCandidates[ipip];
                            const ParticleCandidate &pim = pimCandidates[ipim];

                            const DetectorInfo &pipDet = pipDets[ipip];
                            const DetectorInfo &pimDet = pimDets[ipim];

                            TLorentzVector pip4(
                                pip.px,
                                pip.py,
                                pip.pz,
                                sqrt(pip.p * pip.p + M_pion * M_pion)
                            );

                            TLorentzVector pim4(
                                pim.px,
                                pim.py,
                                pim.pz,
                                sqrt(pim.p * pim.p + M_pion * M_pion)
                            );

                            TLorentzVector pair4 = pip4 + pim4;
                            TVector3 pair3 = pair4.Vect();

                            // ---------- Shared event level information ----------
                            di_run = evInfo.run;
                            di_event = evInfo.event;
                            di_event_file_index = evInfo.event_file_index;

                            di_input_file_index = evInfo.input_file_index;
                            di_input_file_number = evInfo.input_file_number;

                            di_run_config = evInfo.run;
                            di_event_config = evInfo.event;

                            di_trigger = evInfo.trigger;
                            di_timestamp = evInfo.timestamp;
                            di_unixtime = evInfo.unixtime;

                            di_helicity = evInfo.helicity;
                            di_helicity_valid = evInfo.helicity_valid;

                            // ---------- Target/tensor information ----------
                            di_target_state = evInfo.target_state;

                            di_target_pol = evInfo.target_pol;
                            di_target_pol_err = evInfo.target_pol_err;

                            di_tensor_pol = evInfo.tensor_pol;
                            di_tensor_pol_err = evInfo.tensor_pol_err;

                            di_vector_pol = evInfo.vector_pol;
                            di_vector_pol_err = evInfo.vector_pol_err;

                            copyTargetInfoToBranches(evInfo, di_target_extra);

                            // ---------- Pair/multiplicity information ----------
                            di_pair_id = pair_id_di;

                            di_n_good_e = n_good_e;
                            di_n_good_pip = n_good_pip;
                            di_n_good_pim = n_good_pim;
                            di_n_dihadron_pairs = n_dihadron_pairs_this_event;

                            // ---------- Particle indices and status ----------
                            di_e_pindex = bestElectron.pindex;
                            di_pip_pindex = pip.pindex;
                            di_pim_pindex = pim.pindex;

                            di_pid_e = 11;
                            di_pid_pip = 211;
                            di_pid_pim = -211;

                            di_e_status = bestElectron.status;
                            di_pip_status = pip.status;
                            di_pim_status = pim.status;

                            di_e_det_status = bestElectron.det_status;
                            di_pip_det_status = pip.det_status;
                            di_pim_det_status = pim.det_status;

                            // ---------- Electron kinematics ----------
                            di_pxe = bestElectron.px;
                            di_pye = bestElectron.py;
                            di_pze = bestElectron.pz;

                            di_vxe = bestElectron.vx;
                            di_vye = bestElectron.vy;
                            di_vze = bestElectron.vz;
                            di_vte = bestElectron.vt;

                            di_pe = bestElectron.p;
                            di_thetae = bestElectron.theta;
                            di_phie = bestElectron.phi;

                            // ---------- Electron detector/PID ----------
                            di_e_charge = bestElectronDet.charge;

                            di_e_beta = bestElectronDet.beta;
                            di_e_chi2pid = bestElectronDet.chi2pid;

                            di_e_has_htcc = bestElectronDet.has_htcc;
                            di_e_htcc_nphe = bestElectronDet.htcc_nphe;

                            di_e_has_cal = bestElectronDet.has_cal;
                            di_e_cal_sector = bestElectronDet.cal_sector;

                            di_e_pcal_e = bestElectronDet.pcal_e;
                            di_e_ecin_e = bestElectronDet.ecin_e;
                            di_e_ecout_e = bestElectronDet.ecout_e;
                            di_e_cal_e = bestElectronDet.cal_e;
                            di_e_sampling_fraction = bestElectronDet.sampling_fraction;

                            di_e_pcal_x = bestElectronDet.pcal_x;
                            di_e_pcal_y = bestElectronDet.pcal_y;

                            di_e_lu = bestElectronDet.lu;
                            di_e_lv = bestElectronDet.lv;
                            di_e_lw = bestElectronDet.lw;

                            // ---------- pi+ kinematics ----------
                            di_px_pip = pip.px;
                            di_py_pip = pip.py;
                            di_pz_pip = pip.pz;

                            di_vx_pip = pip.vx;
                            di_vy_pip = pip.vy;
                            di_vz_pip = pip.vz;
                            di_vt_pip = pip.vt;

                            di_p_pip = pip.p;
                            di_theta_pip = pip.theta;
                            di_phi_pip = pip.phi;

                            // ---------- pi+ detector/PID ----------
                            di_pip_charge = pipDet.charge;

                            di_pip_beta = pipDet.beta;
                            di_pip_chi2pid = pipDet.chi2pid;

                            di_pip_has_htcc = pipDet.has_htcc;
                            di_pip_htcc_nphe = pipDet.htcc_nphe;

                            di_pip_has_cal = pipDet.has_cal;
                            di_pip_cal_sector = pipDet.cal_sector;

                            di_pip_pcal_e = pipDet.pcal_e;
                            di_pip_ecin_e = pipDet.ecin_e;
                            di_pip_ecout_e = pipDet.ecout_e;
                            di_pip_cal_e = pipDet.cal_e;
                            di_pip_sampling_fraction = pipDet.sampling_fraction;

                            di_pip_pcal_x = pipDet.pcal_x;
                            di_pip_pcal_y = pipDet.pcal_y;

                            di_pip_lu = pipDet.lu;
                            di_pip_lv = pipDet.lv;
                            di_pip_lw = pipDet.lw;

                            // ---------- pi- kinematics ----------
                            di_px_pim = pim.px;
                            di_py_pim = pim.py;
                            di_pz_pim = pim.pz;

                            di_vx_pim = pim.vx;
                            di_vy_pim = pim.vy;
                            di_vz_pim = pim.vz;
                            di_vt_pim = pim.vt;

                            di_p_pim = pim.p;
                            di_theta_pim = pim.theta;
                            di_phi_pim = pim.phi;

                            // ---------- pi- detector/PID ----------
                            di_pim_charge = pimDet.charge;

                            di_pim_beta = pimDet.beta;
                            di_pim_chi2pid = pimDet.chi2pid;

                            di_pim_has_htcc = pimDet.has_htcc;
                            di_pim_htcc_nphe = pimDet.htcc_nphe;

                            di_pim_has_cal = pimDet.has_cal;
                            di_pim_cal_sector = pimDet.cal_sector;

                            di_pim_pcal_e = pimDet.pcal_e;
                            di_pim_ecin_e = pimDet.ecin_e;
                            di_pim_ecout_e = pimDet.ecout_e;
                            di_pim_cal_e = pimDet.cal_e;
                            di_pim_sampling_fraction = pimDet.sampling_fraction;

                            di_pim_pcal_x = pimDet.pcal_x;
                            di_pim_pcal_y = pimDet.pcal_y;

                            di_pim_lu = pimDet.lu;
                            di_pim_lv = pimDet.lv;
                            di_pim_lw = pimDet.lw;

                            // ---------- DIS kinematics ----------
                            di_Q2 = diKin.Q2;
                            di_nu = diKin.nu;

                            di_W2 = diKin.W2;
                            di_W = diKin.W;

                            di_Wd2 = diKin.Wd2;
                            di_Wd = diKin.Wd;

                            di_xB = diKin.xB;
                            di_xD = diKin.xD;
                            di_y = diKin.y;

                            di_is_dis = diKin.is_dis;
                            di_is_resonance = diKin.is_resonance;

                            // ---------- Dihadron pair variables ----------
                            di_pair_px = pair4.Px();
                            di_pair_py = pair4.Py();
                            di_pair_pz = pair4.Pz();
                            di_pair_E = pair4.E();
                            di_pair_p = pair3.Mag();

                            di_pair_theta = pair3.Theta();
                            di_pair_phi = pair3.Phi();

                            di_Mh = pair4.M();

                            if(fabs(di_nu) > 1.0e-12){
                                di_z_pip = pip4.E() / di_nu;
                                di_z_pim = pim4.E() / di_nu;
                                di_z_pair = di_z_pip + di_z_pim;
                            }else{
                                di_z_pip = -9999;
                                di_z_pim = -9999;
                                di_z_pair = -9999;
                            }

                            // Pair transverse momentum relative to q.
                            if(q_vec_di.Mag() > 1.0e-12){
                                TVector3 qhat = q_vec_di.Unit();
                                TVector3 pairT = pair3 - pair3.Dot(qhat) * qhat;
                                di_PhT = pairT.Mag();
                            }else{
                                di_PhT = -9999;
                            }

                            // Leave these as placeholders for this first skim.
                            // We will define phi_h and phi_R carefully in the
                            // analysis stage, using the Trento convention.
                            di_phi_h = -9999;
                            di_phi_R = -9999;

                            if(writeRootOutput){
                                dihadronTree->Fill();
                                runSummary.n_dihadron_rows_written++;                
                            }

                            dihadron_pairs_written++;
                            pair_id_di++;
                        }
                    }

                    if(writeHipoOutput && !hipo_event_written_this_event){
                        hipoWriter.addEvent(event);
                        hipo_events_written++;
                        hipo_events_written_this_file++;
                        hipo_event_written_this_event = true;
                    }

                    counter++;
                    total_events++;

                    continue;
                }

                // ---------- Store all good electron candidates ----------
                // The sidis tree uses only bestElectron.
                // This electronTree stores every good trigger electron candidate in the event.
                // A good electron here means: pid = 11  det_status = -2
                // is_best_e = 1 marks the electron used in the sidis tree.
                // is_best_e = 0 marks the other electron candidates.

                int e_candidate_id = 0;

                for(int erow = 0; erow < nrows; erow++){

                    int pid_now = particles.getInt("pid", erow);
                    short raw_status_now = particles.getShort("status", erow);
                    int det_status_now = raw_status_now / 1000;

                    if(pid_now != 11) continue;
                    if(det_status_now != -2) continue;

                    float px_e = particles.getFloat("px", erow);
                    float py_e = particles.getFloat("py", erow);
                    float pz_e = particles.getFloat("pz", erow);

                    float p_e = sqrt(px_e*px_e + py_e*py_e + pz_e*pz_e);
                    TVector3 eMom(px_e, py_e, pz_e);

                    // Shared event level information.
                    // run/event come from RUN::config if available.
                    // event_file_index preserves the old local counter.
                    elec_run = evInfo.run;
                    elec_event = evInfo.event;
                    elec_event_file_index = evInfo.event_file_index;

                    elec_input_file_index = evInfo.input_file_index;
                    elec_input_file_number = evInfo.input_file_number;

                    elec_run_config = evInfo.run;
                    elec_event_config = evInfo.event;

                    elec_trigger = evInfo.trigger;
                    elec_timestamp = evInfo.timestamp;
                    elec_unixtime = evInfo.unixtime;

                    elec_helicity = evInfo.helicity;
                    elec_helicity_valid = evInfo.helicity_valid;

                    // Target/tensor placeholders.
                    // These will be filled later from a target state map.
                    elec_target_state = evInfo.target_state;

                    elec_target_pol = evInfo.target_pol;
                    elec_target_pol_err = evInfo.target_pol_err;

                    elec_tensor_pol = evInfo.tensor_pol;
                    elec_tensor_pol_err = evInfo.tensor_pol_err;

                    elec_vector_pol = evInfo.vector_pol;
                    elec_vector_pol_err = evInfo.vector_pol_err;

                    copyTargetInfoToBranches(evInfo, elec_target_extra);

                    elec_candidate_id = e_candidate_id;
                    elec_is_best = (erow == bestElectron.pindex) ? 1 : 0;

                    elec_n_good_e = n_good_e;
                    elec_n_good_h = n_good_h;

                    elec_pindex = erow;
                    elec_pid = pid_now;
                    elec_status = raw_status_now;
                    elec_det_status = det_status_now;

                    elec_px = px_e;
                    elec_py = py_e;
                    elec_pz = pz_e;

                    elec_vx = particles.getFloat("vx", erow);
                    elec_vy = particles.getFloat("vy", erow);
                    elec_vz = particles.getFloat("vz", erow);
                    elec_vt = particles.getFloat("vt", erow);

                    elec_p = p_e;
                    elec_theta = eMom.Theta();
                    elec_phi = eMom.Phi();

                    DetectorInfo elecDet = getDetectorInfo(
                        erow,
                        particles,
                        calorimeter,
                        cherenkov,
                        Calomap,
                        Htccmap,
                        p_e
                    );

                    elec_charge = elecDet.charge;

                    elec_beta = elecDet.beta;
                    elec_chi2pid = elecDet.chi2pid;

                    elec_has_htcc = elecDet.has_htcc;
                    elec_htcc_nphe = elecDet.htcc_nphe;

                    elec_has_cal = elecDet.has_cal;
                    elec_cal_sector = elecDet.cal_sector;

                    elec_pcal_e = elecDet.pcal_e;
                    elec_ecin_e = elecDet.ecin_e;
                    elec_ecout_e = elecDet.ecout_e;
                    elec_cal_e = elecDet.cal_e;
                    elec_sampling_fraction = elecDet.sampling_fraction;

                    elec_pcal_x = elecDet.pcal_x;
                    elec_pcal_y = elecDet.pcal_y;

                    elec_lu = elecDet.lu;
                    elec_lv = elecDet.lv;
                    elec_lw = elecDet.lw;

                    if(writeRootOutput && writeElectronTree == 1){
                        electronTree->Fill();
                    }

                    e_candidate_id++;
                }

                // pair_id counts how many selected hadrons were paired with the best electron inside this accepted event.
                int pair_id = 0;

                //printf("---------- PARTICLE BANK CONTENT and Calorimeter energy -------\n");
                for(int row = 0; row < nrows; row++){
                    int   pid = particles.getInt("pid",row);
                    float  px = particles.getFloat("px",row);
                    float  py = particles.getFloat("py",row);
                    float  pz = particles.getFloat("pz",row);
                    float  vx = particles.getFloat("vx",row);
                    float  vy = particles.getFloat("vy",row);
                    float  vz = particles.getFloat("vz",row);
                    float  vt = particles.getFloat("vt",row);
                    short status = particles.getShort("status",row);
                    //cout<<pid<<"\t"<< setw(8)<<px<<"\t"<< setw(8)<<py<<"\t"<< setw(8)<<pz<<"\t"<< setw(8)<<"\n";

                    // Keep the raw status before the existing histogram code divides it by 1000.
                    short raw_status = status;
                    int det_status = raw_status / 1000;

                    float p_mag = sqrt(px*px + py*py + pz*pz);
                    TVector3 hMom(px, py, pz);

                    // ---------- Fill one row for each selected hadron ----------
                    // This is the main skim logic.
                    // Example: selected_hadron_pid = 211
                    // If this event has one trigger electron and three pi+ candidates, then skimTree->Fill() will be called three times.
                    // This avoids variable names like pxpip, pxpip1, pxpip2.
                    // Instead, each pi+ gets its own row with the same electron repeated.

                    if(pid == selected_hadron_pid && det_status == 2){

                        // Shared event level information.
                        // run/event come from RUN::config if available.
                        // event_file_index preserves the old local counter.
                        out_run = evInfo.run;
                        out_event = evInfo.event;
                        out_event_file_index = evInfo.event_file_index;

                        out_input_file_index = evInfo.input_file_index;
                        out_input_file_number = evInfo.input_file_number;

                        out_run_config = evInfo.run;
                        out_event_config = evInfo.event;

                        out_trigger = evInfo.trigger;
                        out_timestamp = evInfo.timestamp;
                        out_unixtime = evInfo.unixtime;

                        out_helicity = evInfo.helicity;
                        out_helicity_valid = evInfo.helicity_valid;

                        // Target/tensor placeholders.
                        // These will be filled later from a target state map.
                        out_target_state = evInfo.target_state;

                        out_target_pol = evInfo.target_pol;
                        out_target_pol_err = evInfo.target_pol_err;

                        out_tensor_pol = evInfo.tensor_pol;
                        out_tensor_pol_err = evInfo.tensor_pol_err;

                        out_vector_pol = evInfo.vector_pol;
                        out_vector_pol_err = evInfo.vector_pol_err;

                        copyTargetInfoToBranches(evInfo, out_target_extra);

                        out_pair_id = pair_id;

                        out_n_good_e = n_good_e;
                        out_n_good_h = n_good_h;

                        out_e_pindex = bestElectron.pindex;
                        out_h_pindex = row;

                        out_pid_e = 11;
                        out_pid_h = selected_hadron_pid;

                        out_e_status = bestElectron.status;
                        out_h_status = raw_status;

                       out_e_det_status = bestElectron.det_status;
                        out_h_det_status = det_status;

                        DetectorInfo hDet = getDetectorInfo(
                            row,
                            particles,
                            calorimeter,
                            cherenkov,
                            Calomap,
                            Htccmap,
                            p_mag
                        );

                        if(applyDetPidCut == 1){
                            bool good_det_pid = passDetPidCuts(bestElectronDet, hDet, selected_hadron_pid);
                            if(!good_det_pid){
                                continue;
                            }
                        }

                        // Electron detector/PID information
                        e_charge = bestElectronDet.charge;

                        e_beta = bestElectronDet.beta;
                        e_chi2pid = bestElectronDet.chi2pid;

                        e_has_htcc = bestElectronDet.has_htcc;
                        e_htcc_nphe = bestElectronDet.htcc_nphe;

                        e_has_cal = bestElectronDet.has_cal;
                        e_cal_sector = bestElectronDet.cal_sector;

                        e_pcal_e = bestElectronDet.pcal_e;
                        e_ecin_e = bestElectronDet.ecin_e;
                        e_ecout_e = bestElectronDet.ecout_e;
                        e_cal_e = bestElectronDet.cal_e;
                        e_sampling_fraction = bestElectronDet.sampling_fraction;

                        e_pcal_x = bestElectronDet.pcal_x;
                        e_pcal_y = bestElectronDet.pcal_y;

                        e_lu = bestElectronDet.lu;
                        e_lv = bestElectronDet.lv;
                        e_lw = bestElectronDet.lw;

                        // Hadron detector/PID information
                        h_charge = hDet.charge;

                        h_beta = hDet.beta;
                        h_chi2pid = hDet.chi2pid;

                        h_has_htcc = hDet.has_htcc;
                        h_htcc_nphe = hDet.htcc_nphe;

                        h_has_cal = hDet.has_cal;
                        h_cal_sector = hDet.cal_sector;

                        h_pcal_e = hDet.pcal_e;
                        h_ecin_e = hDet.ecin_e;
                        h_ecout_e = hDet.ecout_e;
                        h_cal_e = hDet.cal_e;
                        h_sampling_fraction = hDet.sampling_fraction;

                        h_pcal_x = hDet.pcal_x;
                        h_pcal_y = hDet.pcal_y;

                        h_lu = hDet.lu;
                        h_lv = hDet.lv;
                        h_lw = hDet.lw;

                        // Electron information
                        pxe = bestElectron.px;
                        pye = bestElectron.py;
                        pze = bestElectron.pz;

                        vxe = bestElectron.vx;
                        vye = bestElectron.vy;
                        vze = bestElectron.vz;
                        vte = bestElectron.vt;

                        pe = bestElectron.p;
                        thetae = bestElectron.theta;
                        phie = bestElectron.phi;

                        // Hadron information
                        pxh = px;
                        pyh = py;
                        pzh = pz;

                        vxh = vx;
                        vyh = vy;
                        vzh = vz;
                        vth = vt;

                        ph = p_mag;
                        thetah = hMom.Theta();
                        phih = hMom.Phi();

                        if(writeRootOutput){
                            skimTree->Fill();
                            runSummary.n_sidis_rows_written++; 
                        }

                        pair_id++;
                    }
                  
                } //end of particle bank for loop
       
                //printf("---------- END OF PARTICLE BANK and Calorimeter energy -------\n");

                if(writeHipoOutput && pair_id > 0 && !hipo_event_written_this_event){
                    hipoWriter.addEvent(event);
                    hipo_events_written++;
                    hipo_events_written_this_file++;
                    hipo_event_written_this_event = true;
                }

                counter++;
                total_events++;

            } //end of while loop over events
            printf("processed events = %d\n",counter);

            // ---------- Write the skim tree for this input HIPO file ----------
            // This produces one ROOT file per HIPO file.
            // The tree is named "sidis".
            // Another tree name is "electrons".

            if(writeRootOutput){
                fout->cd();

                // Write scaler information first.
                // These trees exist for every mode and make FCup normalized checks possible downstream.
                scalerTree->Write();

                runSummaryTree->Fill();
                runSummaryTree->Write();

                if(skimMode == "sidis"){
                    skimTree->Write();
                    if(writeElectronTree == 1){
                        electronTree->Write();
                    }
                }else if(skimMode == "inclusive"){
                    inclusiveTree->Write();
                }else if(skimMode == "dihadron"){
                    dihadronTree->Write();
                }

                fout->Close();

                cout << "Wrote ROOT skim file: " << outRootName << endl;
            }

            if(writeHipoOutput){
                hipoWriter.close();
                cout << "Wrote HIPO skim file: " << outHipoName << endl;
                cout << "HIPO events written in this file = " << hipo_events_written_this_file << endl;
            }

            // The ROOT file will contain two TTrees:
            // sidis tree
            //     one row = best electron + one pi+ candidate
            // electrons tree
            //     one row = one good electron candidate in an accepted e pi+ event
            // Example: event 100 has 2 good electrons and 3 good pi+
            // sidis tree gets 3 rows:
            // best electron + pi+ 0
            // best electron + pi+ 1
            // best electron + pi+ 2
            // electrons tree gets 2 rows:
            // electron candidate 0
            // electron candidate 1
            // The sidis rows will have n_good_e = 2 and n_good_h = 3.
            // The electron tree will also have n_good_e = 2 and n_good_h = 3
            // and one of them will have is_best_e = 1.

            
        } //end of argv loop
    
    } //end of argc if

        
    printf("Total sidis events = %d\n", sidis_events);

    printf("Total inclusive events seen before region filter = %d\n", inclusive_events_seen);
    printf("Total inclusive DIS candidates W > 2 = %d\n", inclusive_dis_events);
    printf("Total inclusive resonance candidates 0 < W < 2 = %d\n", inclusive_res_events);
    printf("Total inclusive QE candidates broad xB/Q2/y cut = %d\n", inclusive_qe_events);
    printf("Total inclusive events rejected by region filter = %d\n", inclusive_region_rejected);
    printf("Total inclusive events written = %d\n", inclusive_events);

    printf("Total dihadron preselected events = %d\n", dihadron_events);
    printf("Total dihadron pairs written = %lld\n", dihadron_pairs_written);

    printf("Total processed/output events = %d\n", total_events);

    printf("Target map entries = %d\n", target_map_entries);
    printf("Target matched input events = %lld\n", target_matched_input_events);
    printf("Target unmatched input events = %lld\n", target_unmatched_input_events);
    printf("Target species p input events = %lld\n", target_species_p_input_events);
    printf("Target species d input events = %lld\n", target_species_d_input_events);

}
// End of main().

void fillScalerTree(
    hipo::bank bankHelScaler,
    const EventInfo &evInfo,
    TTree *scalerTree,
    ScalerBranchVars &scalerVars,
    RunSummaryVars &runSummary
){
    if(!scalerTree){
        return;
    }

    const int nScalerRows = bankHelScaler.getRows();

    for(int irow = 0; irow < nScalerRows; irow++){

        scalerVars.run = evInfo.run;
        scalerVars.event = evInfo.event;
        scalerVars.event_file_index = evInfo.event_file_index;

        scalerVars.input_file_index = evInfo.input_file_index;
        scalerVars.input_file_number = evInfo.input_file_number;

        // Match FCup_reading_summer2026_V1.cc:
        // use helicity from HEL::scaler, not REC::Event.
        scalerVars.helicity = bankHelScaler.getByte("helicity", irow);

        if(scalerVars.helicity == -1 || scalerVars.helicity == 1){
            scalerVars.helicity_valid = 1;
        }else{
            scalerVars.helicity_valid = 0;
        }

        // Match FCup_reading_summer2026_V1.cc:
        // sum HEL::scaler.fcupgated directly by scaler helicity.
        scalerVars.fcup = bankHelScaler.getFloat("fcup", irow);
        scalerVars.fcupgated = bankHelScaler.getFloat("fcupgated", irow);

        // Keep these names for compatibility with rgc_analysis.cc.
        // For HEL::scaler, these are not differences.
        // They are direct per row HEL::scaler charge values.
        scalerVars.fcup_delta = scalerVars.fcup;
        scalerVars.fcupgated_delta = scalerVars.fcupgated;

        scalerTree->Fill();

        runSummary.n_scaler_rows++;

        if(scalerVars.helicity == 1){
            runSummary.fcup_plus += scalerVars.fcup_delta;
            runSummary.fcupgated_plus += scalerVars.fcupgated_delta;
        }else if(scalerVars.helicity == -1){
            runSummary.fcup_minus += scalerVars.fcup_delta;
            runSummary.fcupgated_minus += scalerVars.fcupgated_delta;
        }else{
            runSummary.fcup_zero += scalerVars.fcup_delta;
            runSummary.fcupgated_zero += scalerVars.fcupgated_delta;
        }
    }
}

// Read event level information from REC::Event and RUN::config.
// Notes:
//   RUN::config is preferred for true run/event identifiers.
//   counter is kept as event_file_index for debugging and backwards checks.
//   helicity is read from REC::Event.helicity.
//   target/tensor polarization values are placeholders for now.
// If a bank has zero rows, fallback values are used.

EventInfo getEventInfo(
    hipo::bank bankEvent,
    hipo::bank bankRunConfig,
    int event_file_index,
    int fallback_run,
    int input_file_index,
    int input_file_number){

    EventInfo info;

    info.event_file_index = event_file_index;
    info.input_file_index = input_file_index;
    info.input_file_number = input_file_number;

    // Fallbacks.
    // These are used if RUN::config is missing or empty.
    info.run = fallback_run;
    info.event = event_file_index;

    // ---------- RUN::config ----------
    // Standard metadata. This should fix the old issue where run was sometimes 0
    // because it was parsed from the filename using a fragile sscanf pattern.

    if(bankRunConfig.getRows() > 0){

        int run_from_config = bankRunConfig.getInt("run", 0);
        int event_from_config = bankRunConfig.getInt("event", 0);

        if(run_from_config > 0){
            info.run = run_from_config;
        }else{
            info.run = fallback_run;
        }

        if(event_from_config > 0){
            info.event = event_from_config;
        }else{
            info.event = event_file_index;
        }

        info.trigger = bankRunConfig.getLong("trigger", 0);
        info.timestamp = bankRunConfig.getLong("timestamp", 0);

        info.unixtime = bankRunConfig.getInt("unixtime", 0);

    }else{

        info.run = fallback_run;
        info.event = event_file_index;
    }

    // ---------- REC::Event ----------
    // Beam helicity.
    // Valid physics helicity states should normally be -1 and +1.
    // Other values are kept but marked invalid.

    if(bankEvent.getRows() > 0){

        info.helicity = bankEvent.getByte("helicity", 0);

        if(info.helicity == -1 || info.helicity == 1){
            info.helicity_valid = 1;
        }else{
            info.helicity_valid = 0;
        }
    }

    return info;
}

string trimString(
    const string &input){

    size_t first = input.find_first_not_of(" \t\r\n");

    if(first == string::npos){
        return "";
    }

    size_t last = input.find_last_not_of(" \t\r\n");

    return input.substr(first, last - first + 1);
}

vector<string> splitCSVLine(
    const string &line){

    vector<string> tokens;
    string token;
    stringstream ss(line);

    while(getline(ss, token, ',')){
        tokens.push_back(trimString(token));
    }

    return tokens;
}

int getTargetSpeciesID(
    const string &species){

    if(species == "p" || species == "P"){
        return 1;
    }

    if(species == "d" || species == "D"){
        return 2;
    }

    return -9999;
}

bool loadTargetMap(
    const string &targetMapFile,
    unordered_map<int, TargetInfo> &targetMap){

    targetMap.clear();

    ifstream fin(targetMapFile.c_str());

    if(!fin.is_open()){
        cout << "Error: cannot open target map file:" << endl;
        cout << "  " << targetMapFile << endl;
        return false;
    }

    string line;
    int n_loaded = 0;
    int n_skipped = 0;

    while(getline(fin, line)){

        line = trimString(line);

        if(line == ""){
            continue;
        }

        // Header line starts with #run in your CSV.
        if(line[0] == '#'){
            continue;
        }

        vector<string> cols = splitCSVLine(line);

        if(cols.size() < 9){
            n_skipped++;
            continue;
        }

        TargetInfo info;

        info.run = atoi(cols[0].c_str());

        info.start_time = cols[1];
        info.stop_time = cols[2];

        info.species = cols[3];
        info.species_id = getTargetSpeciesID(info.species);

        info.cell = atoi(cols[4].c_str());

        info.pol_online = atof(cols[5].c_str());
        info.pol_offline = atof(cols[6].c_str());
        info.pol_offline_err = atof(cols[7].c_str());

        info.run_dose = atof(cols[8].c_str());

        info.target_match = 1;

        if(info.run > 0){
            targetMap[info.run] = info;
            n_loaded++;
        }else{
            n_skipped++;
        }
    }

    fin.close();

    cout << "Target map load summary:" << endl;
    cout << "  loaded  = " << n_loaded << endl;
    cout << "  skipped = " << n_skipped << endl;

    return true;
}

TargetInfo findTargetInfo(
    int run,
    const unordered_map<int, TargetInfo> &targetMap){

    TargetInfo emptyInfo;

    auto it = targetMap.find(run);

    if(it == targetMap.end()){
        return emptyInfo;
    }

    TargetInfo info = it->second;
    info.target_match = 1;

    return info;
}

void applyTargetInfoToEventInfo(
    EventInfo &evInfo,
    const TargetInfo &targetInfo,
    const string &polSource){

    if(targetInfo.target_match != 1){
        return;
    }

    evInfo.target_match = 1;

    evInfo.target_species = targetInfo.species;
    evInfo.target_species_id = targetInfo.species_id;

    evInfo.target_cell = targetInfo.cell;

    evInfo.target_pol_online = targetInfo.pol_online;
    evInfo.target_pol_offline = targetInfo.pol_offline;
    evInfo.target_pol_offline_err = targetInfo.pol_offline_err;

    evInfo.target_run_dose = targetInfo.run_dose;

    evInfo.target_start_time = targetInfo.start_time;
    evInfo.target_stop_time = targetInfo.stop_time;

    float selectedPol = -9999;
    float selectedPolErr = -9999;

    if(polSource == "online"){

        selectedPol = targetInfo.pol_online;

        // The CSV does not provide online polarization uncertainty.
        selectedPolErr = -9999;

    }else{

        selectedPol = targetInfo.pol_offline;
        selectedPolErr = targetInfo.pol_offline_err;
    }

    evInfo.target_pol = selectedPol;
    evInfo.target_pol_err = selectedPolErr;

    // Treat the selected NMR polarization as vector polarization.
    // For deuteron target, calculate tensor polarization from vector polarization.
    evInfo.vector_pol = selectedPol;
    evInfo.vector_pol_err = selectedPolErr;

    evInfo.tensor_pol = -9999;
    evInfo.tensor_pol_err = -9999;

    if(evInfo.target_species_id == 2 && selectedPol > -999.0){

        double vectorPolForCalc = selectedPol;
        double tensorPolScale = 1.0;

        // The formula expects vector polarization as a fraction.
        // If the map value is stored in percent, convert for the calculation.
        // The stored tensor polarization is then converted back to the same unit.
        if(fabs(vectorPolForCalc) > 1.0){
            vectorPolForCalc = vectorPolForCalc / 100.0;
            tensorPolScale = 100.0;
        }

        double tensorArg = 4.0 - 3.0 * vectorPolForCalc * vectorPolForCalc;

        if(tensorArg >= 0.0){

            evInfo.tensor_pol = tensorPolScale * (
                2.0 - sqrt(tensorArg)
            );

            if(selectedPolErr > 0.0 && tensorArg > 0.0){

                double dTensorDVector =
                    3.0 * vectorPolForCalc / sqrt(tensorArg);

                evInfo.tensor_pol_err = fabs(dTensorDVector) * selectedPolErr;
            }
        }
    }

    if(selectedPol > 0.0){
        evInfo.target_state = 1;
    }else if(selectedPol < 0.0){
        evInfo.target_state = -1;
    }else{
        evInfo.target_state = 0;
    }
}

void addTargetBranches(
    TTree *tree,
    TargetBranchVars &targetVars){

    tree->Branch("target_match", &targetVars.target_match, "target_match/I");

    tree->Branch("target_species", &targetVars.target_species);
    tree->Branch("target_species_id", &targetVars.target_species_id, "target_species_id/I");

    tree->Branch("target_cell", &targetVars.target_cell, "target_cell/I");

    tree->Branch("target_pol_online", &targetVars.target_pol_online, "target_pol_online/F");
    tree->Branch("target_pol_offline", &targetVars.target_pol_offline, "target_pol_offline/F");
    tree->Branch("target_pol_offline_err", &targetVars.target_pol_offline_err, "target_pol_offline_err/F");

    tree->Branch("target_run_dose", &targetVars.target_run_dose, "target_run_dose/F");

    tree->Branch("target_start_time", &targetVars.target_start_time);
    tree->Branch("target_stop_time", &targetVars.target_stop_time);
}

void copyTargetInfoToBranches(
    const EventInfo &evInfo,
    TargetBranchVars &targetVars){

    targetVars.target_match = evInfo.target_match;

    targetVars.target_species = evInfo.target_species;
    targetVars.target_species_id = evInfo.target_species_id;

    targetVars.target_cell = evInfo.target_cell;

    targetVars.target_pol_online = evInfo.target_pol_online;
    targetVars.target_pol_offline = evInfo.target_pol_offline;
    targetVars.target_pol_offline_err = evInfo.target_pol_offline_err;

    targetVars.target_run_dose = evInfo.target_run_dose;

    targetVars.target_start_time = evInfo.target_start_time;
    targetVars.target_stop_time = evInfo.target_stop_time;
}

//function definition for the mapping (link different REC_branches)
vectorMap loadMapIndex(hipo::bank bankName, int detectorID){
    vectorMap idxmap;
    idxmap.clear();
    for(int nRow=0; nRow<bankName.getRows(); ++nRow){
        int rValue = bankName.getShort("pindex",nRow);
	    int detector = bankName.getByte("detector",nRow);
	    if(detector==detectorID){	
	        idxmap[rValue].push_back(nRow);
	    }
    }
    return idxmap;
}

// Collect detector/PID information for one REC::Particle row.
// Inputs:
//   pindex          = REC::Particle row number
//   bankParticle    = REC::Particle bank
//   bankCalorimeter = REC::Calorimeter bank
//   bankCherenkov   = REC::Cherenkov bank
//   Calomap         = map from REC::Particle row to REC::Calorimeter rows
//   Htccmap         = map from REC::Particle row to REC::Cherenkov rows
//   momentum        = particle momentum magnitude
// Output:
//   DetectorInfo with scalar quantities that are written to ROOT branches.

DetectorInfo getDetectorInfo(
    int pindex,
    hipo::bank bankParticle,
    hipo::bank bankCalorimeter,
    hipo::bank bankCherenkov,
    const vectorMap &Calomap,
    const vectorMap &Htccmap,
    float momentum){

    DetectorInfo info;

    // ---------- REC::Particle information ----------
    info.charge = bankParticle.getByte("charge", pindex);
    info.beta = bankParticle.getFloat("beta", pindex);
    info.chi2pid = bankParticle.getFloat("chi2pid", pindex);

    // ---------- REC::Cherenkov / HTCC information ----------
    // detector == 15 is used in the older SIDIS code for HTCC.
    // nphe is summed in case more than one Cherenkov row is linked to the same particle.

    auto Htcc_pIdxmap = Htccmap.find(pindex);

    if(Htcc_pIdxmap != Htccmap.end()){
        info.has_htcc = 1;
        info.htcc_nphe = 0.0;
        vector<int> Htcc_mapColns = (*Htcc_pIdxmap).second;
        for(auto HtccColNum : Htcc_mapColns){
            info.htcc_nphe += bankCherenkov.getFloat("nphe", HtccColNum);
        }
    }

    // ---------- REC::Calorimeter information ----------
    // detector == 7 is used for the ECAL/PCAL detector.
    // layer == 1  -> PCAL
    // layer == 4  -> ECIN
    // layer == 7  -> ECOUT
    // cal_e is the total calorimeter energy sum over all linked calorimeter rows.
    // sampling_fraction = cal_e / momentum.

    auto Calo_pIdxmap = Calomap.find(pindex);

    if(Calo_pIdxmap != Calomap.end()){

        info.has_cal = 1;

        vector<int> Calo_mapColns = (*Calo_pIdxmap).second;

        for(auto CaloColNum : Calo_mapColns){

            int cal_layer = bankCalorimeter.getByte("layer", CaloColNum);
            info.cal_sector = bankCalorimeter.getByte("sector", CaloColNum);

            float energy = bankCalorimeter.getFloat("energy", CaloColNum);
            info.cal_e += energy;

            if(cal_layer == 1){

                info.pcal_e += energy;

                info.pcal_x = bankCalorimeter.getFloat("x", CaloColNum);
                info.pcal_y = bankCalorimeter.getFloat("y", CaloColNum);

                info.lu = bankCalorimeter.getFloat("lu", CaloColNum);
                info.lv = bankCalorimeter.getFloat("lv", CaloColNum);
                info.lw = bankCalorimeter.getFloat("lw", CaloColNum);

            }else if(cal_layer == 4){

                info.ecin_e += energy;

            }else if(cal_layer == 7){

                info.ecout_e += energy;
            }
        }

        if(momentum > 0.0){
            info.sampling_fraction = info.cal_e / momentum;
        }
    }

    return info;
}

// Read DC fiducial edge information from REC::Traj for one REC::Particle row.
// This reads the DC edge quantities used in the dedicated inclusive DIS cut.

DCEdgeInfo getDCEdgeInfo(
    int pindex,
    hipo::bank bankTraj
){
    DCEdgeInfo info;

    for(int irow = 0; irow < bankTraj.getRows(); irow++){

        if(bankTraj.getInt("pindex", irow) != pindex) continue;
        if(bankTraj.getInt("detector", irow) != DEDICATED_DIS_DC_DETECTOR) continue;

        int layer = bankTraj.getInt("layer", irow);
        float edge = bankTraj.getFloat("edge", irow);

        if(layer == DEDICATED_DIS_DC_LAYER_R1){
            info.edge_r1 = edge;
            info.has_r1 = 1;
        }

        if(layer == DEDICATED_DIS_DC_LAYER_R2){
            info.edge_r2 = edge;
            info.has_r2 = 1;
        }

        if(layer == DEDICATED_DIS_DC_LAYER_R3){
            info.edge_r3 = edge;
            info.has_r3 = 1;
        }
    }

    return info;
}

// Sector dependent sampling fraction cut used by the dedicated inclusive DIS cut.

bool passInclusiveDISSamplingFraction(
    int sector,
    double p,
    double sampling_fraction
){
    if(sector < 1 || sector > 6) return false;
    if(p <= 0.0) return false;

    int idx = sector - 1;

    double sf_low =
        DEDICATED_DIS_SF_LOW[idx][0] +
        DEDICATED_DIS_SF_LOW[idx][1] * p +
        DEDICATED_DIS_SF_LOW[idx][2] * p * p;

    double sf_high =
        DEDICATED_DIS_SF_HIGH[idx][0] +
        DEDICATED_DIS_SF_HIGH[idx][1] * p +
        DEDICATED_DIS_SF_HIGH[idx][2] * p * p;

    if(sampling_fraction <= sf_low) return false;
    if(sampling_fraction >= sf_high) return false;

    return true;
}

// PCAL fiducial cut used by the dedicated inclusive DIS cut.
// Sector 1 uses a tighter lv/lw minimum than sectors 2 through 6.

bool passInclusiveDISPCALFiducial(
    int sector,
    double lv,
    double lw
){
    if(sector < 1 || sector > 6) return false;

    double min_lvw = DEDICATED_DIS_PCAL_LVW_S26;

    if(sector == 1){
        min_lvw = DEDICATED_DIS_PCAL_LVW_S1;
    }

    if(lv <= min_lvw) return false;
    if(lw <= min_lvw) return false;

    return true;
}

// DC fiducial edge cut used by the dedicated inclusive DIS cut.

bool passInclusiveDISDCFiducial(
    const DCEdgeInfo &dcInfo
){
    if(dcInfo.has_r1 != 1) return false;
    if(dcInfo.has_r2 != 1) return false;
    if(dcInfo.has_r3 != 1) return false;

    if(dcInfo.edge_r1 <= DEDICATED_DIS_DC_EDGE_R1) return false;
    if(dcInfo.edge_r2 <= DEDICATED_DIS_DC_EDGE_R2) return false;
    if(dcInfo.edge_r3 <= DEDICATED_DIS_DC_EDGE_R3) return false;

    return true;
}

// Dedicated inclusive DIS detector cut.
// This should be used only for inclusive DIS when detpidcut is selected.

bool passInclusiveDISDetectorCuts(
    const ParticleCandidate &electron,
    const DetectorInfo &eDet,
    const DCEdgeInfo &dcInfo
){
    if(eDet.has_htcc != 1) return false;
    if(eDet.htcc_nphe < DEDICATED_DIS_NPHE_MIN) return false;

    if(eDet.has_cal != 1) return false;
    if(eDet.cal_e <= 0.0) return false;
    if(eDet.cal_sector < 1 || eDet.cal_sector > 6) return false;

    if(eDet.pcal_e <= DEDICATED_DIS_PCAL_E_MIN) return false;

    if(!passInclusiveDISSamplingFraction(
        eDet.cal_sector,
        electron.p,
        eDet.sampling_fraction
    )){
        return false;
    }

    if(!passInclusiveDISPCALFiducial(
        eDet.cal_sector,
        eDet.lv,
        eDet.lw
    )){
        return false;
    }

    if(!passInclusiveDISDCFiducial(dcInfo)){
        return false;
    }

    return true;
}

// Dedicated inclusive DIS kinematic cut.
// This should be used only for inclusive DIS when kincut is selected.

bool passInclusiveDISKinematicCuts(
    const ParticleCandidate &electron,
    const InclusiveKin &kin
){
    double theta_deg = electron.theta * 180.0 / TMath::Pi();

    if(kin.W <= DEDICATED_DIS_W_MIN) return false;
    if(kin.Q2 <= DEDICATED_DIS_Q2_MIN) return false;
    if(electron.p <= DEDICATED_DIS_EP_MIN) return false;

    if(theta_deg <= DEDICATED_DIS_THETA_MIN_DEG) return false;
    if(theta_deg >= DEDICATED_DIS_THETA_MAX_DEG) return false;

    if(electron.vz < DEDICATED_DIS_VZ_MIN) return false;
    if(electron.vz > DEDICATED_DIS_VZ_MAX) return false;

    return true;
}

// Dedicated inclusive QE detector cut.
// This should be used only for inclusive QE when detpidcut is selected.

bool passInclusiveQEDetectorCuts(
    const ParticleCandidate &electron,
    const DetectorInfo &eDet
){
    if(electron.p <= 0.0) return false;

    if(fabs(eDet.chi2pid) >= DEDICATED_QE_ABS_CHI2PID_MAX) return false;

    if(eDet.has_cal != 1) return false;

    if(eDet.lv <= DEDICATED_QE_LV_MIN) return false;
    if(eDet.lw <= DEDICATED_QE_LW_MIN) return false;

    if(eDet.pcal_e <= DEDICATED_QE_PCAL_E_MIN) return false;

    double sf_qe = (eDet.pcal_e + eDet.ecin_e) / electron.p;

    if(sf_qe >= DEDICATED_QE_SF_MAX) return false;

    double pcal_over_p = eDet.pcal_e / electron.p;
    double ecin_over_p = eDet.ecin_e / electron.p;

    double ecin_min =
        DEDICATED_QE_ECIN_PCAL_SLOPE * pcal_over_p +
        DEDICATED_QE_ECIN_PCAL_INTERCEPT;

    if(ecin_over_p < ecin_min) return false;

    return true;
}

// Dedicated inclusive QE kinematic cut.
// This should be used only for inclusive QE when kincut is selected.

bool passInclusiveQEKinematicCuts(
    const ParticleCandidate &electron,
    const InclusiveKin &kin
){
    double theta_deg = electron.theta * 180.0 / TMath::Pi();

    if(theta_deg <= DEDICATED_QE_THETA_MIN_DEG) return false;
    if(theta_deg >= DEDICATED_QE_THETA_MAX_DEG) return false;

    if(electron.vz <= DEDICATED_QE_VZ_MIN) return false;
    if(electron.vz >= DEDICATED_QE_VZ_MAX) return false;

    if(electron.p <= DEDICATED_QE_P_MIN) return false;

    if(kin.Q2 <= DEDICATED_QE_Q2_MIN) return false;
    if(kin.Q2 >= DEDICATED_QE_Q2_MAX) return false;

    if(kin.W <= DEDICATED_QE_W_MIN) return false;
    if(kin.W >= DEDICATED_QE_W_MAX) return false;

    return true;
}

// Loose detector/PID cuts for inclusive electron only mode.
// This is the electron part of passDetPidCuts(), separated so inclusive mode
// does not need a dummy hadron object.

bool passElectronDetPidCuts(const DetectorInfo &eDet){

    if(eDet.charge != -1) return false;

    if(eDet.has_htcc != 1) return false;
    if(eDet.htcc_nphe <= 2.0) return false;

    if(eDet.has_cal != 1) return false;
    if(eDet.pcal_e <= 0.0) return false;
    if(eDet.cal_e <= 0.0) return false;

    if(eDet.sampling_fraction <= 0.10) return false;

    if(eDet.lu <= 0.0) return false;
    if(eDet.lv <= 0.0) return false;
    if(eDet.lw <= 0.0) return false;

    return true;
}

// Loose detector/PID cuts used only when --detpidcut 1 is selected.
// These cuts are intentionally conservative.
// They remove obvious bad objects, but avoid SIDIS kinematic cuts.
// No Q2, W, x, z, pT, or missing mass cuts are applied here.
// You can tighten these numbers later after checking QA plots from the skim.

bool passDetPidCuts(
    const DetectorInfo &eDet,
    const DetectorInfo &hDet,
    int selected_hadron_pid){

    if(!passElectronDetPidCuts(eDet)) return false;

    // ---------- Hadron detector/PID sanity cuts ----------
    // The skim already required:
    //   pid_h == selected_hadron_pid
    //   h_det_status == 2
    // Here we only add very loose sanity cuts.

    if(hDet.beta <= 0.0) return false;

    // Check charge consistency with selected PID.
    // pi+, K+, proton should have positive charge.
    // pi-, K- should have negative charge.
    if(selected_hadron_pid > 0 && hDet.charge <= 0) return false;
    if(selected_hadron_pid < 0 && hDet.charge >= 0) return false;

    return true;
}

// Generic electron hadron SIDIS event test.
// This replaces the pi only event test for the skim.
// It checks whether the event has:
//   1. at least one trigger electron: pid = 11 and det_status = -2
//   2. at least one selected hadron:  pid = selected_hadron_pid and det_status = 2
// This allows the same code to skim: epi+, epi-, eK+, eK-, e proton, etc, without writing separate .cc files.

void eh_sidis_event_test(hipo::bank bankParticle, int selected_hadron_pid, bool &eh_sidis){

    int n_ele = 0;
    int n_hadron = 0;

    auto n_particles = bankParticle.getRows();

    for(int n_part = 0; n_part < n_particles; ++n_part){

        int pid_n = bankParticle.getInt("pid", n_part);
        short part_status = bankParticle.getShort("status", n_part);
        int det_status = part_status / 1000;

        if(pid_n == 11 && det_status == -2){
            ++n_ele;
        }

        if(pid_n == selected_hadron_pid && det_status == 2){
            ++n_hadron;
        }
    }

    if(n_ele >= 1 && n_hadron >= 1){
        eh_sidis = true;
    }else{
        eh_sidis = false;
    }

    return;
}

void epippim_dihadron_event_test(
    hipo::bank bankParticle,
    bool &dihadron_event
){
    int n_ele = 0;
    int n_pip = 0;
    int n_pim = 0;

    auto n_particles = bankParticle.getRows();

    for(int n_part = 0; n_part < n_particles; ++n_part){

        int pid_n = bankParticle.getInt("pid", n_part);
        short part_status = bankParticle.getShort("status", n_part);
        int det_status = part_status / 1000;

        if(pid_n == 11 && det_status == -2){
            n_ele++;
        }

        if(pid_n == 211 && det_status == 2){
            n_pip++;
        }

        if(pid_n == -211 && det_status == 2){
            n_pim++;
        }
    }

    dihadron_event = (n_ele >= 1 && n_pip >= 1 && n_pim >= 1);
}

// Convert the selected hadron PID into a short readable tag.
// The tag is used in output ROOT file names.

string getHadronTag(int pid){

    if(pid == 211)  return "pip";
    if(pid == -211) return "pim";

    if(pid == 321)  return "kp";
    if(pid == -321) return "km";

    if(pid == 2212) return "p";

    return "pid" + std::to_string(pid);
}


bool stringIsAllDigits(const string &s){
    if(s.empty()){
        return false;
    }

    for(char c : s){
        if(!isdigit(static_cast<unsigned char>(c))){
            return false;
        }
    }

    return true;
}

int extractInputFileNumber(const string &baseName){
    string stem = baseName;

    size_t hipoPos = stem.rfind(".hipo");
    if(hipoPos != string::npos){
        stem = stem.substr(0, hipoPos);
    }

    size_t lastDot = stem.find_last_of('.');
    if(lastDot == string::npos){
        return -9999;
    }

    string token = stem.substr(lastDot + 1);

    if(!stringIsAllDigits(token)){
        return -9999;
    }

    return atoi(token.c_str());
}

// Calculate inclusive electron scattering kinematics.
// Important convention:
//   W and W2 are nucleon level quantities.
//   These are used for DIS/resonance separation:
//      DIS        W > 2 GeV
//      resonance  0 < W < 2 GeV
//   Wd and Wd2 are deuteron system quantities.
//   These are kept separately for tensor deuteron studies.
// Definitions:
//   q  = k - k'
//   Q2 = -q^2
//   nu = q0
//   W_N^2 = M_N^2 + 2 M_N nu - Q2
//   W_D^2 = M_D^2 + 2 M_D nu - Q2
//   xB = Q2 / (2 M_N nu)
//   xD = Q2 / (2 M_D nu)
//   y  = nu / Ebeam

InclusiveKin calculateInclusiveKinematics(
    const ParticleCandidate &electron){

    InclusiveKin kin;

    TLorentzVector e_in4vect(0.0, 0.0, beam_enrg, beam_enrg);

    double e_energy = sqrt(
        electron.p * electron.p + M_electron * M_electron
    );

    TLorentzVector e_out4vect(
        electron.px,
        electron.py,
        electron.pz,
        e_energy
    );

    TLorentzVector q_4vect = e_in4vect - e_out4vect;

    kin.Q2 = -q_4vect.M2();
    kin.nu = q_4vect.E();

    // Use proton mass as the default nucleon mass.
    // For ND3/deuteron inclusive analyses, this is the standard practical W
    // used to separate DIS and nucleon resonance regions.
    const double M_N = M_p;

    // nucleon level W.
    kin.W2 = M_N * M_N + 2.0 * M_N * kin.nu - kin.Q2;

    if(kin.W2 > 0.0){
        kin.W = sqrt(kin.W2);
    }else{
        kin.W = -sqrt(fabs(kin.W2));
    }

    // deuteron system W.
    kin.Wd2 = M_d * M_d + 2.0 * M_d * kin.nu - kin.Q2;

    if(kin.Wd2 > 0.0){
        kin.Wd = sqrt(kin.Wd2);
    }else{
        kin.Wd = -sqrt(fabs(kin.Wd2));
    }

    // Bjorken variables.
    if(fabs(2.0 * M_N * kin.nu) > 1.0e-12){
        kin.xB = kin.Q2 / (2.0 * M_N * kin.nu);
    }

    if(fabs(2.0 * M_d * kin.nu) > 1.0e-12){
        kin.xD = kin.Q2 / (2.0 * M_d * kin.nu);
    }

    if(fabs(beam_enrg) > 1.0e-12){
        kin.y = kin.nu / beam_enrg;
    }

    // Broad inclusive quasi elastic candidate selection.
    // This is deliberately loose for the skim.
    // It keeps the QE peak and xB > 1 region for later tensor analysis.
    // Do not put narrow Hall B note cuts here. Those belong in rgc_analysis.cc.
    if(kin.Q2 > 0.2 &&
    kin.xB > 0.5 &&
    kin.xB < 2.1 &&
    kin.y > 0.0 &&
    kin.y < 1.0){
        kin.is_qe = 1;
    }else{
        kin.is_qe = 0;
    }

    // Region flags based on nucleon W.
    if(kin.W > 2.0){
        kin.is_dis = 1;
        kin.is_resonance = 0;
    }else if(kin.W > 0.0 && kin.W < 2.0){
        kin.is_dis = 0;
        kin.is_resonance = 1;
    }else{
        kin.is_dis = 0;
        kin.is_resonance = 0;
    }

    return kin;
}

// inclusive region filter. This function decides whether an inclusive electron event should be written.
// Allowed regions:
//   all = keep all events
//   dis = keep W > 2
//   res = keep 0 < W < 2
//   qe  = keep broad quasi elastic candidates based on xB, Q2, and y
// The actual W flags are calculated in calculateInclusiveKinematics().

bool passInclusiveRegion(
    const InclusiveKin &kin,
    const string &inclusiveRegion){

    if(inclusiveRegion == "all"){
        return true;
    }

    if(inclusiveRegion == "dis"){
        return kin.is_dis == 1;
    }

    if(inclusiveRegion == "res"){
        return kin.is_resonance == 1;
    }

    if(inclusiveRegion == "qe"){
        return kin.is_qe == 1;
    }

    return false;
}


