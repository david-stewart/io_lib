#ifndef ioJetMatcher__h
#define ioJetMatcher__h

#include <vector>
#include <string>

#include "TH1D.h"
#include "TH2D.h"

#include "RooUnfoldResponse.h"
#include "ioClass.h"
#include "ioCfnc.h"
#include "TRandom3.h"

using std::vector;
using std::sort;

double* ioEdges_pAuJet_prelim_13bins(); // 14 edges for 13 bins

struct ioJetMatcher_float {
    // For filling vectors for Raghav matching algorithm 
    // will have comparison methods in order to std::sort
    float eta;
    float phi;
    float pT;
    bool  is_matched;
	double operator()(const ioJetMatcher_float&, const float jetR2=0.16);
    ioJetMatcher_float(float,float,float);
    ioJetMatcher_float();
};

struct ioJetMatcher_index {
    // for filling vectors for CGAL matching algorithm
    // will have comparison methods in order to std::sort
    int i_MC;   // index MC jet
    int i_reco; // index reco jet
    float dR;   // distance between the matches
};

struct ioJetMatcher {
    // like above, but uses it's own ioXsec
    public:
    string name;
    ioXsec& Xsec;
    RooUnfoldResponse response_noweight;
    RooUnfoldResponse response;
    RooUnfoldResponse response_cut;
    TH2D hg2_cut;
    TRandom3 _rand;
    ioInBounds in_reco_bounds;
    ioInBounds in_meas_bounds;

    TH1D       hg_pthb_cnt;
    vector<TH1D> v_miss{};
    vector<TH1D> v_T{};
    vector<TH1D> v_M{};
    vector<TH1D> v_fake{};
    vector<TH2D> v_match{};

    vector<TH1D> A_miss{};
    vector<TH1D> A_fake{};
    vector<TH2D> A_match{};

    vector<TH1D> B_miss{};
    vector<TH1D> B_fake{};
    vector<TH2D> B_match{};

    double fake_limit;
    ioXYbounder out_of_match_bounds {};
    bool b_ptht_Mlimit{false};
    array<ioXYbounder,9> pthb_Mlimit; // pT-hat-bin Measured Limit

    /* int* eventid; */ 
    /* int* runid; */
    map<int,int> outlier_ids; // runid->eventid

    ioJetMatcher (const char* name, ioXsec& _Xsec, 
            const char* bin_file, 
            const char* tag_M, 
            const char* tag_T, 
            ioOptMap options={}, ioBinVec _hg2ptbins={81,-0.5,80.5}, 
            ioOptMap dict={{
                "make_AB",1,
                "ratio_AtoB",0.5,
                "Xsec_vs_T",1,
                "Xsec_vs_M",1,
                "Xsec_vs_match",1,
                "Xsec_vs_fake",1,
                "Xsec_vs_miss",1,
                "hg1_R2_match",1,
                "R",0.4
            }}
    );
    ~ioJetMatcher(){};
    double pt_fakes, pt_misses;

    /* bool b_make_AB       {false}; */
    bool b_Xsec_vs_miss  {false};
    bool b_Xsec_vs_T     {false};
    bool b_Xsec_vs_M     {false};
    bool b_Xsec_vs_fake  {false};
    bool b_Xsec_vs_match {false};
    bool b_hg1_R2_match  {false};
    /* bool switch_AB{true}; */

    RooUnfoldResponse* response_A; // 1/2 of the data
    RooUnfoldResponse* response_B; // other 1/2 of data
    TH1D* hg1_R2_match{nullptr};


    void addjet_MC(float eta, float phi, float pT);
    void addjet_reco(float eta, float phi, float pT);
    bool do_matching(int pthatbin); // return true if successful matching

    double ratio_AtoB;

    void write();

    private: //internal data to do the matching
    float jet_R2; // jet_R * jet_R   
	vector<ioJetMatcher_float> data_MC;
	vector<ioJetMatcher_float> data_reco;
};


#endif
