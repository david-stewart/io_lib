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
    ioXYbounder out_of_match_bounds;

    /* vector<double> fakes; */
    /* vector<double> missed; */
    /* vector<pair<double,double>> matches; */
    /* ioInBounds bounds_T; */
    /* ioInBounds bounds_M; */
    /* bool cut_outliers{false}; */

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
                "hg2_Xsec_vs_T",1,
                "hg2_Xsec_vs_M",1,
                "hg2_Xsec_vs_match",1,
                "hg2_Xsec_vs_fake",1,
                "hg1_R2_match",1,
                "R",0.4,
                "cut_outliers",0
            }}
    );
    double pt_fakes, pt_misses;

    bool b_make_AB;
    bool b_hg2_Xsec_vs_T;
    bool b_hg2_Xsec_vs_M;
    bool b_hg2_Xsec_vs_match;
    bool b_hg2_Xsec_vs_fake;
    bool b_hg1_R2_match;
    /* bool switch_AB{true}; */

    RooUnfoldResponse* response_A; // 1/2 of the data
    RooUnfoldResponse* response_B; // other 1/2 of data
    TH2D* hg2_Xsec_vs_T{nullptr};
    TH2D* hg2_Xsec_vs_M{nullptr};
    TH2D* hg2_Xsec_vs_match{nullptr};
    TH2D* hg2_Xsec_vs_fake{nullptr};
    TH1D* hg1_R2_match{nullptr};

    void addjet_MC(float eta, float phi, float pT);
    void addjet_reco(float eta, float phi, float pT);
    bool do_matching(int pthatbin); // return true if successful matching

    TRandom3 _rand;
    double ratio_AtoB;

    void write();

    private: //internal data to do the matching
    float jet_R2; // jet_R * jet_R   
	vector<ioJetMatcher_float> data_MC;
	vector<ioJetMatcher_float> data_reco;
};

#endif
