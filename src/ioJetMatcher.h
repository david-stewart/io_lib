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

/* class ioJetMatcher { */
/*     public: */
/*     TH1D* R2_distr_matches; */
/*     RooUnfoldResponse response; */
/*     RooUnfoldResponse* response_A; // 1/2 of the data */
/*     RooUnfoldResponse* response_B; // other 1/2 of data */
   
/*     void addjet_MC(float eta, float phi, float pT); */
/*     void addjet_reco(float eta, float phi, float pT); */
/*     void do_matching_highfirst(double weight=1.);// Raghav's algorithm */
/*     void do_matching(double weight=1.); // my algorithm */

/*     ioJetMatcher( RooUnfoldResponse response, float _jet_R=0.4 ); */
/*     // other constructors use the io_fnc iMakeRooUnfoldResponse */
/*     ioJetMatcher( int nbins, double lo_bin, double hi_bin, */ 
/*                   const char* tag="", const char* title="", */
/*                   float _jet_R=0.4 ); */
/*     ioJetMatcher( int nbins, double* edges, */ 
/*                   const char* tag="", const char* title="", */
/*                   float _jet_R=0.4 ); */
/*     ioJetMatcher( int nb_meas,  double lo_meas,  double hi_mea, */
/*                   int nb_truth, double lo_truth, double hi_truth, */
/*                   const char* tag="", const char* title="", */
/*                   float _jet_R=0.4 ); */
/*     ioJetMatcher( int nb_meas,  double* edge_meas, */
/*                   int nb_truth, double* edge_truth, */
/*                   const char* tag="", const char* title="", */
/*                   float _jet_R=0.4 ); */
/*     // final *and best* constructor gives file, and tag names for binnings */
/*     ioJetMatcher( const char* edge_file, */ 
/*                   const char* meas_tag, */
/*                   const char* truth_tag, */
/*                   const char* name_tag="", */
/*                   const char* title="", */
/*                   float _jet_R=0.4 */ 
/*     ); */
/*     void init(float _jet_R); */
/*     bool fill_A{true}; */

/*     void reset(); */
/*     void write(); */
/*     /1* void write(bool with_miss_fakes=true, *1/ */ 
/*                /1* bool scale_by_bin_width=true, *1/ */
/*                /1* bool write_unified2D=false); *1/ */
/*     std::string tag; */
    
/*     private: //internal data to do the matching */
/*     float jet_R2; // jet_R * jet_R */   
/* 	vector<ioJetMatcher_float> data_MC; */
/* 	vector<ioJetMatcher_float> data_reco; */

/* }; */

struct ioJetMatcher_outlier {
    ioJetMatcher_outlier();
    
    //---------------------------------
    // optionally initialize and use
    //---------------------------------
    void init( 
        map<int,double> boundaries,
        double  _pt_fakes,
        double  _pt_misses,
        ioXsec* _Xsec
    );
    ioXsec* Xsec;
    map<int,double> pt_limits{}; // entries are <pthatbin, limit>
    double pt_fakes, pt_misses; // intiali
    bool is_set{false}; // if not is set, simply always ignore

    // use
    bool operator()(int, double); // check if it is an outlier

    // collect the values of out-of bound values
    map<int, vector<double>> M_matches{};
    map<int, vector<double>> T_matches{};
    map<int, vector<double>> misses{}; // plotted on y-axis at x=pt-misses
    map<int, vector<double>> fakes{};  // plotted on x-axis at y=pt-fakes

    // how to collect above values
    void add_fake(int pthatbin, double pt);
    void add_miss(int pthatbin, double pt);
    void add_match(int pthatbin, double M_pt, double T_pt);

    // write out the 2-D TGraph of matches outliers
    // including 1-D outliers of any misses and fakes
    void write_TGraph(
        int pthatbin, 
        ioOptMap options={}, 
        ioOptMap dict={{
        "MarkerStyle", kFullCircle,
        "MarkerColor", kBlack,
        "MarkerFakeMiss", kOpenCircle,
        "name","outlier_"}}
    );
    void write_TGraph(
        int pthatbin, const char* name, 
        int markerstyle, int markercolor,
        ioOptMap options={},
        ioOptMap dict={{"MarkerFakeMiss",0}}
    );
};

struct ioJetMatcherX {
    // like above, but uses it's own ioXsec
    public:
    string name;
    ioXsec& Xsec;
    RooUnfoldResponse response_noweight;
    RooUnfoldResponse response;
    ioInBounds bounds_T;
    ioInBounds bounds_M;

    bool cut_outliers{false};

    /* int* eventid; */ 
    /* int* runid; */
    map<int,int> outlier_ids; // runid->eventid

    ioJetMatcherX (const char* name, ioXsec& _Xsec, 
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

    ioJetMatcher_outlier& outlier(char); // 'F' 'M' or 'T'
    void set_outlier(char C, map<int, double> boundaries);
    ioJetMatcher_outlier fake_outliers{};
    ioJetMatcher_outlier M_outliers{};
    ioJetMatcher_outlier T_outliers{};
    void set_outlier(ioJetMatcher_outlier&, map<int,double> pt_limits);

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
    void do_matching(int pthatbin);
    void do_matching(pair<double,double> pthatrange);

    TRandom3 _rand;
    double ratio_AtoB;

    void reset();
    void write();

    private: //internal data to do the matching
    float jet_R2; // jet_R * jet_R   
	vector<ioJetMatcher_float> data_MC;
	vector<ioJetMatcher_float> data_reco;
};

#endif
