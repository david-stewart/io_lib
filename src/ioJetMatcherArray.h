#ifndef ioJetMatcherArray__h
#define ioJetMatcherArray__h

#include <vector>
#include <string>

#include "TH1D.h"
#include "TH2D.h"

#include "RooUnfoldResponse.h"
#include "ioClass.h"
#include "ioCfnc.h"
#include "TRandom3.h"

struct ioJetMatcherArray {
    // like above, but uses it's own ioXsec
    public:
    string name;
    bool write_9 { false };
    int  cull_n  { 10 };
    double cut_high_sigma {0};
    double cut_high_sigma_offset{8.};
    bool apply_Mlimit { true };

    vector<string> v_names;
    ioXsec& Xsec;

    void cull_add_array(array<TH2D*,9>&, string which, const char* tag="");
    void cull_add_array(array<TH1D*,9>&, string which, const char* tag="");
    void write_response(TH2D* match, TH1D* miss, string which, const char* tag);

    vector<array<TH2D*,9>> v_response;
    vector<array<TH1D*,9>> v_truth;

    vector<array<TH2D*,9>> v_response_A;
    vector<array<TH1D*,9>> v_truth_A;

    vector<array<TH2D*,9>> v_response_B;
    vector<array<TH1D*,9>> v_truth_B;

    /* vector<RooUnfoldResponse> v_response; // */ 
    /* vector<RooUnfoldResponse> v_response_A; // */ 
    /* vector<RooUnfoldResponse> v_response_B; // */ 
    TRandom3 _rand;

    TH1D       hg_pthb_cnt;
    array<ioXYbounder,9> pthb_Mlimit; // pT-hat-bin Measured Limit

    ioJetMatcherArray (
            const char* _name,
            ioXsec& _Xsec, 
            const char* bin_file, 
            vector<string> bin_names,
            vector<string> bin_tags_M,
            vector<string> bin_tags_T,
            const char* pthb_Mlimit_file,
            double ratio_AtoB=0.3,
            bool debug = true
    );
    ~ioJetMatcherArray(){};

    void addjet_MC(float eta, float phi, float pT);
    void addjet_reco(float eta, float phi, float pT);
    bool do_matching(int pthatbin);
    ioXYbounder out_of_match_bounds {};
    double ratio_AtoB;
    void write();
    bool debug;

    private: //internal data to do the matching
    float jet_R2; // jet_R * jet_R   
	vector<ioJetMatcher_float> data_MC;
	vector<ioJetMatcher_float> data_reco;
};

#endif
