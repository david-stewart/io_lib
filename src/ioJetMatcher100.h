#ifndef ioJetMatcher100__h
#define ioJetMatcher100__h

#include <vector>
#include <string>

#include "TH1D.h"
#include "TH2D.h"

#include "RooUnfoldResponse.h"
#include "ioClass.h"
#include "ioCfnc.h"
#include "TRandom3.h"

struct ioJetMatcher100 {
    // like above, but uses it's own ioXsec
    public:

    /* vector<string> v_names; */
    ioXsec* Xsec{nullptr};
    TH1D hg_pthb_cnt;
    TH1D* sigma;

    array<TH2D*,9> v_response, v_response_A, v_response_B;
    array<TH1D*,9> v_truth,    v_truth_A,    v_truth_B;
    TH2D *hg2_response, *hg2_response_A, *hg2_response_B;
    TH1D *hg1_truth,    *hg1_truth_A,    *hg1_truth_B; 
    TH1D *hg1_measured;

    TRandom3 _rand;
    ioJetMatcher100 (
            ioXsec* _Xsec, 
            double ratio_AtoB=0.3,
            int _ncull=20,
            double high_sig_cut=5.,
            double high_sig_off=8,
            bool debug = false
    );
    ioJetMatcher100 ( const char* file );
    ~ioJetMatcher100(){};

    RooUnfoldResponse* make_ruu(ioBinVec bins_M, ioBinVec bins_T, const char* name, bool write=true);

    void addjet_MC(float eta, float phi, float pT);
    void addjet_reco(float eta, float phi, float pT);
    bool do_matching(int pthatbin);
    ioXYbounder out_of_match_bounds {};
    double ratio_AtoB;
    void write();
    bool debug;
    int    cull_n                { 20 };
    double cut_high_sigma        { 4. };
    double cut_high_sigma_offset { 8. };

    TH2D* rebin(TH2D*, ioBinVec bins_M, ioBinVec bins_T);
    TH1D* rebin(TH1D*, ioBinVec bins_T);
    void process_arrays(array<TH2D*,9>&, TH2D*&, array<TH1D*,9>&, TH1D*&, const char* tag="");

    private: //internal data to do the matching
    float jet_R2; // jet_R * jet_R   
	vector<ioJetMatcher_float> data_MC;
	vector<ioJetMatcher_float> data_reco;
};

#endif
