#ifndef ioJetMatcherGoodBins__h
#define ioJetMatcherGoodBins__h

#include <vector>
#include <string>

#include "TH1D.h"
#include "TH2D.h"

#include "RooUnfoldResponse.h"
#include "ioJetMatcher.h"
#include "ioClass.h"
#include "ioCfnc.h"
#include "TRandom3.h"

struct ioJetMatcherGoodBins {
    // like above, but uses it's own ioXsec
    public:
    string name;
    TRandom3 _rand;
    ioOptMap options;
    int nXsec;
    double jet_R2; // jet_R * jet_R   
    TH1D  hg_pthb_cnt;
    bool debug;
    double ratio_AtoB;
    TH1D hg1_R2_match;
    bool limit_to_goodbins;

    ioInBounds in_T_bounds;
    ioInBounds in_M_bounds;

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

    TH2D* hg_goodmatch {nullptr};
    TH1D* hg_goodmiss {nullptr};
    TH1D* hg_goodfake {nullptr};
    vector<ioIntSet> goodmatchbins;
    vector<ioIntSet> goodmissbins;
    vector<ioIntSet> goodfakebins;

    ioJetMatcherGoodBins (
            const char* name, 
            const char* bin_file, 
            const char* tag_M, 
            const char* tag_T, 
            const char* goodbin_file="",
            ioOptMap opt={},
            ioOptMap dict={{ 
                "R",0.4, 
                "tag_goodbins_M", "bins_M",
                "tag_goodbins_T", "bins_T",
                "nXsec", 9,
                "ratio_AtoB",0.5,
                "debug",false
                
            }}
    );
    ~ioJetMatcherGoodBins(){};

    void addjet_MC(float eta, float phi, float pT);
    void addjet_reco(float eta, float phi, float pT);
    bool do_matching(int pthatbin); // return true if cut for not in a good bin
    bool clear_MCreco();


    void write();

    private: //internal data to do the matching
	vector<ioJetMatcher_float> data_MC;
	vector<ioJetMatcher_float> data_reco;
};

#endif
