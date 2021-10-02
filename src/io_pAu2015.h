#ifndef io_pAu2015__h
#define io_pAu2015__h

#include "ioClass.h"
#include "io_fnc.h"
#include "RooUnfoldResponse.h"
#include "TH1D.h"

struct io_Unfolder_Set {
    RooUnfoldResponse* ruu;
    RooUnfoldResponse* ruu_A;
    RooUnfoldResponse* ruu_B;
    RooUnfoldResponse* ruu_TS;
    RooUnfoldResponse* ruu_TU;
    RooUnfoldResponse* ruu_HC50;

    // these three histograms keep the uncertainty of the closure of A with B
    TH1D* truth_A;
    TH1D* closure_A;
    TH1D* rat_A;
    ioPtrDbl ptr_closure;

    TH1D* hg_IP2; // unfolded spectra
    TH1D* rat_IP2;// uncertainty relative to main unfolded spectra
    ioPtrDbl ptr_IP2;

    TH1D* hg_IP4; // main unfolded spectra

    TH1D* hg_IP6; // unfolded spectra
    TH1D* rat_IP6;// uncertainty relative to main unfolded spectra
    ioPtrDbl ptr_IP6;// uncertainty relative to main unfolded spectra

    TH1D* hg_TS; // unfolded spectra
    TH1D* rat_TS;// uncertainty relative to main unfolded spectra
    ioPtrDbl ptr_TS;// uncertainty relative to main unfolded spectra

    TH1D* hg_TU; // unfolded spectra
    TH1D* rat_TU;// uncertainty relative to main unfolded spectra
    ioPtrDbl ptr_TU;// uncertainty relative to main unfolded spectra

    TH1D* hg_HC50; // unfolded spectra
    TH1D* rat_HC50;// uncertainty relative to main unfolded spectra
    ioPtrDbl ptr_HC50;// uncertainty relative to main unfolded spectra

    ioPtrDbl unf_unc(TH1D* raw, TH1D*& _hg, TH1D*& _rat, RooUnfoldResponse*& _ruu, int n_iter=2);


    io_Unfolder_Set( 
            const char* tag="one",
            const char* file_name="sys_err.root");
    TH1D* unfold(TH1D* hg);
};





// uncertainties to add: (from Isaac Mooney's jet analysis note)
// 1. IP2(6)Iteration: iterations in Bayesian unfoling instead of 4 (fine, no problem! -- not done here)
// 2. TS - tower scale uncertainty of 3.8%
// 3. TU - tracking efficiency uncertainty 4%
// 4. HC50 -- hadronic subtraction uncertainty, use 50% instead of 100%
// 5. GS -- generator smearing of response matrix (10%?)
 
#endif
