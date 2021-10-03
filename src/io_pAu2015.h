#ifndef io_pAu2015__h
#define io_pAu2015__h

#include "ioClass.h"
#include "io_fnc.h"
#include "RooUnfoldResponse.h"
#include "TH1D.h"

void io_draw_boxes(TH1D* x_axis, ioPtrDbl err, int color, double alpha, 
        double i_left=2, double i_right=3, double i_total=5) {
};

struct io_Unfolder_Set {
    bool scale_by_binW {true};

    TH1D* hg_base; // main unfolded spectra
    ioPtrDbl pts_base;
    ioPtrDbl pts_sumsq;

    int n_iter;
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
    ioPtrDbl pts_rat_A;
    ioPtrDbl pts_closure_A;

    TH1D*     hg_base_A;

    TH1D* hg_IPm2; // unfolded spectra
    TH1D* rat_IPm2;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_IPm2;

    TH1D* hg_IPp2; // unfolded spectra
    TH1D* rat_IPp2;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_IPp2;// uncertainty relative to main unfolded spectra

    TH1D* hg_TS; // unfolded spectra
    TH1D* rat_TS;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_TS;// uncertainty relative to main unfolded spectra

    TH1D* hg_TU; // unfolded spectra
    TH1D* rat_TU;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_TU;// uncertainty relative to main unfolded spectra

    TH1D* hg_HC50; // unfolded spectra
    TH1D* rat_HC50;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_HC50;// uncertainty relative to main unfolded spectra

    TH1D* hg_closure;
    TH1D* pts_closure;

    ioPtrDbl unf_unc(TH1D* raw, TH1D*& _hg, TH1D*& _rat, RooUnfoldResponse*& _ruu, int n_iter);


    io_Unfolder_Set( 
            const char* tag="one",
            const char* file_name="sys_err.root",
            int _n_iter=7);
    TH1D* unfold(TH1D* hg);
};


struct io_Unfolder_Set_Ratio {
    io_Unfolder_Set_Ratio( io_Unfolder_Set& num, io_Unfolder_Set& den);
     
    TH1D*    hg_base;
    ioPtrDbl pts_base;

    TH1D* hg_IPm2;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_IPm2;

    TH1D* hg_IPp2;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_IPp2;

    TH1D* hg_TS;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_TS;// uncertainty relative to main unfolded spectra

    TH1D* hg_TU;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_TU;// uncertainty relative to main unfolded spectra

    TH1D* hg_HC50;// uncertainty relative to main unfolded spectra
    ioPtrDbl pts_HC50;// uncertainty relative to main unfolded spectra

    TH1D*    hg_base_A;
    ioPtrDbl pts_base_A;

    ioPtrDbl pts_sumsq;

};



/* // uncertainties to add: (from Isaac Mooney's jet analysis note) */
/* // 1. IP2(6)Iteration: iterations in Bayesian unfoling instead of 4 (fine, no problem! -- not done here) */
/* // 2. TS - tower scale uncertainty of 3.8% */
/* // 3. TU - tracking efficiency uncertainty 4% */
/* // 4. HC50 -- hadronic subtraction uncertainty, use 50% instead of 100% */
/* // 5. GS -- generator smearing of response matrix (10%?) */
 
#endif
