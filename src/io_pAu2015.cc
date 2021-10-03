#include "io_pAu2015.h"

ioPtrDbl io_Unfolder_Set::unf_unc(TH1D* raw, TH1D*& _hg, 
        TH1D*& _rat, RooUnfoldResponse*& _ruu, int _n_iter) {
    _hg =  io_BayesUnfold(raw, _ruu, _n_iter);
    if (scale_by_binW) io_scaleByBinWidth(_hg);
    _rat =  ioDivideTH1(_hg, hg_base);
    ioPtrDbl pts { _hg };
    /* pts -= hg_ref; */
    /* pts -= 1.; */
    /* cout << " pts: " << pts << endl; */
    return pts;
};
io_Unfolder_Set::io_Unfolder_Set( 
        const char* tag,
        const char* file_name,
        int _n_iter)  : n_iter { _n_iter }
{

    ioGetter got{};
    ruu    =   (RooUnfoldResponse*) got( file_name, Form("%s",tag  ) );
    ruu_A  =   (RooUnfoldResponse*) got( file_name, Form("%s_A",tag) );
    ruu_B  =   (RooUnfoldResponse*) got( file_name, Form("%s_B",tag) );
    ruu_TS =   (RooUnfoldResponse*) got( file_name, Form("TS_%s",tag) );
    ruu_TU =   (RooUnfoldResponse*) got( file_name, Form("TU_%s",tag) );
    ruu_HC50 = (RooUnfoldResponse*) got( file_name, Form("HC50_%s",tag) );

    // do the closure test
    truth_A      = (TH1D*) ruu_A->Htruth();
    TH1D* meas_A = (TH1D*) ruu_A->Hmeasured();
    closure_A    = (TH1D*) io_BayesUnfold(meas_A, ruu_B,4);
    if (scale_by_binW) {
        io_scaleByBinWidth(truth_A);
        io_scaleByBinWidth(meas_A);
        io_scaleByBinWidth(closure_A);
    }
    rat_A        = ioDivideTH1(truth_A,closure_A);

};
TH1D* io_Unfolder_Set::unfold(TH1D* hg) {
    hg_base     = io_BayesUnfold(hg, ruu, n_iter);
    if (scale_by_binW) io_scaleByBinWidth(hg_base);
    pts_base    = hg_base;
    pts_rat_A   = rat_A; 
    cout << " pts_rat_A " << pts_rat_A << endl;
    pts_closure_A = rat_A; 
    pts_closure_A *= pts_base;
    hg_base_A   = (TH1D*) hg_base->Clone(ioUniqueName());
    for (int i{1}; i<=hg_base_A->GetXaxis()->GetNbins(); ++i) {
        hg_base_A->SetBinContent(i,pts_closure_A[i-1]);
    }

    pts_IPm2  = unf_unc(hg, hg_IPm2, rat_IPm2, ruu, n_iter-2);
    pts_IPp2  = unf_unc(hg, hg_IPp2, rat_IPp2, ruu, n_iter+2);
    pts_TS    = unf_unc(hg, hg_TS, rat_TS, ruu_TS,    n_iter);
    pts_TU    = unf_unc(hg, hg_TU, rat_TU, ruu_TU,    n_iter);
    pts_HC50  = unf_unc(hg, hg_HC50, rat_HC50, ruu_HC50,  n_iter);
    pts_sumsq = io_calc_quadrature( {pts_IPm2, pts_IPp2, pts_TS, pts_TU, pts_HC50, pts_closure_A}, hg_base );
    for (int i=0; i<pts_IPm2.size; ++i) {
        cout << Form(" x[%4.1f] entry %2i Closure: %5.5f IPm2: %5.5f IPp2: %5.5f TS: %5.5f TU: %5.5f HC50: %5.5f ",
                hg_base->GetXaxis()->GetBinCenter(i),i,pts_closure_A[i], pts_IPm2[i], pts_IPp2[i], pts_TS[i], pts_TU[i], pts_HC50[i]) << endl;
    };
    return hg_base;
};

io_Unfolder_Set_Ratio::io_Unfolder_Set_Ratio(
        io_Unfolder_Set& nSet,
        io_Unfolder_Set& dSet ) {

    hg_base  = ioDivideTH1(nSet.hg_base, dSet.hg_base);
    pts_base = hg_base;

    hg_IPm2  = ioDivideTH1(nSet.hg_IPm2, dSet.hg_IPm2);
    pts_IPm2 = hg_IPm2;

    hg_IPp2  = ioDivideTH1(nSet.hg_IPp2, dSet.hg_IPp2);
    pts_IPp2 = hg_IPp2;

    hg_TS  = ioDivideTH1(nSet.hg_TS, dSet.hg_TS);
    pts_TS = hg_TS;

    hg_TU  = ioDivideTH1(nSet.hg_TU, dSet.hg_TU);
    pts_TU = hg_TU;

    hg_HC50  = ioDivideTH1(nSet.hg_HC50, dSet.hg_HC50);
    pts_HC50 = hg_HC50;

    pts_base_A  = hg_base;
    cout << " RAT pts_rat_A " << nSet.pts_rat_A << endl;
    pts_base_A *= nSet.pts_rat_A;
    hg_base_A = (TH1D*) hg_base->Clone(ioUniqueName());
    for (int i=0; i<hg_base_A->GetXaxis()->GetNbins(); ++i) {
        hg_base_A->SetBinContent(i+1, pts_base_A[i]);
    }

    pts_sumsq = io_calc_quadrature( {pts_IPm2, pts_IPp2, pts_TS, 
            pts_TU, pts_HC50, pts_base_A}, pts_base );
};
        
