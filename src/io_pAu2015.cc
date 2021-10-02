#include "io_pAu2015.h"

ioPtrDbl io_Unfolder_Set::unf_unc(TH1D* raw, TH1D*& _hg, 
        TH1D*& _rat, RooUnfoldResponse*& _ruu, int n_iter) {
    _hg =  io_BayesUnfold(raw, _ruu, n_iter);
    /* cout << " is _hg null? " << (_hg==nullptr) << endl; */
    /* cout << " name: " << _hg->GetName() << " " << _hg->Integral() << endl; */
    _rat =  ioDivideTH1(_hg, hg_IP4);
    ioPtrDbl pts { _rat };
    pts -= 1.;
    /* cout << " pts: " << pts << endl; */
    return pts;
};
io_Unfolder_Set::io_Unfolder_Set( 
        const char* tag,
        const char* file_name) 
{

    ioGetter got{};
    ruu    =   (RooUnfoldResponse*) got( file_name, Form("all_%s",tag  ) );
    ruu_A  =   (RooUnfoldResponse*) got( file_name, Form("all_%s_A",tag) );
    ruu_B  =   (RooUnfoldResponse*) got( file_name, Form("all_%s_B",tag) );
    ruu_TS =   (RooUnfoldResponse*) got( file_name, Form("TS_%s",tag) );
    ruu_TU =   (RooUnfoldResponse*) got( file_name, Form("TU_%s",tag) );
    ruu_HC50 = (RooUnfoldResponse*) got( file_name, Form("HC50_%s",tag) );

    // do the closure test
    truth_A      = (TH1D*) ruu_A->Htruth();
    TH1D* meas_A = (TH1D*) ruu_A->Hmeasured();
    closure_A    = (TH1D*) io_BayesUnfold(meas_A, ruu_B,4);
    rat_A        = ioDivideTH1(truth_A,closure_A);
    ptr_closure = rat_A;
    ptr_closure -= 1.;

};
TH1D* io_Unfolder_Set::unfold(TH1D* hg) {
    array<ioPtrDbl, 10> iter_array;
    for (int i=0; i<10; ++i) {
        iter_array[i] = io_BayesUnfold(hg,ruu,i);
    };
    int z = 0;
    for (auto ia : iter_array){
        cout << " tiger errors i_iter: " << z++ << " ";
        for (int k = 8; k<=10; ++k)
            cout << Form(" %10.8g",ia[k]);
        cout << endl;
    }
    hg_IP4 = io_BayesUnfold(hg, ruu, 60);
    ptr_IP2 = unf_unc(hg, hg_IP2, rat_IP2, ruu, 50);
    ptr_IP6 = unf_unc(hg, hg_IP6, rat_IP6, ruu, 70);
    ptr_TS = unf_unc(hg, hg_TS, rat_TS, ruu_TS, 50);
    ptr_TU = unf_unc(hg, hg_TU, rat_TU, ruu_TU, 50);
    ptr_HC50 = unf_unc(hg, hg_HC50, rat_HC50, ruu_HC50, 50);

    for (auto i : vector<int>{8,9,10}) {
        cout << Form(" entry %2i Closure: %5.5f IP2: %5.5f IP6: %5.5f TS: %5.5f TU: %5.5f HC50: %5.5f ",
                i,ptr_closure[i], ptr_IP2[i], ptr_IP6[i], ptr_TS[i], ptr_TU[i], ptr_HC50[i]) << endl;
    };
    return hg_IP4;
};
