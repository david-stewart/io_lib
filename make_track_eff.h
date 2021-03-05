#include "TF1.h"
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TBranch.h"
#include "TString.h"
#include "src/good_types.h"
#include "src/apply_prior.h"
#include "src/fn_abund.h"
#include "src/divide_TH1D.h"
#include "RooUnfold/RooUnfoldResponse.h"
#include "RooUnfold/RooUnfoldBayes.h"
#include "RooUnfold/RooUnfoldBayes.h"


void make_track_eff(
        int ibin_bbc0, int ibin_bbc1, 
        int ibin_eta0, int ibin_eta1, 
        const char* input_file="out/response_THnSparse.root",
        const char* weighting="dAu",
        string suffix="", // use '*' for no suffix, 
                          // '' for auto-generate, anything else for specific
        vector<string> p_types= vector<string>())
{
// MT THnSparse:
//      0  AB_bins (for A and B randomly separated populations)
//      1  eta-bins
//      2  bbc-bins
//      3  pt-Truth
//      4  pt-Measured
// T  THnSparse:
//      0  AB_bins (for A and B randomly separated populations)
//      1  eta-bins
//      2  bbc-bins
//      3  pt-Truth

    // within the local TFile, write the following:
    // TH2D*  pT {M,T}
    // TH1D*  M
    // TH1D*  T
    // TH1D*  single_bin_eff (i.e. ratio of M to T)
    // TH1D*  unf_eff        (ratio of M to unfolded)

    // Input:
    //  0: bbc bin to start projection in THnSparse input
    //  1: bbc bin to end   projection in THnSparse input
    //  2-3: same as 0-1 but for eta
    //  4: input file with the THnSparse
    //  5: which set to use for the weighting of truth abundancies ("pp" or "dAu")
    //     if "" then don't weight the truth pT spectra
    //  6: vector of particle types to use. If blank, will use all six types.

    
    // Save the current output file
    TFile* s_current = gDirectory->GetFile();
    TH2D *MT_cum;
    TH2D *MT_cumA;
    TH2D *MT_cumB;
    TH1D* T_cum;
    TH1D* T_cumA;
    TH1D* T_cumB;

    // initialize p_types input
    if (p_types.size()==0) {
        p_types.push_back("p");
        p_types.push_back("pbar");
        p_types.push_back("pi");
        p_types.push_back("antipi");
        p_types.push_back("K");
        p_types.push_back("antiK");
    };
    for (unsigned int i=0;i<p_types.size(); ++i) {
        if (!is_good_type(p_types[i].c_str())) {
            cout << "Fatal: input type " << p_types[i] << " is not a good type." << endl;
            cout << "Must be one of p, pbar, pi, antipi, K, antiK" << endl;
            exit(2);
        }
    }

    // Make the output suffix (for histogram names())
    if (suffix == "*") {
        suffix = "";
    } else if (suffix == "") {
        suffix = Form("_%i_%i_bbc__%i_%i_eta",ibin_bbc0, ibin_bbc1, ibin_eta0, ibin_eta1);
    }

    // Get and weight the MT and T matrices and add them together
    TFile* f_input = new TFile(input_file,"read");
    f_input->cd();

    for (unsigned int i=0;i<p_types.size(); ++i) {
        THnSparse* sparse_MT;
        THnSparse* sparse_T;
        f_input->GetObject(Form("sparse_MT_%s", p_types[i].c_str()), sparse_MT);
        f_input->GetObject(Form("sparse_T_%s",  p_types[i].c_str()), sparse_T);

        sparse_MT->GetAxis(1)->SetRange(ibin_eta0, ibin_eta1);
        sparse_T ->GetAxis(1)->SetRange(ibin_eta0, ibin_eta1);

        sparse_MT->GetAxis(2)->SetRange(ibin_bbc0, ibin_bbc1);
        sparse_T ->GetAxis(2)->SetRange(ibin_bbc0, ibin_bbc1);

        sparse_MT->GetAxis(0)->SetRange(1,2);
        sparse_T ->GetAxis(0)->SetRange(1,2);

        TH2D* MT = sparse_MT->Projection(3,4);
        MT->SetName(Form("MT%s%s",suffix.c_str(),p_types[i].c_str()));
        TH1D* T  = sparse_T ->Projection(3);
        T->SetName(Form("T%s%s",suffix.c_str(),p_types[i].c_str()));
        /* cout << " entries: " << T->GetEntries() << endl; */

        //-------get the projections from the A part only to T_ and M_
        sparse_MT->GetAxis(0)->SetRange(1,1);
        sparse_T ->GetAxis(0)->SetRange(1,1);
        TH2D* MT_A = sparse_MT->Projection(3,4);
        MT_A->SetName(Form("MT_A%s%s",suffix.c_str(),p_types[i].c_str()));
        TH1D* T_A  = sparse_T ->Projection(3);
        T_A->SetName(Form("T_A%s%s",suffix.c_str(),p_types[i].c_str()));
        
        //-------get the projections from the B part only for T, M, and MT
        sparse_MT->GetAxis(0)->SetRange(2,2);
        sparse_T ->GetAxis(0)->SetRange(2,2);
        TH2D* MT_B = sparse_MT->Projection(3,4);
        MT_B->SetName(Form("MT_B%s%s",suffix.c_str(),p_types[i].c_str()));
        TH1D* T_B  = sparse_T ->Projection(3);
        T_B->SetName(Form("T_B%s%s",suffix.c_str(),p_types[i].c_str()));

        if ( (!strcmp(weighting,"pp")) || (!strcmp(weighting,"dAu")) ) {
           apply_prior(fn_abund(weighting,  p_types[i]), MT,   T, true);
           apply_prior(fn_abund(weighting,  p_types[i]), MT_A, T_A, true);
           apply_prior(fn_abund(weighting,  p_types[i]), MT_B, T_B, true);
        } else if ( strcmp(weighting,"") ) {
            cout << " Fatal error: input weighting " << weighting << " must be blank, pp, or dAu" << endl;
            exit(2);
        }

        if (i==0) {
            s_current->cd();
            MT_cum = (TH2D*) MT->Clone(Form("MT%s",suffix.c_str()));
            T_cum  = (TH1D*) T->Clone( Form("T%s",suffix.c_str()));
            /* M_cum  = (TH1D*) M->Clone( Form("M%s",suffix.c_str())); */

            MT_cumA = (TH2D*) MT_A->Clone( Form("MT_A%s",suffix.c_str()));
            T_cumA  = (TH1D*) T_A ->Clone( Form("T_A%s", suffix.c_str()));
            /* M_cumA  = (TH1D*) M_A ->Clone( Form("M_A%s", suffix.c_str())); */

            MT_cumB = (TH2D*) MT_B->Clone( Form("MT_B%s",suffix.c_str()));
            T_cumB  = (TH1D*) T_B ->Clone( Form("T_B%s", suffix.c_str()));
            /* M_cumB  = (TH1D*) M_B ->Clone( Form("M_B%s", suffix.c_str())); */

            double eta_lo = sparse_MT->GetAxis(1)->GetBinLowEdge(ibin_eta0);
            double eta_hi = sparse_MT->GetAxis(1)->GetBinUpEdge(ibin_eta1);

            double bbc_lo = sparse_MT->GetAxis(2)->GetBinLowEdge(ibin_bbc0);
            double bbc_hi = sparse_MT->GetAxis(2)->GetBinUpEdge(ibin_bbc1);
            cout << " ibin_bbc0 " << ibin_bbc0 << "   " << bbc_lo << endl;
            cout << " ibin_bbc1 " << ibin_bbc1 << "   " << bbc_hi << endl;
            /* for (int k=1;k<11;++k) cout << " " << k << " " << sparse_MT->GetAxis(2)->GetBinLowEdge(k) << endl; */

            MT_cum->SetTitle(Form("Events cut for #eta#in[%.1f,%.1f) bbc#in[%.0f,%.0f)",eta_lo, eta_hi, bbc_lo,bbc_hi));
            T_cum->SetTitle(Form("Events cut for #eta#in[%.1f,%.1f) bbc#in[%.0f,%.0f)",eta_lo, eta_hi, bbc_lo,bbc_hi));
            /* M_cum->SetTitle(Form("Events cut for #eta#in[%.1f,%.1f) bbc#in[%.0f,%.0f)",eta_lo, eta_hi, bbc_lo,bbc_hi)); */
            f_input->cd();
        } else {
            MT_cum->Add(MT);
            T_cum ->Add(T);

            MT_cumA->Add(MT_A);
            T_cumA ->Add(T_A);

            MT_cumB->Add(MT_B);
            T_cumB ->Add(T_B);
            /* M_cum ->Add(M); */
        }
    }
    // now that MT{,_A,_B} are weighted, project to get M
    TH1D* M_cum  = MT_cum ->ProjectionX(Form("M%s",suffix.c_str()));
    TH1D* M_cumA = MT_cumA->ProjectionX(Form("M_A%s",suffix.c_str()));
    TH1D* M_cumB = MT_cumB->ProjectionX(Form("M_B%s",suffix.c_str()));

    s_current->cd();
    MT_cum->Write(); 
    T_cum->Write(); 
    M_cum->Write(); 

    // get the single bin ratio
    TH1D* ratio = divide_TH1D(M_cum, T_cum, Form("eff_s_bin%s",suffix.c_str()));
    ratio->Write();

    // unfold M to T (should, of course, come close to T given that T & M were used to make the MT matrix)
    // this is not super-useful, but it at least shows that the unfolding is doing what it ought to
    f_input->cd();
    RooUnfoldResponse *rooUnfRes = new RooUnfoldResponse (M_cum,T_cum,MT_cum,Form("unfold_response%s",suffix.c_str()));
    RooUnfoldBayes* bayes = new RooUnfoldBayes(rooUnfRes,M_cum,3);

    cout << " a0 " << endl;
    cout << " suffix: " << suffix << endl;
    TH1D* unfolded = (TH1D*) bayes->Hreco();
    unfolded->SetName(Form("unfoldedAB%s",suffix.c_str()));
    delete rooUnfRes;
    delete bayes;

    s_current->cd();
    TH1D* unf_ratio = divide_TH1D(M_cum, unfolded, Form("eff_unf%s",suffix.c_str()));
    unf_ratio->SetTitle("Ratio of M to T where T is M unf. w/TM");
    unf_ratio->Write();

    // make the closure test (get MT_B and unfold M_A to get unf_M_A and compare to T_A
    f_input->cd();
    RooUnfoldResponse *rooUnfResA = new RooUnfoldResponse (M_cumB,T_cumB,MT_cumB,Form("unfold_response%sB",suffix.c_str()));
    RooUnfoldBayes* bayesA = new RooUnfoldBayes(rooUnfResA,M_cumA,3);
    cout << " a1 " << endl;
    TH1D* unfoldedA = (TH1D*) bayesA->Hreco();
    unfoldedA->SetName(Form("unfoldedA%s",suffix.c_str()));
    delete rooUnfResA;
    delete bayesA;
    s_current->cd();
    TH1D* closure = divide_TH1D(unfoldedA, T_cumA, Form("closure%s",suffix.c_str()));
    closure->SetTitle(Form("Ratio (meausured-A unfolded with MT-B) vs (truth-A) %s",suffix.c_str()));
    closure->Write();
    
    f_input->Close();
    delete f_input;
}
