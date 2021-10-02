#include "ioJetMatcher100.h"
#include <iostream>
#include "io_fnc.h"
#include "ioClass.h"
#include "TH1D.h"
#include "TH2D.h"

// ioJetMatcher100 ::
ioJetMatcher100::ioJetMatcher100 ( const char* file, const char* endtag ) {
    ioGetter got{};
    hg2_response   = (TH2D*) got(file, Form("%sresponse",endtag));
    hg1_truth      = (TH1D*) got(file, Form("%struth",endtag));
    hg1_measured   = (TH1D*) got(file, Form("%smeasured",endtag));
};
RooUnfoldResponse* ioJetMatcher100::make_ruu(const char* fname, const char* M_tag, 
        const char* T_tag, const char* name, bool write) 
{
    return make_ruu( ioBinVec{fname, M_tag}, ioBinVec{fname, T_tag}, name, write);
};
RooUnfoldResponse* ioJetMatcher100::make_ruu(ioBinVec bins_M, ioBinVec bins_T, const char* name, bool write) {
    TH2D* response = rebin(hg2_response, bins_M, bins_T);
    TH1D* truth    = rebin(hg1_truth,    bins_T);
    TH1D* measured = (TH1D*) response->ProjectionX(ioUniqueName());
    RooUnfoldResponse* ruu = new RooUnfoldResponse(measured, truth, response, name);
    if (write) ruu->Write();
    return ruu;
};
ioJetMatcher100::ioJetMatcher100 (
       ioXsec* _Xsec, 
       const char* _tag,
        double _ratioAtoB,
        int cull_n,
        double high_sig_cut,
        double high_sig_offset,
        bool _debug
    ) :
    Xsec{_Xsec}, 
    tag{_tag},
    hg_pthb_cnt { Form("%spthg_cnt",_tag),"Counter;#hat{#it{p}}_{T}-bin;N_{events}",
        Xsec->nbins_pthat, -0.5, Xsec->nbins_pthat-0.5 },
    ratio_AtoB {_ratioAtoB},
    cut_high_sigma { high_sig_cut },
    cut_high_sigma_offset { high_sig_offset },
    debug{_debug}
{
    ioBinVec _bins { Xsec->pthatbins };
    sigma = new TH1D(Form("%ssigma",tag.c_str()),"bins_weights", _bins, _bins);
    for (int i=0;i<9;++i) sigma->SetBinContent(i+1, Xsec->Xsec(i,1));
    ioBinVec bins { {0.,0.,100,100. }};
    for (int k{0}; k<9; ++k) {
        v_response[k]   = new TH2D( Form("%sR_%i",tag.c_str(),k),   ";Measured;Truth", bins,bins,bins, bins );
        v_response_A[k] = new TH2D( Form("%sR_%i_A",tag.c_str(),k), ";Measured;Truth", bins,bins,bins, bins);
        v_response_B[k] = new TH2D( Form("%sR_%i_B",tag.c_str(),k), ";Measured;Truth", bins,bins,bins, bins);

        v_truth[k]      = new TH1D( Form("%sT_%i",tag.c_str(),k), ";Truth;", bins, bins );
        v_truth_A[k]    = new TH1D( Form("%sT_%i_A",tag.c_str(),k), ";Truth;", bins, bins );
        v_truth_B[k]    = new TH1D( Form("%sT_%i_B",tag.c_str(),k), ";Truth;", bins, bins );

        v_response[k]->Sumw2();
        v_response_A[k]->Sumw2();
        v_response_B[k]->Sumw2();

        v_truth[k]->Sumw2();
        v_truth_A[k]->Sumw2();
        v_truth_B[k]->Sumw2();
    }

    jet_R2 = 0.16;
    ioBinVec pthatbins { Xsec->pthatbins };
}

bool ioJetMatcher100::do_matching(int pthatbin) {
    vector<double> misses;
    vector<pair<double,double>> matches;
    vector<double> v_delta_R2;
    /* ioXYbounder& Mlimits { pthb_Mlimit[pthatbin] }; */
    for (auto& MC : data_MC) {
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (reco.is_matched) continue;
            double delta_R2 { MC(reco,jet_R2) };
            if (delta_R2 != 0) {
                /* if (apply_Mlimit && Mlimits(MC.pT,reco.pT)) continue; */
                found_match = true;
                reco.is_matched = true;
                matches.push_back({reco.pT, MC.pT});
                break;
            }
        }
        if (!found_match) { misses.push_back(MC.pT); }
    }

    // now that is it not an outlier event, fill the data
    hg_pthb_cnt.Fill((double)pthatbin);
    bool fillA = _rand.Uniform() < ratio_AtoB;
    int k = pthatbin;

    for (auto& match : matches) {
        v_response[k]->Fill(match.first,match.second); 
        if (fillA)
            v_response_A[k]->Fill(match.first,match.second); 
        else
            v_response_B[k]->Fill(match.first,match.second); 
    }
    for (auto& miss : misses) {
        v_truth[k]->Fill(miss); 
        if (fillA)
            v_truth_A[k]->Fill(miss); 
        else
            v_truth_B[k]->Fill(miss); 
    }
    data_MC.clear();
    data_reco.clear();
    return false;
};

/* void ioJetMatcher100::cull_add_array(array<TH2D*,9>& data, string which, const char* tag) { */
/*     for (int i{0};i<9;++i) { */
/*         io_cullsmallbins(data[i],cull_n); */
/*         if (cut_high_sigma != 0.) io_cut_high_sigmaX(data[i],cut_high_sigma_offset); */
/*         if (write_9) { */
/*             data[i]->SetName(Form("res_%s_%s%s__n%i",name.c_str(),which.c_str(),tag,i)); */
/*             data[i]->Write(); */
/*         } */
/*         data[i]->Scale(Xsec->Xsec(i,(int)hg_pthb_cnt.GetBinContent(i+1))); */
/*     } */
/*     for (int i{1};i<9;++i) { */
/*         data[0]->Add(data[i]); */
/*     } */
/* }; */

void ioJetMatcher100::process_arrays(
        array<TH2D*,9>& arr_resp,  TH2D*& response, 
        array<TH1D*,9>& arr_truth, TH1D*& truth, const char* endtag) 
{
    for (int i{0};i<9;++i) {
        io_cullsmallbins(arr_resp[i], cull_n);
        io_cullsmallbins(arr_truth[i],cull_n);

        cout << " tag: " << tag << "  which " << endtag << " set: " << i << endl;
        cout << " cut_high_sigma " << cut_high_sigma << " " << cut_high_sigma_offset <<  endl;
        if (cut_high_sigma != 0) io_cut_high_sigmaX(arr_resp[i],cut_high_sigma, cut_high_sigma_offset);
        if (cut_high_sigma != 0) io_cut_high_sigmaX(arr_truth[i],cut_high_sigma, cut_high_sigma_offset);

        TH1D* _truth = (TH1D*)arr_resp[i]->ProjectionY(ioUniqueName());
        arr_truth[i]->Add(_truth);
        delete _truth;

        arr_resp [i]->Scale(sigma->GetBinContent(i+1)/hg_pthb_cnt.GetBinContent(i+1));
        arr_truth[i]->Scale(sigma->GetBinContent(i+1)/hg_pthb_cnt.GetBinContent(i+1));
    }
    response = (TH2D*) arr_resp[0] ->Clone(Form("%sresponse%s",tag.c_str(),endtag));
    truth    = (TH1D*) arr_truth[0]->Clone(Form("%struth%s",tag.c_str(),endtag));

    /* response ->Sumw2(); */
    /* truth    ->Sumw2(); */
    
    for (int i{1};i<9;++i) {
        response->Add(arr_resp[i], 1.);
        truth   ->Add(arr_truth[i],1.);
    }
    response->Write();
    truth->Write();
    for (int i{0};i<9;++i) {
        arr_resp [i]->Write();
        arr_truth[i]->Write();
    }
    // write the roounfold 
    TH1D* measured = (TH1D*) response->ProjectionX(Form("%smeasured%s",tag.c_str(),endtag));
    RooUnfoldResponse* ruu = new RooUnfoldResponse( measured, truth, response, Form("%sruu%s",tag.c_str(),endtag) );
    ruu->Write();
};
/* void ioJetMatcher100::write_response(TH2D* h2, TH1D* truth, string which, const char* posttag) { */
/*     TH1D* truth_add = (TH1D*) h2->ProjectionY(ioUniqueName()); */
/*     truth->Add(truth_add); */
/*     truth->SetName(Form("%s_Truth_%s%s",name.c_str(),which.c_str(),posttag)); */
/*     TH1D* meas = (TH1D*) h2->ProjectionX(ioUniqueName()); */
/*     meas->SetName(Form("%s_Measured_%s%s",name.c_str(),which.c_str(),posttag)); */
/*     h2->SetName(Form("%s_Matched_%s%s",name.c_str(),which.c_str(),posttag)); */
/*     RooUnfoldResponse* ruu = new RooUnfoldResponse(meas, truth, */ 
/*             h2, Form("%s_%s%s", name.c_str(), which.c_str(), posttag)); */
/*     ruu->Write(); */
/* }; */

TH2D* ioJetMatcher100::rebin(TH2D* hg2, ioBinVec bins_M, ioBinVec bins_T) {
    TH2D* r_hg = new TH2D(ioUniqueName(), ";Measured;Truth", bins_M, bins_M, bins_T, bins_T);

    TAxis* ax_X = hg2->GetXaxis();
    TAxis* ax_Y = hg2->GetYaxis();

    for (int y=0; y<bins_T; ++y) {
        int y0 = ax_Y->FindBin(bins_T[y]+0.5);
        int y1 = ax_Y->FindBin(bins_T[y+1]-0.5);
        TH1D* proj = (TH1D*) hg2->ProjectionY(ioUniqueName(),y0,y1);
        TH1D* rebin = (TH1D*) hg2->Rebin(bins_M,ioUniqueName(),bins_M);
        rebin->SetBinContent(0,0.);
        rebin->SetBinContent(bins_T+1,0.);
        rebin->SetBinError(0,0.);
        rebin->SetBinError(bins_T+1,0.);
        for (int x=0;x<bins_M;++x) {
            r_hg->SetBinContent(x+1,y+1,rebin->GetBinContent(x+1));
        }
    };
    return r_hg;
};
TH1D* ioJetMatcher100::rebin(TH1D* hg1, ioBinVec bins_T) {
    /* TH1D* r_hg = new TH1D(name, Form("%s;Measured;",name), bins_T, bins_T); */
    TH1D* rebin = (TH1D*) hg1->Rebin(bins_T,ioUniqueName(),bins_T);
    rebin->SetBinContent(0,0.);
    rebin->SetBinContent(bins_T+1,0.);
    rebin->SetBinError(0,0.);
    rebin->SetBinError(bins_T+1,0.);
    return rebin;
};

void ioJetMatcher100::write() {
    process_arrays(v_response, hg2_response, v_truth, hg1_truth, "");
    hg1_measured = (TH1D*) hg2_response->ProjectionX(Form("%smeasured",tag.c_str()));
    hg1_measured->Write();
    process_arrays(v_response_A, hg2_response_A, v_truth_A, hg1_truth_A, "_A");
    process_arrays(v_response_B, hg2_response_B, v_truth_B, hg1_truth_B, "_B");
    hg_pthb_cnt.Write();
    sigma->Write();
};

void ioJetMatcher100::addjet_MC(float eta, float phi, float pT) {
    data_MC.push_back({eta,phi,pT});
};

void ioJetMatcher100::addjet_reco(float eta, float phi, float pT) {
    data_reco.push_back({eta,phi,pT});
};

