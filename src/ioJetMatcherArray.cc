#include "ioJetMatcherArray.h"
#include <iostream>
#include "io_fnc.h"
#include "TH1D.h"
#include "TH2D.h"

// ioJetMatcherArray ::
ioJetMatcherArray::ioJetMatcherArray (
       const char* _name,
       ioXsec& _Xsec, 
       const char* bin_file, 
        vector<string> bin_names,
        vector<string> bin_tags_M,
        vector<string> bin_tags_T,
        const char* pthb_Mlimit_file,
        double _ratioAtoB,
        bool _debug
    ) :
    name{_name},
    Xsec{_Xsec}, 
    v_names{bin_names},
    hg_pthb_cnt { Form("pthg_cnt_%s",_name),"Counter;#hat{#it{p}}_{T}-bin;N_{events}",
        Xsec.nbins_pthat, -0.5, Xsec.nbins_pthat-0.5 },
    ratio_AtoB {_ratioAtoB},
    debug{_debug}
{
    for (int i=0;i<(int)bin_names.size();++i) {
        ioBinVec bins_M { bin_file, bin_tags_M[i].c_str() };
        ioBinVec bins_T { bin_file, bin_tags_T[i].c_str() };
        TH2D form_2D {ioUniqueName(), ";Measured;Truth", bins_M, bins_M, bins_T, bins_T}; 
        TH1D form_1D {ioUniqueName(), ";Measured;Truth", bins_T, bins_T };

        v_response.push_back({});
        v_response_A.push_back({});
        v_response_B.push_back({});

        v_truth.push_back({});
        v_truth_A.push_back({});
        v_truth_B.push_back({});

        for (int k{0}; k<9; ++k) {
            v_response[i][k] = (TH2D*) form_2D.Clone(ioUniqueName());
            v_response_A[i][k] = (TH2D*) form_2D.Clone(ioUniqueName());
            v_response_B[i][k] = (TH2D*) form_2D.Clone(ioUniqueName());

            v_truth[i][k] = (TH1D*) form_1D.Clone(ioUniqueName());
            v_truth_A[i][k] = (TH1D*) form_1D.Clone(ioUniqueName());
            v_truth_B[i][k] = (TH1D*) form_1D.Clone(ioUniqueName());
        }
    }

    jet_R2 = 0.16;
    ioBinVec pthatbins { Xsec.pthatbins };
    auto pt_true = ioReadValVec(pthb_Mlimit_file,"pt_jet_bin_upbound");
    for (int ibin{0}; ibin<Xsec.nbins_pthat; ++ibin) {
        vector<double> pt_meas_limit;
        auto JES = ioReadValVec(pthb_Mlimit_file,Form("JES__%i",ibin));
        auto JER = ioReadValVec(pthb_Mlimit_file,Form("JER__%i",ibin));
        for (unsigned int i{0}; i<JES.size(); ++i) {
            pt_meas_limit.push_back(pt_true[i]+JES[i]+5.*JER[i]*8.);
        }
        pthb_Mlimit[ibin] = {pt_true,pt_meas_limit};
    }
}

bool ioJetMatcherArray::do_matching(int pthatbin) {
    vector<double> misses;
    vector<pair<double,double>> matches;
    vector<double> v_delta_R2;
    ioXYbounder& Mlimits { pthb_Mlimit[pthatbin] };
    for (auto& MC : data_MC) {
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (reco.is_matched) continue;
            double delta_R2 { MC(reco,jet_R2) };
            if (delta_R2 != 0) {
                if (apply_Mlimit && Mlimits(MC.pT,reco.pT)) continue;
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
        for (auto& rep : v_response) rep[k]->Fill(match.first,match.second); 
        if (fillA)
            for (auto& rep : v_response_A) rep[k]->Fill(match.first,match.second); 
        else
            for (auto& rep : v_response_B) rep[k]->Fill(match.first,match.second); 
    }
    for (auto& miss : misses) {
        for (auto& rep : v_truth) rep[k]->Fill(miss); 
        if (fillA)
            for (auto& rep : v_truth_A) rep[k]->Fill(miss); 
        else
            for (auto& rep : v_truth_B) rep[k]->Fill(miss); 
    }
    data_MC.clear();
    data_reco.clear();
    return false;
};

void ioJetMatcherArray::cull_add_array(array<TH2D*,9>& data, string which, const char* tag) {
    for (int i{0};i<9;++i) {
        io_cullsmallbins(data[i],cull_n);
        if (cut_high_sigma != 0.) io_cut_high_sigmaX(data[i],cut_high_sigma_offset);
        if (write_9) {
            data[i]->SetName(Form("res_%s_%s%s__n%i",name.c_str(),which.c_str(),tag,i));
            data[i]->Write();
        }
        data[i]->Scale(Xsec.Xsec(i,(int)hg_pthb_cnt.GetBinContent(i+1)));
    }
    for (int i{1};i<9;++i) {
        data[0]->Add(data[i]);
    }
};
void ioJetMatcherArray::cull_add_array(array<TH1D*,9>& data, string which, const char* tag) {
    for (int i{0};i<9;++i) {
        io_cullsmallbins(data[i],cull_n);
        if (write_9) {
            data[i]->SetName(Form("miss_%s_%s%s__n%i",name.c_str(),which.c_str(),tag,i));
            data[i]->Write();
        }
        data[i]->Scale(Xsec.Xsec(i,(int)hg_pthb_cnt.GetBinContent(i+1)));
    }
    for (int i{1};i<9;++i) {
        data[0]->Add(data[i]);
    }
};

void ioJetMatcherArray::write_response(TH2D* h2, TH1D* truth, string which, const char* posttag) {
    TH1D* truth_add = (TH1D*) h2->ProjectionY(ioUniqueName());
    truth->Add(truth_add);
    truth->SetName(Form("%s_Truth_%s%s",name.c_str(),which.c_str(),posttag));
    TH1D* meas = (TH1D*) h2->ProjectionX(ioUniqueName());
    meas->SetName(Form("%s_Measured_%s%s",name.c_str(),which.c_str(),posttag));
    h2->SetName(Form("%s_Matched_%s%s",name.c_str(),which.c_str(),posttag));
    RooUnfoldResponse* ruu = new RooUnfoldResponse(meas, truth, 
            h2, Form("%s_%s%s", name.c_str(), which.c_str(), posttag));
    ruu->Write();
};

void ioJetMatcherArray::write() {
    // cull data with bins < 10 entries
    for (unsigned int i=0; i<v_response.size(); ++i) {
        cull_add_array(v_response[i],    v_names[i],"");
        cull_add_array(v_response_A[i],  v_names[i],"_A");
        cull_add_array(v_response_B[i],  v_names[i],"_B");

        cull_add_array(v_truth[i],    v_names[i],"");
        cull_add_array(v_truth_A[i],  v_names[i],"_A");
        cull_add_array(v_truth_B[i],  v_names[i],"_B");

        write_response(v_response[i][0], v_truth[i][0],     v_names[i], "");
        write_response(v_response_A[i][0], v_truth_A[i][0], v_names[i], "_A");
        write_response(v_response_B[i][0], v_truth_B[i][0], v_names[i], "_B");
    }
    hg_pthb_cnt.Write();
};

void ioJetMatcherArray::addjet_MC(float eta, float phi, float pT) {
    data_MC.push_back({eta,phi,pT});
};

void ioJetMatcherArray::addjet_reco(float eta, float phi, float pT) {
    data_reco.push_back({eta,phi,pT});
};

