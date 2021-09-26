#include "ioJetMatcher.h"
#include <iostream>
#include "io_fnc.h"
#include "TH1D.h"
#include "TH2D.h"

double* ioEdges_pAuJet_prelim_13bins(){
    return  new double[14]{0.,1.,2.,3.,4.,5.,6.,8.,10.,15.,20.,25.,36.,50.};
};

// jet_match_floats
ioJetMatcher_float::ioJetMatcher_float(float _eta, float _phi, float _pT) :
    eta{_eta}, phi{_phi}, pT{_pT}, is_matched{false} {};
ioJetMatcher_float::ioJetMatcher_float() :
    eta{0.}, phi{0.}, pT{0.}, is_matched{false} {};

double ioJetMatcher_float::operator()(const ioJetMatcher_float& B, const float jetR2) {
    const float deta {eta-B.eta};
    const float dphi {io_dphi(phi, B.phi)};
    const float dist2 = deta*deta + dphi*dphi;
    if (dist2 == 0.) return 0.0001;
    else if (dist2 > jetR2) return 0.;
    else return dist2;
};

bool operator==(const ioJetMatcher_float& L, const ioJetMatcher_float R) 
    { return L.pT == R.pT; };

bool operator<(const ioJetMatcher_float& L, const ioJetMatcher_float R) 
    { return L.pT < R.pT; };

bool operator>(const ioJetMatcher_float& L, const ioJetMatcher_float R) 
    { return L.pT > R.pT; };

ioJetMatcher::ioJetMatcher (const char* _name, ioXsec& _Xsec, 
            const char* bin_file, const char* tag_M, const char* tag_T, 
            ioOptMap options, ioBinVec _hg2ptbins, ioOptMap dict
    ) :
    name {_name}, 
    Xsec{_Xsec}, 
    response{ioMakeRooUnfoldResponse(_name,bin_file,tag_M,tag_T)},
    response_noweight{ioMakeRooUnfoldResponse(Form("%s_noweight",_name),bin_file,tag_M,tag_T)},
    response_cut{ioMakeRooUnfoldResponse(Form("%s_cut",_name),bin_file,tag_M,tag_T)},
    hg2_cut{Form("hg2_cut_%s",_name),"title;measured;truth",80,-10.,70.,80,-10.,70.},
    _rand{},
    in_reco_bounds {bin_file,tag_M},
    in_meas_bounds {bin_file,tag_T},
    hg_pthb_cnt { Form("pthg_cnt_%s",_name),"Counter;#hat{#it{p}}_{T}-bin;N_{events}",
        Xsec.nbins_pthat, -0.5, Xsec.nbins_pthat-0.5 }
{
    dict += options;
    bool debug = dict.has("debug");

    if (dict.has("JES_JER_file") && dict.has("JER_limit") && dict.has("match_fake_limit")) {
        // set a limit of measured jets for matched true pt
        auto pt_true = ioReadValVec(dict["JES_JER_file"],dict["pt_true"].c_str());
        auto JES =     ioReadValVec(dict["JES_JER_file"],"JES");
        auto JER =     ioReadValVec(dict["JES_JER_file"],"JER");

        auto match_fake_limit = dict["match_fake_limit"]();
        auto JER_limit  = dict["JER_limit"]();

        vector<double> pt_meas_limit;
        for (unsigned int i{0}; i<JES.size(); ++i) {
            pt_meas_limit.push_back(pt_true[i]+JES[i]+JER_limit*JER[i]+match_fake_limit);
            if (debug) cout << Form(
             " debug : Truth %4.1f  JES %5.2f JER %4.2f nJERlim %3.1f fakept_lim %4.2f -> M_lim %4.2f",
             pt_true[i], JES[i], JER[i], JER_limit, match_fake_limit, 
                            (pt_true[i]+JES[i]+JER_limit*JER[i]+match_fake_limit)) << endl;
            
        }
        out_of_match_bounds = {pt_true,pt_meas_limit};
    } 
    fake_limit = dict("fake_limit",in_meas_bounds.hi_bound);
    if (debug) {
        cout << " debug: fake_limit " << fake_limit << endl;
        cout << " out_of_match_bounds.size " << out_of_match_bounds.size << endl;
    }

    if (dict.has("pthb_Mlimit_file") && dict.has("JER_limit") 
            && dict.has("match_fake_limit")) {
        b_ptht_Mlimit = true;
        auto pt_true = ioReadValVec(dict["pthb_Mlimit_file"],"pt_jet_bin_upbound");
        auto match_fake_limit = dict["match_fake_limit"]();
        auto JER_limit  = dict["JER_limit"]();
        for (int ibin{0}; ibin<Xsec.nbins_pthat; ++ibin) {
            vector<double> pt_meas_limit;
            auto JES = ioReadValVec(dict["pthb_Mlimit_file"],Form("JES__%i",ibin));
            auto JER = ioReadValVec(dict["pthb_Mlimit_file"],Form("JER__%i",ibin));
            for (unsigned int i{0}; i<JES.size(); ++i) {
                pt_meas_limit.push_back(pt_true[i]+JES[i]+JER_limit*JER[i]+match_fake_limit);
                if (debug) cout << Form(
                 " debug : bin %i Truth %4.1f  JES %5.2f JER %4.2f nJERlim %3.1f fakept_lim %4.2f -> M_lim %4.2f",
                 ibin,pt_true[i], JES[i], JER[i], JER_limit, match_fake_limit, 
                                (pt_true[i]+JES[i]+JER_limit*JER[i]+match_fake_limit)) << endl;
                
            }
            pthb_Mlimit[ibin] = {pt_true,pt_meas_limit};
        }
    }

    ratio_AtoB = dict["ratio_AtoB"]();

    ioBinVec binsM { bin_file, tag_M };
    pt_misses = binsM[0];

    ioBinVec binsT { bin_file, tag_T };
    /* pt_fakes = binsT[0]; */

    int nbins = Xsec.nbins_pthat;
    ioBinVec bins_M { bin_file, tag_M };
    ioBinVec bins_T { bin_file, tag_T };
    if (dict["Xsec_vs_miss"]==1) {
        b_Xsec_vs_miss = true;
        for (int i{0}; i<nbins; ++i) {
            v_miss.push_back({Form("miss_%i__%s",i,_name),
                Form("Misses for #hat{#it{p}}_{T}-bin %i;#it{p}_{T}-miss;N",i),
                bins_T, bins_T});
            A_miss.push_back({Form("A_miss_%i__%s",i,_name),
                Form("Misses for #hat{#it{p}}_{T}-bin %i;#it{p}_{T}-miss;N",i),
                bins_T, bins_T});
            B_miss.push_back({Form("B_miss_%i__%s",i,_name),
                Form("Misses for #hat{#it{p}}_{T}-bin %i;#it{p}_{T}-miss;N",i),
                bins_T, bins_T});
        }
    }
    if (dict["Xsec_vs_T"]==1) {
        b_Xsec_vs_T = true;
        for (int i{0}; i<nbins; ++i) {
            v_T.push_back({Form("T_%i__%s",i,_name),
                Form("Truth jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T}-true;N",i),
                bins_T, bins_T});
        }
    }
    if (dict["Xsec_vs_M"]==1) {
        b_Xsec_vs_M = true;
        for (int i{0}; i<nbins; ++i) {
            v_M.push_back({Form("M_%i__%s",i,_name),
                Form("Measured jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T};N",i),
                bins_M, bins_M});
        }
    }
    if (dict["Xsec_vs_fake"]==1) {
        b_Xsec_vs_fake = true;
        for (int i{0}; i<nbins; ++i) {
            v_fake.push_back({Form("fake_%i__%s",i,_name),
                Form("Fake jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T};N",i),
                bins_M, bins_M});
            A_fake.push_back({Form("A_fake_%i__%s",i,_name),
                Form("Fake jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T};N",i),
                bins_M, bins_M});
            B_fake.push_back({Form("B_fake_%i__%s",i,_name),
                Form("Fake jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T};N",i),
                bins_M, bins_M});
        }
    }
    if (dict["Xsec_vs_match"]==1) {
        b_Xsec_vs_match = true;
        for (int i{0}; i<nbins; ++i) {
            v_match.push_back({Form("match_%i__%s",i,_name),
                Form("Matched jets for #hat{#it{p}}_{T}-bin %i;"
                    "#it{p}_{T}-measured;#it{p}_{T}-truth;",i),
                bins_M, bins_M, bins_T, bins_T});
            A_match.push_back({Form("A_match_%i__%s",i,_name),
                Form("Matched jets for #hat{#it{p}}_{T}-bin %i;"
                    "#it{p}_{T}-measured;#it{p}_{T}-truth;",i),
                bins_M, bins_M, bins_T, bins_T});
            B_match.push_back({Form("B_match_%i__%s",i,_name),
                Form("Matched jets for #hat{#it{p}}_{T}-bin %i;"
                    "#it{p}_{T}-measured;#it{p}_{T}-truth;",i),
                bins_M, bins_M, bins_T, bins_T});
        }
    }
    if (dict["hg1_R2_match"]==1) {
        b_hg1_R2_match = true;
        hg1_R2_match = new TH1D(
        Form("hg1_R2_match_%s",_name), 
        "Matched jets; #sqrt((#Delta#phi)^2+(#delta#eta)^2)",
        100, 0., 0.2 );
    }
    
    jet_R2 = dict["R"]()*dict["R"]();

    response_A = (RooUnfoldResponse*)
        response.Clone(Form("%s_A",response.GetName()));
    response_B = (RooUnfoldResponse*)
        response.Clone(Form("%s_B",response.GetName()));
    ioBinVec pthatbins { Xsec.pthatbins };
}

bool ioJetMatcher::do_matching(int pthatbin, double weight) {
    double W { Xsec.Xsec(pthatbin) };
    if (weight) W *= weight;
    /* vector<double> fakes; */ // Using fakes is incorrect here.
                                // leftover reconstructed jets are actual jets from 
                                // the embedded event
    vector<double> misses;
    vector<pair<double,double>> matches;
    vector<double> v_delta_R2;

    // Find fakes, misses, and matches first, so that if an outlier event
    // is found, it can be filled without filling any other event.
    ioXYbounder& Mlimits { b_ptht_Mlimit ? pthb_Mlimit[pthatbin] : out_of_match_bounds };

    bool out_bounds_cut {false};
    for (auto& MC : data_MC) {
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (reco.is_matched) continue;
            double delta_R2 { MC(reco,jet_R2) };
            if (delta_R2 != 0) {
                // cut the event if matched pair is out of bounds
                if (Mlimits(MC.pT,reco.pT)) {
                    out_bounds_cut = true;
                    response_cut.Fill(reco.pT, MC.pT,W);
                    hg2_cut.Fill(reco.pT, MC.pT);
                    break;
                } 
                found_match = true;
                reco.is_matched = true;
                matches.push_back({reco.pT, MC.pT});
                v_delta_R2.push_back(delta_R2);
                break;
            }
        }
        if (!found_match) { misses.push_back(MC.pT); }
    }
    /* for (auto& reco : data_reco) { */
    /*     if (!reco.is_matched) { */
    /*         if (reco.pT > fake_limit) { */
    /*             out_bounds_cut = true; */
    /*         } else { */
    /*             fakes.push_back(reco.pT); */
    /*         } */
    /*     } */
    /* } */
    if (out_bounds_cut) {
        data_MC.clear();
        data_reco.clear();
        return true;
    }

    // now that is it not an outlier event, fill the data
    hg_pthb_cnt.Fill((double)pthatbin);
    double pthat_val { (0.5)*(Xsec.pthatbins[pthatbin]+
                       Xsec.pthatbins[pthatbin+1]) };
    bool fillA = _rand.Uniform() < ratio_AtoB;

    int n{0};
    for (auto& match : matches) {
        response.Fill(match.first,match.second,W); 
        response_noweight.Fill(match.first,match.second); 
        // tag split
        if (fillA) response_A->Fill(match.first, match.second,W);
        else       response_B->Fill(match.first, match.second,W);


        if (b_Xsec_vs_match) {
            v_match[pthatbin].Fill(match.first, match.second);
            if (fillA) A_match[pthatbin].Fill(match.first, match.second);
            else       B_match[pthatbin].Fill(match.first, match.second);
        }
        if (b_Xsec_vs_T)     v_T[pthatbin].Fill(match.second);
        if (b_Xsec_vs_M)     v_M[pthatbin].Fill(match.first);

        if (b_hg1_R2_match) hg1_R2_match->Fill(v_delta_R2[n]);
        ++n;
    }
    for (auto& miss : misses) {
        response.Miss(miss,W);
        response_noweight.Miss(miss);
        // tag split
        if (fillA) response_A->Miss(miss, W);
        else       response_B->Miss(miss, W);
        if (b_Xsec_vs_miss) {
            v_miss[pthatbin].Fill(miss);
            if (fillA) A_miss[pthatbin].Fill(miss);
            else       B_miss[pthatbin].Fill(miss);
        }
        if (b_Xsec_vs_T)    v_T[pthatbin].Fill(miss);
    }
    /* for (auto& fake : fakes) { */
    /*     response.Fake(fake,W); */
    /*     response_noweight.Fake(fake); */

    /*     if (b_Xsec_vs_fake) { */
    /*         v_fake[pthatbin].Fill(fake); */
    /*         if (fillA) A_fake[pthatbin].Fill(fake); */
    /*         else       B_fake[pthatbin].Fill(fake); */
    /*     } */
    /*     if (b_Xsec_vs_M)    v_M[pthatbin].Fill(fake); */

    /*     if (fillA) response_A->Fake(fake, W); */
    /*     else       response_B->Fake(fake, W); */
    /* }; */
    data_MC.clear();
    data_reco.clear();
    return false;
};

void ioJetMatcher::write() {
    /* R2_distr_matches->Write(); */
    hg_pthb_cnt.Write();
    hg2_cut.Write();
    response.Write();
    response_noweight.Write();
    response_cut.Write();
    response_A->Write();
    response_B->Write();

    for (auto& h : v_T)      h.Write();
    for (auto& h : v_M)      h.Write();
    for (auto& h : v_fake)   h.Write();
    for (auto& h : v_miss)   h.Write();
    for (auto& h : v_match)  h.Write();

    for (auto& h : A_fake)   h.Write();
    for (auto& h : A_miss)   h.Write();
    for (auto& h : A_match)  h.Write();

    for (auto& h : B_fake)   h.Write();
    for (auto& h : B_miss)   h.Write();
    for (auto& h : B_match)  h.Write();

    if (b_hg1_R2_match)  hg1_R2_match->Write();
};

void ioJetMatcher::addjet_MC(float eta, float phi, float pT) {
    /* if (bounds_T(pT)) data_MC.push_back({eta,phi,pT}); */
    if (in_meas_bounds(pT)) data_MC.push_back({eta,phi,pT});
    /* else cout << " out of bounds! " << pT << "  " << in_meas_bounds << endl; */
};

void ioJetMatcher::addjet_reco(float eta, float phi, float pT) {
    /* if (bounds_M(pT)) data_reco.push_back({eta,phi,pT}); */
    if (in_reco_bounds(pT)) data_reco.push_back({eta,phi,pT});
    /* else cout << " out of bounds! " << pT << "  " << in_reco_bounds << endl; */
};

// ioJetMatcherArray ::
ioJetMatcherArray::ioJetMatcherArray (
       const char* _name,
       ioXsec& _Xsec, 
       const char* bin_file, 
        vector<string> bin_names,
        vector<string> bin_tags,
        const char* pthb_Mlimit_file,
        double _ratioAtoB,
        bool _debug
    ) :
    Xsec{_Xsec}, 
    v_names{bin_names},
    /* response{ioMakeRooUnfoldResponse(_name,bin_file,tag_M,tag_T)}, */
    /* in_reco_bounds {bin_file,tag_M}, */
    /* in_meas_bounds {bin_file,tag_T}, */
    hg_pthb_cnt { Form("pthg_cnt_%s",_name),"Counter;#hat{#it{p}}_{T}-bin;N_{events}",
        Xsec.nbins_pthat, -0.5, Xsec.nbins_pthat-0.5 },
    ratio_AtoB {_ratioAtoB},
    debug{_debug}
{
    for (int i=0;i<(int)bin_names.size();++i) {
        ioBinVec bins { bin_file, bin_tags[i].c_str() };
        TH2D form_2D {ioUniqueName(), ";Measured;Truth", bins, bins, bins, bins}; 
        TH2D form_1D {ioUniqueName(), ";Measured;Truth", bins, bins, bins, bins}; 

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
    auto pt_true = ioReadValVec("pthb_Mlimit_file","pt_jet_bin_upbound");
    for (int ibin{0}; ibin<Xsec.nbins_pthat; ++ibin) {
        vector<double> pt_meas_limit;
        auto JES = ioReadValVec(pthb_Mlimit_file,Form("JES__%i",ibin));
        auto JER = ioReadValVec(pthb_Mlimit_file,Form("JER__%i",ibin));
        for (unsigned int i{0}; i<JES.size(); ++i) {
            pt_meas_limit.push_back(pt_true[i]+JES[i]+5.*JER[i]+8.);
        }
        pthb_Mlimit[ibin] = {pt_true,pt_meas_limit};
    }
}

bool ioJetMatcherArray::do_matching(int pthatbin, double weight) {
    double W { Xsec.Xsec(pthatbin) };
    if (weight) W *= weight;
    /* vector<double> fakes; */ // Using fakes is incorrect here.
                                // leftover reconstructed jets are actual jets from 
                                // the embedded event
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
                if (Mlimits(MC.pT,reco.pT)) continue;
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
    double pthat_val { (0.5)*(Xsec.pthatbins[pthatbin]+
                       Xsec.pthatbins[pthatbin+1]) };
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

void ioJetMatcherArray::cull_add_array(array<TH2D*,9>& data) {
    for (int i{0};i<9;++i) {
        io_cullsmallbins(data[i],10.);
        data[i]->Scale(Xsec.Xsec(i));
    }
    for (int i{1};i<9;++i) data[0]->Add(data[1]);
};
void ioJetMatcherArray::cull_add_array(array<TH1D*,9>& data) {
    for (int i{0};i<9;++i) {
        io_cullsmallbins(data[i],10.);
        data[i]->Scale(Xsec.Xsec(i));
    }
    for (int i{1};i<9;++i) data[0]->Add(data[1]);
};

void write_response(TH2D* h2, TH1D* truth, const char* name, const char* posttag) {
    TH1D* truth_add = (TH1D*) h2->ProjectionY(ioUniqueName());
    truth->Add(truth_add);
    TH1D* meas = (TH1D*) h2->ProjectionX(ioUniqueName());
    RooUnfoldResponse* ruu = new RooUnfoldResponse(meas, truth, h2, Form("%s%s", name, posttag));
    ruu->Write();
};

void ioJetMatcherArray::write() {
    // cull data with bins < 10 entries
    for (auto& rep : v_response)   cull_add_array(rep);
    for (auto& rep : v_response_A) cull_add_array(rep);
    for (auto& rep : v_response_B) cull_add_array(rep);

    for (unsigned int i=0; i<v_response.size(); ++i) {
        write_response(v_response[i][0], v_truth[i][0],     v_names[i].c_str(), "");
        write_response(v_response_A[i][0], v_truth_A[i][0], v_names[i].c_str(), "_A");
        write_response(v_response_B[i][0], v_truth_B[i][0], v_names[i].c_str(), "_B");
    }
    hg_pthb_cnt.Write();
};

void ioJetMatcherArray::addjet_MC(float eta, float phi, float pT) {
    data_MC.push_back({eta,phi,pT});
};

void ioJetMatcherArray::addjet_reco(float eta, float phi, float pT) {
    data_reco.push_back({eta,phi,pT});
};

