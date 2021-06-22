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
    in_meas_bounds {bin_file,tag_T}
{
    dict += options;
    fake_limit = dict("fake_limit",0.);
    if (dict.has("match_bounds_file")) {
        out_of_match_bounds = 
        {dict["match_bounds_file"], dict["pt_true"], dict["pt_measured"] };
    }
    ratio_AtoB = dict["ratio_AtoB"]();

    ioBinVec binsM { bin_file, tag_M };
    pt_misses = binsM[0];

    ioBinVec binsT { bin_file, tag_T };
    pt_fakes = binsT[0];

    b_make_AB       = (dict["make_AB"]==1);
    b_hg2_Xsec_vs_M = (dict["hg2_Xsec_vs_M"]==1);
    b_hg2_Xsec_vs_T = (dict["hg2_Xsec_vs_T"]==1);
    b_hg2_Xsec_vs_match = (dict["hg2_Xsec_vs_match"]==1);
    b_hg2_Xsec_vs_fake  = (dict["hg2_Xsec_vs_fake"]==1);
    b_hg1_R2_match  = (dict["hg1_R2_match"]==1);

    jet_R2 = dict["R"]()*dict["R"]();

    if (b_make_AB) {
        response_A = (RooUnfoldResponse*)
            response.Clone(Form("%s_A",response.GetName()));
        response_B = (RooUnfoldResponse*)
            response.Clone(Form("%s_B",response.GetName()));
    }
    ioBinVec pthatbins { Xsec.pthatbins };
    if (b_hg2_Xsec_vs_T) hg2_Xsec_vs_T = new TH2D(
        Form("hg2_Xsec_vs_T_%s",_name),
        ";#it{p}_{T} Truth;#vec{#it{p}}_{T}",
        _hg2ptbins, _hg2ptbins, pthatbins, pthatbins);
    if (b_hg2_Xsec_vs_M) hg2_Xsec_vs_M = new TH2D(
        Form("hg2_Xsec_vs_M_%s",_name),
        ";#it{p}_{T} Measured;#vec{#it{p}}_{T}",
        _hg2ptbins, _hg2ptbins, pthatbins, pthatbins);
    if (b_hg2_Xsec_vs_match) hg2_Xsec_vs_match = new TH2D(
        Form("hg2_Xsec_vs_match_%s",_name),
        ";#it{p}_{T} Truth matched;#vec{#it{p}}_{T}",
        _hg2ptbins, _hg2ptbins, pthatbins, pthatbins);
    if (b_hg2_Xsec_vs_fake) hg2_Xsec_vs_fake = new TH2D(
        Form("hg2_Xsec_vs_fake_%s",_name),
        ";#it{p}_{T} Measured fakes;#vec{#it{p}}_{T}",
        _hg2ptbins, _hg2ptbins, pthatbins, pthatbins);
    if (b_hg1_R2_match) hg1_R2_match = new TH1D(
        Form("hg1_R2_match_%s",_name), 
        "Matched jets; #sqrt((#Delta#phi)^2+(#delta#eta)^2)",
        100, 0., 0.2 );
}

bool ioJetMatcher::do_matching(int pthatbin) {
    vector<double> fakes;
    vector<double> misses;
    vector<pair<double,double>> matches;
    vector<double> v_delta_R2;

    // Find fakes, misses, and matches first, so that if an outlier event
    // is found, it can be filled without filling any other event.

    bool out_bounds_cut {false};
    for (auto& MC : data_MC) {
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (reco.is_matched) continue;
            double delta_R2 { MC(reco,jet_R2) };
            if (delta_R2 != 0) {
                // cut the event if matched pair is out of bounds
                if (out_of_match_bounds(MC.pT,reco.pT)) {
                    out_bounds_cut = true;
                    response_cut.Fill(reco.pT, MC.pT);
                    hg2_cut.Fill(reco.pT, MC.pT);
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
    for (auto& reco : data_reco) {
        if (!reco.is_matched) {
            if (fake_limit == 0 || (reco.pT <= fake_limit)) fakes.push_back(reco.pT);
        }
    }
    if (out_bounds_cut) {
        data_MC.clear();
        data_reco.clear();
        return true;
    }

    // now that is it not an outlier event, fill the data
    double W { Xsec.Xsec(pthatbin) };
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

        if (hg2_Xsec_vs_T)       hg2_Xsec_vs_T    ->Fill(match.second, pthat_val);
        if (hg2_Xsec_vs_M)       hg2_Xsec_vs_M    ->Fill(match.first,  pthat_val);
        if (b_hg2_Xsec_vs_match) hg2_Xsec_vs_match->Fill(match.second, pthat_val);
        if (b_hg1_R2_match) hg1_R2_match->Fill(v_delta_R2[n]);
        ++n;
    }
    for (auto& miss : misses) {
        response.Miss(miss,W);
        response_noweight.Miss(miss);
        // tag split
        if (fillA) response_A->Miss(miss, W);
        else       response_B->Miss(miss, W);
        if (hg2_Xsec_vs_T)       hg2_Xsec_vs_T    ->Fill(miss, pthat_val);
    }
    for (auto& fake : fakes) {
        response.Fake(fake,W);
        response_noweight.Fake(fake);
        if (hg2_Xsec_vs_M)       hg2_Xsec_vs_M   ->Fill(fake, pthat_val);
        if (b_hg2_Xsec_vs_fake)  hg2_Xsec_vs_fake->Fill(fake, pthat_val);
        if (fillA) response_A->Fake(fake, W);
        else       response_B->Fake(fake, W);
    };
    data_MC.clear();
    data_reco.clear();
    return false;
};

void ioJetMatcher::write() {
    /* R2_distr_matches->Write(); */
    hg2_cut.Write();
    response.Write();
    response_noweight.Write();
    response_cut.Write();
    if (b_make_AB) {
        response_A->Write();
        response_B->Write();
    }
    if (b_hg2_Xsec_vs_M) hg2_Xsec_vs_M->Write();
    if (b_hg2_Xsec_vs_T) hg2_Xsec_vs_T->Write();
    if (b_hg2_Xsec_vs_match) hg2_Xsec_vs_match->Write();
    if (b_hg2_Xsec_vs_fake) hg2_Xsec_vs_fake->Write();
    if (b_hg1_R2_match) hg1_R2_match->Write();
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
