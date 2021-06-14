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

// constructors: either provide a RooUnfoldResponse to copy, 
// or provide enough information for them to be constructed
/* void ioJetMatcher::init(float _jet_R) { */
/* 	jet_R2 = _jet_R*_jet_R; */ 
/*     R2_distr_matches = new TH1D(Form("R2_distr_matches_%s",response.GetName()), */
/*             ";R-matches;N",100,0.,1.); */
/*     /1* cout << jet_R2 << endl; *1/ */
/*     data_MC.clear(); */ 
/*     data_reco.clear(); */
/*     response_A = (RooUnfoldResponse*) response.Clone(Form("%s_A",response.GetName())); */
/*     response_B = (RooUnfoldResponse*) response.Clone(Form("%s_B",response.GetName())); */
/* } */
/* ioJetMatcher::ioJetMatcher( RooUnfoldResponse _response, float _jet_R) : */
/*     response{_response } */
/* { init(_jet_R); }; */
/* // construct with new RooUnfoldResponse */
/* ioJetMatcher::ioJetMatcher( int nbins, double lo_bin, double hi_bin, */ 
/*     const char* tag, const char* title, float _jet_R) : */
/*     response{ ioMakeRooUnfoldResponse(nbins, lo_bin, hi_bin, tag, title) } */
/* { init(_jet_R); }; */
/* // construct with new RooUnfoldResponse */
/* ioJetMatcher::ioJetMatcher( int nbins, double* edges, */
/*     const char* tag, const char* title, float _jet_R) : */
/*     response{ ioMakeRooUnfoldResponse(nbins, edges, tag, title) } */
/* { init(_jet_R); }; */
/* // construct with new RooUnfoldResponse */
/* ioJetMatcher::ioJetMatcher( */
/*         int nb_measured, double lo_measured, double hi_measured, */
/*         int nb_truth, double lo_truth, double hi_truth, */ 
/*         const char* tag, const char* title, float _jet_R) : */
/*     response{ ioMakeRooUnfoldResponse(nb_measured, lo_measured, hi_measured, */
/*             nb_truth, lo_truth, hi_truth, tag, title) } */
/* { init(_jet_R); }; */
/* // construct with new RooUnfoldResponse */
/* ioJetMatcher::ioJetMatcher( */
/*         int nb_measured, double* edge_measured, */
/*         int nb_truth, double* edge_truth, */
/*         const char* tag, const char* title, float _jet_R) : */
/*     response{ ioMakeRooUnfoldResponse(nb_measured, edge_measured, */
/*             nb_truth, edge_truth, tag, title) } */
/* { init(_jet_R); }; */
/* // construct with new RooUnfoldResponse */
/* ioJetMatcher::ioJetMatcher( */ 
/*         const char* edge_file, */ 
/*         const char* meas_tag, */
/*         const char* truth_tag, */
/*         const char* name_tag, */
/*         const char* title, */
/*         float _jet_R */
/* ) { */
/*     pair<int,double*> meas_bins  = ioReadValsPtr(edge_file, {{"tag",meas_tag}}); */
/*     pair<int,double*> truth_bins = ioReadValsPtr(edge_file, {{"tag",truth_tag}}); */
/*     response = ioMakeRooUnfoldResponse( */
/*             meas_bins.first-1, meas_bins.second, */
/*             truth_bins.first-1, truth_bins.second, */
/*             name_tag, title); */
/*     init(_jet_R); */
/* }; */

/* void ioJetMatcher::addjet_MC(float eta, float phi, float pT) { */
/*     data_MC.push_back({eta,phi,pT}); */
/* }; */

/* void ioJetMatcher::addjet_reco(float eta, float phi, float pT) { */
/*     data_reco.push_back({eta,phi,pT}); */
/* }; */

/* void ioJetMatcher::do_matching_highfirst(double W) { */
/*     //sort into high to low pT order */
/*     /1* std::sort(data_MC.begin(), data_MC.end(), std::greater<ioJetMatcher_float>()); *1/ */
/*     /1* std::sort(data_reco.begin(), data_reco.end(), std::greater<ioJetMatcher_float>()); *1/ */

/*     fill_A = !fill_A; */
/*     for (auto& MC : data_MC) { */
/*         // cut out out-of-bounds jet */
/*         bool found_match {false}; */
/*         for (auto& reco : data_reco) { */
/*             if (reco.is_matched) continue; */
/*             double delta_R2 { MC(reco,jet_R2) }; */
/*             /1* if (MC(reco,jet_R2)!=0.) { *1/ */
/*             if (delta_R2 != 0) { */
/*                 R2_distr_matches->Fill(delta_R2); */
/*                 reco.is_matched = true; */
/*                 response.Fill(reco.pT, MC.pT, W); */
/*                 if (fill_A) response_A->Fill(reco.pT, MC.pT, W); */
/*                 else             response_B->Fill(reco.pT, MC.pT, W); */
/*                 found_match = true; */
/*                 break; */
/*             } */
/*         } */
/*         if (!found_match) { */
/*             response.Miss(MC.pT,W); */
/*             if (fill_A) response_A->Miss(MC.pT, W); */
/*             else             response_B->Miss(MC.pT, W); */
/*         } */
/*     } */
/*     for (auto& reco : data_reco) { */
/*         if (!reco.is_matched) { */
/*             response.Fake(reco.pT,W); */
/*             if (fill_A) response_A->Fake(reco.pT, W); */
/*             else             response_B->Fake(reco.pT, W); */
/*         }// enter a fake */
/*     } */
/*     /1* cout << " entries: " << response.Hresponse()->GetEntries() << " " *1/ */ 
/*                          /1* << response_A->Hresponse()->GetEntries() << " " *1/ */ 
/*                          /1* << response_B->Hresponse()->GetEntries(); *1/ */
/*     reset(); */
/*     // algorithm: */
/*     // starting with the highest to lowest MC jet -- try and match to the highest to lowest reco jets */
/*     // and keep first match */

/* }; */
/* void ioJetMatcher::reset() { */
/*     data_MC.clear(); */
/*     data_reco.clear(); */
/* }; */
/* void ioJetMatcher::write() { */
/*     R2_distr_matches->Write(); */
/*     response.Write(); */
/*     response_A->Write(); */
/*     response_B->Write(); */
/* }; */

ioJetMatcher_outlier::ioJetMatcher_outlier() {};
void ioJetMatcher_outlier::init(
        map<int,double> boundaries,
        double _pt_fakes,
        double _pt_misses,
        ioXsec* _Xsec
    ) {
    pt_limits = boundaries;
    pt_fakes = _pt_fakes;
    pt_misses = _pt_misses;
    Xsec = _Xsec;
    is_set = true;
};

/* ioJetMatcher_outlier::operator bool () { return is_set; }; */
bool ioJetMatcher_outlier::operator() (int pthatbin, double pt) {
    if (!is_set) return false;
    if (pt_limits.count(pthatbin)==0) return false;
    return (pt > pt_limits[pthatbin]);
};

void ioJetMatcher_outlier::add_match(int pthatbin, double M, double T) {
    M_matches[pthatbin].push_back(M);
    T_matches[pthatbin].push_back(T);
};

void ioJetMatcher_outlier::add_fake(int pthatbin, double M) {
    fakes[pthatbin].push_back(M); };

void ioJetMatcher_outlier::add_miss(int pthatbin, double T) {
    misses[pthatbin].push_back(T); };

void ioJetMatcher_outlier::write_TGraph(
    int pthatbin, const char* name, 
    int markerstyle, int markercolor,
    ioOptMap options,
    ioOptMap dict
) {
    dict += options;
    if (dict["MarkerFakeMiss"]==0) 
        dict["MarkerFakeMiss"]=markerstyle;
    dict["MarkerStyle"]=markerstyle;
    dict["MarkerColor"]=markercolor;
    dict["name"]=name;
    write_TGraph(pthatbin, dict);
};
void ioJetMatcher_outlier::write_TGraph(
      int pthatbin, 
      ioOptMap options,
      ioOptMap dict
) {
    double pthatval = Xsec->pthatbin_center(pthatbin);
    dict += options;
    /* cout << "name " << dict["name"].str() << "  markercolor " << dict["MarkerColor"] << endl; */
    if (pt_limits.count(pthatbin)==0) {
        throw std::runtime_error(
        Form(
            "fatal error in ioJetMatcher_outlier: "
            " limit for requested pthatbin (%i) not set",
            pthatbin)
        );
    }

    // draw TGraph the single point for the boundary
    {
        double *x = new double[1];
        double *y = new double[1];
        x[0] = pt_limits[pthatbin];
        y[0] = pthatval;
        TGraph gr { 1, x, y };
        io_fmt(&gr, dict);
        gr.SetName(Form("%s_limit_%i",dict["name"].c_str(),pthatbin));
        gr.Write();
    }
    
    // draw the matches TGraph
    if (M_matches[pthatbin].size()>0) {
        ioBinVec x { M_matches[pthatbin], false };
        ioBinVec y { T_matches[pthatbin], false };
        TGraph gr { x.size, x, y };
        io_fmt( &gr, options );
        gr.Write( Form("%s_match_%i",dict["name"].c_str(),pthatbin) );
    }
    // draw the miss-fakes TGraph
    dict["MarkerStyle"] = dict["MarkerFakeMiss"];

    ioBinVec x_fakes { fakes[pthatbin], false };
    ioBinVec y_fakes { vector<double>(fakes[pthatbin].size(), pt_fakes), false };

    ioBinVec x_misses { misses[pthatbin], false };
    ioBinVec y_misses { vector<double>(misses[pthatbin].size(), pt_fakes), false };

    ioBinVec x {{x_fakes.vec, x_misses}};
    ioBinVec y {{y_fakes.vec, y_misses}};

    if (x.size > 0) {
        TGraph gr { x.size, x, y };
        io_fmt( &gr, options );
        gr.Write( Form("%s_missfakes_%i",dict["name"].c_str(),pthatbin));
    }
};

ioJetMatcherX::ioJetMatcherX (const char* _name, ioXsec& _Xsec, 
            const char* bin_file, const char* tag_M, const char* tag_T, 
            ioOptMap options, ioBinVec _hg2ptbins, ioOptMap dict
    ) :
    name {_name}, 
    Xsec{_Xsec}, 
    response{ioMakeRooUnfoldResponse(_name,bin_file,tag_M,tag_T)},
    bounds_T{bin_file,tag_T},
    bounds_M{bin_file,tag_M},
    _rand{}
{
    dict += options;
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

ioJetMatcher_outlier& ioJetMatcherX::outlier(char C) {
    if (C == 'F') return fake_outliers;
    if (C == 'M') return M_outliers;
    if (C == 'T') return T_outliers;
    throw std::runtime_error(
    "Fatal call to ioJetMatcherX::outlier, can only use 'F', 'M', or 'T'");
};
void ioJetMatcherX::set_outlier(char C, map<int,double> boundaries) {
    outlier(C).init(boundaries, pt_fakes, pt_misses, &Xsec);
};

void ioJetMatcherX::do_matching(pair<double,double> pthatrange) {
    do_matching(Xsec.pthatbin(pthatrange));
};

void ioJetMatcherX::do_matching(int pthatbin) {
    double W { Xsec.Xsec(pthatbin) };
    double pthat_val { (0.5)*(Xsec.pthatbins[pthatbin]+
                       Xsec.pthatbins[pthatbin+1]) };
    bool b_FillA = _rand.Uniform() < ratio_AtoB;
    /* switch_AB = !switch_AB; */

    for (auto& MC : data_MC) {
        if (cut_outliers && T_outliers(pthatbin,MC.pT)) continue;
        if (b_hg2_Xsec_vs_T) 
            hg2_Xsec_vs_T->Fill(MC.pT,pthat_val);
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (cut_outliers && M_outliers(pthatbin,reco.pT)) continue;
            if (reco.is_matched) continue;
            double delta_R2 { MC(reco,jet_R2) };
            if (delta_R2 != 0) {
                found_match = true;
                reco.is_matched = true;
                if (b_hg2_Xsec_vs_match) 
                    hg2_Xsec_vs_match->Fill(MC.pT,pthat_val);
                if (b_hg1_R2_match) hg1_R2_match->Fill(delta_R2);

                if (T_outliers(pthatbin, MC.pT))   T_outliers.add_match(pthatbin, reco.pT, MC.pT);
                if (M_outliers(pthatbin, reco.pT)) M_outliers.add_match(pthatbin, reco.pT, MC.pT);
                response.Fill(reco.pT, MC.pT, W);
                if (b_FillA)  response_A->Fill(reco.pT, MC.pT, W);
                else            response_B->Fill(reco.pT, MC.pT, W);
                break;
            }
        }
        if (!found_match) {
            response.Miss(MC.pT,W);
            if (T_outliers(pthatbin, MC.pT))   T_outliers.add_miss(pthatbin, MC.pT);
            if (b_FillA) response_A->Miss(MC.pT, W);
            else           response_B->Miss(MC.pT, W);
        }
    }
    for (auto& reco : data_reco) {
        if (hg2_Xsec_vs_M) hg2_Xsec_vs_M->Fill(reco.pT,pthat_val);
        if (!reco.is_matched) {
            if (cut_outliers && fake_outliers(pthatbin,reco.pT)) continue;
            if (hg2_Xsec_vs_fake) hg2_Xsec_vs_fake->Fill(reco.pT,pthat_val);
            if (M_outliers(pthatbin, reco.pT)) M_outliers.add_fake(pthatbin, reco.pT);
            if (fake_outliers(pthatbin, reco.pT)) fake_outliers.add_fake(pthatbin, reco.pT);
            response.Fake(reco.pT,W);
            if (b_FillA) response_A->Fake(reco.pT, W);
            else           response_B->Fake(reco.pT, W);
        }// enter a fake
    }
    reset();
};
void ioJetMatcherX::reset() {
    data_MC.clear();
    data_reco.clear();
};
void ioJetMatcherX::write() {
    /* R2_distr_matches->Write(); */
    response.Write();
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

void ioJetMatcherX::addjet_MC(float eta, float phi, float pT) {
    /* if (bounds_T(pT)) data_MC.push_back({eta,phi,pT}); */
    data_MC.push_back({eta,phi,pT});
};

void ioJetMatcherX::addjet_reco(float eta, float phi, float pT) {
    /* if (bounds_M(pT)) data_reco.push_back({eta,phi,pT}); */
    data_reco.push_back({eta,phi,pT});
};
