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
void ioJetMatcher::init(float _jet_R) {
	jet_R2 = _jet_R*_jet_R; 
    data_MC.clear(); 
    data_reco.clear();
    cout << " name: " << response.GetName() << endl; 
    response_A = (RooUnfoldResponse*) response.Clone(Form("%s_A",response.GetName()));
    response_B = (RooUnfoldResponse*) response.Clone(Form("%s_B",response.GetName()));
}
ioJetMatcher::ioJetMatcher( RooUnfoldResponse _response, float _jet_R) :
    response{_response }
{ init(_jet_R); };
// construct with new RooUnfoldResponse
ioJetMatcher::ioJetMatcher( int nbins, double lo_bin, double hi_bin, 
    const char* tag, const char* title, float _jet_R) :
    response{ ioMakeRooUnfoldResponse(nbins, lo_bin, hi_bin, tag, title) }
{ init(_jet_R); };
// construct with new RooUnfoldResponse
ioJetMatcher::ioJetMatcher( int nbins, double* edges,
    const char* tag, const char* title, float _jet_R) :
    response{ ioMakeRooUnfoldResponse(nbins, edges, tag, title) }
{ init(_jet_R); };
// construct with new RooUnfoldResponse
ioJetMatcher::ioJetMatcher(
        int nb_measured, double lo_measured, double hi_measured,
        int nb_truth, double lo_truth, double hi_truth, 
        const char* tag, const char* title, float _jet_R) :
    response{ ioMakeRooUnfoldResponse(nb_measured, lo_measured, hi_measured,
            nb_truth, lo_truth, hi_truth, tag, title) }
{ init(_jet_R); };
// construct with new RooUnfoldResponse
ioJetMatcher::ioJetMatcher(
        int nb_measured, double* edge_measured,
        int nb_truth, double* edge_truth,
        const char* tag, const char* title, float _jet_R) :
    response{ ioMakeRooUnfoldResponse(nb_measured, edge_measured,
            nb_truth, edge_truth, tag, title) }
{ init(_jet_R); };
// construct with new RooUnfoldResponse
ioJetMatcher::ioJetMatcher( 
        const char* edge_file, 
        const char* meas_tag,
        const char* truth_tag,
        const char* name_tag,
        const char* title,
        float _jet_R
) {
    pair<int,double*> meas_bins  = ioReadValsPtr(edge_file, {{"tag",meas_tag}});
    pair<int,double*> truth_bins = ioReadValsPtr(edge_file, {{"tag",truth_tag}});
    response = ioMakeRooUnfoldResponse(
            meas_bins.first-1, meas_bins.second,
            truth_bins.first-1, truth_bins.second,
            name_tag, title);
    init(_jet_R);
};

void ioJetMatcher::addjet_MC(float eta, float phi, float pT) {
    data_MC.push_back({eta,phi,pT});
};

void ioJetMatcher::addjet_reco(float eta, float phi, float pT) {
    data_reco.push_back({eta,phi,pT});
};

void ioJetMatcher::do_matching_highfirst(double W) {
    //sort into high to low pT order
    /* std::sort(data_MC.begin(), data_MC.end(), std::greater<ioJetMatcher_float>()); */
    /* std::sort(data_reco.begin(), data_reco.end(), std::greater<ioJetMatcher_float>()); */

    half_switch = !half_switch;
    for (auto& MC : data_MC) {
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (reco.is_matched) continue;
            if (MC(reco,jet_R2)) {
                reco.is_matched = true;
                response.Fill(reco.pT, MC.pT, W);
                if (half_switch) response_A->Fill(reco.pT, MC.pT, W);
                else             response_B->Fill(reco.pT, MC.pT, W);
                found_match = true;
                break;
            }
        }
        if (!found_match) {
            response.Miss(MC.pT,W);
            if (half_switch) response_A->Miss(MC.pT, W);
            else             response_B->Miss(MC.pT, W);
        }
    }
    for (auto& reco : data_reco) {
        if (!reco.is_matched) {
            response.Fake(reco.pT,W);
            if (half_switch) response_A->Fake(reco.pT, W);
            else             response_B->Fake(reco.pT, W);
        }// enter a fake
    }
    reset();
    // algorithm:
    // starting with the highest to lowest MC jet -- try and match to the highest to lowest reco jets
    // and keep first match

};
void ioJetMatcher::reset() {
    data_MC.clear();
    data_reco.clear();
};
void ioJetMatcher::write() {
    response.Write();
    response_A->Write();
    response_B->Write();
};
    /* response.Write(); */
    /* TH1D* hg_truth    = (TH1D*) response.Htruth(); */
    /* TH1D* hg_measured = (TH1D*) response.Hmeasured(); */
    /* TH2D* hg_response = (TH2D*) response.Hresponse(); */
    /* /1* cout << " z0 " << endl; *1/ */
    /* if (scale_by_bin_width) { */
    /*     io_scaleByBinWidth(hg_truth); */
    /*     io_scaleByBinWidth(hg_measured); */
    /*     io_scaleByBinWidth(hg_response); */
    /* } */
    /* hg_truth->Write(); */
    /* hg_measured->Write(); */
    /* hg_response->Write(); */

    /* if (with_miss_fakes) { */
    /*     TH1D* miss = (TH1D*) hg_response->ProjectionY(Form("%s_miss",tag.c_str())); */
    /*     miss->Scale(-1.); */
    /*     miss->Add(hg_truth); */
    /*     miss->Write(); */

    /*     TH1D* fakes = (TH1D*) hg_response->ProjectionX(Form("%s_fakes",tag.c_str())); */
    /*     fakes->Scale(-1.); */
    /*     fakes->Add(hg_measured); */
    /*     fakes->Write(); */

    /* /1* cout << " z10 " << endl; *1/ */
    /*     if (add_unified2D) { // make a TH2D histogram is misses and face in a new, */ 
    /*                          // dummy bin to the left and bottem of the TH2D. */
    /*         TAxis *x_axis = hg_response->GetXaxis(); */
    /*         TAxis *y_axis = hg_response->GetYaxis(); */

    /*         int x_nbins = x_axis->GetNbins(); */
    /*         int y_nbins = y_axis->GetNbins(); */

    /*         double *x_edges = new double [ x_nbins+2 ]; */
    /*         double *y_edges = new double [ y_nbins+2 ]; */
    /* /1* TH1D* __1unified = new TH1D( "__1name__","title;x;y", 13, 0., 50.); *1/ */

    /* /1* cout << " z11 " << endl; *1/ */
    /*         x_edges[0] = x_axis->GetBinLowEdge(1) - x_axis->GetBinWidth(1); */
    /*         for (int i{1}; i<=x_nbins; ++i) x_edges[i] = x_axis->GetBinLowEdge(i); */
    /*         x_edges[x_nbins+1] = x_axis->GetBinUpEdge(x_nbins); */

    /*         y_edges[0] = y_axis->GetBinLowEdge(1) - y_axis->GetBinWidth(1); */
    /*         for (int i{1}; i<=y_nbins; ++i) y_edges[i] = y_axis->GetBinLowEdge(i); */
    /*         y_edges[y_nbins+1] = y_axis->GetBinUpEdge(y_nbins); */

    /*         TH2D unified { Form("%s_responseFakeMiss",tag.c_str()), */ 
    /*             Form("%s;Reconstructed (first row is fake jets);Pythia (first column is misses)", */
    /*                     tag.c_str()), x_nbins+1, x_edges, y_nbins+1, y_edges }; */
    /*     /1* TH1D* __unified = new TH1D( "name__","title;x;y", 13, 0., 50.); *1/ */

    /*         for (int ix{1}; ix <=x_nbins; ++ix) { */
    /*             unified.SetBinContent(ix,1, fakes->GetBinContent(ix)); */
    /*             unified.SetBinError  (ix,1, fakes->GetBinError(ix)); */
    /*         } */
    /*         for (int iy{1}; iy <=y_nbins; ++iy) { */
    /*             unified.SetBinContent(1,iy, miss->GetBinContent(iy)); */
    /*             unified.SetBinError  (1,iy, miss->GetBinError(iy)); */
    /*         } */
    /*         for (int ix{1}; ix <=x_nbins; ++ix) */
    /*         for (int iy{1}; iy <=y_nbins; ++iy) { */
    /*             unified.SetBinContent(ix+1,iy+1, hg_response->GetBinContent(ix,iy)); */
    /*             unified.SetBinError  (ix+1,iy+1, hg_response->GetBinError  (ix,iy)); */
    /*         } */
    /*         unified.Write(); */
    /*     } */
    /* } */
/* }; */


