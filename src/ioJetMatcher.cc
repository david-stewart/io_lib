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

    fill_A = !fill_A;
    for (auto& MC : data_MC) {
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (reco.is_matched) continue;
            if (MC(reco,jet_R2)) {
                reco.is_matched = true;
                response.Fill(reco.pT, MC.pT, W);
                if (fill_A) response_A->Fill(reco.pT, MC.pT, W);
                else             response_B->Fill(reco.pT, MC.pT, W);
                found_match = true;
                break;
            }
        }
        if (!found_match) {
            response.Miss(MC.pT,W);
            if (fill_A) response_A->Miss(MC.pT, W);
            else             response_B->Miss(MC.pT, W);
        }
    }
    for (auto& reco : data_reco) {
        if (!reco.is_matched) {
            response.Fake(reco.pT,W);
            if (fill_A) response_A->Fake(reco.pT, W);
            else             response_B->Fake(reco.pT, W);
        }// enter a fake
    }
    /* cout << " entries: " << response.Hresponse()->GetEntries() << " " */ 
                         /* << response_A->Hresponse()->GetEntries() << " " */ 
                         /* << response_B->Hresponse()->GetEntries(); */
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

