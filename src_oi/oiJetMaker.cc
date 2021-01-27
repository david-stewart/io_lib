#include "oiJetMaker.h"

oiJetMaker::oiJetMaker(ioOptMap _opt, ioOptMap _defaults) :
    opt    { _defaults+_opt },
    jet_R  { opt["jet_R"]() },
    jetrap { opt["jetrap"]() == -1 ? (1.-jet_R) : opt["jetrap"]() },
    jet_def { opt["jet_def"] == "antikt"    ? fastjet::antikt_algorithm
            : opt["jet_def"] == "kt"        ? fastjet::kt_algorithm
            : opt["jet_def"] == "cambridge" ? fastjet::cambridge_algorithm 
            : fastjet::antikt_algorithm, jet_R
    },
    jet_selection { !fastjet::SelectorIsPureGhost() && fastjet::SelectorAbsEtaMax(jetrap)},
    min_jet_pt { opt["min_jet_pt"]() }
{};

void oiJetMaker::add_particle(double pt, double eta, double phi){
    particles.push_back(PseudoJet());
    particles[n_particles++].reset_PtYPhiM( pt, eta, phi, pi0mass );
};

void oiJetMaker::reset() {
    particles.clear();
    jets.clear();
    pseudo_jets.clear();
    n_particles = 0;
    njets = 0;
    n_next = -1;
};

bool oiJetMaker::next() {
    ++n_next;
    if (n_next < njets) return true;
    else {
        n_next = -1;  
        return false; // automatically reset for use in a while-loop
    }
};

int oiJetMaker::cluster_jets() {

    fastjet::ClusterSequence clustSeq(particles, jet_def);
    pseudo_jets = sorted_by_pt( jet_selection( clustSeq.inclusive_jets(min_jet_pt) ));
    njets = pseudo_jets.size();
    while (next()) {
        jets.push_back(
            {pseudo_jets[n_next].perp(), 
             pseudo_jets[n_next].eta(), 
             pseudo_jets[n_next].phi()});
    };
    return njets;
};
double oiJetMaker::pT(int i) {
    if (i == -1) return jets[n_next].pT;
    return jets[i].pT;
};
double oiJetMaker::eta(int i) {
    if (i == -1) return jets[n_next].eta;
    return jets[i].eta;
};
double oiJetMaker::phi(int i) {
    if (i == -1) return jets[n_next].phi;
    return jets[i].phi;
};
/* oiJetMaker::~oiJetMaker(){}; */

__oiJetMaker_Jet::__oiJetMaker_Jet(double _pT,double _eta, double _phi) :
    pT{_pT}, eta{_eta}, phi{_phi} {};

