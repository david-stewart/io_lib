#include "oiJetMaker.h"

oiJetMaker::oiJetMaker(ioOptMap _opt, ioOptMap _defaults) :
    opt    { _defaults+_opt },
    jet_R  { opt["jet_R"]() },
    jetrap { opt["jetrap"]() == -1 ? (1.-jet_R) : opt["jetrap"]() },
    calc_areas { opt["calc_areas"]() == 1 },
    ghost_max_rap { opt["ghost_max_rap"]() },
    ghost_R       { opt["ghost_R"]() },
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
    particles[n_particles++].reset_PtYPhiM( pt, eta, phi, PIPLUS_MASS );
};

void oiJetMaker::reset() {
    particles.clear();
    jets.clear();
    pseudo_jets.clear();
    n_particles = 0;
    njets = 0;
    n_next = -1;
};
void oiJetMaker::reset_jets() {
    jets.clear();
    pseudo_jets.clear();
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

int oiJetMaker::size() { return njets; };

__oiJetMaker_Jet& oiJetMaker::operator[](int i) { return jets[i]; };

int oiJetMaker::cluster_jets() {
    if (calc_areas) {
        fastjet::AreaDefinition area_def( 
            fastjet::active_area_explicit_ghosts, 
            fastjet::GhostedAreaSpec(ghost_max_rap, 1, ghost_R)
        );
        fastjet::ClusterSequenceArea clustSeq(particles, jet_def, area_def);
        pseudo_jets = sorted_by_pt( jet_selection( clustSeq.inclusive_jets(min_jet_pt) ));
        while (next()) {
            jets.push_back(
                    {pseudo_jets[n_next].perp(), 
                    pseudo_jets[n_next].eta(), 
                    pseudo_jets[n_next].phi(),
                    pseudo_jets[n_next].area()});
        };
    } else {
        fastjet::ClusterSequence clustSeq(particles, jet_def);
        pseudo_jets = sorted_by_pt( jet_selection( clustSeq.inclusive_jets(min_jet_pt) ));
        while (next()) {
            jets.push_back(
                    {pseudo_jets[n_next].perp(), 
                    pseudo_jets[n_next].eta(), 
                    pseudo_jets[n_next].phi()});
        };
    }
    njets = jets.size();
    return njets;
};
double oiJetMaker::pt(int i) {
    if (i == -1) return jets[n_next].pt;
    return jets[i].pt;
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

__oiJetMaker_Jet::__oiJetMaker_Jet(double _pt,double _eta, double _phi) :
    pt{_pt}, eta{_eta}, phi{_phi}, area{0.} {};
__oiJetMaker_Jet::__oiJetMaker_Jet(double _pt,double _eta, double _phi, double _area) :
    pt{_pt}, eta{_eta}, phi{_phi}, area{_area} {};

