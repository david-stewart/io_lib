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
void oiJetMaker::add_particle(double pt, double eta, double phi, int index, bool is_neutral){
    particles.push_back(PseudoJet());
    particles[n_particles].reset_PtYPhiM( pt, eta, phi, PIPLUS_MASS );
    if (is_neutral) index = -index-2; // must preserve -1 for ghost particles
    particles[n_particles].set_user_index(index);
    n_particles++;
};

vector<int> oiJetMaker::get_indices(PseudoJet& jet) {
    // get the indicess of all the sub-jets
    vector<int> indices;
    for (auto& P : jet.constituents()) {
        if (P.user_index()==-1) continue;
        /* if (P.is_pure_ghost()) continue; */
        indices.push_back(P.user_index()); 
    }
    return indices;
};

vector<int> oiJetMaker::get_pos_indices(PseudoJet& jet) {
    // get the indicess of all the sub-jets
    vector<int> indices;
    for (auto& P : jet.constituents()) {
        /* if (P.is_pure_ghost()) continue; */
        int index{P.user_index()};
        if (index>=0) indices.push_back(index); 
    }
    return indices;
};

vector<int> oiJetMaker::get_neg_indices(PseudoJet& jet) {
    /* cout << " b0 " << endl; */
    // get the indicess of all the sub-jets
    vector<int> indices;
    // /*junk*/ cout << " b1 " << endl;
    for (auto& P : jet.constituents()) {
        /* cout << " index: " << jet.user_index() << endl; */
    // /*junk*/ cout << " b2 " << endl;
        /* if (P.user_index()==-2) continue; */
    // /*junk*/ cout << " b3 " << endl;
        int index{P.user_index()};
        if (index ==-1) continue;
    // /*junk*/ cout << " b4 " << endl;
        if (index<0) {
            index = -(index+2);
            indices.push_back(index); 
        }
    }
    return indices;
};

void oiJetMaker::reset() {
    particles.clear();
    jets.clear();
    pseudojets.clear();
    n_particles = 0;
    njets = 0;
    n_next = -1;
    /* cout << " cseq == nullptr " << (cseq==nullptr) <<  endl; */
    /* if (cseq!=nullptr) delete cseq; */
    if (calc_areas) {
        if (cseqarea!=nullptr) delete cseqarea;
        /* cout << " cseqarea == nullptr " << (cseqarea==nullptr) <<  endl; */
    }
};
void oiJetMaker::reset_jets() {
    jets.clear();
    pseudojets.clear();
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
        cseqarea = new fastjet::ClusterSequenceArea(particles, jet_def, area_def);
        /* fastjet::ClusterSequenceArea clustSeq(particles, jet_def, area_def); */
        pseudojets = sorted_by_pt( jet_selection( cseqarea->inclusive_jets(min_jet_pt) ));
        for (auto& jet : pseudojets) {
        /* while (next()) { */
            jets.push_back(
                    {jet.perp(), 
                    jet.eta(), 
                    jet.phi(),
                    jet.area()});
        };
    } else {
        if (cseq!=nullptr) delete cseq;
        cseq = new fastjet::ClusterSequence(particles, jet_def);
        pseudojets = sorted_by_pt( jet_selection( cseq->inclusive_jets(min_jet_pt) ));
        for (auto& jet : pseudojets) {
            jets.push_back(
                    {jet.perp(), 
                    jet.eta(), 
                    jet.phi()});
        };
    }
    njets = jets.size();
    return njets;
};
double oiJetMaker::pt(int i) {
    if (i == -1) return jets[n_next].pt;
    return jets[i].pt;
};
double oiJetMaker::area(int i) {
    if (i == -1) return jets[n_next].area;
    return jets[i].area;
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

bool oiJetMaker::has_N_index(PseudoJet& jet, int index) {
    for (auto& i : oiJetMaker::get_neg_indices(jet)) {
        if (i == index) return true;
    }
    return false;
};

bool oiJetMaker::has_C_index(PseudoJet& jet, int index) {
    for (auto& i : oiJetMaker::get_pos_indices(jet)) {
        if (i == index) return true;
    }
    return false;
};

__oiJetMaker_Jet::__oiJetMaker_Jet(double _pt,double _eta, double _phi) :
    pt{_pt}, eta{_eta}, phi{_phi}, area{0.} {};
__oiJetMaker_Jet::__oiJetMaker_Jet(double _pt,double _eta, double _phi,  double _area) :
    pt{_pt}, eta{_eta}, phi{_phi},  area{_area} {};


