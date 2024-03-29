#include "oiJetMaker.h"

ostream& operator<<(ostream& os, oiJetMaker& out) {
    cout << " Contents of oiJetMaker, ("<<out.pseudojets.size() <<") PseudoJets :" << endl;
    cout << Form("%-5s %-9s %-9s %-9s","index","pt","phi","eta") << endl;
    int i{0};
    for (auto& jet : out.pseudojets) cout << Form("%-5i %-9.5f %-9.5f %-9.5f",i++,jet.perp(),jet.phi(),jet.eta()) << endl;
    return os;
};
oiJetMaker::oiJetMaker(ioOptMap _opt, ioOptMap _defaults) :
    opt    { _defaults+_opt },
    jet_R  { opt["jet_R"]() },
    jetrap { opt["jetrap"]() == -1 ? (1.-jet_R) : opt["jetrap"]() },
    calc_areas { opt["calc_areas"]() == 1 },
    ghost_max_rap { opt["ghost_max_rap"]() },
    ghost_R       { opt["ghost_R"]() },
    jet_def { opt["jet_def"].str() == "antikt"    ? fastjet::antikt_algorithm
            : opt["jet_def"].str() == "kt"        ? fastjet::kt_algorithm
            : opt["jet_def"].str() == "cambridge" ? fastjet::cambridge_algorithm 
            : fastjet::antikt_algorithm, jet_R
    },
    /* jet_selection {  fastjet::SelectorAbsEtaMax(jetrap)}, */
    jet_selection { !fastjet::SelectorIsPureGhost() && fastjet::SelectorAbsEtaMax(jetrap)},
    min_jet_pt { opt["min_jet_pt"]() }
{
    /* if (jet_def.jet_algorithm() == fastjet::antikt_algorithm) cout << " antikt! " << endl; */
    /* if (jet_def.jet_algorithm() == fastjet::kt_algorithm) cout << " kT! " << endl; */
    /* if (jet_def.jet_algorithm() == fastjet::cambridge_algorithm) cout << " cambrdige! " << endl; */
    /* cout << " opt[] " << opt["jet_def"] << "  jet_R " << jet_R << " " << jet_def.R() << endl; */
    /* cout << opt << " end OPT " << endl; */
    /* string JD { opt["jet_def"].str() == "antikt"    ? "fastjet::antikt_algorithm" */
            /* : opt["jet_def"].str() == "kt"        ? "fastjet::kt_algorithm" */
            /* : opt["jet_def"].str() == "cambridge" ? "fastjet::cambridge_algorithm" */
            /* : " DEFAULT"}; */
    /* cout << "JD: " << JD << endl; */
    /* cout << " _opt " << _opt << "||" << endl; */
 /* cout << Form("JetR : %f, jetrap : %f, min_jet_pt: %f", */
        /* jet_R, jetrap, min_jet_pt) << endl; */
};

void oiJetMaker::remass(double mass, int index) {
    if (index == -1) index = n_particles-1;
    auto& p= particles[index];
    p.reset_PtYPhiM( p.perp(), p.eta(), p.phi(), mass);
};
void oiJetMaker::add_particle(PseudoJet jet, double mass) {
    ++n_particles;
    particles.push_back(jet);
    remass(mass);
};
void oiJetMaker::add_particle(PseudoJet jet) {
    ++n_particles;
    particles.push_back(jet);
};
void oiJetMaker::add_particle(double pt, double eta, double phi, double mass){
    ++n_particles;
    particles.push_back(PseudoJet());
    particles[n_particles-1].reset_PtYPhiM( pt, eta, phi, mass );
};
/* void oiJetMaker::add_particle(double pt, double eta, double phi){ */
/*     particles.push_back(PseudoJet()); */
/*     particles[n_particles-1].reset_PtYPhiM( pt, eta, phi, PIPLUS_MASS ); */
/* }; */
void oiJetMaker::add_particle(double pt, double eta, double phi, int index, bool is_neutral, double mass){
    ++n_particles;
    particles.push_back(PseudoJet());
    particles[n_particles-1].reset_PtYPhiM( pt, eta, phi, mass );
    if (is_neutral) index = -index-2; // must preserve -1 for ghost particles
    particles[n_particles-1].set_user_index(index);
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
PseudoJet& oiJetMaker::operator()(int i) { return pseudojets[i]; };

int oiJetMaker::cluster_jets() {
    /* cout << " calc_areas: " << calc_areas << " jet_R " << jet_R << "  min_jet_pt " << min_jet_pt << endl; */
    /* cout << " PAR: " << particles.size() << endl; */
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
    /* cout << " njets:: " << njets << endl; */
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


