#ifndef oiJetMaker__h
#define oiJetMaker__h

#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include <vector>
#include "ioOptMap.h"

using fastjet::PseudoJet;
using std::vector;

class oiJetMaker {
    static constexpr double pi0mass {0.13498};
    bool   make_areas {false};
    double jet_R;
    double ghost_R;
    double ghost_max_rap;
    const  string jet_algo;
    double max_abs_eta_jet = -1.;  // if negative, defaults to 1.-jet_R

    public:
    vector<PseudoJet> in_particles {};
    int n_in{0};
    oiJetMaker();
    ~oiJetMaker();
    void add_particle(double pt, double eta, double phi);
};


#endif
