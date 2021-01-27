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

struct __oiJetMaker_Jet {
    __oiJetMaker_Jet(double _pT,double _eta, double _phi);
    double pT, eta, phi;
};

class oiJetMaker {
    static constexpr double pi0mass {0.13498};

    public:
    ioOptMap opt;
    const double jet_R;
    const double                 jetrap;
    const fastjet::JetDefinition jet_def; // {antikt,kt,cambridge}_algorithm
    const fastjet::Selector      jet_selection; // default to not_pure_ghost && jetrap
    const double min_jet_pt;

    // set the values for const members in the constructor
    oiJetMaker(ioOptMap options={},
            ioOptMap defaults={{
            "jet_R",  0.4, 
            "jetrap",  -1, // if -1. will default to 1 - jetR
            "jet_def", "antikt",
            "min_jet_pt", 0.2
        }}
    );

    // add options in the future to jet jet areas with ghost_particles

    vector<PseudoJet> particles {};
    void add_particle(double pT, double eta, double phi);
    int  cluster_jets(); // returns how many jets
    void reset();

    int    n_particles{0};
    int    njets{0};
    int    n_next{-1};
    bool   next(); // used to loop through all data using n_next and njets
    double pT(int i=-1);  // get current jet pT
    double eta(int i=-1); // get current jet eta
    double phi(int i=-1); // get current jet phi

    vector<PseudoJet> pseudo_jets {};
    vector<__oiJetMaker_Jet> jets {}; // filled from pseudo_jets

    /* ~oiJetMaker(); */
};

#endif
