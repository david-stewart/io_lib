#ifndef oiJetMaker__h
#define oiJetMaker__h

#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include <vector>
#include "ioOptMap.h"

using fastjet::PseudoJet;
using fastjet::ClusterSequence;
using fastjet::ClusterSequenceArea;
using std::vector;

struct __oiJetMaker_Jet {
    __oiJetMaker_Jet(double _pt,double _eta, double _phi);
    __oiJetMaker_Jet(double _pt,double _eta, double _phi, double _area);
    double pt, eta, phi, area;
};

class oiJetMaker {
    static constexpr double PI0MASS     {0.13498};
    static constexpr double PIPLUS_MASS {0.136957};

    public:
    ioOptMap opt;
    double jet_R;
    double                 jetrap;
    fastjet::JetDefinition jet_def; // {antikt,kt,cambridge}_algorithm
    fastjet::Selector      jet_selection; // default to not_pure_ghost && jetrap
    double min_jet_pt;
    bool   calc_areas;
    double ghost_max_rap;
    double ghost_R;

    friend ostream& operator<<(ostream& os, oiJetMaker& out);
    // set the values for const members in the constructor
    oiJetMaker(ioOptMap options={},
            ioOptMap defaults={{
            "jet_R",  0.4, 
            "jetrap",  -1, // if -1. will default to 1 - jetR
            "jet_def", "antikt",
            "min_jet_pt", 0.2,
            "calc_areas", 0,
            "ghost_max_rap", 4.,   // 1 for true, 0 for false
            "ghost_R",  0.01
        }}
    );

    // add options in the future to jet jet areas with ghost_particles

    vector<PseudoJet> particles {};
    void remass(double mass, int index=-1);
    void add_particle(PseudoJet);
    void add_particle(PseudoJet, double mass);
    void add_particle(double pt, double eta, double phi, double mass=0.);
    /* void add_particle(double pt, double eta, double phi); */
    void add_particle(double pt, double eta, double phi, int index, bool is_neutral=false, double mass=0.);

    // get the indices in the jets:
    static bool has_N_index(PseudoJet& jet, int index); // see if part. is in jet
    static bool has_C_index(PseudoJet& jet, int index); // " ditto 
    static vector<int> get_indices(PseudoJet&); // return a vector of indicess of a givne PseudoJet
    static vector<int> get_pos_indices(PseudoJet&); // only positive indices
    static vector<int> get_neg_indices(PseudoJet&); // return each negative index (i)-> -i-1


    int  cluster_jets(); // returns how many jets
    void reset();      // reset both jets and particles
    void reset_jets(); // reset only jets (in order to add more particles)

    int    n_particles{0};
    int    njets{0};
    int    size();
    int    n_next{-1};
    bool   next(); // used to loop through all data using n_next and njets
    double pt(int i=-1);  // get current jet pt
    double eta(int i=-1); // get current jet eta
    double phi(int i=-1); // get current jet phi
    double area(int i=-1); // get current jet phi

    vector<PseudoJet> pseudojets {};

    ClusterSequence*     cseq     {nullptr};
    ClusterSequenceArea* cseqarea {nullptr};
    
    vector<__oiJetMaker_Jet> jets {}; // filled from pseudojets
    __oiJetMaker_Jet& operator[](int);
    PseudoJet& operator()(int);

    /* ~oiJetMaker(); */
};

#endif
