#ifndef ioJetMatcher__h
#define ioJetMatcher__h

#include <vector>
#include <string>

#include "TH1D.h"
#include "TH2D.h"

#include "RooUnfoldResponse.h"

using std::vector;
using std::sort;

double* ioEdges_pAuJet_prelim_13bins(); // 14 edges for 13 bins

struct ioJetMatcher_float {
    // For filling vectors for Raghav matching algorithm 
    // will have comparison methods in order to std::sort
    float eta;
    float phi;
    float pT;
    bool  is_matched;
	double operator()(const ioJetMatcher_float&, const float jetR2=0.16);
    ioJetMatcher_float(float,float,float);
    ioJetMatcher_float();
};

struct ioJetMatcher_index {
    // for filling vectors for CGAL matching algorithm
    // will have comparison methods in order to std::sort
    int i_MC;   // index MC jet
    int i_reco; // index reco jet
    float dR;   // distance between the matches
};

class ioJetMatcher {
    public:
    RooUnfoldResponse response;
   
    void addjet_MC(float eta, float phi, float pT);
    void addjet_reco(float eta, float phi, float pT);
    void do_matching_highfirst(double weight=1.);// Raghav's algorithm
    void do_matching(double weight=1.); // my algorithm

    ioJetMatcher( RooUnfoldResponse& response, float _jet_R=0.4 );
    void reset();
    void write();
    /* void write(bool with_miss_fakes=true, */ 
               /* bool scale_by_bin_width=true, */
               /* bool write_unified2D=false); */
    std::string tag;
    
    private: //internal data to do the matching
    float jet_R2; // jet_R * jet_R   
	vector<ioJetMatcher_float> data_MC;
	vector<ioJetMatcher_float> data_reco;

};

#endif
