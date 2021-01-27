#include "oiJetMaker.h"

oiJetMaker::oiJetMaker() {};
void oiJetMaker::add_particle(double pt, double eta, double phi){
    in_particles.push_back(PseudoJet());
    in_particles[n_in++].reset_PtYPhiM( pt, eta, phi, pi0mass );
};

