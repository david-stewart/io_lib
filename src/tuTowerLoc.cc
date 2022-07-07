#include "tuTowerLoc.h"
#include "tu_fnc.h"

float tuTowerEta(int i_tower) { return __tuTowerEta[i_tower]; };
float tuTowerPhi(int i_tower) { return __tuTowerPhi[i_tower]; };

/* float tuTowerEta_lobound(int i_tower) { */
    /* float val = __tuTowerEta[i_tower]; */
    /* return val; */
/* }; */
void tuDrawTowerLocBoundaryTLines(tuOptMap opt, tuOptMap def){
    opt += def;
    for (int i{1}; i<41; ++i) tuDrawTLine(
            __tuTowerEtaBounds[i],-M_PI,__tuTowerEtaBounds[i],M_PI,opt);
    for (int i{1}; i<121; ++i) tuDrawTLine(
            -1,__tuTowerPhiBounds[i],1.,__tuTowerPhiBounds[i],opt);
};

tuTowerLocator::tuTowerLocator() {
    for (int i{0}; i<41;  ++i) _eta_bounds[i] = __tuTowerEtaBounds[i];
    for (int i{0}; i<121; ++i) _phi_bounds[i] = __tuTowerPhiBounds[i];
};
pair<float,float> tuTowerLocator::eta_bounds(int i_tower) {
    int which = (lower_bound(_eta_bounds.begin(), _eta_bounds.end(), 
                     tuTowerEta(i_tower))-_eta_bounds.begin())-1;
    return {_eta_bounds[which],_eta_bounds[which+1]};
};
pair<float,float> tuTowerLocator::phi_bounds(int i_tower) {
    int which = (lower_bound(_phi_bounds.begin(), _phi_bounds.end(), 
                     tuTowerPhi(i_tower))-_phi_bounds.begin())-1;
    return {_phi_bounds[which],_phi_bounds[which+1]};
};
float tuTowerLocator::eta(int i_tower) { return tuTowerEta(i_tower); };
float tuTowerLocator::phi(int i_tower) { return tuTowerPhi(i_tower); };
float tuTowerLocator::eta(int i_tower, double ratio) {
    auto bounds = eta_bounds(i_tower);
    return bounds.first + (bounds.second-bounds.first)*ratio;
};
float tuTowerLocator::phi(int i_tower, double ratio) {
    auto bounds = phi_bounds(i_tower);
    return bounds.first + (bounds.second-bounds.first)*ratio;
};

TGraph* tuTowerLoc_TGraph(vector<int> i_towers, double x_offset, double y_offset) {
    int npts { (int) i_towers.size() };
    double* phi_pts = new double[npts];
    double* eta_pts = new double[npts];
    int n{0};
    if (x_offset!=0.5 || y_offset!=0.5) {
        tuTowerLocator loc{};
        for (auto& i : i_towers) {
            eta_pts[n] = loc.eta(i,x_offset);
            phi_pts[n] = loc.phi(i,y_offset);
            ++n;
        };
        return new TGraph(npts, eta_pts, phi_pts);
    } else {
        for (auto& i : i_towers) {
            eta_pts[n] = tuTowerEta(i);
            phi_pts[n] = tuTowerPhi(i);
            ++n;
        };
        return new TGraph(npts, eta_pts, phi_pts);
    }
};

TH2D* tuTowerLoc_TH2D(){
    TH2D* hg = new TH2D(tuUniqueName(),";#eta;#phi",40,__tuTowerEtaBounds,120,__tuTowerPhiBounds);
    hg->GetXaxis()->SetTickLength(0.);
    hg->GetYaxis()->SetTickLength(0.);
    hg->SetStats(0);
    return hg;
};

