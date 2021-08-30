#include "ioTowerLoc.h"
#include "io_fnc.h"

float ioTowerEta(int i_tower) { return __ioTowerEta[i_tower]; };
float ioTowerPhi(int i_tower) { return __ioTowerPhi[i_tower]; };

/* float ioTowerEta_lobound(int i_tower) { */
    /* float val = __ioTowerEta[i_tower]; */
    /* return val; */
/* }; */
void ioDrawTowerLocBoundaryTLines(ioOptMap opt, ioOptMap def){
    opt += def;
    for (int i{1}; i<41; ++i) ioDrawTLine(
            __ioTowerEtaBounds[i],-IO_pi,__ioTowerEtaBounds[i],IO_pi,opt);
    for (int i{1}; i<121; ++i) ioDrawTLine(
            -1,__ioTowerPhiBounds[i],1.,__ioTowerPhiBounds[i],opt);
};

ioTowerLocator::ioTowerLocator() {
    for (int i{0}; i<41;  ++i) _eta_bounds[i] = __ioTowerEtaBounds[i];
    for (int i{0}; i<121; ++i) _phi_bounds[i] = __ioTowerPhiBounds[i];
};
pair<float,float> ioTowerLocator::eta_bounds(int i_tower) {
    int which = (lower_bound(_eta_bounds.begin(), _eta_bounds.end(), 
                     ioTowerEta(i_tower))-_eta_bounds.begin())-1;
    return {_eta_bounds[which],_eta_bounds[which+1]};
};
pair<float,float> ioTowerLocator::phi_bounds(int i_tower) {
    int which = (lower_bound(_phi_bounds.begin(), _phi_bounds.end(), 
                     ioTowerPhi(i_tower))-_phi_bounds.begin())-1;
    return {_phi_bounds[which],_phi_bounds[which+1]};
};
float ioTowerLocator::eta(int i_tower) { return ioTowerEta(i_tower); };
float ioTowerLocator::phi(int i_tower) { return ioTowerPhi(i_tower); };
float ioTowerLocator::eta(int i_tower, double ratio) {
    auto bounds = eta_bounds(i_tower);
    return bounds.first + (bounds.second-bounds.first)*ratio;
};
float ioTowerLocator::phi(int i_tower, double ratio) {
    auto bounds = phi_bounds(i_tower);
    return bounds.first + (bounds.second-bounds.first)*ratio;
};

TGraph* ioTowerLoc_TGraph(vector<int> i_towers, double x_offset, double y_offset) {
    int npts { (int) i_towers.size() };
    double* phi_pts = new double[npts];
    double* eta_pts = new double[npts];
    int n{0};
    if (x_offset!=0.5 || y_offset!=0.5) {
        ioTowerLocator loc{};
        for (auto& i : i_towers) {
            eta_pts[n] = loc.eta(i,x_offset);
            phi_pts[n] = loc.phi(i,y_offset);
            ++n;
        };
        return new TGraph(npts, eta_pts, phi_pts);
    } else {
        for (auto& i : i_towers) {
            eta_pts[n] = ioTowerEta(i);
            phi_pts[n] = ioTowerPhi(i);
            ++n;
        };
        return new TGraph(npts, eta_pts, phi_pts);
    }
};

TH2D* ioTowerLoc_TH2D(){
    TH2D* hg = new TH2D(ioUniqueName(),";#eta;#phi",40,__ioTowerEtaBounds,120,__ioTowerPhiBounds);
    hg->GetXaxis()->SetTickLength(0.);
    hg->GetYaxis()->SetTickLength(0.);
    hg->SetStats(0);
    return hg;
};

