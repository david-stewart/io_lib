#include "ioTowerLoc.h"
#include "io_fnc.h"

float ioTowerEta(int i_tower) { return __ioTowerEta[i_tower]; };
float ioTowerPhi(int i_tower) { return __ioTowerPhi[i_tower]; };

TGraph* ioTowerLoc_TGraph(vector<int> i_towers) {
    int npts { (int) i_towers.size() };
    double* phi_pts = new double[npts];
    double* eta_pts = new double[npts];
    int n{0};
    for (auto& i : i_towers) {
        phi_pts[n] = ioTowerPhi(i);
        eta_pts[n] = ioTowerEta(i);
        ++n;
    };
    return new TGraph(npts, eta_pts, phi_pts);
};

TH2D* ioTowerLoc_TH2D(){
    TH2D* hg = new TH2D(ioUniqueName(),";#eta;#phi",40,__ioTowerEtaBounds,120,__ioTowerPhiBounds);
    hg->SetStats(0);
    return hg;
};

