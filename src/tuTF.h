#ifndef tu_TF_h
#define tu_TF_h

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "tuOptMap.h"
/* #include "tu_fnc.h" */

TF1*        tu_TsallisFit            (double m0, double A, double T, double n,  double x_min=0.1, double x_max=10.);
TF1*        tu_dAu_200GeV_TsallisFit (const char* name, double x_min=0.1, double x_max=10.);
TF1*        tu_pp_200GeV_TsallisFit  (const char* name, double x_min=0.1, double x_max=10.);
void        tu_apply_prior           (TF1*, TH1D*); // weight TH1D* by intergral of TF1*
void        tu_apply_prior           (TF1*, TH2D*, TH1D*, bool weight_both=false); 
vector<double> tuQuantilesTF(TH1D* hg, vector<double> percents);
pair<TF1*, tuOptMap> tuFitJESJER(TH1D* hg, double pt_jet, 
        double quant_lo=0.3, double quant_hi=0.99, 
        const char* tag="");
const char* tuUniqueNameTF  (int i=0); // return a unique name to the directory

#endif
