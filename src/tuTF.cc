#include "tuTF.h"
#include "TMath.h"
#include <numeric>
#include <iostream>
#include <fstream>
#include "TDirectory.h"

const char* tuUniqueNameTF(int i) {
    while (gDirectory->FindObjectAny(Form("unique_nameTF__%i",i))!=nullptr) ++i;
    return Form("unique_nameTF__%i",i);
};

TF1* tu_pp_200GeV_TsallisFit(const char* name, double x_min, double x_max) {
    if (!strcmp(name,"pi")) return tu_TsallisFit( 0.135, 5.503706201866,
        0.1277746601471, 9.759524221237, x_min, x_max);
    else if (!strcmp(name,"antipi")) return tu_TsallisFit( 0.135, 5.58650216135,
        0.1270067076772, 9.73258646095,  x_min, x_max);
    else if (!strcmp(name,"p")) return tu_TsallisFit( 0.938272, 0.07499819855665,
        0.1758814595453, 10.54524945929, x_min, x_max);
    else if (!strcmp(name,"pbar")) return tu_TsallisFit( 0.938272, 0.06427254263927,
        0.1691919676508, 10.06359691238, x_min, x_max);
    else if (!strcmp(name,"K")) return tu_TsallisFit( 0.493677, 0.06934945866744,
        0.1929677543032, 11.82478291302, x_min, x_max);
    else if (!strcmp(name,"antiK")) return tu_TsallisFit( 0.493677, 0.06934945866744,
        0.1929677543032, 11.82478291302, x_min, x_max);
    else {
        cout << "fatal error: particle name \"" << name << "\" not valid." << endl;
        throw std::runtime_error("Particle name must be pi, antipi, p, pbar, K, antiK");
        return nullptr;
    }
};
//*
void tu_apply_prior(TF1* fn, TH1D* T) {
    // respone: x-axis=Meas. y-axis=True
    // Note: x-axis underflow may contain all the misses
    /* TH1D* T = response->ProjectionY(); */
    TAxis* Y_ax = T->GetXaxis();
    /* TAxis* X_ax = MT->GetXaxis(); */
    for (int y=0; y<Y_ax->GetNbins()+2; ++y){
        if (T->GetBinContent(y) == 0) continue;
        double w = fn->Integral(Y_ax->GetBinLowEdge(y), Y_ax->GetBinUpEdge(y)) /
                T->GetBinContent(y);
            double val = T->GetBinContent(y)*w;
            double err = T->GetBinError(y)*w;
            T->SetBinContent(y,val);
            T->SetBinError(y,err);
    };
};
//*
void tu_apply_prior(TF1* fn, TH2D* MT, TH1D* T, bool both) {
    // Note: this only applies the prior to MT, NOT to T
    // respone: x-axis=Meas. y-axis=True
    // NB: assumed that the x-axis underflow bin contains all the misses
    /* TH1D* T = response->ProjectionY(); */
    TAxis* Y_ax = MT->GetYaxis();
    TAxis* X_ax = MT->GetXaxis();

    for (int y=0; y<Y_ax->GetNbins()+2; ++y){
        if (T->GetBinContent(y) == 0) continue;
        double w = fn->Integral(Y_ax->GetBinLowEdge(y), Y_ax->GetBinUpEdge(y)) /
                T->GetBinContent(y);
        for (int x = 0; x<X_ax->GetNbins()+2; ++x) {
            double val = MT->GetBinContent(x,y) * w;
            double err = MT->GetBinError(x,y) * w;
            MT->SetBinContent(x,y,val);
            MT->SetBinError(x,y,err);
        }
        if (both) {
            double val = T->GetBinContent(y)*w;
            double err = T->GetBinError(y)*w;
            T->SetBinContent(y,val);
            T->SetBinError(y,err);
        }
    };
};
TF1*  tu_TsallisFit(double m0, double A, double T, double n, double x_min, double x_max) {
    int i = 0;
    TF1* fun = new TF1(tuUniqueNameTF(),"[0] / TMath::Power((1.+"
            "(TMath::Sqrt(x*x+[1]*[1])-[1])/([2]*[3])), [3])", x_min, x_max);
    fun->SetParameter(0, A);
    fun->SetParameter(1, m0);
    fun->SetParameter(2, T);
    fun->SetParameter(3, n);
    return fun;
};
//*
TF1* tu_dAu_200GeV_TsallisFit(const char* name, double x_min, double x_max) {
    if ("pi") return tu_TsallisFit( 0.135, 13.62511532585, 
            0.1505198389696, 10.0949606259, x_min, x_max);
    else if (!strcmp(name,"antipi")) return tu_TsallisFit( 0.135, 13.62511532585, 
                0.1505198389696, 10.0949606259, x_min, x_max);
    else if (!strcmp(name,"p")) return tu_TsallisFit( 0.938272, 0.2873036847272, 
            0.2114254602693, 11.44109908352, x_min, x_max);
    else if (!strcmp(name,"pbar")) return tu_TsallisFit( 0.938272, 0.2257531137243, 
            0.2195209346999, 13.02603055528, x_min, x_max);
    else if (!strcmp(name,"K")) return tu_TsallisFit( 0.493677, 0.461550522549, 
            0.2215729988205, 11.61375838522, x_min, x_max);
    else if (!strcmp(name,"antiK")) return tu_TsallisFit( 0.493677, 0.4675391506256, 
            0.2192846302987, 11.44758460939, x_min, x_max);
    else {
        cout << "fatal error: particle name \"" << name << "\" not valid." << endl;
        throw std::runtime_error("Particle name must be pi, antipi, p, pbar, K, antiK");
        return nullptr;
    }
};
/* vector<double> tuQuantilesTF(TH1D* hg, vector<double> percents){ */
/*     cout << hg->GetName() << endl; */
/*     int n {(int)percents.size()}; */
/*     double *x = new double[n]; */
/*     double *q = new double[n]; */
/*     for (int i{0}; i<n; ++i) x[i] = percents[i]; */
/*     hg->GetQuantiles(n,q,x); */
/*     vector<double> vec; */
/*     for (int i{0}; i<n; ++i) vec.push_back(q[i]); */
/*     delete[] x; */
/*     delete[] q; */
/*     return vec; */
/* }; */
/* pair<TF1*, tuOptMap> tuFitJESJER(TH1D* hg, double pt_jet, */ 
/*         double quant_lo, double quant_hi, */ 
/*         const char* tag) */
/* { */
/*     string my_tag = (strcmp(tag,"")==0) ? tuUniqueNameTF() : tag; */
/*     double center = hg->GetXaxis()->GetBinCenter(hg->GetMaximumBin()); */
/*     auto quant = tuQuantilesTF(hg,{quant_lo,quant_hi}); */
/*     if (quant[0] >= center || quant[1] <= center) { */
/*         cout << "Can't fit in tuFitJESJER : the max value bin at " */
/*              << Form("%f doesn't lie between quantiles (%f,%f) at (%f,%f)", */
/*                      center, quant_lo, quant_hi, quant[0], quant[1]) << endl; */
/*         return {nullptr, {{ */
/*             "a0",0., */
/*             "a1",0., */
/*             "a2",0., */
/*             "pt_jet", pt_jet, */
/*             "JES", 0., */
/*             "JES_err",0., */
/*             "JER", 0., */
/*             "JER_err", 0., */
/*             "bound_lo", 0., */
/*             "bound_hi", 0., */
/*             "quant_lo", 0., */
/*             "quant_hi", 0.}} */
/*         }; */
/*     } */
/*     TF1* f1 = new TF1(my_tag.c_str(),"gaus(0)",quant[0],quant[1]); */
/*     f1->SetParameter(0, hg->GetMaximum()); */
/*     f1->SetParameter(1, center); */
/*     f1->SetParameter(2, center*0.1); */
/*     hg->Fit(my_tag.c_str(),"","",quant[0],quant[1]); */
/*     return {f1, {{ */
/*             "a0",f1->GetParameter(0), */
/*             "a1",f1->GetParameter(1), */
/*             "a2",f1->GetParameter(2), */
/*             "pt_jet", pt_jet, */
/*             "JES", (f1->GetParameter(1) - pt_jet), */
/*             "JES_err", (f1->GetParError(1)), */
/*             "JER", (f1->GetParameter(2) ), */
/*             "JER_err", (f1->GetParError(2)), */
/*             "bound_lo", quant[0], */
/*             "bound_hi", quant[1], */
/*             "quant_lo", quant_lo, */
/*             "quant_hi", quant_hi}}}; */
/* }; */
