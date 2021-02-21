#ifndef io_func__h
#define io_func__h

#include "TDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include <iostream>
#include "TLatex.h"
#include "TPaveText.h"
#include "io_fmt.h"
#include "TH2D.h"
#include "TPRegexp.h"
#include "TF1.h"

#include "ioOptMap.h"

#define IO_pi      3.14159265
#define IO_twopi   6.28318531
#define IO_piless1 2.14159265

using std::pair;

TH1D* ioDivideTH1(TH1* num, TH1* den, bool norm=false);

const char* ioUniqueName(int i=0); // return a unique name to the directory

void ioWaitPrimitive(int i=0);

vector<int> ioColorVec(int n_colors=2, int palette=kCMYK, bool print=false);

void ioDrawTLatex(const char* msg, double x, double y, 
        ioOptMap options={}, ioOptMap dict= {{
        "TextColor",kBlack, "TextColorAlpha",1., "TextSize",22, "TextFont",43,
        "TextAlign",12, "TextAngle",0., }});

void ioDrawTPaveText(double x0, double y0, double x1, double y1,
        vector<const char*> msg, ioOptMap options={}, ioOptMap dict={{
            "TextColor",kBlack, "TextColorAlpha",1., "FillColor",kWhite, "FillColorAlpha",0.,
            "TextSize",22, "TextFont",43, "TextAlign",12, "TextAngle",0., "BorderOpt","NBndc",
            "LineWidth",0, "FillStyle",4000
        }});

void ioDrawTLine(double x0, double y0, double x1, double y1, 
        ioOptMap options={}, ioOptMap dict = {{
            "LineColor",kBlack, "LineColorAlpha",1., 
            "LineStyle",1,
            "LineWidth",1 }});

void ioDrawTLineBox(double x0, double y0, double x1, double y1, 
        ioOptMap options={});

double io_get_box_integral(TH2D* hg, vector<double>p, vector<double>q);

vector<double> io_vec_BinContent(TH1* hg, bool under_over_flow=false);

vector<double> io_vec_BinError(TH1* hg, bool under_over_flow=false);

vector<double> io_vec_BinCenter(TH1* hg, bool under_over_flow=false);

pair<int, double*> io_vec_BinEdge(TH1* hg) ;

// split strings like: "first name|| second name|| third name|| fourth name"
vector<string> io_split_string(string str) ;

// scale each bin by the 1./width 
void io_scaleByBinWidth(TH1D* hg, double scale_factor=1.);

// scale each bin by the 1./widthx/widthy
void io_scaleByBinWidth(TH2D* hg, double scale_factor=1., bool byXwidth=true, bool byYwidth=true);

// functions to get vector<double> of bin content, error, and edges
// using TH1D*
vector<double> io_vecBinContent(TH1D* hg, bool under_over_flow=false);
vector<double> io_vecBinError  (TH1D* hg, bool under_over_flow=false);
// using TProfile*
vector<double> io_vecBinContent(TProfile* hg, bool under_over_flow=false);
vector<double> io_vecBinError  (TProfile* hg, bool under_over_flow=false);
vector<double> io_vecBinEntries(TProfile* hg, bool under_over_flow=false);
vector<double> io_vecAxisBinCenter (TAxis* axis, bool under_over_flow=false);
vector<double> io_vecAxisBinEdges  (TAxis* axis, bool under_over_flow=false);


// map pi-K-p numbers to 0,1,2,3,4,5 : pi, pi-, K, K-, p, pbar
int    io_geant05(int geantid);
const char* io_geant05_ascii(int geantid);
const char* io_geant05_greek(int geantid);

// get dAu_200_Tsallis and pp_200_Tsallic Fn fits
// Fit for 
TF1*  io_TsallisFit(double m0, double A, double T, double n,  double x_min=0.1, double x_max=10.);
TF1*  io_dAu_200GeV_TsallisFit(const char* name, double x_min=0.1, double x_max=10.);
TF1*  io_pp_200GeV_TsallisFit (const char* name, double x_min=0.1, double x_max=10.);

void io_apply_prior(TF1*, TH1D*); // weight TH1D* by intergral of TF1*
// weight TH2D* by (y-axis) by TH1D*/THF* (per bin)
void io_apply_prior(TF1*, TH2D*, TH1D*, bool weight_both=false); 

TH1D* io_BayesUnfold(TH1D* data, TH1D* T, TH2D* R, int iRepUnfold=3, TH1D* M=nullptr);

TLegend* ioNewTLegend();

float io_dphi(float phi0, float phi1);
bool io_AbsTransPhi(float phi0, float phi1, float lo_bound=1., float hi_bound=IO_piless1);
float io_02pi(float &phi);
float io_02pi(float  phi);

#endif
