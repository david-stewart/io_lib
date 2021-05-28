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
#include "TProfile2D.h"

#include "ioOptMap.h"
#include "RooUnfoldResponse.h"

#define IO_pi      3.14159265
#define IO_twopi   6.28318531
#define IO_piless1 2.14159265

using std::pair;

TH1D* ioDivideTH1(TH1* num, TH1* den, bool norm=false);

// Scale a TH2* histogram by columns (or rows) by values in TH1D*
// dimensions must add up
TH2* ioDivideTH2byTH1(TH2* num, TH1* den, bool scale_by_cols=true); //

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

double ioPadxRat(double x_in); // get x-coordinatio of ratio x_in
double ioPadyRat(double y_in); // get y-coordinatio of ratio y_in

/* double io_yPadPerc(double); */
/* double io_xPadPerc(double); */

void ioDrawTLineHorizontal(double y, ioOptMap options={});
void ioDrawTLineVertical(double x, ioOptMap options={});

double io_get_box_integral(TH2D* hg, vector<double>p, vector<double>q);
double io_get_box_integral(TProfile* hg, pair<double,double>p={0.,0.}, pair<double,double>q={0.,0.});
double io_get_box_mean(TProfile2D* hg, pair<double,double>p={0.,0.}, pair<double,double>q={0.,0.});

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
TH1D* io_BayesUnfold(TH1D* data, RooUnfoldResponse* response, int iRepUnfold=3);

TLegend* ioNewTLegend();

float io_dphi(float phi0, float phi1); // phi1 - phi0, returns in range (-pi,pi)
bool  io_isAbsTransPhi(float phi0, float phi1, float lo_bound=1., float hi_bound=IO_piless1); // checks if is transverse to pi
float io_02pi(float &phi); // puts phi in range [0,2pi]
float io_02pi(float  phi); // same as above

vector<double> ioQuantiles(TH1D* hg, vector<double> percents);
string  ioStringVec(vector<double>, const char* name="vec", const char* fmt="5.2f");

int io_count_digits(int i, int min_val=1); // returns number of digits in int; for 12->2; for 102311->6
int ioSum(const vector<int>);
int ioSum(const vector<bool>);

vector<int> ioReadIntVec(const char* file, int col=0, bool sort=true, bool strip_commas=true);
// notes: (1) all lines that start with a non-numeric word are treated as comments
//        (2) if col == -1, then read all values

/* vector<double> ioReadFloatVec(const char* file, int col=0, bool sort=false, bool strip_commas=true); */
// notes: (1) all lines that start with a non-numeric word are treated as comments
//        (2) if col == -1, then read all values

//----------------
vector<double> ioReadValVec(const char* file, ioOptMap options={}, 
        ioOptMap dict= {{"tag","none","sort",false, "strip_commas",true,"column","all"}});
pair<int,double*> ioReadValsPtr(const char* file, ioOptMap ptions={}, 
        ioOptMap dict= {{"begin_index",0,"end_index",-1,"tag","none",
        "sort",false, "strip_commas",true,"column","all"}});

// same as ioReadValVec, but makes a new array[double] (optionally starting and ending 
// offset in the vector) and returns a pointer to the beginning
double* ioRangeValsPtr(int nbins, double bin_lo, double bin_hi);
 // |-> return pointer to an evenly spaced range of values
double* ioSetValsPtr(vector<double>);
 // |->  use all values in vector; if there is a repeating integer then use that as number
 // of bins leading up to the next number
 //   example:
 //         0, 5, 5, 1. 2. 3. = 0 .2 .4 .6 .8 1.0 2. 3.


// notes:  will look for <tag> ... data ... </tag>
//         will skip all values after non-float or non-tag on each line
//         will read all values in "... data ..."
//         can optionally sort data
//         tag=="none" will read all values in file (outside of comments)
//         any word that is not <tag> </tag> or double will terminate the line
//         column==val will only pick out that val per line
//         skip_commas will remove all commas from all lines

TGraph* ioMakeTGraph(vector<double>& x, vector<double>& y);
TGraph* ioMakeTGraph(vector<double> x, vector<double> y);

// find the bin in a vector
int iowhichbin0(double val, vector<double>&); // return which bin (starting from 0) the data is in: lower bound <= val < upper bound
int iowhichbin0(double val, TH1D*); // remember that the first bin is zero-indexed
int iowhichbin1(double val, vector<double>&); // return which bin (starting from 0) the data is in: lower bound <= val < upper bound
int iowhichbin1(double val, TH1D*); // remember that the first bin is zero-indexed

double* ax_doubleptr(vector<int> vec);
/* int iowhichbin(double val, int nBec, double*); */
/* int iowhichbin(double val, vector<double>&, vector<int> remap); // return which bin (starting from 0) the data is in: lower bound <= val < upper bound */
/* int iowhichbin(double val, int, double*,   vector<int> remap); */

double io_R(double x0,double y0,double x1,double y1);

// return the ratio of a circle of radius R outside of a line distance d away 
double ioRatCircleOverLine (double R, double d);
double ioRatCircleInTwoParallelLines (const double d0,const double d1, double C, double R);

RooUnfoldResponse ioMakeRooUnfoldResponse(int nbins, double lo_bin, double hi_bin, 
            const char* tag="", const char* title="");
RooUnfoldResponse ioMakeRooUnfoldResponse(int nbins, double* edges, 
        const char* tag="", const char* title="");
// non-symmetric in truth and measured
RooUnfoldResponse ioMakeRooUnfoldResponse(
        int nb_measured, double lo_measured, double hi_measured,
        int nb_truth, double lo_truth, double hi_truth, 
        const char* tag="", const char* title="");
RooUnfoldResponse ioMakeRooUnfoldResponse(
        int nb_measured, double* edges_measured, 
        int nb_truth,    double* edges_truth,
        const char* tag="", const char* title="");

double ioPolyP6_a0_a1x_a2xx_a3y_a4yy_a5xy(double* x, double *p);

/* int ioCheckWordTag(string word, string tag); // if tag=="name", then */
bool ioWordIsTag   (string word, string tag); // if tag=="name", then
bool ioWordIsEndTag(string word, string tag); 
bool ioWordIsTag   (TString word, string tag); // if tag=="name", then
bool ioWordIsEndTag(TString word, string tag); 
bool ioIsAnyTag    (string word); // does it match <\\S*> or </\\S*> ?
bool ioIsAnyTag    (TString word);
// if word==<name> return 1 for start, if word==</name> return 2 for end; else return 0

#endif