#ifndef tu_func__h
#define tu_func__h

#include "TDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include <iostream>
#include "TLatex.h"
#include "TPaveText.h"
#include "TH2D.h"
#include "TPRegexp.h"
/* #include "TF1.h" */
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "tuOptMap.h"
/* #include "RooUnfold/RooUnfoldResponse.h" */
#include "tu_fmt.h"
#include "tuConst.h"
/* #include "tu_IOS.h" */
/* #include "tu_enum.h" */
#include "tuOptMap.h"

/* #define TU_pi      3.14159266 // use M_PI */
/* #define TU_2PI   6.28318532 */
/* #define TU_piless1 2.14159265 */
/* #define TU_3to8pi  1.17809725 */
/* #define TU_5to8pi  1.96349541 */
/* #define TU_halfpi  1.5707963267948966 */

using std::pair;

int         tuOpenShape  (int i, int cycle=5); //5 is rotate after kOpenStar
int         tuFullShape  (int i, int cycle=5); //5 is rotate after kFullStar

TH1D*       tuHgBlank     (TH1* hg, tuOptMap opt={}, tuOptMap dict={});
TH1D*       tuHgBlank     (int nbins, double* edges, tuOptMap opt={}, tuOptMap dict={{"MarkerStyle",kDot,"MarkerColor",kWhite,"MarkerAlpha",0.,"LineAlpha",0.}});
void        tu_fmt_ax     (TH1* to, TH1* from);

void        tuPause       (int i=0);
const char* tuUniqueName  (int i=0, const char* prefix="unique_name"); // return a unique name to the directory
const char* tuUniqueName  (int i, string prefix); // return a unique name to the directory
const char* tuUniqueName  (string prefix, int i=0); // return a unique name to the directory
/* const char* tuUniqueName  (const char*,   int i=0); // return a unique name to the directory */
                                                                  //
vector<int> tuColorVec    (int n_colors=2, int palette=kCMYK, bool print=false);
array<double,3> tuPercentAround(TH1*, double);
TH1*        tuDivide      (TH1* num, TH1* den, tuOptMap opt={}, tuOptMap dict={}); // possible numbers: norm, and style-den 
TH1*        tuMultiply    (TH1* num, TH1* den, tuOptMap opt={}, tuOptMap dict={}); // possible numbers: norm, and style-den 
void        tuDrawTLatex  (
                            const char* msg, double x, double y, 
                            tuOptMap options={}, tuOptMap dict= {{
                            "TextColor",kBlack, "TextColorAlpha",1., "TextSize",22, "TextFont",43,
                            "TextAlign",12, "TextAngle",0., }});
void        tuDrawTLine   (
                             double x0, double y0, double x1, double y1, 
                             tuOptMap options={}, tuOptMap dict = {{
                             "LineColor",kBlack, "LineColorAlpha",1., 
                             "LineStyle",1, "LineWidth",1 }});
void        tuDrawTLineBox ( double x0, double y0, double x1, double y1, tuOptMap options={});
void        tuDrawBoxErrors( TGraphAsymmErrors* tgas,
                             tuOptMap options={}, tuOptMap dict = {{
                             "LineColor",kBlack, "LineColorAlpha",1., 
                             "LineStyle",1,
                             "LineWidth",1 }});
double      tuPadxRat      ( double x_in); // get x-coordinatio of ratio x_in
double      tuPadyRat      ( double y_in); // get y-coordinatio of ratio y_in
void        tuDrawTLineHorizontal (double y, tuOptMap options={});
void        tuDrawTLineVertical   (double x, tuOptMap options={});
void        tu_scaleByBinWidth    (TH1* hg, double scale_factor=1.);
void        tu_scaleByBinWidth    (TH2* hg, double scale_factor=1., bool byXwidth=true, bool byYwidth=true);
vector<string> tu_split_string    (string str) ;
void        tu_print(TH1*, const char* tag="");
TLegend*    tuNewTLegend();

// map pi-K-p numbers to 0,1,2,3,4,5 : pi, pi-, K, K-, p, pbar
int         tu_geant05               (int geantid);
const char* tu_geant05_ascii         (int geantid);
const char* tu_geant05_greek         (int geantid);
/* TH1D*       tu_BayesUnfold           (TH1D* data, TH1D* T, TH2D* R, int iRepUnfold=3, TH1D* M=nullptr); */
TH1D*       tuNorm                   (TH1D*, const char which='o'); // 0 for no, 1 for yes, 2 for variable bin-width
// xTH2*  tuDivideTH2byTH1(TH2* num, TH1* den, bool scale_by_cols=true); //
// xTH1D* tuRebin(TH1D*, int, double* =nullptr);
double      tu_get_box_integral      (TH2D* hg, vector<double>p, vector<double>q);
double      tu_get_box_integral      (TProfile* hg, pair<double,double>p={0.,0.}, pair<double,double>q={0.,0.});
double      tu_get_box_mean          (TProfile2D* hg, pair<double,double>p={0.,0.}, pair<double,double>q={0.,0.});
vector<double> tu_vec_BinContent(TH1* hg, bool under_over_flow=false);
vector<double> tu_vec_BinError(TH1* hg, bool under_over_flow=false);
vector<double> tu_vec_BinCenter(TH1* hg, bool under_over_flow=false);
pair<int, double*> tu_vecBinEdge(TH1* hg) ;

// split strings like: "first name|| second name|| third name|| fourth name"

// scale each bin by the 1./width 

// scale each bin by the 1./widthx/widthy

// functions to get vector<double> of bin content, error, and edges
// using TH1D*
/* vector<double> tu_vecBinContent(TH1* hg, bool under_over_flow=false); */
/* vector<double> tu_vecBinError  (TH1D* hg, bool under_over_flow=false); */
/* // using TProfile* */
/* vector<double> tu_vecBinContent(TProfile* hg, bool under_over_flow=false); */
/* vector<double> tu_vecBinError  (TProfile* hg, bool under_over_flow=false); */
/* vector<double> tu_vecBinEntries(TProfile* hg, bool under_over_flow=false); */
/* vector<double> tu_vecAxisBinCenter (TAxis* axis, bool under_over_flow=false); */
/* vector<double> tu_vecAxisBinEdges  (TAxis* axis, bool under_over_flow=false); */


/* array<int,6> tu_geant05_colors { 1179, 1230, 1281, 1332, 1383, 1433 }; */
// get dAu_200_Tsallis and pp_200_Tsallic Fn fits
// Fit for 
// weight TH2D* by (y-axis) by TH1D*/THF* (per bin)
/* TH1D* tu_BayesUnfold(TH1D* data, RooUnfoldResponse* response, int iRepUnfold=3); */
/* RooUnfoldResponse* rebinRooUnfoldResponse(TH1D* hg_bins, RooUnfoldResponse* response); */

float tu_dphi(float phi0, float phi1); // phi1 - phi0, returns in range (0,2pi)
float tu_absDphi(float phi0, float phi1); // phi1 - phi0, returns in range (-pi,pi)
/* bool  tu_isAbsTransPhi(float phi0, float phi1, float lo_bound=1., float hi_bound=tu_piless1); // checks if is transverse to pi */
/* bool  tu_isAbsTrans358pi(float phi0, float phi1, float lo_bound=tu_3to8pi, float hi_bound=tu_5to8pi); // checks if is transverse to pi */
float tu_02pi(float &phi); // puts phi in range [0,2pi]
float tu_02pi(float  phi); // same as above

vector<double> tuQuantiles(TH1D* hg, vector<double> ratios);
vector<int>    tuBinQuantiles(TH1D* hg, vector<double> ratios);
string  tuStringVec(vector<double>, const char* name="vec", const char* fmt="5.2f");

int tu_count_digits(int i, int min_val=1); // returns number of digits in int; for 12->2; for 102311->6
int tuSum(const vector<int>);
int tuSum(const vector<bool>);

vector<int> tuReadIntVec(const char* file, int col=0, bool sort=true, bool strip_commas=true);
// notes: (1) all lines that start with a non-numeric word are treated as comments
//        (2) if col == -1, then read all values

/* vector<double> tuReadFloatVec(const char* file, int col=0, bool sort=false, bool strip_commas=true); */
// notes: (1) all lines that start with a non-numeric word are treated as comments
//        (2) if col == -1, then read all values

//----------------

map<int,string>    tuReadIntStrMap(const char* file, tuOptMap options={},
        tuOptMap dict= {{"tag","none","sort",false, "strip_commas",false,"column","all"}});
map<string,string> tu_VecStrToMapStrStr(vector<string>);
map<string,string> tuReadMapStrStr(const char* file, tuOptMap options={}, 
        tuOptMap dict= {{"tag","none","sort",false, "strip_commas",false,"column","all"}});
map<string,string> tuReadMapStrStr(const char* file, const char* tag, tuOptMap options={}, 
        tuOptMap dict= {{"sort",false, "strip_commas",false,"column","all"}});
vector<string> tuReadStrVec(const char* file, tuOptMap options={}, 
        tuOptMap dict= {{"tag","none","sort",false, "strip_commas",false,"column","all"}});
vector<string> tuReadStrVec(const char* file, const char* tag, tuOptMap options={}, 
        tuOptMap dict= {{"sort",false, "strip_commas",false,"column","all"}});
vector<double> tuReadValVec(const char* file, tuOptMap options={}, 
        tuOptMap dict= {{"tag","none","sort",false, "strip_commas",true,"column","all"}});
vector<double> tuReadValVec(const char* file, const char* tag, tuOptMap options={}, 
        tuOptMap dict= {{"sort",false, "strip_commas",true,"column","all"}});
pair<int,double*> tuReadValsPtr(const char* file, tuOptMap ptions={}, 
        tuOptMap dict= {{"begin_index",0,"end_index",-1,"tag","none",
        "sort",false, "strip_commas",true,"column","all"}});

// same as tuReadValVec, but makes a new array[double] (optionally starting and ending 
// offset in the vector) and returns a pointer to the beginning
double* tuRangeValsPtr(int nbins, double bin_lo, double bin_hi);
 // |-> return pointer to an evenly spaced range of values
double* tuSetValsPtr(vector<double>);
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

TGraph* tuMakeTGraph(vector<double> x, vector<double> y, tuOptMap dict={});
TGraph* tuMakeTGraph(TH1D* hg, bool invert_XY=false, bool skip_zeros=true, bool normalize=false);
TGraph* tuMakeTGraph(TProfile* pr, bool invert_XY=false, bool skip_zeros=true, bool normalize=false);

TGraphErrors* tuMakeTGraphErrors(vector<double> x, vector<double> y, vector<double> y_err, 
                                 vector<double> x_err={}, vector<bool> use={});
TGraphErrors* tuMakeTGraphErrors(TH1D* hg, bool invert_XY=false, bool skip_zeros=true, bool normalize=false );
/* TGraph* tuMakeTGraph(vector<double> x, vector<double> y); */

// find the bin in a vector
int tuwhichbin0(double val, vector<double>&); // return which bin (starting from 0) the data is in: lower bound <= val < upper bound
int tuwhichbin0(double val, TH1D*); // remember that the first bin is zero-indexed
int tuwhichbin1(double val, vector<double>&); // return which bin (starting from 0) the data is in: lower bound <= val < upper bound
int tuwhichbin1(double val, TH1D*); // remember that the first bin is zero-indexed

void tu_normalize_per_row(TH2D*, double overall_factor=1.); // re-weight all entries to 1 per row
void tu_normalize_per_col(TH2D*, double overall_factor=1.); // re-weight all entries to 1 per row

// return the ratio of a circle of radius R outside of a line distance d away 
double tuRatCircleOverLine (double R, double d);
double tuRatCircleInTwoParallelLines (const double d0,const double d1, double C, double R);

/* RooUnfoldResponse tuMakeRooUnfoldResponse(int nbins, double lo_bin, double hi_bin, */ 
            /* const char* tag="", const char* title=""); */
/* RooUnfoldResponse tuMakeRooUnfoldResponse(int nbins, double* edges, */ 
        /* const char* tag="", const char* title=""); */
/* // non-symmetric in truth and measured */
/* RooUnfoldResponse tuMakeRooUnfoldResponse( */
        /* int nb_measured, double lo_measured, double hi_measured, */
        /* int nb_truth, double lo_truth, double hi_truth, */ 
        /* const char* tag="", const char* title=""); */
/* RooUnfoldResponse tuMakeRooUnfoldResponse( */
        /* int nb_measured, double* edges_measured, */ 
        /* int nb_truth,    double* edges_truth, */
        /* const char* tag="", const char* title=""); */

double tuPolyP6_a0_a1x_a2xx_a3y_a4yy_a5xy(double* x, double *p);

/* int tuCheckWordTag(string word, string tag); // if tag=="name", then */
bool tuWordIsTag   (string word, string tag); // if tag=="name", then
bool tuWordIsEndTag(string word, string tag); 
bool tuWordIsTag   (TString word, string tag); // if tag=="name", then
bool tuWordIsEndTag(TString word, string tag); 
bool tuIsAnyTag    (string word); // does it match <\\S*> or </\\S*> ?
bool tuIsAnyTag    (TString word);
// if word==<name> return 1 for start, if word==</name> return 2 for end; else return 0

void tu_normByRow(TH2D* hg, double factor=1.0, bool use_max_val=false);

vector<double> tu_print_first_blank(TH2D*);


const char* tu_cutdiff(int a, int b, const char* fmt = "%6i(%6i,%5.2f)");
// return the fit, along with a dictionary of:
// a0 : height of Gaus fit
// a1 : mean of Gaus fit
// a2 : width of Gaus fit
// JES : a1 - pt_jet
// JER : a2 / pt_jet
// pt_jet : copy of pt_jet in
// bound_lo : quantile from quant_lo
// bound_hi : quantile from quant_hi

/* tuS_keepcut_stats tu_keepcutbin(TH1*, int ibin, double minval, */ 
/*         double zero_val=0., double zero_err=0.); */

/* tuS_keepcut_stats tu_cullsmallbins(TH1*, double min_val, tu area=tu::in, */ 
/*                                    bool remove_underoverflow=true); */
/* vector<int> tu_binvec(TH1* _h, tu=tu::in); */

TH2D* tu_cut_high_sigmaX(TH2D* hg, double n_sigma=5., double offset=0., bool print=true);
TH1D* tu_cut_high_sigmaX(TH1D* hg, double n_sigma=5., double offset=0., bool print=true);

// set the bin to the new value, return the sum of change of values
double tu_setbinzero(TH1* hg, int bin, double val=0, double err=0.);

TH1* tuSetCntErrors(TH1* hg);
TH1* tuAddBinCnt(TH1* hg_to, TH1* hg_from, 
        bool set_bin_errors=false,
        bool rm_under_overflow=true); // returs the hg_to

TH1D* tu_build_CDF(TH1* hg, int first_bin=1, int last_bin=0, double alpha=1.);

// Scrubbers for messy histograms
double tuScrubBlock(TH2* hg, int x0, int x1, int y0, int y1);
double tuScrubBlock(TH1* hg, int x0, int x1);
double tuScrubBins (TH1* hg, int min_bins);
double tuScrubIslands(TH2* hg,  bool isX=true, int nblank=1, double max_scrub_rat = 0.05);
double tuScrubNsigs(TH2D* hg, double nsig, bool isX=true, double q0=0., double q1=1.);
double tuScrubNsigs(TH1D* hg, double nsig, double q0=0., double q1=1.);

vector<int> tuVecScrubNsig (TH2D* hg, double nsig=2., double q0=0, double q1=1., int which=kRight, bool scrub=true);
vector<int> tuVecScrubQuant(TH2D* hg, double quantile=0, int which=kRight, bool scrub=true);
/* vector<int> tuVecScrubIsland(TH2D* hg, double quantile=0, int which=kRight, bool scrub=true); */

int tuVecScrubNsig (TH1D* hg, double nsig=2., double q0=0, double q1=1., int which=kRight, bool scrub=true);
int tuScrubQuant(TH1D* hg, double quantile=0, int which=kRight, bool scrub=true);

double tuVecScrub(TH2D* hg, vector<int> index, int which=kRight);


int   tuZeroCopy(TH2D* h_data, TH2D* h_filter); // return number of bins zero'd

tuOptMap tuCalcRowStats(TH1* hg, double q0=0., double q1=1., double nSig=3., bool cut=false);

TH2D* tuNaiveRebin2D (TH2D* hg, vector<double> x_bins, vector<double> y_bins, string name="", bool in_place=false);
TH1D* tuNaiveRebin1D (TH1D* hg, vector<double> x_bins, string name="", bool in_place=false);

/* void tuAdjustRTM(TH2D* resp_before, TH2D* resp_after, TH1D* M, TH1D* T */

void tuInflate(TH1* hg);
void tuSqrtErr(TH1* hg, bool print_if_wrong=false); // match up bin contents with Sqrt of contents
void tuUpdateTMfromR(TH2D* resp, TH2D* mod_resp, TH1D* truth, TH1D* meas); // update truth and measured from changes (cuts) in resp
                                                                           //
string tuFileName(string input, vector<string> more_name={});
/* string tuFileName(vector<string> input); */

/* void tuInflate(TH2* hg); */
// Add content of hg_from to bins in hg_to
//

// function to trim all bins < n_min from a TH2D or TH1D
/* void tuTrimSmallBins(TH2D* hg, int Nmin, bool cut_underover_flow=true); */
/* void tuTrimSmallBins(TH1D* hg, int Nmin, bool cut_underover_flow=true); */

#endif
