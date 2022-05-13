# ifndef tu_fmt__h
# define tu_fmt__h
// test a comment line

#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "tuOptMap.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

// a few default dictionaries
tuOptMap tu_fmt__hg_dict();
tuOptMap tu_fmt__leg_dict();

// default dictionary options for formating TH1D

TLine* tu_fmt(TLine* line, tuOptMap options={}, tuOptMap dict={{
        "LineColor",kBlack, "LineColorAlpha",1., "LineStyle",1,
        "LineWidth",1 }});

TH1* tu_fmt (TH1* hg, tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict());
/* TH2D* tu_fmt (TH2D* hg, tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict()); */
/* TProfile* tu_fmt (TProfile* hg, tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict()); */

/* void tu_fmt_ranges(vector<TH1D*> hgs, tuOptMap dict={}); */
/* void tu_fmt_ranges(vector<TH1D*> hgs, vector<TH1D*> hgs_2, tuOptMap dict={}); */

TCanvas*  tu_fmt (TCanvas* canv, tuOptMap dict={});
TPad*     tu_fmt (TPad* pad,     tuOptMap dict={});
TLegend*  tu_fmt (TLegend* leg,  tuOptMap _override={}, tuOptMap dict=tu_fmt__leg_dict());

/* TH1D* tu_fmt (TH1D* hg, TPad* pad, tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict()); */
/* THStack* tu_fmt (THStack* hg,   tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict()); */
//repeat of the above, but for TGraph
TGraph*            tu_fmt (TGraph* hg,            tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict());
/* TGraphErrors*      tu_fmt (TGraphErrors* hg,      tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict()); */
/* TGraphAsymmErrors* tu_fmt (TGraphAsymmErrors* hg, tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict()); */
TMultiGraph*       tu_fmt (TMultiGraph* hg,       tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict());
TBox* tu_fmt(TBox*, tuOptMap _opts);
/* TGraph*            tu_fmt (TGraph* hg, TPad* pad, tuOptMap _override={}, tuOptMap dict=tu_fmt__hg_dict()); */

# endif
