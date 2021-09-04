# ifndef io_fmt__h
# define io_fmt__h
// test a comment line

#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "ioOptMap.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

// default dictionary options for formating TH1D

TLine* io_fmt(TLine* line, ioOptMap options={}, ioOptMap dict={{
        "LineColor",kBlack, "LineColorAlpha",1., "LineStyle",1,
        "LineWidth",1 }});

ioOptMap io_fmt__hg_dict();

void io_fmt_ranges(vector<TH1D*> hgs, ioOptMap dict={});

void io_fmt_ranges(vector<TH1D*> hgs, vector<TH1D*> hgs_2, ioOptMap dict={});

TCanvas* io_fmt (TCanvas* canv, ioOptMap dict={});
TPad*    io_fmt (TPad* pad, ioOptMap dict={});
TLegend* io_fmt (TLegend* leg, ioOptMap dict={});

TProfile* io_fmt (TProfile* hg, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());

TH1D* io_fmt (TH1D* hg,   ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());
TH1D* io_fmt (TH1D* hg, TPad* pad, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());

TH2D* io_fmt (TH2D* hg,   ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());

THStack* io_fmt (THStack* hg,   ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());

//repeat of the above, but for TGraph
TGraph* io_fmt (TGraph* hg, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());
TGraphErrors* io_fmt (TGraphErrors* hg, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());
TGraphAsymmErrors* io_fmt (TGraphAsymmErrors* hg, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());
TMultiGraph* io_fmt (TMultiGraph* hg, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());
TGraph* io_fmt (TGraph* hg, TPad* pad, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());

# endif
