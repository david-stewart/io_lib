# ifndef io_fmt__h
# define io_fmt__h

#include "TH1D.h"
#include "ioOptMap.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"

// default dictionary options for formating TH1D

void io_fmt(TLine* line, ioOptMap options={}, ioOptMap dict={{
        "LineColor",kBlack, "LineColorAlpha",1., "LineStyle",1,
        "LineWidth",1 }});

ioOptMap io_fmt__hg_dict();

void io_fmt_ranges(vector<TH1D*> hgs, ioOptMap dict={});

void io_fmt_ranges(vector<TH1D*> hgs, vector<TH1D*> hgs_2, ioOptMap dict={});

void io_fmt (TCanvas* canv, ioOptMap dict={});
void io_fmt (TPad* pad, ioOptMap dict={});
void io_fmt (TLegend* leg, ioOptMap dict={});

void io_fmt (TProfile* hg, ioOptMap _override={}, 
        ioOptMap dict=io_fmt__hg_dict());

void io_fmt (TH1D* hg, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());
void io_fmt (TH1D* hg, TPad* pad, ioOptMap _override={}, ioOptMap dict=io_fmt__hg_dict());

# endif
