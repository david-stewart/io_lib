#ifndef ioPads__h
#define ioPads__h
/*
   readme.md :
# ioPads
    This is a year 2021 refresh on the myPads class.
    It is used to make sets of TPads on a TCanvas.

    It is not as yet full-featured as myPads.cc, but is considerably cleaner.
    For now, use in ROOT6 with `.L ioPads.cc+`
 
 */
#include <vector>
#include "TPad.h"
#include "TCanvas.h"
#include <iostream>

#include "io_fnc.h"
#include "io_fmt.h"

using namespace std;

struct ioPadDim {
    // a structure to contain the four coordinates requisite for a TPad (in x or y):
    //    low   : left/bottom edge of TPad (outer edge of margin)
    //    p_low : left/bottom edge of the plot area
    //    p_up  : right/top edge of the plot area
    //    up    : right/top edge of TPad
    //
    //    If initialized with only two arguments (low, up), then set p_low to low and p_up to up
    //    If initialized with three arguments (low, plow, up) set p_up to up

    double low;
    double p_low;
    double p_up;
    double up;

    void check_input() {
        if (   low   < 0. || low   > 1. 
            || p_low < 0. || p_low > 1. 
            || p_up  < 0. || p_up  > 1. 
            || up    < 0. || up    > 1. ) {
            cout << " Fatal error: input coordinates for ioPadDim for pads must all "
                " be in range [0,1] " << endl;
            print();
            exit (2);
        } else if ( low > p_low || p_low > p_up || p_up > up ) {
            cout << " Fatal error: input coordinates must monotonically increase " << endl;
            print();
            exit(2);
        }
    };

    ioPadDim( double _low, double _p_low, double _p_up, double _up ) :
        low{_low}, p_low{_p_low}, p_up{_p_up}, up{_up} { check_input(); };
    ioPadDim( double _low, double _up ) : 
        low{_low}, p_low{_low}, p_up{_up}, up{_up} { check_input(); };
    ioPadDim( double _low, double _p_low, double _up ) : 
        low{_low}, p_low{_p_low}, p_up{_up}, up{_up} { check_input(); };
    ioPadDim( ) :
        low{0.}, p_low{0.}, p_up{1.}, up{1.} { check_input(); };

    void print() const {
            cout << Form(" Four points are: (%.2f, %.2f %.2f, %.2f)",low,p_low,p_up,up) << endl;
    };

    double low_margin () const {
        // use to get set the lower margin
        double margin { (p_low - low) / (up - low) };
        if (margin < 0) margin = 0;
        return margin;
    };
    double up_margin () const {
        // use to get set the lower margin
        double margin { (up - p_up) / (up - low) };
        if (margin < 0) margin = 0;
        return margin;
    };
};


struct ioPads {
    //internal data
    TCanvas* canvas = nullptr;
    vector<pair<ioPadDim,ioPadDim>> pad_dimensions;
    vector<ioPadDim> y_fourPoint;
    vector<TPad*> pads; // all the generated smaller pads
    TPad* canvas_pad;   // a single big pad the size of the canvas

    // constructor: vector of coordinates to make each of the pads (8 coordinates each)
    string prefix;
    ioPads ( vector<pair<ioPadDim, ioPadDim>> pad_dimensions={{{},{}}}, int canvas_width=0, int canvas_height=0, const char* prefix="");

    // must initialize separate from constructor (so that the user has a chance to initialize all the 
    // required options for)
    void init();

    void add_pad(pair<ioPadDim,ioPadDim>&);
    void add_pad(vector<pair<ioPadDim,ioPadDim>> input);

    // TText Preference
    int canvas_width         { 1200 };
    int canvas_height        {  800 };

};

ioPads::ioPads ( vector<pair<ioPadDim, ioPadDim>> _pad_dimensions, int _canvas_width, int _canvas_height, const char* _prefix) 
    : pad_dimensions{ _pad_dimensions }, prefix{_prefix}
{
    if (_canvas_width)  canvas_width  = _canvas_width;
    if (_canvas_height) canvas_height = _canvas_height;
};



vector<pair<ioPadDim,ioPadDim>> ioPadDimGrid( vector<ioPadDim> x_coord, vector<ioPadDim> y_coord, 
        bool by_rows = true) {
    // make a grid of inputs (pairs of  ioPadDim) from x_coord[x0, x1, ... xn] and y_coord[y0, y1, ... yn] to
    // <x0,y0>, <x1,y0> ... <xn,y0>
    // <x0,y1>, <x1,y1> ... <xn,y1>
    //       ...
    // <x0,yn>, <x1,yn> ... <xn,yn>
    vector<pair<ioPadDim, ioPadDim>> vec_fourPoint{};
    if (by_rows) { 
        // go across by rows and then down by columns
        for (const auto& y : y_coord) {
            for (const auto& x : x_coord) {
                vec_fourPoint.push_back({x,y});
            }
        }
    } else {
        // go down first column, then next, etc...
        for (const auto& x : x_coord) {
            for (const auto& y : y_coord) {
                vec_fourPoint.push_back({x,y});
            }
        }
    }
    return vec_fourPoint;
};

void ioPads::init() {
    // make and stylize the TCanvas and pads currently in the list
    int i{0};
    canvas = new TCanvas(ioUniqueName(),"",canvas_width, canvas_height);
    canvas->Draw();

    cout << " name " << canvas->GetName() << endl;
    io_fmt(this->canvas);
    /* stylize(this->canvas); */

    i=0;
    canvas_pad = new TPad(ioUniqueName(),"",0.,0.,1.,1.);
    io_fmt(canvas_pad);
    canvas_pad->Draw();

    // add all pads
    add_pad(pad_dimensions);
};


void ioPads::add_pad(pair<ioPadDim,ioPadDim>& coord){
    canvas->cd();

    const ioPadDim x { coord.first };
    const ioPadDim y { coord.second };
    int i{0};
    while (gDirectory->FindObjectAny(Form("%s_%i",canvas_pad->GetName(),i))) { ++i; }
    TPad* p = new TPad(Form("%s_%i",canvas_pad->GetName(),i),"",x.low,y.low,x.up,y.up);
    
    // set the boundaries left(l), right(r), top(t), bottom(b)
    p->SetLeftMargin(x.low_margin());
    p->SetRightMargin(x.up_margin());
    p->SetBottomMargin(y.low_margin());
    p->SetTopMargin(y.up_margin());

    io_fmt(p);
    p->Draw();
    pads.push_back(p);
};

#endif
