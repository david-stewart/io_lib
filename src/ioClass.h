#ifndef ioClass__h
#define ioClass__h

#include <iostream>
#include "TFile.h"
#include <vector>
#include "TPad.h"
#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include <map>
#include <string>

#include "io_fnc.h"
#include "io_fmt.h"

class ioGetter{
    public:
    map<string, vector<TObject*> > got {};
    map <string, TFile*> files    {};
    string path;
    int n_objects;

    ioGetter(string _="");
    TFile* get_file(string _f_name);

    TObject* operator()(string f_name, string object_name);
};

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

    void check_input();
    ioPadDim( double _low, double _p_low, double _p_up, double _up );
    ioPadDim( double _low, double _up );
    ioPadDim( double _low, double _p_low, double _up );
    ioPadDim( );
    void print() const;
    double low_margin () const;
    double up_margin () const;
    bool operator==(ioPadDim& b) const; 
};


struct ioPads {
    //internal data
    TCanvas* canvas = nullptr;
    vector<pair<ioPadDim,ioPadDim>> pad_dimensions;
    vector<TPad*> pads; // all the generated smaller pads
    TPad* canvas_pad;   // a single big pad the size of the canvas

    // constructor: vector of coordinates to make each of the pads (8 coordinates each)
    ioPads ( vector<pair<ioPadDim, ioPadDim>> pad_dimensions={{{},{}}}, int canvas_width=0, int canvas_height=0);

    //FIXME
    // * Initialize with either a single vector of ioPadDim (which will go as x0, y0, x1, 
    //   y1, etc...) or two vectors of ioPadDim, which will go as {x0,x1,...} {y0,y1,...}
    // * Add operator()(int,int=0) for accessing y,x pad
    ioPads ( vector<ioPadDim>, int canvas_width, int canv_heigth );
    ioPads ( vector<ioPadDim>y_dim, vector<ioPadDim>x_dim={}, int canvas_width=0, int canv_heigth=0 );
    ioPads ( int nPads=1, int c_wide=0, int c_height=0 ); // default of 1 and 2 pad TPad sets
    TPad*  operator()(int col=0, int row=0);

    int nCol{1};

    // must initialize separate from constructor (so that the user has a chance to initialize all the 
    // required options for)
    void init();

    void add_pad(pair<ioPadDim,ioPadDim>&);
    void add_pad(vector<pair<ioPadDim,ioPadDim>> input);

    int canvas_width         { 1200 };
    int canvas_height        {  800 };

    // To do here:

};

struct ioIntList {
    vector<int> list;
    bool operator()(int); // check if argument is in the list
    bool has(int);
    bool has_not(int);
    ioIntList(const char* in_file, ofstream& log, bool print=true);
    ioIntList(const char* in_file, bool print=true);
    private:
    string make(const char* in_file, bool print);
};

/* vector<pair<ioPadDim,ioPadDim>> ioPadDimGrid( */ 
        /* vector<ioPadDim> x_coord, vector<ioPadDim> y_coord, bool by_rows = true); */

#endif
