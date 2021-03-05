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
#include "TLine.h"
#include "TGraph.h"
#include "TTree.h"

#include "io_fnc.h"
#include "io_fmt.h"

#include <iostream>

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

// meant to map run-id's
struct ioIntMap {
    /* Reads a file and stores a map of Int -> Int
     * Has the option to test for existence of map (like ioIntList above)
     * Intended to be used primarily as runId -> id on histogram maps.
     */
    ioIntMap(const char* file_name, 
            int index_column=0, 
            int data_column =1,
            bool echo_print=true,
            vector<int> skip_vals={}
    );
    ioIntMap(const char* file_name, 
            int index_column=0, 
            int data_column =1,
            bool echo_print=true,
            ofstream& log,
            vector<int> skip_vals={}
    );
    private:
    ioIntMap_constructor(const char* file_name, 
            int index_column, 
            int data_column,
            bool echo_print,
            vector<int> skip_vals={}
    );

    public:

    /* map<int,double> map_id;  // runid -> map_id */
    map<int,int>    data_map; // runid -> duration
    bool   has_key(int);
    int    operator()(int key); // get data_map[key]; if false, returns 0. as default
    double operator[](int key); // same as (double) operator()(int key)
    int    last_val; // last val returned
    int    size();   // size of map
};

/* vector<pair<ioPadDim,ioPadDim>> ioPadDimGrid( */ 
        /* vector<ioPadDim> x_coord, vector<ioPadDim> y_coord, bool by_rows = true); */

struct ioRunListId {
    public:
    ioRunListId(const char* file_name, bool skip_127053_138064=true);
    map<int,double> map_id;  // runid -> map_id
    map<int,int>    run_sec; // runid -> duration
    bool has_run(int);
    double set_id(int run_id);
    double id; // last set id
    int  size();
    double operator()(int);  // return map_id[runid], unless not present, then -1.
};

struct ioHgStats {
    // used to get basic statistics and cuts
    double mean;
    double stddev;
    void calc_stats(); // just the 

    TAxis* axis;

    /* vector<bool>   cut ; // points that are cut */
    vector<double> vals;
    vector<double> errs;
    vector<double> weight;
    int nbins;

    TLine* get_horizontal_TLine(double cut, bool cut_times_sigma=false);

    TGraph* points_above(double cut, bool cut_times_sigma=false);
    TGraph* points_below(double cut, bool cut_times_sigma=false);
    TGraph* points_between(double cut0, double cut1, bool cut_times_sigma=false);
    TGraph* points();

    double mean_Xsigma(double X); // return StdDev * X + mean

    ioHgStats& cut_above(double cut, bool cut_times_sigma=false);
    ioHgStats& cut_below(double cut, bool cut_times_sigma=false);
    ioHgStats& cut_zeros();
    ioHgStats& cut_to_range(double cut0, double cut1, bool cut_times_sigma=false);
    ioHgStats& restore_points();

    ioHgStats(TH1D* hg, bool cut_zeros=false);
    ioHgStats(TProfile* hg, bool cut_zeros=false, bool weight_by_entries=false);

    int count_points();
    vector<int> bin_indices();
};

// ioMsgTree:
// Used to write messages to a MsgTree in the local file, and read from it
class ioMsgTree {
    private:
        string b_msg;
        TTree tree;
    public:
        ioMsgTree(bool set_echo=true);
        static void read_messages(const char* f_name);
        void msg(string msg);                  // write a message
        void msg(vector<string> messages);     // *ditto*
        void write(); // write to tree
        void dash();  // write dashes to tree
        bool echo_to_cout {false};
};



#endif
