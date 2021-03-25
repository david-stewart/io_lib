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
#include "io_operators.h"

#include <iostream>

// --------------------------------------------------------------------------------------
// ioGetter:
//   Use: 
//     Easy way to get TObject* from TFiles:
//   Example:
//     initialize:  $:   ioGetter got{};
//     get objects: $:   auto hg = (TH1D*) got("input_file.root","histogram_name");
// --------------------------------------------------------------------------------------
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

// --------------------------------------------------------------------------------------
// ioPadDim:
//   Use:
//      Contains four numbers, used for low, low-margin,high,high-margin
//      for TPad*, for dimensions from left-to-right and bottom-to-top
//      See notes below.
// --------------------------------------------------------------------------------------
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


// --------------------------------------------------------------------------------------
// ioPads:
//   Use:
//      Make a TCanvas with various TPads
//   How:
//      Initialize:
//        Default single canvas:
//          $: ioPads pads{1};
//        Default double canvas:
//          $: ioPads pads{2};
//        Grid of pads (4 rows x 7 cols):
//          $: auto pads = ioPads{ 4, 8 };
//        Grid of pads with custom edges:
//                          {{y dim}  { ydim }       },{{x dim}      {x dim}          width height
//          $: ioPads pads{ {{.5,.98},{0,0.1,0.5,0.5}},{{0.,0.1,0.4},{0.4,.4,.9,.99}},1200,800}
//      Use pad <i>:
//          $: pads(i); // this will activate the given pad
// --------------------------------------------------------------------------------------
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
    /* ioPads ( int nPads=1, int c_wide=0, int c_height=0 ); // default of 1 and 2 pad TPad sets */
    ioPads ( int nYpads=1, int nXpads=1, double y_margin=0.15, double x_margin=0.17, 
             int c_wide=0, int c_height=0 ); // default of 1 and 2 pad TPad sets
    TPad*  operator()(int col=0, int row=0);

    int nRow{1};

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
    int  operator[](int); // return location of arg in list (!Warning: does not check for existence)
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
            int index_column,
            int data_column,
            bool echo_print,
            ofstream& log,
            vector<int> skip_vals={}
    );
    private:
    string ioIntMap_constructor(const char* file_name, 
            int index_column, 
            int data_column,
            bool echo_print,
            vector<int> skip_vals={}
    );
    public:
    map<int,int> data_map; // runid -> duration
    bool         has(int);
    vector<int>  keys(); // returns a sorted vector of the keys present
    int&         operator[](int key); // get data_map[key]; if false, returns 0. as default
    int          size();   // size of map
};


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
    vector<double> unmasked_vals();
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
    ioHgStats& mask(const vector<bool>  mask); // mask all points in vector which are TRUE
    /* ioHgStats& mask(vector<int>&  mask); */
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
        void slurp_file(const char* which_file); // write a given file to input
};


// ioIntVec
// NOTE: This is not currently implemented, as it is currently lower priority.
// -------------------------------
//   Data is a table with
//            tag1   tag2   tag3
//   runid.0  true   false  true
//   runid.1  true   true   true
//   runid.2  false  false  true
//   runid.3  true   false  false
// -------------------------------
//
//
/* struct ioIntVec { */
/*     // constructors and ostream */
/*     ioIntVec( const char* file_name, bool echo_print=true, vector<string>_tags={} ); */
/*     ioIntVec( const char* file_name, ofstream& log, bool echo_print=true, vector<string>_tags={}); */
/*     string ioIntVec_constructor( const char* file_name, bool echo_print, vector<string>_tags={} ); */


/*     // data members */
/*     /1* map<int,vector<int>> data_map {}; // runid -> vector<bool> *1/ */

/*     vector<int> keys; */
/*     vector<vector<int>> data; */

/*     vector<string> tags {};   // names of all the columns */

/*     // access data and manipulate data */
/*     int             i_tag(string tag); // returns column value of tag */
/*     vector<int>     tag_cols(vector<string> _tags, const char* err_msg_name=""); */
/*     bool            has_tag(string tag); */
/*     ioIntVec&       swap_tags(string tag0, string tag1); // swap column locations */
/*     ioIntVec&       subset(vector<string> tag); */
/*     ioIntVec&       rm_tags(vector<string> tag); */
/*     ioIntVec&       add_tag(string tag, int default_val=-1); */ 
/*     void            rename_tag(string, string); */
/*     bool            has_key(int key);     // checks if it has key */
/*     vector<int>&    operator[](int key);   // returns the vector at entry */
/*     int             size();   // size of map */

/*     vector<int>     vals(string tag, vector<bool> mask={}, bool keep_on_true=true); // keep all values */

/*     vector<int>     vals(string tag, string mask); // keep all values */
/*     vector<int>     vals(string tag, vector<string> mask); // TODO [ ] new */
/*     vector<int>     vals(string tag, vector<bool>); // TODO [ ] new */

/*     vector<int>     get_keys(vector<bool> mask, bool keep_on_true=true); // keep all values */
/*                     // from mask */
/*     vector<bool>    is_any(vector<pair<string,bool>> mask_keep_true); // per row, see if any match bool */
/*     vector<bool>    is_all(vector<pair<string,bool>> mask_keep_true); // per row, see if all match bool */

/*     friend ostream& operator<<(ostream& os, ioIntVec& dt); */
/*     void            write_to_file(const char* file_name, vector<string> comments={}); */
/* }; */
struct ioIntVec {
    // constructors and ostream
    ioIntVec( const char* file_name, bool echo_print=true, vector<string>_tags={} );
    ioIntVec( const char* file_name, ofstream& log, bool echo_print=true, vector<string>_tags={});
    string ioIntVec_constructor( const char* file_name, bool echo_print, vector<string>_tags={} );


    // data members
    /* map<int,vector<int>> data_map {}; // runid -> vector<bool> */

    vector<int> keys;
    vector<vector<int>> data;

    vector<string> tags {};   // names of all the columns

    // access data and manipulate data
    int             i_tag(string tag); // returns column value of tag
    vector<int>     tag_cols(vector<string> _tags, const char* err_msg_name="");
    bool            has_tag(string tag);
    ioIntVec&       swap_tags(string tag0, string tag1); // swap column locations
    ioIntVec&       subset(vector<string> tag);
    ioIntVec&       rm_tags(vector<string> tag);
    ioIntVec&       add_tag(string tag, int default_val=-1); 
    void            rename_tag(string, string);
    bool            has_key(int key);     // checks if it has key
    vector<int>&    operator[](int key);   // returns the vector at entry
    int             size();   // size of map
    vector<int>     vals(string tag); // column of data under "tag"
    vector<int>     vals(string tag, vector<bool> mask);// column of data under "tag" where mask is True
    vector<int>     get_keys(vector<bool> mask, bool keep_on_true=true); // keep all values
    vector<bool>    mask(string tag,const char* module="mask");
    vector<bool>    is_any(vector<string> tag);
    vector<bool>    is_all(vector<string> tag);

    private:
    bool flip_tag(string& tag); // tag starts with "!" then strip it and return false, else return true;
    void assert_tag(string tag, const char* routine=""); // fail is tag not in tags

    public:
    friend ostream& operator<<(ostream& os, ioIntVec& dt);
    void            write_to_file(const char* file_name, vector<string> comments={});
};

struct ioMinMax {
    ioMinMax(string _name="");

    ioMinMax(double*,       int& _index, string _name="");
    ioMinMax(int*,          int& _index, string _name="");
    ioMinMax(unsigned int*, int& _index, string _name="");
    ioMinMax(short*,        int& _index, string _name="");

    ioMinMax(double&,       string _name="");
    ioMinMax(int&,          string _name="");
    ioMinMax(unsigned int&, string _name="");
    ioMinMax(short&,        string _name="");

    int fill_option{-1};
    string name;
    long long int n_entries{0};
    double min{0.};
    double max{0.};
    int    nbins(); // return (int)(max-min)+1
    long long int operator()(double);
    long long int operator()();

    double* ptr_double{nullptr}; // option 0
    int*    ptr_int{nullptr};    // option 1
    unsigned int* ptr_uint{nullptr};   // option 2
    short*  ptr_short{nullptr};  // option 3
    int* index{nullptr};

    friend ostream& operator<<(ostream& os, ioMinMax& self);
};
#endif
