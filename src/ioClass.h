#ifndef ioClass__h
#define ioClass__h

#include <iostream>
#include "TFile.h"
#include <vector>
#include "TPad.h"
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
#include "TF1.h"
#include <set>
#include "TRandom3.h"
#include "RooUnfoldResponse.h"
#include "fastjet/PseudoJet.hh"
#include "TH1D.h"
#include "TH2D.h"
using fastjet::PseudoJet;
using std::array;

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

struct ioRanger {
    // used to give relative location in range
    double lo_range;
    double hi_range;
    double lo_out;
    double hi_out;
    bool has_min_out = false;
    bool has_max_out = false;
    double min_out{0.};
    double max_out{0.};
    void set_min_out(double);
    void set_max_out(double);
    double operator()(double x); // returns ratio between lo_out and hi_out of X in range
    ioRanger(double in_range_lo, double in_range_hi, 
            double out_range_lo=0., 
            double out_range_hi=1.);
    //lo and hi out default to the range
};



struct ioBinVec {
    ioBinVec(int nbins, double lo, double hi);
    ioBinVec(vector<vector<double>>);
    vector<double> bin_centers();
    ioBinVec(TAxis* ax);
    /* ioBinVec operator+=(const ioBinVec& _); */
    /* friend ioBinVec operator+ (const ioBinVec lhs, const ioBinVec& rhs); */

    // two constructors either enter a vector or read one and then use init
    ioBinVec(vector<double>, bool range_double=true);
    ioBinVec(const char* file, ioOptMap options={}, bool use_binspacer=true); 
    ioBinVec(const char* file, const char* tag, ioOptMap options={}, bool use_binspacer=true); 
    ioBinVec(TH1*, const char axis='x');
    ioBinVec(const ioBinVec&);
    ioBinVec& operator= (const ioBinVec& _);
    void init(vector<double>, bool range_double=true);
    /* void update(); */
    // all constructors use build_ptr
    void build_ptr();
    // read the vec from a file with ioReadValVec
    ~ioBinVec();
    double*        ptr;
    bool           has_ptr{false};
    int            size {0}; // number of edges
    vector<double> vec {};
    void set_val(int i, double val);
    void print(ostream& os=std::cout);

    /* int nbins(); // return size_ptr-1 */
    /* operator int (); */
    operator int () const;
    operator double* () const;
    /* operator const double* (); */
    operator vector<double> ();
    double operator[](int); 
    double bin_underflow(); 
    double bin_overflow();
    friend ostream& operator<<(ostream& os, ioBinVec& val);

    vector<double>::const_iterator begin() const;
    vector<double>::const_iterator end() const;
};
// small class used to see if things are in bounds
struct ioInBounds {
    double lo_bound;
    double hi_bound;
    bool operator()(double); // return if in bounds
    void init(ioBinVec);
    ioInBounds(ioBinVec);
    ioInBounds(const char* file, const char* tag="");
    ioInBounds(double, double);
    friend ostream& operator<<(ostream& os, ioInBounds& rhs);
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
struct ioPadDimSet {
    /* ioPadDimSet(int _nPads, // make negative is wish to flip directions */
            /* vector<double> _lefts={}, */ 
            /* vector<double> _rights={}); */
    ioPadDimSet(
            vector<double> _lefts={},  // number of pads is first entry in _lefts if >= 1.
            vector<double> _rights={});
    vector<double> lefts;
    vector<double> rights;
    int            nPads;
    vector<ioPadDim> calc_pads();
    static ioPadDim make_pad(double left, 
            double left_margin, double pad_widht, 
            double right_margin); // 
    /* operator vector<ioPadDim> (); */
};

/* vector<ioPadDim> ioPadDimSet_(int nPads, double left_margin, double right_margin=0.005, */
/*         double overall_left=0.01, double overall_right=0.01); */
/* vector<ioPadDim> ioPadDimSet( */
/*         int nPads, */ 
/*         double left=0.15,  // first left margin */
/*         bool reverse=false, // swap the ordering (for y it's nice to go from top to bottom) */
/*         double right=0.01, // last right margin */
/*         double left_margin=0.05, // how far in on the canvas pad */
/*         double right_margin=0.05, // how far in on the canvas pad */
/*         double left_in=0.001, // inner left margins */
/*         double right_in=0.001 // inner right margins */
/*         ); */
// make a set of nPad diminsions, each with left_margin, right_margin, with
// overall_left on the left, and overall_right on the right, with all remaining space
// in the plots
// example:   ioPadDimSet(2,0.1) results 2 pads, both with left margin of 0.1, right margin
// of 0.01; there is 0.01 to the left and right of both of these; all remaining space is in
// the pads


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

    //FIXME
    /* ioPads ( vector<ioPadDim>, int canvas_width, int canv_heigth ); */
    // * Initialize with either a single vector of ioPadDim (which will go as x0, y0, x1, 
    //   y1, etc...) or two vectors of ioPadDim, which will go as {x0,x1,...} {y0,y1,...}
    // * Add operator()(int,int=0) for accessing y,x pad
    ioPads ( vector<pair<ioPadDim, ioPadDim>> pad_dimensions={{{},{}}}, 
            int canvas_width=0, int canvas_height=0);
    ioPads ( vector<ioPadDim>y_dim, vector<ioPadDim>x_dim={}, int canvas_width=0, int canv_heigth=0 );
    /* ioPads ( int nPads=1, int c_wide=0, int c_height=0 ); // default of 1 and 2 pad TPad sets */
    /* ioPads ( int nYpads=1, int nXpads=1, double y_margin=0.15, double x_margin=0.17, */ 
    /*          int c_wide=0, int c_height=0 ); // default of 1 and 2 pad TPad sets */
    ioPads ( int nYpads=1, int nXpads=1, int c_side=0, int c_height=0, 
            ioPadDimSet Y={{0.15}}, ioPadDimSet X={{0.15}} );
    TPad*  operator()(int col=0, int row=0);

    int nRow{1};
    int nCol{1};

    // must initialize separate from constructor (so that the user has a chance to initialize all the 
    // required options for)
    void init();

    void add_pad(pair<ioPadDim,ioPadDim>&);
    void add_pad(vector<pair<ioPadDim,ioPadDim>> input);

    int canvas_width         { 1200 };
    int canvas_height        {  800 };
    void stamp(const char*, ioOptMap opt={}, 
        ioOptMap dict = {{ 
        "TextColor", (kGray+2), 
        "TextSize", 12,
        "x-loc", .05,
        "y-loc", .05}} );

    // To do here:

};

struct ioIntSet {
    /* ioIntSet(); */
    /* void sort(); */
    vector<int> list {};
    bool operator()(int); // check if argument is in the list
    /* bool add_data(const ioIntSet&); // union with a second set */
    bool has(int) const;
    ioIntSet(vector<int> list); // added for the sake of consts
    ioIntSet(const char* in_file, const char* tag);
    ioIntSet(const char* in_file="", int col=0, bool print=true, bool strip_commas=true);
    ioIntSet(const char* in_file, std::ofstream& log, int col=0, bool print=true, bool strip_commas=true);
    std::ostringstream read_file(const char* in_file, int col=0, bool print=true, bool strip_commas=true);
    int  operator[](int); // return location of arg in list (!Warning: does not check for existence)
                          // warning: may not be meaningful with sort, and existence
    ioIntSet& operator+=(const ioIntSet& rhs);
    ioIntSet& operator-=(const ioIntSet& rhs); // remove union
    ioIntSet& operator*=(const ioIntSet& rhs); // get the union
    friend ostream& operator<<(ostream& os, const ioIntSet& dt);
    int size();
    void clear();
    void write_to_file(const char* file_name, vector<string> comments={});
};

struct ioIntList {
    vector<int> list;
    bool operator()(int); // check if argument is in the list
    bool has(int);
    bool has_not(int);
    ioIntList(const char* in_file, std::ofstream& log, bool print=true);
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
            std::ofstream& log,
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
    bool         has(int) const;
    vector<int>  keys() const; // returns a sorted vector of the keys present
    int&         operator[](int key); // get data_map[key]; if false, returns 0. as default
    int          size() const;   // size of map
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
    ioHgStats(vector<double> x_vals, vector<double> vals, vector<double> errs, bool cut_zeros = false);

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

struct ioIntVec {
    // A class that has at it's heart a vector<vector<int>> of data
    // The first row of data is are the keys of the data, for a kind of data frame table
    // This is a poor man's implementation
    ioIntVec( const char* file_name, bool echo_print=true, vector<string>_tags={} );
    ioIntVec( const char* file_name, std::ofstream& log, bool echo_print=true, vector<string>_tags={});
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

struct ioMinMaxPtr {
    ioMinMaxPtr(){};
    double min_val;
    double max_val;
    bool has_data{false};
    void* min_ptr{nullptr};
    void* max_ptr{nullptr};
    void operator()(double, void* ptr=nullptr);
};

struct ioMinMax {
    ioMinMax(string _name="");

    ioMinMax(double*,       string _name="");
    ioMinMax(int*,          string _name="");
    ioMinMax(unsigned int*, string _name="");
    ioMinMax(short*,        string _name="");
    ioMinMax(char*,         string _name="");

    ioMinMax(double*,       int& _index, string _name="");
    ioMinMax(int*,          int& _index, string _name="");
    ioMinMax(unsigned int*, int& _index, string _name="");
    ioMinMax(short*,        int& _index, string _name="");
    ioMinMax(char*,        int& _index, string _name="");

    ioMinMax(double&,       string _name="");
    ioMinMax(int&,          string _name="");
    ioMinMax(unsigned int&, string _name="");
    ioMinMax(short&,        string _name="");
    ioMinMax(char&,        string _name="");

    int fill_option{-1};
    string name;
    long long int n_entries{0};
    double min{0.};
    double max{0.};
    int    nbins(); // return (int)(max-min)+1
    long long int operator()();
    long long int operator()(int); // to fill an array then
    long long int fill(double val);

    double* ptr_double{nullptr}; // option 0
    int*    ptr_int{nullptr};    // option 1
    unsigned int* ptr_uint{nullptr};   // option 2
    short*  ptr_short{nullptr};  // option 3
    char*   ptr_char{nullptr};  // option 3
    int*    index{nullptr};

    friend ostream& operator<<(ostream& os, ioMinMax& self);
};

struct ioIntBinCnt {
    // used to keep count for inputs are discreet integer locations
    TH1D* hg1;
    TH2D* hg2;
    bool is2D{false};
    void fill(int i, double weight=1.);
    void fill(int i, int j, double weight=1.);
    static double getcnt(TH1D* hg, int i);
    static double getcnt(TH2D* hg, int i, int j);
    void write();
    ioIntBinCnt(const char* name, vector<int> x_dim, const char* title);
    ioIntBinCnt(const char* name, vector<int> x_dim, vector<int> y_dim={}, const char* title="");
};

struct ioFnCaller{
    ioFnCaller(const char* file_data, double(&_fn)(double*,double*));
    double (&fn) (double*, double*);
    double operator()(double,double);
    double *p; // fit parameters, read from data_file
    double *x;
};

struct ioXsec{
    // generalization of src/ioXsec_pAu2015.h
    // ptHat bins, X-sections, and nEvents are read in from file

    // data
    ioXsec (
        const char* tag_file,
        const char* Xsection_tag="Xsection",
        const char* pthatbin_tag="pThat",
        const char* nEvents_tag="nEvents"
    );
    vector<double> Xsection;
    vector<double> pthatbins;
    double pthatbin_center(int);
    int nbins_pthat;
    vector<int> Nevents;
    vector<int> Ncollected;
    long int n_collected_total{0};

    // functions
    int pthatbin(pair<double,double> pthat_bounds);
    
    double Xsec(pair<double,double>, int number_of_events=0); 
    double Xsec(int pthatbin, int number_of_events=0); 
    double operator()(int pthatbin, int num_of_events=0) 
        { return Xsec(pthatbin, num_of_events); };
    double operator()(pair<double,double> pthatbounds, int num_of_events=0) 
        { return Xsec(pthatbounds, num_of_events); };

    void collect (int pthatbin);
    void collect (pair<double,double> pthatbounds);
    void check_pthatbin(int bin);

    TH1D* hg_collected;
};

struct ioIntStrFunctor {
    ioIntStrFunctor ( const char* file, ioOptMap options={}, 
        ioOptMap dict= {{"tag","none","sort",false, "strip_commas",false,"column","all"}});
    const char* operator()(int);
    map<int,string> data;
};
struct ioStrStrFunctor {
    ioStrStrFunctor ( const char* file, ioOptMap options={}, 
        ioOptMap dict= {{"tag","none","sort",false, "strip_commas",false,"column","all"}});
    const char* operator()(const char*);
    map<string,string> data;
};

struct ioFirst {
    bool is_first{true};
    ioFirst(){};
    operator bool();
    bool operator()();
};

struct ioCycleTrue {
    int period;
    int cnt;
    ioCycleTrue(int period_in);
    bool operator()();
    operator bool();
    void reset();
};

struct ioCycleSpacer {
    ioCycleTrue cycle;
    string spacer;
    int n_width;
    void reset();
    ioCycleSpacer(int period, int n_width=0, const char* def_spacer=" ");
    friend ostream& operator<<(ostream& os, ioCycleSpacer& self);
};

struct ioXYbounder {
    // put in two vectors: (sorted) X & (not-sorted) Y
    // then when asking about input x, it will find bin-X and return bound-Y
    vector<double> X;
    vector<double> Y;
    int size;
    double lodef;  // default value if x not found within X
    double hidef;  //
    bool operator()(double x,double y); // returns if val Y is out of bounds at X
    double operator()(double x); // returns boundary value y at x
    ioXYbounder(const char* file, const char* tagX, const char* tagY, 
            ioOptMap options={}); // options are hidef for default-hi and default-lo
                                 // if not set will just use the first and last values in 
                                 // Y for under and overflow
    ioXYbounder(vector<double> x, vector <double> y, ioOptMap options={}); 
    // options are hidef for default-hi and default-lo
    ioXYbounder(); // default to always returning false for no bounds set
};



/* TGraphAsymmErrors* ioMakeTGASE(TH1D* hg, bool invert, array<double,4> x_center, bool skip_zeros, bool normalize ) { */
/*     vector<double> x, x_err, y, y_err; */
/*     TAxis* axis = hg->GetXaxis(); */
/*     if (normalize) hg->Scale(1./hg->Integral()); */
/*     for (int i{1}; i<=hg->GetXaxis()->GetNbins(); ++i) { */
/*         if (hg->GetBinContent(i) == 0 && hg->GetBinError(i) == 0 && skip_zeros) continue; */
/*         x.push_back(axis->GetBinCenter(i)); */
/*         x_err.push_back(0.); */
/*         y.push_back(hg->GetBinContent(i)); */
/*         y_err.push_back(hg->GetBinError(i)); */
        
/*     } */
/*     TGraphErrors* gr; */
/*     double lo = hg->GetXaxis()->GetBinLowEdge(1); */
/*     double hi = hg->GetXaxis()->GetBinUpEdge( hg->GetXaxis()->GetNbins()); */
/*     if (invert) { */
/*         gr = ioMakeTGraphErrors(y,x,x_err,y_err); */
/*         gr->SetMinimum(lo); */
/*         gr->SetMaximum(hi); */
/*     } else { */
/*         gr = ioMakeTGraphErrors(x,y,y_err,x_err); */
/*         gr->GetXaxis()->SetLimits(lo,hi); */
/*     } */
/*     return gr; */
/* }; */

/* TGraphAsymmErrors* ioMakeTGASE(int n, ioBinVec x, ioBinVec y, ioBinVec exleft, ioBinVec exright, */
/*          ioBinVec eyleft, ioBinVec eyright={}); */

/* TGraphAsymmErrors* ioMakeTGASE(TH1D* hg, bool invert_XY=false, */ 
/*         array<double,4> x_center={0.,1.,-1.}, // left/right/center ratio of bin; */ 
/*                                               // center=-1 or < left defaults to middle of left */ 
/*                                               // and right */
/*         bool skip_zeros=true); */
/* /1* TGraph* ioMakeTGraph(vector<double> x, vector<double> y); *1/ */

struct ioPtrDbl {
    ioPtrDbl(TAxis* ax, double bin_loc=0.5, bool get_widths=false);  //<-
    ioPtrDbl(vector<double>); // <-
    ioPtrDbl(TH1*, bool get_errors=false); // also Errors // <-
    ioPtrDbl(const ioPtrDbl&); // <-
    /* ioPtrDbl(ioPtrDbl); // <- */
    ioPtrDbl(const char* file, const char* tag); // <-
    ioPtrDbl(int i); // <-
    ioPtrDbl(){};

    // internal
    vector<double> vec{};
    double*        ptr{nullptr};
    int            size{0};

    /* void init(vector<double>, bool range_double=true); */
    void build_ptr();
    ioPtrDbl& update();
    // read the vec from a file with ioReadValVec
    ~ioPtrDbl();
    /* void set_val(int i, double val); */

    /* int nbins(); // return size_ptr-1 */
    operator int ();
    operator double* ();
    operator vector<double> ();
    double& operator[](int); 
    friend ostream& operator<<(ostream& os, ioPtrDbl& val);
    int throw_error(const char* msg="");

    ioPtrDbl& operator +=(const ioPtrDbl& rhs);
    ioPtrDbl& operator -=(const ioPtrDbl& rhs);
    ioPtrDbl& operator *=(const ioPtrDbl& rhs);
    ioPtrDbl& operator /=(const ioPtrDbl& rhs);
    ioPtrDbl& operator /=(const double rhs);
    ioPtrDbl& operator *=(const double rhs);
    ioPtrDbl& operator +=(const double rhs);
    ioPtrDbl& operator -=(const double rhs);

    ioPtrDbl& make_min(const ioPtrDbl& rhs); // adjust to have just the min values
    ioPtrDbl& make_max(const ioPtrDbl& rhs); // adjust to have just the max values
    ioPtrDbl& abs_diff(const ioPtrDbl& rhs); // adjust to have just the max values
    ioPtrDbl& square_diff(const ioPtrDbl& rhs); // adjust to have just the max values
    ioPtrDbl& abs();
    ioPtrDbl& zero_negatives();
    ioPtrDbl& sqrt();

    vector<double>::iterator begin();
    vector<double>::iterator end();
    ioPtrDbl& operator=(const ioPtrDbl&);

};
ioPtrDbl  operator+(const ioPtrDbl& lhs, const ioPtrDbl& rhs);
ioPtrDbl  operator-(const ioPtrDbl& lhs, const ioPtrDbl& rhs);
ioPtrDbl  io_calc_quadrature(vector<ioPtrDbl>);
ioPtrDbl  io_calc_quadrature(vector<ioPtrDbl>, ioPtrDbl);
ioPtrDbl  io_calc_mean      (vector<ioPtrDbl>);
ioPtrDbl  io_calc_max_bound (vector<ioPtrDbl>);
ioPtrDbl  io_calc_min_bound (vector<ioPtrDbl>);
ioPtrDbl  io_calc_max_berr  (vector<ioPtrDbl>, ioPtrDbl);
ioPtrDbl  io_calc_min_berr  (vector<ioPtrDbl>, ioPtrDbl);
pair<ioPtrDbl,ioPtrDbl>  io_calc_bounds    (vector<ioPtrDbl>);

// some drawing options
TGraphAsymmErrors* io_draw_error_boxes(TH1D* mean, ioPtrDbl, ioOptMap opts, 
        array<double,4> = {-1,-1}, double x_range_lo=0., double x_range_hi=0.);
TGraphAsymmErrors* io_draw_error_boxes(TH1D* mean, array<ioPtrDbl,2>, 
        ioOptMap opts, array<double,4> = {-1,-1}, double x_range_lo=0., double x_range_hi=0.);

struct ioSysErrors {
    /* ioSysErrors(); */
    ioSysErrors(const ioSysErrors&);
    ioSysErrors(TH1*, array<double,4> x_rat={-1,-1,0.5}, array<ioPtrDbl,2>_err={});
    ioSysErrors(TGraphAsymmErrors*, array<double,4> x_rat={-1,-1,0.5});
    ioSysErrors& set_x_range(double, double);
    ioSysErrors& swap_xy ();

    ioSysErrors& setYlow(ioPtrDbl&);
    ioSysErrors& setYhigh(ioPtrDbl&);
    ioSysErrors& setYhilo(ioPtrDbl&); // set Y symmetric

    /* void draw_boxes(ioOptMap dict={}); */

    /* ioSysErrors divide(const ioSysErrors& other); */

    /* ioSysErrors& set (TH1*); // sets x,y, and errors */
    /* ioSysErrors& set (TGraphAsymmErrors*); // sets x,y, and errors */

    int size{0};
    TGraphAsymmErrors* tgase {nullptr};
    ioPtrDbl getY();
    ioPtrDbl getYlow();
    ioPtrDbl getYhigh();

    /* vector<ioPtrDbl> vec_data   {}; // to add in quadrature */
    /* ioSysErrors& add_data   (ioPtrDbl); */
    /* ioSysErrors& add_data   (vector<ioPtrDbl>); */

    /* ioSysErrors& calc_mean(vector<ioPtrDbl> data={}); */
    /* ioSysErrors& calc_quadrature(vector<ioPtrDbl> data ={}); // calculate quadratue relative to the mean */
    /* ioSysErrors& calc_bounds(vector<ioPtrDbl> data ={});     // calculate bounds relative to mean */
    /* ioSysErrors& calc_symmetric_bounds(vector<ioPtrDbl> data ={}); */

    ioSysErrors& set_rat_xbins(array<double,4> rat_rel);
    /* ioSysErrors& set_center_xbins(double rat_center); */
    operator TGraphAsymmErrors* (); // case it to TGraphAsymmErrors
    TGraphAsymmErrors* operator-> ();
};

struct ioJetMatcher_float {
    // For filling vectors for Raghav matching algorithm 
    // will have comparison methods in order to std::sort
    float eta;
    float phi;
    float pT;
    bool  is_matched;
	double operator()(const ioJetMatcher_float&, const float jetR2=0.16);
    ioJetMatcher_float(float,float,float);
    ioJetMatcher_float();
};
bool operator==(const ioJetMatcher_float& L, const ioJetMatcher_float R);

bool operator<(const ioJetMatcher_float& L, const ioJetMatcher_float R) ;

bool operator>(const ioJetMatcher_float& L, const ioJetMatcher_float R) ;

struct io_TF1fitter {
    string name;
    string fn_str;
    TF1* fn;
    unsigned int n_par;
    double lo, hi;

    void fix_match_params(io_TF1fitter& f1_to_match, std::set<int> which_pars);
    // fix parameters to the values found in f1_to_match, set others free
    
    vector<double> fval;
    vector<double> ferr;
    vector<double>& fit(TH1D*, double lo=0, double hi=0);

    vector<double>&     operator()(TH1D*, double lo=0, double hi=0);
    vector<double>&     operator()();
    pair<double,double> operator()(int i); // return <fval, ferr> of parameter i
    double              operator[](int i); // return fval of parameter i

    io_TF1fitter(const char* fnc_str, const char* name="",TH1D* hg=nullptr, 
            double lo=0, double hi=0);
    void set_presets(vector<double>);
    friend ostream& operator<<(ostream& os, io_TF1fitter& ft);
};

struct ioParticleThrower {
    TRandom3 r3;
    vector<unsigned int>    whole;
    vector<double> remainder;
    vector<double> l_bound;
    vector<double> u_bound;
    ioParticleThrower( TH1D* _hg_probs, double multiple=1., unsigned int _seed=0 );
    friend ostream& operator<<(ostream& os, ioBinVec& ioParticleThrower);

    bool throw_particle();
    void set_rand(unsigned int);
    bool is_thrown;
    unsigned int i_bin{0};
    unsigned int i_whole{0};
    unsigned int i_size;
    double pt;
    double phi;
    double eta;
};

struct ioPoissonParticleThrower {
    ioPoissonParticleThrower(TRandom3& _0, TH1D& _1, double mult=1.);
    // Give it a distribution of particle pT values
    // It will then integrate for the mean number of particles, and throw
    // that number, distributed from the distribution
    
    // public
    TRandom3& r3;
    TH1D&     dist;
    bool      throw_particle(); // will keep returning true until the last particle has been thrown.
    double pt;
    double phi;
    double eta;
    double mean;

    // internal
    unsigned int  i_to_throw;
    unsigned int  i_thrown   {0};
};

struct ioHopper1D {
    vector<double> edges;
    vector<pair<double,double>> hopper;
    void reset();
    int nbins;
    double total {0.};
    ioHopper1D(vector<double> _edges);
    vector<pair<double,double>>::iterator begin();
    vector<pair<double,double>>::iterator end();
    int fill(double val, double weight=1.);
};

struct ioSimpleJetMatcher {
    vector<unsigned int> misses {};
    vector<unsigned int> fakes {};
    vector<pair<unsigned int, unsigned int>> matches {};
    ioSimpleJetMatcher(vector<PseudoJet>& jets_T, vector<PseudoJet>& jets_M, double R=0.4);
    ioSimpleJetMatcher(vector<PseudoJet>& jets_T, vector<PseudoJet>& jets_M,  double trigger_phi, 
            pair<double,double> delta_phi_range={0, M_PI/3.}, double R=0.4);
    ioSimpleJetMatcher(vector<PseudoJet>& jets_T, vector<PseudoJet>& jets_M,  double trigger_phi, 
            pair<double,double> delta_phi_range, TH1D* _JES_mean, TH1D* _JER, double _nJER, double R=0.4);
};

struct ioRooUnfoldResponseFiller {
    array<RooUnfoldResponse*,9> response, response_A, response_B;
    ioRooUnfoldResponseFiller(TH1D* form_true, TH1D* form_meas, string tag="");
    vector<int> pthatbins {5,7,9,11,15,25,35,45,55,65};
    double min_T, max_T, min_M, max_M;
    void write(bool print=false);
    void fill(int pthatbin, vector<PseudoJet>& T, vector<PseudoJet>& M, bool is_A);
    void fill(int pthatbin, vector<PseudoJet>& T, vector<PseudoJet>& M, bool is_A, double phi_trig,
            pair<double,double> trig_range);
};

// the same as above, but with the ability to add one set with a weight always of 1., the other with a weighting of W
struct ioRooUnfoldResponseFiller_W {
    array<RooUnfoldResponse*,9> response, response_A, response_B;
    array<RooUnfoldResponse*,9> response_W, response_A_W, response_B_W;
    ioRooUnfoldResponseFiller_W(TH1D* form_true, TH1D* form_meas, string tag="");
    vector<int> pthatbins {5,7,9,11,15,25,35,45,55,65};
    double min_T, max_T, min_M, max_M;
    void write(bool print=false);
    /* bool fill(int pthatbin, vector<PseudoJet>& T, vector<PseudoJet>& M, bool is_A, double W); */
    bool fill(int pthatbin, vector<PseudoJet>& T, vector<PseudoJet>& M, bool
            is_A, double phi_trig, pair<double,double> trig_range, double W);
    /* bool fill(double etatrig, int pthatbin, vector<PseudoJet>& T, */
    /*         vector<PseudoJet>& M, bool is_A, double phi_trig, */
    /*         pair<double,double> trig_range, double W); */
};

struct ioRooUnfoldResponseFiller_W_JESJER {
    TH1D* JES_mean;
    TH1D* JER;
    double nJER_limit;
    
    array<RooUnfoldResponse*,9> response, response_A, response_B;
    array<RooUnfoldResponse*,9> response_W, response_A_W, response_B_W;
    ioRooUnfoldResponseFiller_W_JESJER(TH1D* form_true, TH1D* form_meas, 
            TH1D* _JES_mean, TH1D* _JER, const double _nJER_limit, string tag="");
    vector<int> pthatbins {5,7,9,11,15,25,35,45,55,65};
    double min_T, max_T, min_M, max_M;
    void write(bool print=false);
    void fill(int pthatbin, vector<PseudoJet>& T, vector<PseudoJet>& M, bool is_A, double W);
    void fill(int pthatbin, vector<PseudoJet>& T, vector<PseudoJet>& M, bool is_A, double phi_trig,
            pair<double,double> trig_range, double W);
};


#endif
