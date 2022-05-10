#ifndef tuClass__h
#define tuClass__h

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
#include <string>
#include "tuOptMap.h"
/* #include "tu_fnc.h" */
/* #include "tu_fmt.h" */
/* #include "tu_operators.h" */
/* #include "TF1.h" */ // FIXME :: compiler doesn't like the TF1.h on the imac
#include <set>
#include "TRandom3.h"
#include "TAxis.h"

// // --------------------------------------------------------------------------------------
// // tuGetter:
// //   Use: 
// //     Easy way to get TObject* from TFiles:
// //   Example:
// //     initialize:  $:   tuGetter got{};
// //     get objects: $:   auto hg = (TH1D*) got("input_file.root","histogram_name");
// // --------------------------------------------------------------------------------------
class tuGetter{
    public:
    map<string, vector<TObject*> > got {};
    map <string, TFile*> files    {};
    /* string path; */
    int n_objects;

    tuGetter();
    TFile* get_file(string f_name);

    TObject* operator()(string f_name, string object_name);
};
// // --------------------------------------------------------------------------------------
// // tuBinVec:
// //   Use: 
// //     Easy way to get TObject* from TFiles:
// //   Example:
// //     initialize:  $:   tuGetter got{};
// //     get objects: $:   auto hg = (TH1D*) got("input_file.root","histogram_name");
// // --------------------------------------------------------------------------------------
 struct tuBinVec {
     tuBinVec(vector<double>);
     vector<double> bin_centers();
     tuBinVec(TAxis* ax);
     tuBinVec(const char* file, tuOptMap options={});
     tuBinVec(const char* file, const char* tag, tuOptMap options={});
     tuBinVec(TH1*, const char axis='x');
     tuBinVec(const tuBinVec&);
     /* tuBinVec(const tuBinVec&); */
     void init(vector<double>);
     void build_ptr();
     ~tuBinVec();
     double*        ptr;
     int            size;
     vector<double> vec;

     operator int ();
     operator double* ();
     operator vector<double> ();
     double operator[](int); 
     double bin_underflow(); 
     double bin_overflow();
     friend ostream& operator<<(ostream& os, tuBinVec& val);
 
     vector<double>::iterator begin();
     vector<double>::iterator end();
 };
 // small class used to see if things are in bounds
 struct tuInBounds {
     double lo_bound;
     double hi_bound;
     bool operator()(double); // return if in bounds
     void init(tuBinVec);
     tuInBounds(tuBinVec);
     tuInBounds(const char* file, const char* tag="");
     tuInBounds(double, double);
     friend ostream& operator<<(ostream& os, tuInBounds& rhs);
 };

// --------------------------------------------------------------------------------------
// tuPadDim:
//   Use:
//      Contains four numbers, used for low, low-margin,high,high-margin
//      for TPad*, for dimensions from left-to-right and bottom-to-top
//      See notes below.
// --------------------------------------------------------------------------------------



struct tuPadDim {
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
    tuPadDim( double _low, double _p_low, double _p_up, double _up );
    tuPadDim( double _low, double _up );
    tuPadDim( double _low, double _p_low, double _up );
    tuPadDim( );
    void print() const;
    double low_margin () const;
    double up_margin () const;
    bool operator==(tuPadDim& b) const; 
};
struct tuPadDimSet {
    /* tuPadDimSet(int _nPads, // make negative is wish to flip directtuns */
            /* vector<double> _lefts={}, */ 
            /* vector<double> _rights={}); */
    tuPadDimSet(
            vector<double> _lefts={},  // number of pads is first entry in _lefts if >= 1.
            vector<double> _rights={});
    vector<double> lefts;
    vector<double> rights;
    int            nPads;
    vector<tuPadDim> calc_pads();
    static tuPadDim make_pad(double left, 
            double left_margin, double pad_widht, 
            double right_margin); // 
    /* operator vector<tuPadDim> (); */
};

/* vector<tuPadDim> tuPadDimSet_(int nPads, double left_margin, double right_margin=0.005, */
/*         double overall_left=0.01, double overall_right=0.01); */
/* vector<tuPadDim> tuPadDimSet( */
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
// example:   tuPadDimSet(2,0.1) results 2 pads, both with left margin of 0.1, right margin
// of 0.01; there is 0.01 to the left and right of both of these; all remaining space is in
// the pads


// --------------------------------------------------------------------------------------
// tuPads:
//   Use:
//      Make a TCanvas with vartuus TPads
//   How:
//      Initialize:
//        Default single canvas:
//          $: tuPads pads{1};
//        Default double canvas:
//          $: tuPads pads{2};
//        Grid of pads (4 rows x 7 cols):
//          $: auto pads = tuPads{ 4, 8 };
//        Grid of pads with custom edges:
//                          {{y dim}  { ydim }       },{{x dim}      {x dim}          width height
//          $: tuPads pads{ {{.5,.98},{0,0.1,0.5,0.5}},{{0.,0.1,0.4},{0.4,.4,.9,.99}},1200,800}
//      Use pad <i>:
//          $: pads(i); // this will activate the given pad
// --------------------------------------------------------------------------------------




// inputs to these pads:
// vector of doubles: { nCol, nRow, {kLeft,kRight,kBottom,kTop, 1-4 values <1, val<1, val<1, val<1}, width, height }
struct tuPads {
    //internal data
    TCanvas* canvas = nullptr;
    vector<pair<tuPadDim,tuPadDim>> pad_dimensions;
    vector<TPad*> pads; // all the generated smaller pads
    TPad* canvas_pad;   // a single big pad the size of the canvas

    //FIXME
    /* tuPads ( vector<tuPadDim>, int canvas_width, int canv_heigth ); */
    // * Initialize with either a single vector of tuPadDim (which will go as x0, y0, x1, 
    //   y1, etc...) or two vectors of tuPadDim, which will go as {x0,x1,...} {y0,y1,...}
    // * Add operator()(int,int=0) for accessing y,x pad
    tuPads ( vector<pair<tuPadDim, tuPadDim>> pad_dimensions={{{},{}}}, 
            int canvas_width=0, int canvas_height=0);
    tuPads ( vector<tuPadDim>y_dim, vector<tuPadDim>x_dim={}, int canvas_width=0, int canv_heigth=0 );
    /* tuPads ( int nPads=1, int c_wide=0, int c_height=0 ); // default of 1 and 2 pad TPad sets */
    /* tuPads ( int nYpads=1, int nXpads=1, double y_margin=0.15, double x_margin=0.17, */ 
    /*          int c_wide=0, int c_height=0 ); // default of 1 and 2 pad TPad sets */
    /* tuPads ( int nYpads=1, int nXpads=1, int c_side=0, int c_height=0, */ 
    /*         tuPadDimSet Y={{0.15}}, tuPadDimSet X={{0.15}} ); */
    TPad*  operator()(int col=0, int row=0);

    tuPads ( int nYpads=1, vector<double> dimensions={}, int nXpads=1 );

    int nRow{1};
    int nCol{1};

    // must initialize separate from constructor (so that the user has a chance to initialize all the 
    // required options for)
    void init();

    void add_pad(pair<tuPadDim,tuPadDim>&);
    void add_pad(vector<pair<tuPadDim,tuPadDim>> input);

    int canvas_width         { 1200 };
    int canvas_height        {  800 };
    void stamp(const char*, tuOptMap opt={}, 
        tuOptMap dict = {{ 
        "TextColor", (kGray+2), 
        "TextSize", 12,
        "x-loc", .05,
        "y-loc", .05}} );

    // To do here:

};
// 
// struct tuIntSet {
//     /* tuIntSet(); */
//     /* void sort(); */
//     vector<int> list {};
//     bool operator()(int); // check if argument is in the list
//     /* bool add_data(const tuIntSet&); // untun with a second set */
//     bool has(int);
//     tuIntSet(const char* in_file, const char* tag);
//     tuIntSet(const char* in_file="", int col=0, bool print=true, bool strip_commas=true);
//     tuIntSet(const char* in_file, std::ofstream& log, int col=0, bool print=true, bool strip_commas=true);
//     std::ostringstream read_file(const char* in_file, int col=0, bool print=true, bool strip_commas=true);
//     int  operator[](int); // return location of arg in list (!Warning: does not check for existence)
//                           // warning: may not be meaningful with sort, and existence
//     tuIntSet& operator+=(const tuIntSet& rhs);
//     tuIntSet& operator-=(const tuIntSet& rhs); // remove union
//     tuIntSet& operator*=(const tuIntSet& rhs); // get the union
//     friend ostream& operator<<(ostream& os, tuIntSet& dt);
//     int size();
//     void clear();
//     void write_to_file(const char* file_name, vector<string> comments={});
// };
// 
struct tuIntList {
     vector<int> list;
     bool operator()(int); // check if argument is in the list
     bool has(int);
     bool has_not(int);
     tuIntList(const char* in_file, std::ofstream& log, bool print=true);
     tuIntList(const char* in_file, bool print=true);
     int  operator[](int); // return location of arg in list (!Warning: does not check for existence)
     private:
     string make(const char* in_file, bool print);
 };
// 
// // meant to map run-id's
// struct tuIntMap {
//     /* Reads a file and stores a map of Int -> Int
//      * Has the option to test for existence of map (like tuIntList above)
//      * Intended to be used primarily as runId -> id on histogram maps.
//      */
//     tuIntMap(const char* file_name, 
//             int index_column=0, 
//             int data_column =1,
//             bool echo_print=true,
//             vector<int> skip_vals={}
//     );
//     tuIntMap(const char* file_name, 
//             int index_column,
//             int data_column,
//             bool echo_print,
//             std::ofstream& log,
//             vector<int> skip_vals={}
//     );
//     private:
//     string tuIntMap_constructor(const char* file_name, 
//             int index_column, 
//             int data_column,
//             bool echo_print,
//             vector<int> skip_vals={}
//     );
//     public:
//     map<int,int> data_map; // runid -> duration
//     bool         has(int);
//     vector<int>  keys(); // returns a sorted vector of the keys present
//     int&         operator[](int key); // get data_map[key]; if false, returns 0. as default
//     int          size();   // size of map
// };
// 
// 
// struct tuRunListId {
//     public:
//     tuRunListId(const char* file_name, bool skip_127053_138064=true);
//     map<int,double> map_id;  // runid -> map_id
//     map<int,int>    run_sec; // runid -> duration
//     bool has_run(int);
//     double set_id(int run_id);
//     double id; // last set id
//     int  size();
//     double operator()(int);  // return map_id[runid], unless not present, then -1.
// };
// 
// struct tuHgStats {
//     // used to get basic statistics and cuts
//     double mean;
//     double stddev;
//     void calc_stats(); // just the 
// 
//     TAxis* axis;
// 
//     /* vector<bool>   cut ; // points that are cut */
//     vector<double> unmasked_vals();
//     vector<double> vals;
//     vector<double> errs;
//     vector<double> weight;
//     int nbins;
// 
//     TLine* get_horizontal_TLine(double cut, bool cut_times_sigma=false);
// 
//     TGraph* points_above(double cut, bool cut_times_sigma=false);
//     TGraph* points_below(double cut, bool cut_times_sigma=false);
//     TGraph* points_between(double cut0, double cut1, bool cut_times_sigma=false);
//     TGraph* points();
// 
//     double mean_Xsigma(double X); // return StdDev * X + mean
// 
//     tuHgStats& cut_above(double cut, bool cut_times_sigma=false);
//     tuHgStats& cut_below(double cut, bool cut_times_sigma=false);
//     tuHgStats& cut_zeros();
//     tuHgStats& cut_to_range(double cut0, double cut1, bool cut_times_sigma=false);
//     tuHgStats& mask(const vector<bool>  mask); // mask all points in vector which are TRUE
//     /* tuHgStats& mask(vector<int>&  mask); */
//     tuHgStats& restore_points();
// 
//     tuHgStats(TH1D* hg, bool cut_zeros=false);
//     tuHgStats(TProfile* hg, bool cut_zeros=false, bool weight_by_entries=false);
//     tuHgStats(vector<double> x_vals, vector<double> vals, vector<double> errs, bool cut_zeros = false);
// 
//     int count_points();
//     vector<int> bin_indices();
// };
// 
// // tuMsgTree:
// // Used to write messages to a MsgTree in the local file, and read from it
// class tuMsgTree {
//     private:
//         string b_msg;
//         TTree tree;
//     public:
//         tuMsgTree(bool set_echo=true);
//         static void read_messages(const char* f_name);
//         void msg(string msg);                  // write a message
//         void msg(vector<string> messages);     // *ditto*
//         void write(); // write to tree
//         void dash();  // write dashes to tree
//         bool echo_to_cout {false};
//         void slurp_file(const char* which_file); // write a given file to input
// };
// 
// struct tuIntVec {
//     // A class that has at it's heart a vector<vector<int>> of data
//     // The first row of data is are the keys of the data, for a kind of data frame table
//     // This is a poor man's implementation
//     tuIntVec( const char* file_name, bool echo_print=true, vector<string>_tags={} );
//     tuIntVec( const char* file_name, std::ofstream& log, bool echo_print=true, vector<string>_tags={});
//     string tuIntVec_constructor( const char* file_name, bool echo_print, vector<string>_tags={} );
// 
// 
//     // data members
//     /* map<int,vector<int>> data_map {}; // runid -> vector<bool> */
// 
//     vector<int> keys;
//     vector<vector<int>> data;
// 
//     vector<string> tags {};   // names of all the columns
// 
//     // access data and manipulate data
//     int             i_tag(string tag); // returns column value of tag
//     vector<int>     tag_cols(vector<string> _tags, const char* err_msg_name="");
//     bool            has_tag(string tag);
//     tuIntVec&       swap_tags(string tag0, string tag1); // swap column locations
//     tuIntVec&       subset(vector<string> tag);
//     tuIntVec&       rm_tags(vector<string> tag);
//     tuIntVec&       add_tag(string tag, int default_val=-1); 
//     void            rename_tag(string, string);
//     bool            has_key(int key);     // checks if it has key
//     vector<int>&    operator[](int key);   // returns the vector at entry
//     int             size();   // size of map
//     vector<int>     vals(string tag); // column of data under "tag"
//     vector<int>     vals(string tag, vector<bool> mask);// column of data under "tag" where mask is True
//     vector<int>     get_keys(vector<bool> mask, bool keep_on_true=true); // keep all values
//     vector<bool>    mask(string tag,const char* module="mask");
//     vector<bool>    is_any(vector<string> tag);
//     vector<bool>    is_all(vector<string> tag);
// 
//     private:
//     bool flip_tag(string& tag); // tag starts with "!" then strip it and return false, else return true;
//     void assert_tag(string tag, const char* routine=""); // fail is tag not in tags
// 
//     public:
//     friend ostream& operator<<(ostream& os, tuIntVec& dt);
//     void            write_to_file(const char* file_name, vector<string> comments={});
// };
// 
// struct tuIntBinCnt {
//     // used to keep count for inputs are discreet integer locations
//     TH1D* hg1;
//     TH2D* hg2;
//     bool is2D{false};
//     void fill(int i, double weight=1.);
//     void fill(int i, int j, double weight=1.);
//     static double getcnt(TH1D* hg, int i);
//     static double getcnt(TH2D* hg, int i, int j);
//     void write();
//     tuIntBinCnt(const char* name, vector<int> x_dim, const char* title);
//     tuIntBinCnt(const char* name, vector<int> x_dim, vector<int> y_dim={}, const char* title="");
// };
// 
// struct tuFnCaller{
//     tuFnCaller(const char* file_data, double(&_fn)(double*,double*));
//     double (&fn) (double*, double*);
//     double operator()(double,double);
//     double *p; // fit parameters, read from data_file
//     double *x;
// };
// 
// struct tuXsec{
//     // generalization of src/tuXsec_pAu2015.h
//     // ptHat bins, X-sections, and nEvents are read in from file
// 
//     // data
//     tuXsec (
//         const char* tag_file,
//         const char* Xsection_tag="Xsection",
//         const char* pthatbin_tag="pThat",
//         const char* nEvents_tag="nEvents"
//     );
//     vector<double> Xsection;
//     vector<double> pthatbins;
//     double pthatbin_center(int);
//     int nbins_pthat;
//     vector<int> Nevents;
//     vector<int> Ncollected;
//     long int n_collected_total{0};
// 
//     // functions
//     int pthatbin(pair<double,double> pthat_bounds);
//     
//     double Xsec(pair<double,double>, int number_of_events=0); 
//     double Xsec(int pthatbin, int number_of_events=0); 
//     double operator()(int pthatbin, int num_of_events=0) 
//         { return Xsec(pthatbin, num_of_events); };
//     double operator()(pair<double,double> pthatbounds, int num_of_events=0) 
//         { return Xsec(pthatbounds, num_of_events); };
// 
//     void collect (int pthatbin);
//     void collect (pair<double,double> pthatbounds);
//     void check_pthatbin(int bin);
// 
//     TH1D* hg_collected;
// };
// 
// struct tuIntStrFunctor {
//     tuIntStrFunctor ( const char* file, tuOptMap options={}, 
//         tuOptMap dict= {{"tag","none","sort",false, "strip_commas",false,"column","all"}});
//     const char* operator()(int);
//     map<int,string> data;
// };
// struct tuStrStrFunctor {
//     tuStrStrFunctor ( const char* file, tuOptMap options={}, 
//         tuOptMap dict= {{"tag","none","sort",false, "strip_commas",false,"column","all"}});
//     const char* operator()(const char*);
//     map<string,string> data;
// };
// 
// struct tuFirst {
//     bool is_first{true};
//     tuFirst(){};
//     operator bool();
//     bool operator()();
// };
// 
// struct tuCycleTrue {
//     int period;
//     int cnt;
//     tuCycleTrue(int period_in);
//     bool operator()();
//     operator bool();
//     void reset();
// };
// 
// struct tuCycleSpacer {
//     tuCycleTrue cycle;
//     string spacer;
//     int n_width;
//     void reset();
//     tuCycleSpacer(int period, int n_width=0, const char* def_spacer=" ");
//     friend ostream& operator<<(ostream& os, tuCycleSpacer& self);
// };
// 
// struct tuXYbounder {
//     // put in two vectors: (sorted) X & (not-sorted) Y
//     // then when asking about input x, it will find bin-X and return bound-Y
//     vector<double> X;
//     vector<double> Y;
//     int size;
//     double lodef;  // default value if x not found within X
//     double hidef;  //
//     bool operator()(double x,double y); // returns if val Y is out of bounds at X
//     double operator()(double x); // returns boundary value y at x
//     tuXYbounder(const char* file, const char* tagX, const char* tagY, 
//             tuOptMap options={}); // options are hidef for default-hi and default-lo
//                                  // if not set will just use the first and last values in 
//                                  // Y for under and overflow
//     tuXYbounder(vector<double> x, vector <double> y, tuOptMap options={}); 
//     // options are hidef for default-hi and default-lo
//     tuXYbounder(); // default to always returning false for no bounds set
// };
// 
// 
// 
// /* TGraphAsymmErrors* tuMakeTGASE(TH1D* hg, bool invert, array<double,4> x_center, bool skip_zeros, bool normalize ) { */
// /*     vector<double> x, x_err, y, y_err; */
// /*     TAxis* axis = hg->GetXaxis(); */
// /*     if (normalize) hg->Scale(1./hg->Integral()); */
// /*     for (int i{1}; i<=hg->GetXaxis()->GetNbins(); ++i) { */
// /*         if (hg->GetBinContent(i) == 0 && hg->GetBinError(i) == 0 && skip_zeros) continue; */
// /*         x.push_back(axis->GetBinCenter(i)); */
// /*         x_err.push_back(0.); */
// /*         y.push_back(hg->GetBinContent(i)); */
// /*         y_err.push_back(hg->GetBinError(i)); */
//         
// /*     } */
// /*     TGraphErrors* gr; */
// /*     double lo = hg->GetXaxis()->GetBinLowEdge(1); */
// /*     double hi = hg->GetXaxis()->GetBinUpEdge( hg->GetXaxis()->GetNbins()); */
// /*     if (invert) { */
// /*         gr = tuMakeTGraphErrors(y,x,x_err,y_err); */
// /*         gr->SetMinimum(lo); */
// /*         gr->SetMaximum(hi); */
// /*     } else { */
// /*         gr = tuMakeTGraphErrors(x,y,y_err,x_err); */
// /*         gr->GetXaxis()->SetLimits(lo,hi); */
// /*     } */
// /*     return gr; */
// /* }; */
// 
// /* TGraphAsymmErrors* tuMakeTGASE(int n, tuBinVec x, tuBinVec y, tuBinVec exleft, tuBinVec exright, */
// /*          tuBinVec eyleft, tuBinVec eyright={}); */
// 
// /* TGraphAsymmErrors* tuMakeTGASE(TH1D* hg, bool invert_XY=false, */ 
// /*         array<double,4> x_center={0.,1.,-1.}, // left/right/center ratio of bin; */ 
// /*                                               // center=-1 or < left defaults to middle of left */ 
// /*                                               // and right */
// /*         bool skip_zeros=true); */
// /* /1* TGraph* tuMakeTGraph(vector<double> x, vector<double> y); *1/ */
// 
// struct tuPtrDbl {
//     tuPtrDbl(TAxis* ax, double bin_loc=0.5, bool get_widths=false);  //<-
//     tuPtrDbl(vector<double>); // <-
//     tuPtrDbl(TH1*, bool get_errors=false); // also Errors // <-
//     tuPtrDbl(const tuPtrDbl&); // <-
//     /* tuPtrDbl(tuPtrDbl); // <- */
//     tuPtrDbl(const char* file, const char* tag); // <-
//     tuPtrDbl(int i); // <-
//     tuPtrDbl(){};
// 
//     // internal
//     vector<double> vec{};
//     double*        ptr{nullptr};
//     int            size{0};
// 
//     /* void init(vector<double>, bool range_double=true); */
//     void build_ptr();
//     tuPtrDbl& update();
//     // read the vec from a file with tuReadValVec
//     ~tuPtrDbl();
//     /* void set_val(int i, double val); */
// 
//     /* int nbins(); // return size_ptr-1 */
//     operator int ();
//     operator double* ();
//     operator vector<double> ();
//     double& operator[](int); 
//     friend ostream& operator<<(ostream& os, tuPtrDbl& val);
//     int throw_error(const char* msg="");
// 
//     tuPtrDbl& operator +=(const tuPtrDbl& rhs);
//     tuPtrDbl& operator -=(const tuPtrDbl& rhs);
//     tuPtrDbl& operator *=(const tuPtrDbl& rhs);
//     tuPtrDbl& operator /=(const tuPtrDbl& rhs);
//     tuPtrDbl& operator /=(const double rhs);
//     tuPtrDbl& operator *=(const double rhs);
//     tuPtrDbl& operator +=(const double rhs);
//     tuPtrDbl& operator -=(const double rhs);
// 
//     tuPtrDbl& make_min(const tuPtrDbl& rhs); // adjust to have just the min values
//     tuPtrDbl& make_max(const tuPtrDbl& rhs); // adjust to have just the max values
//     tuPtrDbl& abs_diff(const tuPtrDbl& rhs); // adjust to have just the max values
//     tuPtrDbl& square_diff(const tuPtrDbl& rhs); // adjust to have just the max values
//     tuPtrDbl& abs();
//     tuPtrDbl& zero_negatives();
//     tuPtrDbl& sqrt();
// 
//     vector<double>::iterator begin();
//     vector<double>::iterator end();
//     tuPtrDbl& operator=(const tuPtrDbl&);
// 
// };
// tuPtrDbl  operator+(const tuPtrDbl& lhs, const tuPtrDbl& rhs);
// tuPtrDbl  operator-(const tuPtrDbl& lhs, const tuPtrDbl& rhs);
// tuPtrDbl  tu_calc_quadrature(vector<tuPtrDbl>);
// tuPtrDbl  tu_calc_quadrature(vector<tuPtrDbl>, tuPtrDbl);
// tuPtrDbl  tu_calc_mean      (vector<tuPtrDbl>);
// tuPtrDbl  tu_calc_max_bound (vector<tuPtrDbl>);
// tuPtrDbl  tu_calc_min_bound (vector<tuPtrDbl>);
// tuPtrDbl  tu_calc_max_berr  (vector<tuPtrDbl>, tuPtrDbl);
// tuPtrDbl  tu_calc_min_berr  (vector<tuPtrDbl>, tuPtrDbl);
// pair<tuPtrDbl,tuPtrDbl>  tu_calc_bounds    (vector<tuPtrDbl>);
// 
// // some drawing options
// TGraphAsymmErrors* tu_draw_error_boxes(TH1D* mean, tuPtrDbl, tuOptMap opts, 
//         array<double,4> = {-1,-1}, double x_range_lo=0., double x_range_hi=0.);
// TGraphAsymmErrors* tu_draw_error_boxes(TH1D* mean, array<tuPtrDbl,2>, 
//         tuOptMap opts, array<double,4> = {-1,-1}, double x_range_lo=0., double x_range_hi=0.);
// 
// struct tuSysErrors {
//     /* tuSysErrors(); */
//     tuSysErrors(const tuSysErrors&);
//     tuSysErrors(TH1*, array<double,4> x_rat={-1,-1,0.5}, array<tuPtrDbl,2>_err={});
//     tuSysErrors(TGraphAsymmErrors*, array<double,4> x_rat={-1,-1,0.5});
//     tuSysErrors& set_x_range(double, double);
//     tuSysErrors& swap_xy ();
// 
//     tuSysErrors& setYlow(tuPtrDbl&);
//     tuSysErrors& setYhigh(tuPtrDbl&);
//     tuSysErrors& setYhilo(tuPtrDbl&); // set Y symmetric
// 
//     /* void draw_boxes(tuOptMap dict={}); */
// 
//     /* tuSysErrors divide(const tuSysErrors& other); */
// 
//     /* tuSysErrors& set (TH1*); // sets x,y, and errors */
//     /* tuSysErrors& set (TGraphAsymmErrors*); // sets x,y, and errors */
// 
//     int size{0};
//     TGraphAsymmErrors* tgase {nullptr};
//     tuPtrDbl getY();
//     tuPtrDbl getYlow();
//     tuPtrDbl getYhigh();
// 
//     /* vector<tuPtrDbl> vec_data   {}; // to add in quadrature */
//     /* tuSysErrors& add_data   (tuPtrDbl); */
//     /* tuSysErrors& add_data   (vector<tuPtrDbl>); */
// 
//     /* tuSysErrors& calc_mean(vector<tuPtrDbl> data={}); */
//     /* tuSysErrors& calc_quadrature(vector<tuPtrDbl> data ={}); // calculate quadratue relative to the mean */
//     /* tuSysErrors& calc_bounds(vector<tuPtrDbl> data ={});     // calculate bounds relative to mean */
//     /* tuSysErrors& calc_symmetric_bounds(vector<tuPtrDbl> data ={}); */
// 
//     tuSysErrors& set_rat_xbins(array<double,4> rat_rel);
//     /* tuSysErrors& set_center_xbins(double rat_center); */
//     operator TGraphAsymmErrors* (); // case it to TGraphAsymmErrors
//     TGraphAsymmErrors* operator-> ();
// };
// 
// struct tuJetMatcher_float {
//     // For filling vectors for Raghav matching algorithm 
//     // will have comparison methods in order to std::sort
//     float eta;
//     float phi;
//     float pT;
//     bool  is_matched;
// 	double operator()(const tuJetMatcher_float&, const float jetR2=0.16);
//     tuJetMatcher_float(float,float,float);
//     tuJetMatcher_float();
// };
// bool operator==(const tuJetMatcher_float& L, const tuJetMatcher_float R);
// 
// bool operator<(const tuJetMatcher_float& L, const tuJetMatcher_float R) ;
// 
// bool operator>(const tuJetMatcher_float& L, const tuJetMatcher_float R) ;
// 
// struct tu_TF1fitter {
//     string name;
//     string fn_str;
//     TF1* fn;
//     unsigned int n_par;
//     double lo, hi;
// 
//     void fix_match_params(tu_TF1fitter& f1_to_match, std::set<int> which_pars);
//     // fix parameters to the values found in f1_to_match, set others free
//     
//     vector<double> fval;
//     vector<double> ferr;
//     vector<double>& fit(TH1D*, double lo=0, double hi=0);
// 
//     vector<double>&     operator()(TH1D*, double lo=0, double hi=0);
//     vector<double>&     operator()();
//     pair<double,double> operator()(int i); // return <fval, ferr> of parameter i
//     double              operator[](int i); // return fval of parameter i
// 
//     tu_TF1fitter(const char* fnc_str, const char* name="",TH1D* hg=nullptr, 
//             double lo=0, double hi=0);
//     void set_presets(vector<double>);
//     friend ostream& operator<<(ostream& os, tu_TF1fitter& ft);
// };
// 
// struct tuParticleThrower {
//     TRandom3 r3;
//     vector<unsigned int>    whole;
//     vector<double> remainder;
//     vector<double> l_bound;
//     vector<double> u_bound;
//     tuParticleThrower( TH1D* _hg_probs, double multiple=1., unsigned int _seed=0 );
//     friend ostream& operator<<(ostream& os, tuBinVec& tuParticleThrower);
// 
//     bool throw_particle();
//     void set_rand(unsigned int);
//     bool is_thrown;
//     unsigned int i_bin{0};
//     unsigned int i_whole{0};
//     unsigned int i_size;
//     double pt;
//     double phi;
//     double eta;
// };
// 
// struct tuPoissonParticleThrower {
//     tuPoissonParticleThrower(TRandom3& _0, TH1D& _1, double mult=1.);
//     // Give it a distribution of particle pT values
//     // It will then integrate for the mean number of particles, and throw
//     // that number, distributed from the distribution
//     
//     // public
//     TRandom3& r3;
//     TH1D&     dist;
//     bool      throw_particle(); // will keep returning true until the last particle has been thrown.
//     double pt;
//     double phi;
//     double eta;
//     double mean;
// 
//     // internal
//     unsigned int  i_to_throw;
//     unsigned int  i_thrown   {0};
// };
// 
// struct tuHopper1D {
//     vector<double> edges;
//     vector<pair<double,double>> hopper;
//     void reset();
//     int nbins;
//     tuHopper1D(vector<double> _edges);
//     vector<pair<double,double>>::iterator begin();
//     vector<pair<double,double>>::iterator end();
//     int fill(double val, double weight=1.);
// };
// 

#endif
