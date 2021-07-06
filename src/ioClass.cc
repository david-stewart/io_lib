#include "ioClass.h"
#include "io_fmt.h"
#include "io_fnc.h"

#include <sstream>
#include <algorithm>
#include "TString.h" // not used, but present for some header 
#include "TMath.h"
#include <iostream>
#include <string>
#include <set>
#include "TTreeReader.h"

using namespace std;

// ioGetter Class
ioGetter::ioGetter(string _) : path{_}, n_objects{0} {};

TFile* ioGetter::get_file(string _f_name) {
    string f_name = path + _f_name;
    if (f_name.find(".root") == string::npos) f_name += ".root";
    TFile *f = ( files.count(f_name) ? files[f_name] : new TFile(f_name.c_str(), "read") );
    if (!f->IsOpen()) {
        cout << " Fatal error, cannot open file: " << f_name << endl;
        exit (1);
    }
    return f;
};

TObject* ioGetter::operator()(string f_name, string object_name) {
    TFile *f = get_file(f_name);
    TObject* obj;
    f->GetObject(object_name.c_str(), obj);
    if (obj == nullptr) {
        cout << " !fatal: failed to get object " << object_name 
            << " from file " << f_name << endl;
        exit(2);
    }
    ++n_objects;
    /* obj->SetName(unique_name()); */
    if (got.count(object_name) != 0) got[object_name].push_back(obj);
    else got[object_name].push_back({obj});
    return obj;
};

ioRanger::ioRanger(double _lo_range, double _hi_range, double _lo_out, double _hi_out) :
    lo_range {_lo_range},
    hi_range {_hi_range},
    lo_out   {_lo_out},
    hi_out   {_hi_out}
{};
double ioRanger::operator()(double x) {
    return lo_out + (hi_out-lo_out)*(x-lo_range)/(hi_range-lo_range);
};

vector<double> ioBinVec::bin_centers() {
    vector<double> V;
    for (int i{0}; i<(int)vec.size()-1; ++i) V.push_back(0.5*(vec[i]+vec[i+1]));
    return V;
};
ioBinVec::ioBinVec(TAxis* ax) {
    vector<double> build_vec;
    for (int i{1}; i<=ax->GetNbins()+1; ++i) build_vec.push_back(ax->GetBinLowEdge(i));
    init(build_vec);
};
ioBinVec::ioBinVec(int nbins, double lo, double hi) :
    size {nbins} 
{
    double step = (hi-lo)/nbins;
    for (int i{1}; i<=size; ++i) vec.push_back(lo+i*step);
    build_ptr();
};
void ioBinVec::build_ptr() {
    size = vec.size();
    ptr = new double[size];
    for (int i{0}; i<size; ++i) ptr[i] = vec[i];
};
ioBinVec::operator int () { return size-1; };
ioBinVec::operator double* () { return ptr; };
ioBinVec::operator vector<double> () { return vec; };

ioBinVec::ioBinVec(vector<double> V, bool range_setter) { 
    init(V, range_setter);
};
ioBinVec::ioBinVec(const char* file, ioOptMap options, bool nbin_range){
    init( ioReadValVec(file, options), nbin_range );
};
ioBinVec::ioBinVec(const char* file, const char* tag, ioOptMap options, bool nbin_range){
    options["tag"] = tag;
    init( ioReadValVec(file, options), nbin_range );
};
ioBinVec::ioBinVec(vector<vector<double>> V_in) {
    for (auto& VEC : V_in)
    for (auto    v : VEC) 
        vec.push_back(v);
    build_ptr();
};
// copy constructor
ioBinVec::ioBinVec(const ioBinVec& cp) {
    init(cp.vec, true);
};
ioBinVec::ioBinVec(TH1* h, const char axis) {
    TAxis* ax;
    switch (axis) {
        case 'x':
        case 'X':
            ax = h->GetXaxis();
            break;
        case 'y':
        case 'Y':
            ax = h->GetYaxis();
            break;
        case 'z':
        case 'Z':
            ax = h->GetZaxis();
            break;
        default:
            throw std::runtime_error("in ioBinVec initializer, must select axis 'xXyYzZ'");
    }
    vector<double> build_vec;
    for (int i{1}; i<=ax->GetNbins()+1; ++i) build_vec.push_back(ax->GetBinLowEdge(i));
    init(build_vec);
};

void ioBinVec::set_val(int i, double val) {
    if (i >= size) throw std::runtime_error(
    Form("fatal in ioBinVec::set_val(), trying to change entry %i of vector size %i",
        i, size)
    );
    ptr[i] = val;
    vec[i] = val;
};

/* ioBinVec ioBinVec::operator+=(const ioBinVec& _) { */
    /* for (auto v : _.vec) vec.push_back(v); */
    /* delete[] ptr; */
    /* build_ptr(); */
    /* return *this; */
/* }; */
/* ioBinVec operator+ (ioBinVec lhs, const ioBinVec& rhs) { lhs += rhs; return lhs; }; */

void ioBinVec::init(vector<double> V, bool range_repeat) {
    if (!range_repeat || V.size()==0) {
        for (auto v : V) vec.push_back(v);
        build_ptr();
        return;
    }
    // range repeate will add a range leading to the next number
    // it is triggered by a repeat value of the last number, followed by the number
    // of bins
    //   example:
    //         0, 0, 5, 1. 2. 3. = 0 .2 .4 .6 .8 1.0 2. 3.
    vec.push_back(V[0]);
    int S = V.size();
    int i{1}; 
    while (i<(int)S) {
        if ( V[i] == V[i-1] ) {
            if (i>(S-3)) throw std::runtime_error( "fatal in ioBinVec with range_repeat");
            double step = (V[i+2]-V[i])/V[i+1];
            for (int k{1}; k<=V[i+1]; ++k) vec.push_back(V[i]+k*step);
            i+=3;
        } else {
            vec.push_back(V[i]);
            ++i;
        }
    }
    build_ptr();
};
ioBinVec::~ioBinVec() {
    delete[] ptr;
};
/* int ioBinVec::nbins() { return (int) size-1; }; */
vector<double>::iterator ioBinVec::begin() { return vec.begin(); };
vector<double>::iterator ioBinVec::end()   { return vec.end(); };
double ioBinVec::operator[](int i) { return vec[i]; };
double ioBinVec::bin_underflow() { 
    if (vec.size()<2)  return 0.;
    return vec[0]-(vec[1]-vec[0]);
};
double ioBinVec::bin_overflow() { 
    if (vec.size()<2)  return 0.;
    int i { static_cast<int>(vec.size())-1 };
    return vec[i]+(vec[i]-vec[i-1]);
};
ostream& operator<<(ostream& os, ioBinVec& io) {
    for (auto v : io) cout << " " << v;
    cout << endl;
    return os;
};

bool ioInBounds::operator()(double x) {
    return (x >= lo_bound && x <= hi_bound);
};
void ioInBounds::init(ioBinVec bins) {
    if (bins.vec.size()<2) {
        throw std::runtime_error(
        "Fatal: tried to initialize a ioInBounds with "
        "an input ioBinVec with less than 2 entries!"
        );
    };
    lo_bound = bins[0];
    hi_bound = bins[bins.vec.size()-1];
};
ioInBounds::ioInBounds(const char* file, const char* tag){
    init( {file, tag} );
};
ioInBounds::ioInBounds(ioBinVec bins) { init(bins); };
ioInBounds::ioInBounds(double lo, double hi) :
    lo_bound{lo}, hi_bound(hi) {}; 

ostream& operator<<(ostream& os, ioInBounds& rhs) {
    os << " ioInBounds("<<rhs.lo_bound<<","<<rhs.hi_bound<<")"<<endl;
    return os;
};

// ioPads class (with helper class ioPadDim)
void ioPadDim::check_input() {
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

ioPadDim::ioPadDim( double _low, double _p_low, double _p_up, double _up ) :
    low{_low}, p_low{_p_low}, p_up{_p_up}, up{_up} { check_input(); };
ioPadDim::ioPadDim( double _low, double _up ) : 
    low{_low}, p_low{_low}, p_up{_up}, up{_up} { check_input(); };
ioPadDim::ioPadDim( double _low, double _p_low, double _up ) : 
    low{_low}, p_low{_p_low}, p_up{_up}, up{_up} { check_input(); };
ioPadDim::ioPadDim( ) :
    low{0.}, p_low{0.}, p_up{1.}, up{1.} { check_input(); };


void ioPadDim::print() const {
    cout << Form(" Four points are: (%.2f, %.2f %.2f, %.2f)",low,p_low,p_up,up) << endl;
};

double ioPadDim::low_margin () const {
    double margin { (p_low - low) / (up - low) };
    if (margin < 0) margin = 0;
    return margin;
};
double ioPadDim::up_margin () const {
    // use to get set the lower margin
    double margin { (up - p_up) / (up - low) };
    if (margin < 0) margin = 0;
    return margin;
};

bool ioPadDim::operator==(ioPadDim& B) const {
    return low == B.low
        && p_low == B.p_low
        && p_up  == B.p_up
        && up    == B.up;
};
vector<ioPadDim> ioPadDimSet_(int nPads, double left_margin, double right_margin,
        double overall_left, double overall_right){
    vector<ioPadDim> vec;
    double space = {(1.-(left_margin+right_margin)*nPads-overall_left-overall_right)/nPads};
    if (space<=0) throw std::runtime_error(
            "fatal in ioPadDimSet: margins have consumed more than 100% of TCanvas");
    double left = overall_left;
    for (int i{0}; i<nPads; ++i) {
        /* cout << "i: " << i << endl; */
        vec.push_back( {left,left+left_margin, 
                left+left_margin+space, left+left_margin+space+right_margin});
        left += left_margin+space+right_margin;
    }
    return vec;
};

vector<ioPadDim> ioPadDimSet(
        int nPads, 
        double leftM,  // first left margin
        bool b_reverse, // swap the ordering (for y it's nice to go from top to bottom)
        double rightM, // last right margin
        double left_margin, // how far in on the canvas pad
        double right_margin, // how far in on the canvas pad
        double left_in, // inner left margins
        double right_in // inner right margins
        ){
    vector<ioPadDim> vec;
    double space = {(1.-(left_in+right_in)*(nPads-1)
            -left_margin-right_margin-leftM-rightM)/nPads};
    /* double space = {(1.-(left_in+right_in)*(nPads-1) */
    /* -left_margin-right_margin-left_in-right_in)/nPads}; */
    if (space<=0) throw std::runtime_error(
            "fatal in ioPadDimSet: margins have consumed more than 100% of TCanvas");

    if (nPads == 1) {
        vec.push_back( {
                left_margin, 
                left_margin+leftM, 
                left_margin+leftM+space, 
                left_margin+leftM+space+right_margin});
        return vec;
    }
    vec.push_back( {
            left_margin, 
            left_margin+leftM, 
            left_margin+leftM+space, 
            left_margin+leftM+space+right_in});
    double left = left_margin+leftM+space+right_in;
    for (int i{1}; i<nPads-1; ++i) {
        vec.push_back( {left,left+left_in,left+left_in+space,left+left_in+space+right_in} );
        left=left+left_in+space+right_in;
    }
    vec.push_back( {left,left+left_in,left+left_in+space,left+left_in+space+right_margin} );

    if (b_reverse) reverse(vec.begin(), vec.end());
    return vec;
};

ioPads::ioPads ( vector<pair<ioPadDim, ioPadDim>> _pad_dimensions, int
        _canvas_width, int _canvas_height) :
    pad_dimensions{ _pad_dimensions }
{
    if (_canvas_width)  canvas_width  = _canvas_width;
    if (_canvas_height) canvas_height = _canvas_height;
    /* init(); */
};
ioPads::ioPads( vector<ioPadDim> _pad_dim, int c_wide, int c_height) {
    if (_pad_dim.size() % 2 != 0) {
        cout << "Error in constructor ioPads(vector<ioPadDim>...) :" << endl;
        cout << "   vector<ioPadDim> does *not* have even number of elements" << endl;
        throw std::runtime_error(" fatal: Odd number of constructor elements.");
    }
    int i{0};
    while (i< (int)_pad_dim.size()-2) {
        pad_dimensions.push_back({_pad_dim[i],_pad_dim[i+1]});
        i -= 2;
    };
    if (c_wide) canvas_width = c_wide;
    if (c_height) canvas_height = c_height;
    /* init(); */
};
/* ioPads::ioPads( vector<ioPadDim> y_dim, vector<ioPadDim> x_dim, int c_wide, int c_height) { */
ioPads::ioPads ( int nYpads, int nXpads, double y_margin, double x_margin, 
        int c_wide, int c_height ){ // default of 1 and 2 pad TPad sets
    // make a set of pads with y_margin on bottom, x_margin on left, 0.05 on top and right,
    // and solid packed between
    vector<ioPadDim> y_dim{};
    vector<ioPadDim> x_dim{};

    if (nYpads==1) y_dim.push_back({0.,y_margin,0.95,0.99});
    else {
        double space { (1.-0.05-y_margin)/nYpads };
        y_dim.push_back({0.95-space,0.95-space,0.95,0.99});
        for (int i=1;i<nYpads-1;++i) {
            double top { 0.95 - i*space };
            double bottom { top - space };
            y_dim.push_back({bottom,bottom,top,top});
        }
        y_dim.push_back({0.,y_margin,y_margin+space,y_margin+space});
    }

    if (nXpads==1) x_dim.push_back({0.,x_margin,0.95,0.99});
    else {
        double space { (1.-0.05-x_margin)/nXpads };
        x_dim.push_back({0.,x_margin,x_margin+space,x_margin+space});
        for (int i=1;i<nXpads-1;++i) {
            double left  { x_margin + i*space };
            double right { left + space };
            x_dim.push_back( {left,left,right,right} );
        }
        x_dim.push_back({0.95-space, 0.95-space, 0.95, 0.99});
    }
    for (auto& x : x_dim) 
        for (auto& y : y_dim)
            pad_dimensions.push_back({y,x});

    nRow = y_dim.size();
    if (c_wide) canvas_width = c_wide;
    if (c_height) canvas_height = c_height;
    /* init(); */
};
ioPads::ioPads( vector<ioPadDim> y_dim, vector<ioPadDim> x_dim, int c_wide, int c_height) {
    if (x_dim.size()==1) {
        ioPadDim temp_pad;
        if (x_dim[0] == temp_pad) {
            x_dim[0].low = 0.;
            x_dim[0].p_low = 0.15;
            x_dim[0].p_up  = 0.95;
            x_dim[0].up    = 0.99;
        }
    }
    if (y_dim[0].low < 0) {
        y_dim[0].low = TMath::Abs(y_dim[0].low);
        if (y_dim[0].low < 0.01) y_dim[0].low = 0.;
        vector<ioPadDim> temp;
        for (int i{(int)y_dim.size()-1};i>=0;--i) temp.push_back(y_dim[i]);
        y_dim = temp;
    }
    for (auto& x : x_dim) 
        for (auto& y : y_dim)
            pad_dimensions.push_back({y,x});

    nRow = x_dim.size();
    if (c_wide)   canvas_width = c_wide;
    if (c_height) canvas_height = c_height;
    /* init(); */
};
void ioPads::stamp(const char* msg, ioOptMap options, ioOptMap dict) {
    dict += options;
    canvas_pad->cd();
    /* cout << " x: " << dict["x-loc"] << "  " << dict["y-loc"] << endl; */
    ioDrawTLatex(msg,dict["x-loc"](), dict["y-loc"](), dict);
};
/* ioPads::ioPads(int nPads, int c_wide, int c_high){ */
/*     if (nPads==2) { */
/*         pad_dimensions.push_back( {{0.55,0.55,0.95,0.99},{0.,0.15,0.9,0.99}} ); */
/*         pad_dimensions.push_back( {{0.00,0.15,0.55     },{0.,0.15,0.9,0.99}} ); */
/*         canvas_width  = 800; */
/*         canvas_height = 800; */
/*     } else if (nPads ==1) { */
/*         pad_dimensions.push_back( {{0.00,0.15,0.95,0.99},{0.,0.17,0.9,0.99}} ); */
/*         canvas_width  = 800; */
/*         canvas_height = 800; */
/*     } else { */
/*         throw std::runtime_error(" fatal: Called ioPads(int nPads...) with nPads > 2"); */
/*     } */
/*     /1* init(); *1/ */
/* }; */
TPad* ioPads::operator()(int row, int col) {
    if (pads.size() == 0) init();
    int i_pad = row+col*nRow;
    if (i_pad >= (int)pads.size()) {
        cout << " warning! asking for pad " << i_pad << " in vector of " << pads.size() << "!" << endl;
        cout << "   returning pad[0] instead" << endl;
        i_pad = 0;
    }
    pads[i_pad]->cd();
    return pads[i_pad];
};
/* vector<pair<ioPadDim,ioPadDim>> ioPadDimGrid( vector<ioPadDim> x_coord, */
/*         vector<ioPadDim> y_coord, bool by_rows) { */
/*     // make a grid of inputs (pairs of  ioPadDim) from x_coord[x0, x1, ... xn] and y_coord[y0, y1, ... yn] to */
/*     // <x0,y0>, <x1,y0> ... <xn,y0> */
/*     // <x0,y1>, <x1,y1> ... <xn,y1> */
/*     //       ... */
/*     // <x0,yn>, <x1,yn> ... <xn,yn> */
/*     vector<pair<ioPadDim, ioPadDim>> vec_fourPoint{}; */
/*     if (by_rows) { */ 
/*         // go across by rows and then down by columns */
/*         for (const auto& y : y_coord) { */
/*             for (const auto& x : x_coord) { */
/*                 vec_fourPoint.push_back({x,y}); */
/*             } */
/*         } */
/*     } else { */
/*         // go down first column, then next, etc... */
/*         for (const auto& x : x_coord) { */
/*             for (const auto& y : y_coord) { */
/*                 vec_fourPoint.push_back({x,y}); */
/*             } */
/*         } */
/*     } */
/*     return vec_fourPoint; */
/* }; */

void ioPads::init() {
    // make and stylize the TCanvas and pads currently in the list

    const char* t_name = Form("canv_%s",ioUniqueName(101));
    canvas = new TCanvas(t_name, "",canvas_width, canvas_height);
    /* cout << " a0 " << endl; */
    io_fmt(canvas);
    /* cout << " a1 " << endl; */
    canvas->Draw();
    /* cout << " a2 " << endl; */
    canvas->cd();
    /* cout << " a3 " << endl; */


    // add all pads
    add_pad(pad_dimensions);

    canvas->cd();
    /* canvas_pad = new TPad("canvas_pad","",0.,0.,1.,1.); */
    /* io_fmt(canvas_pad); */
    /* canvas_pad->Draw(); */
};


void ioPads::add_pad(pair<ioPadDim,ioPadDim>& coord){
    canvas->cd();

    if (pads.size()==0) {
        canvas_pad = new TPad("canvas_pad","",0.,0.,1.,1.);
        io_fmt(canvas_pad);
        canvas_pad->Draw();
        canvas->cd();
    }

    const ioPadDim x { coord.second };
    const ioPadDim y { coord.first  };
    int i{0};
    while (gDirectory->FindObjectAny(Form("loc_pad_%i",i))) { ++i; }
    TPad* p = new TPad(Form("loc_pad_%i",i),"",x.low,y.low,x.up,y.up);

    // set the boundaries left(l), right(r), top(t), bottom(b)
    p->SetLeftMargin(x.low_margin());
    p->SetRightMargin(x.up_margin());
    p->SetBottomMargin(y.low_margin());
    p->SetTopMargin(y.up_margin());

    io_fmt(p);
    p->Draw();
    pads.push_back(p);
};

void ioPads::add_pad(vector<pair<ioPadDim,ioPadDim>> input) {
    for (auto& inp : input) add_pad(inp);
};

// Implementation of ioIntList
string ioIntList::make(const char* in_file, bool print) {
    ostringstream msg;
    ifstream file;
    file.open(in_file);
    if (!file.is_open()) {
        msg << "Could not open int list \"" << in_file << "\". No entries entered." << endl;
        cout << msg.str() << endl;
        return msg.str();
    }
    string line;
    while (getline(file,line)) {
        line.append(" ");
        stringstream words(line);
        TString word;
        int n_words {0};
        while (words >> word) {
            if (word.BeginsWith("//") || word.BeginsWith("#")) break;
            // skip lines that start with a non-number (assume it is a header)
            if (n_words==0) {
                if (!word.IsAlnum()) break;
            }
            ++n_words;
            list.push_back(word.Atoi());
        }
    }
    sort(list.begin(),list.end());
    file.close();
    msg << " Successfully read in integer list from \"" << in_file << "\"";
    if (print) msg << ". Values: ";
    msg << endl;
    if (print) {
        for (auto& i : list) msg << "  " << i << endl;
    }
    cout << msg.str();
    return msg.str();
};

ioIntList::ioIntList(const char* in_file, ofstream& log, bool print) {
    log << make(in_file, print);
};
ioIntList::ioIntList(const char* in_file, bool print) { make(in_file, print); };

bool ioIntList::operator()(int val) {
    return std::binary_search(list.begin(), list.end(), val);
};
int ioIntList::operator[](int val) {
    return (int)(std::lower_bound(list.begin(), list.end(), val) - list.begin());
};
bool ioIntList::has(int val) { return this->operator()(val); };
bool ioIntList::has_not(int val) { return !(this->operator()(val)); };


bool ioRunListId::has_run(int id) {
    return static_cast<bool>(map_id.count(id));
};

int ioRunListId::size() { return map_id.size(); };

double ioRunListId::operator()(int run_id) {
    set_id(run_id);
    return id;
};

double ioRunListId::set_id(int run_id) {
    if (has_run(run_id)) {
        id = static_cast<double>(map_id[run_id]);
    } else {
        id = -1.;
    }
    return id;
};

// Implementation of ioIntMap
ioIntMap::ioIntMap(const char* file,
        int index_column, 
        int data_column,
        bool echo_print,
        vector<int> skip_vals
        ) {
    ioIntMap_constructor(file, index_column, data_column, 
            echo_print, skip_vals);
};
ioIntMap::ioIntMap(const char* file,
        int index_column, 
        int data_column,
        bool echo_print,
        ofstream& log,
        vector<int> skip_vals
        ) {
    log << ioIntMap_constructor(file, index_column, data_column, 
            echo_print, skip_vals);
};
string ioIntMap::ioIntMap_constructor (
        const char* in_file,
        int index_column, 
        int data_column,
        bool echo_print,
        vector<int> skip_vals
        ) {
    set<int> skip_val_set;
    for (auto& v : skip_vals) skip_val_set.insert(v);
    /* string ioIntMap::make(const char* in_file, bool print) { */
    ostringstream msg;
    ifstream file;
    file.open(in_file);
    if (!file.is_open()) {
        msg << "Could not open int map file \"" << in_file 
            << "\". No entries entered." << endl;
        cout << msg.str() << endl;
        return msg.str();
    } else {
        if (echo_print) msg << " Reading file " << in_file << " for map of col("
            << index_column << ") to col(" << data_column << ")" << endl;
    }
    /* int n_req { index_column > data_column ? index_column : data_column }; */
    string line;
    while (getline(file,line)) {
        line.append(" ");
        stringstream words(line);
        TString word;
        bool has_index { false };
        bool has_data  { false };
        int  val_index{-1};
        int  val_data{-1};
        int i {0};
        bool comment_flag {false};
        int n_words{0};
        while (words >> word) {
            /* cout << " Line: " << line << "  ->(word): " << word << "  IsAlnum: " << word.IsAlnum()<< endl; */
            /* cout << " word: " << word << endl; */
            if (word.BeginsWith("//") || word.BeginsWith("#")) {
                comment_flag = true;
                break;
            }
            // ignore columns that start with a non-number word
            if (n_words == 0 && !word.IsDigit()) { comment_flag = true; break; }
            ++n_words;
            if (i == index_column) {
                has_index = true;
                val_index = word.Atoi();
                if (skip_val_set.count(val_index)) continue;
            } 
            if ( i == data_column) {
                has_data = true;
                val_data = word.Atoi();
            }
            if (has_index && has_data) break;
            ++i;
        }
        if (comment_flag) continue;
        if (!has_index || !has_data) {
            /* ostringstream loc_msg; */
            msg << "fatal error in reading ioIntMap line from file: "
                << "  -> " << in_file << endl;
            if (!has_index) 
                msg << " Couldn't read index column("<<index_column 
                    <<")" <<endl;
            if (!has_data) 
                msg << " Couldn't read data column("<<data_column 
                    <<")" <<endl;
            msg << "  from line:" << endl
                << "  -> " << line << endl;
            msg << "  Therefore skipping entries on this line." << endl;
        } else {
            if (echo_print) msg << "  " << val_index 
                << " -> " << val_data << endl;
            data_map[val_index] = val_data;
        }
    }
    file.close();
    msg << " Done reading columns from file \"" << in_file <<"\"" << endl;
    cout << msg.str();
    return msg.str();
};

bool ioIntMap::has(int key) {
    return (bool) data_map.count(key);
};

int& ioIntMap::operator[](int key) { return data_map[key]; };

vector<int> ioIntMap::keys() {
    vector<int> vec;
    for (auto m : data_map) vec.push_back(m.first);
    sort(vec.begin(), vec.end());
    return vec;
};
int ioIntMap::size() { return data_map.size(); };


// ioHgStats
ioHgStats::ioHgStats(vector<double> ax_vals, vector<double> _vals, vector<double> _errs, bool cut_zeros) {
    if (ax_vals.size() != _vals.size() || _vals.size() != _errs.size())
        throw std::runtime_error(
                "ioHgStats(vec, vec, vec) required vectors of same length");
    if (ax_vals.size() < 2) 
        throw std::runtime_error(
                "ioHgStats(vec, vec, vec) required vectors size > 1");
    vals = _vals;
    errs = _errs;
    nbins = ax_vals.size();

    double *p = new double[nbins+1];
    p[0] = ax_vals[0]-0.5*(ax_vals[1]-ax_vals[0]);
    for (int i{1}; i<nbins; ++i) p[i] = 0.5*(ax_vals[i-1]+ax_vals[i]);
    p[nbins] = ax_vals[nbins-1]+0.5*(ax_vals[nbins-1]-ax_vals[nbins-2]);
    axis = new TAxis(nbins, p);

    if (cut_zeros) {
        for (auto i{0}; i<nbins;++i) {
            weight.push_back( vals[i] == 0 ? 0. : 1.);
        }
    } else {
        for (auto i{0}; i<nbins; ++i) {
            weight.push_back(1.);
        }
    }
    calc_stats();
};
ioHgStats::ioHgStats(TH1D* hg, bool cut_zeros) {
    vals = io_vecBinContent(hg);
    errs = io_vecBinError  (hg);
    axis=hg->GetXaxis();
    nbins = axis->GetNbins();
    if (cut_zeros) {
        for (auto i{0}; i<nbins;++i) {
            weight.push_back( vals[i] == 0 ? 0. : 1.);
        }
    } else {
        for (auto i{0}; i<nbins; ++i) {
            weight.push_back(1.);
        }
    }
    calc_stats();
};
ioHgStats::ioHgStats(TProfile* hg, bool cut_zeros, bool weight_by_entries) {
    vals = io_vecBinContent(hg);
    errs = io_vecBinError  (hg);
    axis=hg->GetXaxis();
    nbins = axis->GetNbins();
    if (weight_by_entries) {
        weight = io_vecBinEntries(hg);
    } else {
        if (cut_zeros) {
            for (auto i{0}; i<nbins;++i) {
                weight.push_back( vals[i] == 0 ? 0. : 1.);
            }
        } else {
            for (auto i{0}; i<nbins; ++i) {
                weight.push_back(1.);
            }
        }
    }
    calc_stats();
};
void ioHgStats::calc_stats() {
    mean   = TMath::Mean(vals.begin(), vals.end(), weight.begin());
    stddev = TMath::StdDev(vals.begin(), vals.end(), weight.begin());
};
double ioHgStats::mean_Xsigma(double X) {
    /* double sumW{0}; */
    /* double sumV{0}; */
    /* double sumE{0}; */
    /* for (auto w : weight) sumW += w; */
    /* for (auto w : vals) sumV += w; */
    /* for (auto w : errs) sumE += w; */
    /* cout << " mean: " << mean << "  stddev " << stddev << "  X " << X << "  val: " << */
    /* mean+X*stddev << " >> npts" << nbins << Form("stats:%f,%f,%f",sumV,sumE,sumW) <<endl; */
    return mean+ X * stddev;
};
TLine* ioHgStats::get_horizontal_TLine(double cut, bool cut_times_sigma) {
    double x0 = axis->GetBinLowEdge(1);
    double x1 = axis->GetBinUpEdge(axis->GetNbins());
    double y = cut_times_sigma ? mean_Xsigma(cut) : cut;
    return new TLine(x0,y,x1,y);
};
TGraph* ioHgStats::points_above(double cut, bool cut_times_sigma) {
    // count how many points are above
    double y = cut_times_sigma ? mean_Xsigma(cut) : cut;
    int n_above {0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && vals[i]>y) ++n_above;
    }
    double* xpts = new double[n_above];
    double* ypts = new double[n_above];
    int n{0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && vals[i]>y) {
            xpts[n] = axis->GetBinCenter(i+1);
            ypts[n] = vals[i];
            ++n;
        }
    }
    return new TGraph(n_above,xpts,ypts);
};
TGraph* ioHgStats::points_below(double cut, bool cut_times_sigma) {
    // count how many points are above
    double y = cut_times_sigma ? mean_Xsigma(cut) : cut;
    int n_pts {0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && vals[i]<y) ++n_pts;
    }
    double* xpts = new double[n_pts];
    double* ypts = new double[n_pts];
    int n{0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && vals[i]<y) {
            xpts[n] = axis->GetBinCenter(i+1);
            ypts[n] = vals[i];
            ++n;
        }
    }
    return new TGraph(n_pts,xpts,ypts);
};
TGraph* ioHgStats::points_between(double locut, double hicut, bool cut_times_sigma) {
    // count how many points are above
    double y_lo = cut_times_sigma ? mean_Xsigma(locut) : locut;
    double y_hi = cut_times_sigma ? mean_Xsigma(hicut) : hicut;
    int n_pts {0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && vals[i]>y_lo && vals[i]<y_hi) ++n_pts;
    }
    double* xpts = new double[n_pts];
    double* ypts = new double[n_pts];
    int n{0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && vals[i]<y_hi && vals[i]>y_lo) {
            xpts[n] = axis->GetBinCenter(i+1);
            ypts[n] = vals[i];
            ++n;
        }
    }
    return new TGraph(n_pts,xpts,ypts);
};
// cuts
ioHgStats& ioHgStats::cut_above(double cut, bool cut_times_sigma) {
    // count how many points are above
    double y = cut_times_sigma ? mean_Xsigma(cut) : cut;
    int n_above {0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && vals[i]>y) weight[i] = 0.;
    }
    calc_stats();
    return *this;
};
ioHgStats& ioHgStats::cut_below(double cut, bool cut_times_sigma) {
    // count how many points are above
    double y = cut_times_sigma ? mean_Xsigma(cut) : cut;
    int n_pts {0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && vals[i]<y) weight[i] = 0.;
    }
    calc_stats();
    return *this;
};
vector<double> ioHgStats::unmasked_vals() {
    vector<double> r_vec;
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0.) r_vec.push_back(vals[i]);
    }
    return r_vec;
};
ioHgStats& ioHgStats::mask(const vector<bool>  mask) {
    if ((int)mask.size() != nbins) {
        cout << " error in ioHgStats::mask " << endl
            << " There are " << nbins << " bins, but only " << mask.size()
            << " points in in put mask " << endl << endl
            << " Therefore mask is not applied." << endl;
        return *this;
    }
    for (int i{0}; i<nbins; ++i) {
        if (!mask[i]) weight[i] = 0.;
        calc_stats();
    }
    return *this;
};
/* ioHgStats& ioHgStats::mask(vector<int>& mask, bool mask_keep_true) { */
/*     if ((int)mask.size() != nbins) { */
/*         cout << " error in ioHgStats::mask " << endl */
/*              << " There are " << nbins << " bins, but only " << mask.size() */
/*              << " points in in put mask " << endl << endl */
/*              << " Therefore mask is not applied." << endl; */
/*         return *this; */
/*     } */
/*     for (int i{0}; i<nbins; ++i) { */
/*         if ((mask[i] == 0)) weight[i] = 0.; */
/*         calc_stats(); */
/*     } */
/*     return *this; */
/* }; */
ioHgStats& ioHgStats::cut_to_range(double locut, double hicut, bool cut_times_sigma) {
    // count how many points are above
    double y_lo = cut_times_sigma ? mean_Xsigma(locut) : locut;
    double y_hi = cut_times_sigma ? mean_Xsigma(hicut) : hicut;
    int n_pts {0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0. && (vals[i]<y_lo || vals[i]>y_hi)) weight[i] = 0.;
    }
    calc_stats();
    return *this;
};
TGraph* ioHgStats::points() {
    int n_pts {0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0.) ++n_pts;
    }
    double* xpts = new double[n_pts];
    double* ypts = new double[n_pts];
    int n{0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0.) {
            xpts[n] = axis->GetBinCenter(i+1);
            ypts[n] = vals[i];
            ++n;
        }
    }
    return new TGraph(n_pts,xpts,ypts);
};
int ioHgStats::count_points() {
    int n_pts {0};
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0.) ++n_pts;
    }
    return n_pts;
};
vector<int> ioHgStats::bin_indices() {
    vector<int> vec;
    for (auto i{0}; i<nbins; ++i) {
        if (weight[i]!=0.) vec.push_back(i+1);
    }
    return vec;
};
ioHgStats& ioHgStats::restore_points() {
    for (auto i{0}; i<nbins; ++i) {
        weight[i] = 1.;
    }
    calc_stats();
    return *this;
};
ioHgStats& ioHgStats::cut_zeros() {
    for (auto i{0}; i<nbins; ++i) {
        if (vals[i]==0) weight[i] =0.;
    }
    calc_stats();
    return *this;
};

// ------------------------------------------------------
// | Implementation of ioMsgTree                        |
// ------------------------------------------------------
ioMsgTree::ioMsgTree(bool set_echo) : 
    b_msg{""}, 
    tree{"Messages", "Tree of messages"}
{
    tree.Branch("messages", &b_msg);
    dash();
    msg(" Start of msg tree ");
    dash();
    echo_to_cout = set_echo;
};
void ioMsgTree::msg(string msg) {
    b_msg = msg;
    if (echo_to_cout) cout << b_msg << endl;
    tree.Fill();
};
void ioMsgTree::msg(vector<string> messages) {
    for (auto& msg : messages) {
        b_msg = msg;
        if (echo_to_cout) cout << b_msg << endl;
        tree.Fill();
    }
};
void ioMsgTree::dash() {
    b_msg = "---------------";
    if (echo_to_cout) cout << b_msg << endl;
    tree.Fill();
};
void ioMsgTree::write(){
    tree.Write();
};
void ioMsgTree::read_messages(const char* f_name){
    cout << " Reading file: " << f_name << endl;
    /* TTree *tree; */
    TFile* fin  = new TFile(f_name, "read");
    if (!fin->IsOpen()) {
        cout << " Input file: " << f_name << " is not open." << endl;
        delete fin;
        return;
    }
    TTreeReader myReader("Messages",fin);
    TTreeReaderValue<string> msg(myReader, "messages");
    cout << "  Contents of TFile(\""<<f_name<<"\") TTree(\"Messages\"):" << endl;
    while (myReader.Next()) cout << " " << *msg << endl;
    fin->Close();
    delete fin;
    /* } */
    };
void ioMsgTree::slurp_file(const char* which_file) {
    // try and read all lines of which_file into the tree
    msg(Form("--Begin contents of file \"%s\"",which_file));
    ifstream f_in {which_file};
    if (!f_in.is_open()) {
        msg(Form("  error: couldn't open file \"%s\"",which_file));
    } else {
        string line;
        while (getline(f_in,line)) msg(line);
    }
    msg(Form("--End contents of file \"%s\"",which_file));
    return;
};

vector<int> ioIntVec::vals(string tag) {
    assert_tag(tag,"vals");
    int col { i_tag(tag) };
    vector<int> r_vec {};
    for (auto i{0}; i<size(); ++i) r_vec.push_back(data[i][col]);
    return r_vec;
};

vector<int> ioIntVec::vals(string tag, vector<bool> mask) {
    assert_tag(tag,"vals");
    int col { i_tag(tag) };

    if (mask.size() != data.size())
        throw std::runtime_error(
                Form("fatal error in ioIntVec::vals: size of mask and data do not match"));

    vector<int> r_vec {};
    for (auto i{0}; i<(int)mask.size(); ++i) if (mask[i]) r_vec.push_back(data[i][col]);
    return r_vec;
};

bool ioIntVec::flip_tag(string& tag) {
    int t_size = tag.size();
    if (t_size==0) return true;
    if (tag.substr(0,1)=="!") {
        tag = tag.substr(1,t_size-1);
        return false;
    }
    return true;
};

vector<bool> ioIntVec::mask(string tag, const char* module) {
    bool on_true = flip_tag(tag);
    assert_tag(tag,module);
    int col { i_tag(tag) };
    vector<bool> r_vec;
    for (auto& vec : data) r_vec.push_back((vec[col]!=0) == on_true);
    return r_vec;
};

vector<bool> ioIntVec::is_any(vector<string> _tags){
    int s = _tags.size();
    if (!s) return {};
    auto r_vec = mask(_tags[0]);
    for (int i{1}; i<s; ++i) r_vec = r_vec || mask(_tags[i],"is_any");
    return r_vec;
};
vector<bool> ioIntVec::is_all(vector<string> _tags){
    int s = _tags.size();
    if (!s) return {};
    auto r_vec = mask(_tags[0]);
    for (int i{1}; i<s; ++i) r_vec = r_vec && mask(_tags[i],"is_all");
    return r_vec;
};


ioIntVec::ioIntVec( const char* file_name, bool echo_print, vector<string>_tags) 
{ ioIntVec_constructor(file_name, echo_print, _tags); };

ioIntVec::ioIntVec( const char* file_name, ofstream& log, 
        bool echo_print, vector<string>_tags) 
{ 
    log << ioIntVec_constructor(file_name, echo_print, _tags); 
    log << *this;
};

string ioIntVec::ioIntVec_constructor(const char* in_file, 
        bool echo_print, vector<string>tags_requested) 
{
    // read until the tag line is read
    // if given input tags, check they are all present, and use them in that order
    // 

    ostringstream msg;
    ifstream file;
    file.open(in_file);
    if (!file.is_open()) {
        msg << "Could not open int to vector<int> map file \"" << in_file 
            << "\"." << endl << "  -> No entries entered." << endl;
        cout << msg.str() << endl;
        return msg.str();
    } else {
        if (echo_print) msg << " Reading file " << in_file << endl;
    }

    // Read until get the column tags (first uncommented line)

    string line;
    /* vector<string> tags_read; */
    while (getline(file,line)) {
        /* cout << " reading line: " << line << endl; */
        line.append(" ");
        stringstream words(line);
        TString word;
        bool is_key_tag {true};
        while (words >> word) {
            if (word.BeginsWith("//") || word.BeginsWith("#")){
                break;
            }
            if (is_key_tag) {
                is_key_tag = false;
            } else {
                tags.push_back(word.Data());
            }
        }
        if (tags.size() != 0) break;
    }
    if (tags_requested.size() == 0) {
        tags_requested = tags;
    } 
    auto read_cols = tag_cols(tags_requested);
    int max_col = 0;
    for (auto val : read_cols) if (val > max_col) max_col = val;
    tags = tags_requested;

    // read all remaining lines into data
    vector<pair<int,vector<int>>> data_in;
    while(getline(file,line)) {
        vector<int> c_vec;
        int key;
        bool has_key {false};

        line.append(" ");
        stringstream words(line);
        TString word;
        while (words >> word) {
            if (word.BeginsWith("//") || word.BeginsWith("#")){
                break;
            }
            if (!has_key) {
                has_key = true;
                key = word.Atoi();
            } else {
                c_vec.push_back(word.Atoi());
            }
        }
        if (has_key) {
            data_in.push_back({key,c_vec});
            if ((int)c_vec.size() < (max_col+1))
                throw std::runtime_error(
                        Form("In ioIntVec needs at least %i entries (+id) in line \"%s\"",
                            max_col+1, line.c_str())
                        );
        }
    }
    std::sort(data_in.begin(), data_in.end());

    for (auto entry : data_in) {
        keys.push_back(entry.first);
        vector<int> vec;
        for (auto i : read_cols) vec.push_back(entry.second[i]);
        data.push_back(vec);
    };

    file.close();
    msg << " Done reading columns from file \"" << in_file <<"\"" << endl;
    if (echo_print) {
        cout << msg.str();
        cout << *this << endl;
    }
    return msg.str();
};

int ioIntVec::i_tag(string tag) {
    auto it = std::find(tags.begin(), tags.end(), tag);
    if (it == tags.end()) return -1;
    else return (int)(it-tags.begin());
};

void ioIntVec::assert_tag(string tag, const char* module) {
    if (!has_tag(tag)) 
        throw std::runtime_error(Form("fatal in ioIntVec::%s, couldn't find tag \"%s\"",module,tag.c_str()));
};
vector<int> ioIntVec::tag_cols(vector<string> _tags, const char* name) {
    vector<int> cols;
    if (_tags.size()==0) {
        for (int i{0}; i<(int)_tags.size(); ++i) cols.push_back(i);
    } else {
        for (auto& tag : _tags) {
            assert_tag(tag,"tag_cols");
            cols.push_back(i_tag(tag));
        }
    }
    return cols;
};
bool ioIntVec::has_tag(string tag) {
    return ( i_tag(tag) != -1 );
};
ioIntVec& ioIntVec::swap_tags(string tag0, string tag1) {
    // find the two tags among tags
    int i0 = i_tag(tag0);
    if (i0 == -1) { 
        cout << "fatal error in ioIntVec::swap_tags "
            << "could not find tag \"" << tag0 << "\"" << endl;
        return *this;
    }

    int i1 = i_tag(tag1);
    if (i1 == -1) {
        cout << "fatal error in ioIntVec::swap_tags "
            << "could not find tag \"" << tag1 << "\"" << endl;
        return *this;
    }

    for (auto& vec : data) {
        auto temp = vec[i0];
        vec[i0] = vec[i1];
        vec[i1] = temp;
    }
    tags[i0] = tag1;
    tags[i1] = tag0;
    return *this;
};
ioIntVec& ioIntVec::subset(vector<string> sub_tags) {
    // cut down the table to only the subset of tags
    auto keep_cols = tag_cols(sub_tags,"subset");
    for (auto& vec : data) {
        vector<int> copy;
        for (auto& i : keep_cols) {
            copy.push_back(vec[i]);
        }
        vec = copy;
    }
    tags = sub_tags;
    return *this;
};
ioIntVec& ioIntVec::rm_tags(vector<string> tags) {
    auto rm_cols = tag_cols(tags,"rm_tags");
    sort(rm_cols.begin(), rm_cols.end());
    vector<string> keep_tags {};
    for (int i{0}; i<(int)tags.size(); ++i) {
        if (find(rm_cols.begin(), rm_cols.end(), i) == rm_cols.end())
            keep_tags.push_back(tags[i]);
    }
    return subset(tags);
};
ioIntVec& ioIntVec::add_tag(string tag, int def_val) {
    int i_col = i_tag(tag);
    if (i_col != -1) {
        cout << " error in ioIntVec::add_tag; added tag " << tag << " already exists! " << endl;
        return *this;
    }
    tags.push_back(tag);
    cout << " Adding tag: " << tag << endl;
    for (auto& vec : data) vec.push_back(def_val);
    /* for (auto key : keys()) data_map[key].push_back(def_val); */
    return *this;
};
void ioIntVec::rename_tag(string tag0, string tag1) {
    assert_tag(tag0,"rename_tag");
    tags[i_tag(tag0)] = tag1;
};
bool ioIntVec::has_key(int key) { 
    return binary_search(keys.begin(), keys.end(), key);
};
vector<int>& ioIntVec::operator[](int key) { 
    if (!has_key(key)) {
        throw std::runtime_error(Form("fatal: no key \"%i\" in ioIntVec",key));
    }
    int i = (int)(std::lower_bound(keys.begin(), keys.end(), key) - keys.begin());
    return data[i];
};
/* vector<int>    ioIntVec::keys() const { */
/*     vector<int> vec; */
/*     for (auto m : data_map) vec.push_back(m.first); */
/*     sort(vec.begin(), vec.end()); */
/*     return vec; */
/* }; */
vector<int> ioIntVec::get_keys(vector<bool> mask, bool keep_on_true) {
    if (mask.size() != data.size())
        throw std::runtime_error(
                Form("fatal error in ioIntVec::keys: size of mask and data do not match"));

    vector<int> r_vec {};
    for (auto i{0}; i<(int)mask.size(); ++i) {
        if (mask[i] == keep_on_true) r_vec.push_back(keys[i]);
    }
    return r_vec;
};
/* vector<bool> ioIntVec::is_any(vector<pair<string,bool>> mask_keep_true) { */
/*     vector<bool>   v_keep_true; */
/*     vector<string> v_tags; */
/*     for (auto mpair : mask_keep_true) { */
/*         v_keep_true.push_back(mpair.second); */
/*         v_tags.push_back(mpair.first); */
/*     } */
/*     auto cols = tag_cols(v_tags, "is_any"); */
/*     int  n {(int) cols.size()}; */
/*     vector<bool> r_vec {}; */
/*     for (auto& vec : data) { */
/*         bool keep {false}; */
/*         for (int i{0}; i<n; ++i) { */
/*             if ( (vec[cols[i]] != 0) == v_keep_true[i]) { */
/*                 keep = true; */
/*                 break; */
/*             } */
/*         } */
/*         r_vec.push_back(keep); */
/*     } */
/*     return r_vec; */
/* }; */
int ioIntVec::size(){ return (int)keys.size(); };
/* vector<bool> ioIntVec::is_all(vector<pair<string,bool>> mask_keep_true) { */
/*     vector<bool>   v_keep_true; */
/*     vector<string> v_tags; */
/*     for (auto mpair : mask_keep_true) { */
/*         v_keep_true.push_back(mpair.second); */
/*         v_tags.push_back(mpair.first); */
/*     } */
/*     auto cols = tag_cols(v_tags, "is_all"); */
/*     int  n {(int) cols.size()}; */
/*     vector<bool> r_vec {}; */
/*     for (auto& vec : data) { */
/*         bool keep {true}; */
/*         for (int i{0}; i<n; ++i) { */
/*             if ( (vec[cols[i]] != 0) != v_keep_true[i]) { */
/*                 keep = false; */
/*                 break; */
/*             } */
/*         } */
/*         r_vec.push_back(keep); */
/*     } */
/*     return r_vec; */
/* }; */
ostream& operator<<(ostream& os, ioIntVec& io) {
    // generate the max charactures needed in each column

    // max length id line
    int max_id = 2;
    for (auto key : io.keys) {
        int n = io_count_digits(key);
        if (n > max_id) max_id = n;
    }
    os << Form(Form(" %%%is",max_id),"id");

    vector<int> max_chars;
    for (auto tag : io.tags) {
        max_chars.push_back(tag.size());
    }
    for (auto& vec : io.data) {
        int i{0};
        for (auto val : vec) {
            int n = io_count_digits(val);
            if (n > max_chars[i]) max_chars[i] = n;
            ++i;
        }
    }

    int i{0};
    for (auto tag : io.tags) {
        os << Form(Form(" %%%is",max_chars[i]),tag.c_str());
        ++i;
    }
    os << endl;

    const char* fmt_key = Form(" %%%ii",max_id);

    vector<const char*> this_fmt;
    for (auto n : max_chars) {
        this_fmt.push_back(Form(" %%%ii",n));
    }

    int z { 0};
    const int n_tags = io.tags.size();
    for (int i{0}; i<(int)io.keys.size(); ++i) {
        os << " " << std::right << std::setw(max_id) << io.keys[i];
        for (int k{0}; k<n_tags; ++k) {
            os << " " << std::setw(max_chars[k]) << io.data[i][k];
        }
        os << endl;
    }
    return os;
};

void ioIntVec::write_to_file(const char* which_file, vector<string>comments) {
    ofstream f_out { which_file };
    if (!f_out.is_open()) {
        cout <<  "  fatal error in ioIntVec::write_to_file : " << endl
            <<  "    Couldn't open file \""<<which_file<<"\""<< endl
            <<  "    -> not file being written to " << endl;
    }
    for (auto comment : comments) f_out << "// " << comment << endl;
    f_out << *this;
    cout << " Write table to file: " << which_file << endl;
    f_out.close();
};

void ioMinMaxPtr::operator()(double _, void* p){
    if (!has_data) {
        has_data = true;
        min_ptr = p;
        max_ptr = p;
        min_val = _;
        max_val = _;
    } else {
        if (_ > max_val) {
            max_val = _;
            max_ptr = p;
        }
        if (_ < min_val) {
            min_val = _;
            min_ptr = p;
        }
    }
};

ioMinMax::ioMinMax(string _name) : name{_name} {};
long long int ioMinMax::fill(double val) {
    if (n_entries == 0) {
        min = val;
        max = val;
    } else {
        if (val < min) min = val;
        if (val > max) max = val;
    }
    ++n_entries;
    return n_entries;
};
ioMinMax::ioMinMax(double& ptr, string _name) : 
    name{_name}, fill_option{5}, ptr_double{&ptr} {};
ioMinMax::ioMinMax(int& ptr, string _name) : 
    name{_name}, fill_option{6}, ptr_int{&ptr} {};
ioMinMax::ioMinMax(unsigned int& ptr, string _name) : 
    name{_name}, fill_option{7}, ptr_uint{&ptr} {};
ioMinMax::ioMinMax(short& ptr, string _name) : 
    name{_name}, fill_option{8}, ptr_short{&ptr} {};
ioMinMax::ioMinMax(char& ptr, string _name) : 
    name{_name}, fill_option{9}, ptr_char{&ptr} {};

ioMinMax::ioMinMax(double* ptr, string _name) : 
    name{_name}, fill_option{10}, ptr_double{ptr} {};
ioMinMax::ioMinMax(int* ptr, string _name) : 
    name{_name}, fill_option{11}, ptr_int{ptr} {};
ioMinMax::ioMinMax(unsigned int* ptr, string _name) : 
    name{_name}, fill_option{12}, ptr_uint{ptr} {};
ioMinMax::ioMinMax(short* ptr, string _name) : 
    name{_name}, fill_option{13}, ptr_short{ptr} {};
ioMinMax::ioMinMax(char* ptr, string _name) : 
    name{_name}, fill_option{14}, ptr_char{ptr} {};

ioMinMax::ioMinMax(double* ptr, int& _index, string _name) : 
    name{_name}, fill_option{0}, ptr_double{ptr}, index{&_index} {};
ioMinMax::ioMinMax(int* ptr, int& _index, string _name) : 
    name{_name}, fill_option{1}, ptr_int{ptr}, index{&_index} {};
ioMinMax::ioMinMax(unsigned int* ptr, int& _index, string _name) : 
    name{_name}, fill_option{2}, ptr_uint{ptr}, index{&_index} {};
ioMinMax::ioMinMax(short* ptr, int& _index, string _name) : 
    name{_name}, fill_option{3}, ptr_short{ptr}, index{&_index} {};
ioMinMax::ioMinMax(char* ptr, int& _index, string _name) : 
    name{_name}, fill_option{4}, ptr_char{ptr}, index{&_index} {};

long long int ioMinMax::operator()() {
    switch (fill_option) {
        case 0:
            return fill(ptr_double[*index]);
        case 1:
            return fill(ptr_int[*index]);
        case 2:
            return fill(ptr_uint[*index]);
        case 3:
            return fill(ptr_short[*index]);
        case 4:
            return fill(ptr_char[*index]);

        case 5:
            return fill(*ptr_double);
        case 6:
            return fill(*ptr_int);
        case 7:
            return fill(*ptr_uint);
        case 8:
            return fill(*ptr_short);
        case 9:
            return fill(*ptr_char);
        default:
            throw std::runtime_error("ioMinMax::operator()() called with no pointer set");
    }
};
long long int ioMinMax::operator()(int index) {
    switch (fill_option) {
        case 10:
            return fill(ptr_double[index]);
        case 11:
            return fill(ptr_int[index]);
        case 12:
            return fill(ptr_uint[index]);
        case 13:
            return fill(ptr_short[index]);
        case 14:
            return fill(ptr_char[index]);
        default:
            throw std::runtime_error("ioMinMax::operator()(int) called with no pointer set");
    }
};
int ioMinMax::nbins() { return (int)(max-min)+1; };
ostream& operator<<(ostream& os, ioMinMax& self) {
    if (self.name != "") cout << self.name << ": ";
    os << self.min << " " << self.max;
    return os;
};

ioIntSet& ioIntSet::operator+=(const ioIntSet& rhs) {
    vector<int> new_vals;
    for (auto val : rhs.list) {
        if (!binary_search(list.begin(),list.end(),val)) {
            new_vals.push_back(val);
        }
    }
    for (auto v : new_vals) list.push_back(v);
    sort(list.begin(),list.end());
    return *this;
};
ioIntSet& ioIntSet::operator*=(const ioIntSet& sec) {
    vector<int> new_list;
    for (auto val : sec.list)
        if (binary_search(list.begin(),list.end(),val)) 
            new_list.push_back(val);
    list = new_list;
    return *this;
};


int ioIntSet::size() { return list.size(); };
void ioIntSet::clear() { list.clear(); };
ostringstream ioIntSet::read_file(const char* in_file, int col, bool print, bool strip_commas) {
    ostringstream msg;
    if (!strcmp(in_file,"")) return msg;
    if (print) {
        msg << " Reading following values from col " << col << " from " << in_file << endl;
    }
    try {
        auto new_data = ioReadIntVec(in_file, col, true, strip_commas);
        if (list.size()>0) {
            vector<int> add_vals{};
            for (int i{0}; i<(int)new_data.size(); ++i) {
                if (i>0 && new_data[i] == new_data[i-1]) continue;
                if (has(new_data[i])) continue;
                if (print) msg << " set add " << new_data[i] << endl;
                add_vals.push_back(new_data[i]);
            }
            for (auto v : add_vals) list.push_back(v);
        } else {
            for (int i{0}; i<(int)new_data.size(); ++i) {
                if (i>0 && new_data[i] == new_data[i-1]) continue;
                list.push_back(new_data[i]);
            }
        }
        sort(list.begin(), list.end());
        if (print) cout << " Done reading col " << col << " from " << in_file << endl;
    }
    catch (std::runtime_error err) {
        cerr << " fatal error in ioIntSet::read_file " << endl;
        cerr << err.what() << endl;
        /* cout << err << endl; */
    }
    return msg;
};
ostream& operator<<(ostream& os, ioIntSet& dt) { 
    for (auto& v : dt.list) cout << v << endl;
    return os;
};
ioIntSet::ioIntSet(const char* in_file, ofstream& log, int col, bool print, bool strip_commas) {
    log << read_file(in_file, col, print, strip_commas).str() << endl;
};
ioIntSet::ioIntSet(const char* in_file, int col, bool print, bool strip_commas) {
    read_file(in_file, col, print, strip_commas);
};
ioIntSet::ioIntSet(const char* file, const char* tag) {
    for (auto val : ioReadValVec(file,tag,{{"sort",true}})) {
        list.push_back((int)val);
    }
};
bool ioIntSet::operator()(int val) { return std::binary_search(list.begin(), list.end(), val); };
bool ioIntSet::has(int i) { return binary_search(list.begin(),list.end(),i); };
int ioIntSet::operator[](int val) {
    return (int)(std::lower_bound(list.begin(), list.end(), val) - list.begin());
};

ioIntBinCnt::ioIntBinCnt(const char* name, vector<int> x_dim, const char* title) :
    ioIntBinCnt{name, x_dim, {}, title} {};
ioIntBinCnt::ioIntBinCnt(const char* name, vector<int> x_dim, 
        vector<int> y_dim, const char* title) {
    if (y_dim.size()>0) {
        is2D = true;
        const char* use_title = (strcmp(title,"")) ? name : title;
        hg2 = new TH2D(name, use_title, x_dim.size(), ax_doubleptr(x_dim),
                y_dim.size(), ax_doubleptr(y_dim));
    } else {
        const char* use_title = (strcmp(title,"")) ? name : title;
        hg1 = new TH1D(name,use_title, x_dim.size(), ax_doubleptr(x_dim));
    }
};
double ioIntBinCnt::getcnt(TH1D* hg, int i) {
    int i_bin = hg->FindBin((double)i);
    return hg->GetBinContent(i_bin);
};
double ioIntBinCnt::getcnt(TH2D* hg, int i, int j) {
    int i_bin = hg->FindBin((double)i,(double)j);
    return hg->GetBinContent(i_bin);
};
void ioIntBinCnt::fill(int i, double weight){
    hg1->Fill( (double)i,  weight);
};
void ioIntBinCnt::fill(int i, int j, double weight){
    if (!is2D) {
        cout << " Error: only 1D counter generated " << endl;
        return;
    }
    hg2->Fill( (double)i, (double)j, weight);
};
void ioIntBinCnt::write() {
    if (is2D) hg2->Write();
    else hg1->Write();
};

ioFnCaller::ioFnCaller(const char* file_data, double(&_fn)(double*,double*)) :
    fn{_fn}, x{new double[2]}
{
    auto p_vals = ioReadValVec(file_data);
    p = new double[p_vals.size()];
    int i{0};
    for (auto val : p_vals) p[i++] = val;
};
double ioFnCaller::operator()(double x0, double x1){
    x[0] = x0;
    x[1] = x1;
    return fn(x,p);
};

//-------------------------------------------------
ioXsec::ioXsec(
     const char* tag_file,
     const char* Xsection_tag,
     const char* pthatbin_tag,
     const char* nEvents_tag
 ) :
    Xsection  { ioReadValVec(tag_file, Xsection_tag) },
    pthatbins { ioReadValVec(tag_file, pthatbin_tag) },
    nbins_pthat { static_cast<int>( Xsection.size()) },
    Nevents {},
    Ncollected ( nbins_pthat, 0 )
{
    for (auto v : ioReadValVec(tag_file, nEvents_tag)) {
        Nevents.push_back(static_cast<int>(v));
    }
    ioBinVec pthat_edges { tag_file, {{"tag",pthatbin_tag}} };
    hg_collected = new TH1D("pthat_bins_collected",
        "Number of events collected in each pthat bin;#hat{p}_{T};n-collected",
        pthat_edges, pthat_edges);
};

double ioXsec::pthatbin_center(int bin) {
    return 0.5*(pthatbins[bin]+pthatbins[bin+1]);
};

int ioXsec::pthatbin(pair<double,double> bounds) {
    double first { bounds.first };
    auto iter = std::lower_bound(pthatbins.begin(), pthatbins.end(), bounds.first);
    if (iter==pthatbins.end() || (*(iter)!=bounds.first)) { 
        throw std::runtime_error(" fatal in ioXsec::pthatbin: pthatbin not found " );
        return -1;
    }
    if (*(iter+1) != bounds.second) {
        cout << " warning in ioXsec::pthatbin: pthatbin upperbound sought " << endl
             << " doesnt match pthatbin lower found." << endl;
    }
    return (iter - pthatbins.begin());
};

double ioXsec::Xsec(pair<double,double> bounds, int numEvents) {
    return Xsec(pthatbin(bounds), numEvents);
};
double ioXsec::Xsec(int pthatbin, int numEvents) {
    check_pthatbin(pthatbin);
    if (numEvents != 0) return Xsection[pthatbin] / numEvents;

    if (n_collected_total != 0) {
        if (Ncollected[pthatbin]==0) {
            cout << " fatal error: no events for bin " << pthatbin 
                << " have been collected." << endl
                << " Cannot generate weighted cross section." << endl;
            throw std::runtime_error("Asking for ptbin with no events in ioXsec::Xsec");
        } else {
            return Xsection[pthatbin] / Ncollected[pthatbin] ;
        }
    }
    return Xsection[pthatbin] / Nevents[pthatbin];
};

void ioXsec::collect (int pthatbin) {
    ++n_collected_total;
    ++Ncollected[pthatbin];
};
void ioXsec::collect(pair<double,double>bounds) { collect(pthatbin(bounds)); };

void ioXsec::check_pthatbin(int bin) {
    if (bin >= nbins_pthat) {
        throw std::runtime_error( Form(
        "fatal in ioXsec : asked for pthatbin %i but there are only %i bins",
        bin, nbins_pthat));
    }
};

/* io_pThatOutliers::io_pThatOutliers( */
/*         map<int,double> _pt_limits, */ 
/*         double _pt_fakes, double _pt_misses */
/* ) : */
/*     pt_limits {_pt_limits}, */ 
/*     pt_fakes {_pt_fakes}, */ 
/*     pt_misses {_pt_misses} */ 
/* { */
/*     for (auto p : pt_limits) { */
/*         M_matches[p.first] = {}; */
/*         T_matches[p.first] = {}; */
/*         misses[p.first] = {}; */
/*         fakes[p.first] = {}; */
/*     } */
/* }; */

/* bool io_pThatOutliers::check_if_outlier(int _pthatbin, double pt) { */
/*     if (pt_limits.count(_pthatbin) == 0) { */
/*         throw std::runtime_error( */
/*             Form("fatal error in io_pThatOutliers: " */
/*             " limit for required pthatbin (%i) not set", */
/*             _pthatbin) */
/*         ); */
/*     } */
/*     is_outlier = (pt >= pt_limits[_pthatbin]); */
/*     pthatbin = _pthatbin; */
/*     return is_outlier; */
/* }; */

/* void io_pThatOutliers::match(double M, double T) { */
/*     if (is_outlier) { */
/*         M_matches[pthatbin].push_back(M); */
/*         T_matches[pthatbin].push_back(T); */
/*     } */
/* }; */
/* void io_pThatOutliers::miss(double T) { */
/*     if (is_outlier) misses[pthatbin].push_back(T); */
/* }; */
/* void io_pThatOutliers::fake(double M) { */
/*     if (is_outlier) fakes[pthatbin].push_back(M); */
/* }; */


/* void io_pThatOutliers::write_TGraph(int _pthatbin, */
/*         ioOptMap options, ioOptMap dict) { */
/*     dict += options; */
/*     if (pt_limits.count(_pthatbin)==0) { */
/*         throw std::runtime_error( */
/*         Form( */
/*             "fatal error in io_pThatOutliers: " */
/*             " limit for required pthatbin (%i) not set", */
/*             _pthatbin) */
/*         ); */
/*     } */

/*     // draw the matches TGraph */
/*     if (M_matches[pthatbin].size()>0) { */
/*         ioBinVec x { M_matches[pthatbin], false }; */
/*         ioBinVec y { T_matches[pthatbin], false }; */
/*         TGraph gr { x.size, x, y }; */
/*         io_fmt( &gr, options ); */
/*         gr.Write( Form("%s_matches",dict["prefix"].c_str()) ); */
/*     } */
/*     dict["MarkerStyle"] = dict["MarkerFakeMiss"]; */
/*     if (fakes[pthatbin].size() > 0) { */
/*         ioBinVec x { fakes[pthatbin], false }; */
/*         ioBinVec y { vector<double>(fakes[pthatbin].size(), */ 
/*                 pt_fakes) }; */
/*         TGraph gr { x.size, x, y }; */
/*         io_fmt( &gr, options ); */
/*         gr.Write( Form("%s_fakes",dict["prefix"].c_str()) ); */
/*     } */
/*     if (misses[pthatbin].size() > 0) { */
/*         ioBinVec x { vector<double>(misses[pthatbin].size(), */ 
/*                 pt_misses) }; */
/*         ioBinVec y { misses[pthatbin], false }; */
/*         TGraph gr { x.size, x, y }; */
/*         io_fmt( &gr, options ); */
/*         gr.Write( Form("%s_misses",dict["prefix"].c_str()) ); */
/*     } */
/* }; */
    
ioIntStrFunctor::ioIntStrFunctor ( const char* file, ioOptMap options, ioOptMap dict)
{
    dict += options;
    data = ioReadIntStrMap(file, dict);
};
const char* ioIntStrFunctor::operator()(int index) {
    try {
        return data[index].c_str();
    }
    catch (...) {
        cout << " Key failure in ioIntStrFunctor with index " << index << endl;
        throw;
    };
};


//StrStr
ioStrStrFunctor::ioStrStrFunctor ( const char* file, ioOptMap options, ioOptMap dict)
{
    dict += options;
    data = ioReadMapStrStr(file, dict);
};
const char* ioStrStrFunctor::operator()(const char* key) {
    try {
        return data[key].c_str();
    }
    catch (...) {
        cout << " Key failure in ioStrStrFunctor with key " << key << endl;
        throw;
    };
};

ioFirst::operator bool() { 
    if (is_first) {
        is_first = false;
        return true;
    }
    return false;
};
bool ioFirst::operator()() {
    if (is_first) {
        is_first = false;
        return true;
    }
    return false;
};

ioCycleTrue::ioCycleTrue(int period_in) :
    period { period_in }, cnt{0}
{};
bool ioCycleTrue::operator()() {
    ++ cnt;
    if (cnt == period) {
        cnt = 0;
        return true;
    } else {
        return false;
    }
};
ioCycleTrue::operator bool() { return this->operator()(); };
void ioCycleTrue::reset() { cnt = 0; };

ioXYbounder::ioXYbounder(vector<double> x, vector <double> y, ioOptMap opt) :
    X{x}, Y{y},
    size { (int) X.size() },
    lodef{ opt.has("default-lo") ? opt["default-lo"]() : 
           size > 0 ? Y[0] : 0.
    },
    hidef{ opt.has("default-hi") ? opt["default-hi"]() : 
           size > 0 ? Y[size-1] : 0.
    }
{};

ioXYbounder::ioXYbounder() : X {}, Y{}, size{0}, lodef{0.}, hidef{0.}
{};

ioXYbounder::ioXYbounder(
    const char* file, const char* tagX, 
    const char* tagY, ioOptMap opt
) :
    X { ioReadValVec(file, {{"tag",tagX,"sort",true}}) },
    Y { ioReadValVec(file, {{"tag",tagY,"sort",false}}) },
    size { (int) X.size() },
    lodef{ opt.has("default-lo") ? opt["default-lo"]() : 
           size > 0 ? Y[0] : 0.
    },
    hidef{ opt.has("default-hi") ? opt["default-hi"]() : 
           size > 0 ? Y[size-1] : 0.
    }
{};

bool ioXYbounder::operator()(double x, double y) {
    if (size == 0) return false;
    int bin = (int)(std::lower_bound(X.begin(), X.end(), x) - X.begin());
    /* cout << " bin: " << bin << endl; */
    if (bin == 0) return y>lodef;
    else if (bin == size) return y>hidef;
    else return y>Y[bin];
};

double ioXYbounder::operator()(double x) {
    if (size == 0) return -1;
    int bin = (int)(std::lower_bound(X.begin(), X.end(), x) - X.begin());
    if (bin == 0) return lodef;
    else if (bin == size) return hidef;
    else return Y[bin];
};

ioCycleSpacer::ioCycleSpacer(int period, int _n_width, const char* _spacer) :
    cycle{period}, spacer{_spacer}, n_width{_n_width} 
{
    cycle.cnt=-1;
};

ostream& operator<<(ostream& os, ioCycleSpacer& cs) {
    if (cs.cycle) os << endl;
    os << cs.spacer;
    if (cs.n_width) os << setw(cs.n_width);
    return os;
};

void ioCycleSpacer::reset() { cycle.cnt=-1; };

