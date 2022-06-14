#include "tuClass.h"
#include "tu_fmt.h"
#include "tu_fnc.h"

#include <sstream>
#include <algorithm>
#include "TString.h" // not used, but present for some header 
#include "TMath.h"
#include <iostream>
#include <string>
#include <set>
#include "TTreeReader.h"
#include "tuConst.h"

using namespace std;

tuGetter::tuGetter() : n_objects{0} {};

TFile* tuGetter::get_file(string f_name) {
    if (f_name.find(".root") == string::npos) f_name += ".root";
    TFile *f = ( files.count(f_name) ? files[f_name] : new TFile(f_name.c_str(), "read") );
    if (!f->IsOpen()) {
        cout << " Fatal error, cannot open file: " << f_name << endl;
        exit (1);
    }
    return f;
};

TObject* tuGetter::operator()(string f_name, string object_name) {
    TFile *s_current = gDirectory->GetFile();
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

    if (s_current!=nullptr) s_current->cd();
    return obj;
};

vector<double> tuBinVec::bin_centers() {
    vector<double> V;
    for (int i{0}; i<(int)vec.size()-1; ++i) V.push_back(0.5*(vec[i]+vec[i+1]));
    return V;
};
tuBinVec::tuBinVec(TH1* hg, const char xyz) {
    TAxis *ax = (xyz == 'x' ? hg->GetXaxis() :
                 xyz == 'y' ? hg->GetYaxis() :
                              hg->GetZaxis());
    vector<double> build_vec;
    for (int i{1}; i<=ax->GetNbins()+1; ++i) build_vec.push_back(ax->GetBinLowEdge(i));
    init(build_vec);
};
tuBinVec::tuBinVec(TAxis* ax) {
    vector<double> build_vec;
    for (int i{1}; i<=ax->GetNbins()+1; ++i) build_vec.push_back(ax->GetBinLowEdge(i));
    init(build_vec);
};
void tuBinVec::build_ptr() {
    size = vec.size();
    ptr = new double[size];
    for (int i{0}; i<size; ++i) ptr[i] = vec[i];
};
tuBinVec::operator int () const { return size-1; };
tuBinVec::operator double* () const { return ptr; };
tuBinVec::operator vector<double> () { return vec; };

tuBinVec::tuBinVec(vector<double> V) { init(V); };
tuBinVec::tuBinVec(const char* file, tuOptMap options) {
    init( tuReadValVec(file, options) );;
};
tuBinVec::tuBinVec(const char* file, const char* tag, tuOptMap options){
    options("tag") = tag;
    init( tuReadValVec(file, options) ) ;
};

int tuBinVec::bin_from_0(double val) {
    auto loc = upper_bound(vec.begin(), vec.end(), val);
    return loc-vec.begin()-1;
};
int tuBinVec::bin_from_1(double val) {
    return bin_from_0(val) + 1;
};
// copy constructor
tuBinVec::tuBinVec(const tuBinVec& cp) { init(cp.vec); };
void tuBinVec::init(vector<double> V) {
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
            if (i>(S-3)) throw std::runtime_error( "fatal in tuBinVec with range_repeat");
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
tuBinVec::~tuBinVec() {
    delete[] ptr;
};
/* int tuBinVec::nbins() { return (int) size-1; }; */
vector<double>::const_iterator tuBinVec::begin() const { return vec.begin(); };
vector<double>::const_iterator tuBinVec::end()   const { return vec.end(); };
double tuBinVec::operator[](int i) { return vec[i]; };
double tuBinVec::bin_underflow() { 
    if (vec.size()<2)  return 0.;
    return vec[0]-(vec[1]-vec[0]);
};
double tuBinVec::bin_overflow() { 
    if (vec.size()<2)  return 0.;
    int i { static_cast<int>(vec.size())-1 };
    return vec[i]+(vec[i]-vec[i-1]);
};
ostream& operator<<(ostream& os, tuBinVec& tu) {
    for (auto v : tu) cout << " " << v;
    return os;
};

bool tuInBounds::operator()(double x) {
    return (x >= lo_bound && x <= hi_bound);
};
void tuInBounds::init(tuBinVec bins) {
    if (bins.vec.size()<2) {
        throw std::runtime_error(
        "Fatal: tried to initialize a tuInBounds with "
        "an input tuBinVec with less than 2 entries!"
        );
    };
    lo_bound = bins[0];
    hi_bound = bins[bins.vec.size()-1];
};
tuInBounds::tuInBounds(const char* file, const char* tag){
    init( {file, tag} );
};
tuInBounds::tuInBounds(tuBinVec bins) { init(bins); };
tuInBounds::tuInBounds(double lo, double hi) :
    lo_bound{lo}, hi_bound(hi) {}; 

/* ostream& operator<<(ostream& os, tuInBounds& rhs) { */
/*     os << " tuInBounds("<<rhs.lo_bound<<","<<rhs.hi_bound<<")"<<endl; */
/*     return os; */
/* }; */

// tuPads class (with helper class tuPadDim)
void tuPadDim::check_input() {
    if (   low   < 0. || low   > 1. 
            || p_low < 0. || p_low > 1. 
            || p_up  < 0. || p_up  > 1. 
            || up    < 0. || up    > 1. ) {
        cout << " Fatal error: input coordinates for tuPadDim for pads must all "
            " be in range [0,1] " << endl;
        print();
        exit (2);
    } else if ( low > p_low || p_low > p_up || p_up > up ) {
        cout << " Fatal error: input coordinates must monotonically increase " << endl;
        print();
        exit(2);
    }
};

tuPadDim::tuPadDim( double _low, double _p_low, double _p_up, double _up ) :
    low{_low}, p_low{_p_low}, p_up{_p_up}, up{_up} { check_input(); };
tuPadDim::tuPadDim( double _low, double _up ) : 
    low{_low}, p_low{_low}, p_up{_up}, up{_up} { check_input(); };
tuPadDim::tuPadDim( double _low, double _p_low, double _up ) : 
    low{_low}, p_low{_p_low}, p_up{_up}, up{_up} { check_input(); };
tuPadDim::tuPadDim( ) :
    low{0.}, p_low{0.}, p_up{1.}, up{1.} { check_input(); };


void tuPadDim::print() const {
    cout << Form(" Four points are: (%.2f, %.2f %.2f, %.2f)",low,p_low,p_up,up) << endl;
};

double tuPadDim::low_margin () const {
    double margin { (p_low - low) / (up - low) };
    if (margin < 0) margin = 0;
    return margin;
};
double tuPadDim::up_margin () const {
    // use to get set the lower margin
    double margin { (up - p_up) / (up - low) };
    if (margin < 0) margin = 0;
    return margin;
};

bool tuPadDim::operator==(tuPadDim& B) const {
    return low == B.low
        && p_low == B.p_low
        && p_up  == B.p_up
        && up    == B.up;
};

tuPadDimSet::tuPadDimSet(vector<double> _lefts, vector<double> _rights ) :
          rights{_rights} 
{
    if (_lefts.size() == 0) nPads = 1;
    else if (_lefts[0] >= 1.) {
        nPads = (int) _lefts[0];
        for (int i{0}; i<(int)_lefts.size()-1; ++i) lefts.push_back(_lefts[i+1]);
    } else {
        nPads = 1;
        lefts = _lefts;
    }
};

tuPadDim tuPadDimSet::make_pad(double left, 
            double left_margin, double pad_width, 
            double right_margin) 
{
    return tuPadDim{ left, 
                     left+left_margin,
                     left+left_margin+pad_width,
                     left+left_margin+pad_width+right_margin };
};

vector<tuPadDim> tuPadDimSet::calc_pads() {
    int npads = nPads;
    bool flip_direction = false;
    if (npads < 0) { 
        npads = -npads; 
        flip_direction=true;
    };
    vector<tuPadDim> pads (npads) ;

    double first_left = (lefts.size() > 0) ? lefts[0] : 0.2;
    double inner_left = (lefts.size() > 1) ? lefts[1] : 0.0001;
    double page_left  = (lefts.size() > 2) ? lefts[2] : 0.01;

    double last_right  = (rights.size() > 0) ? rights[0] : 0.0001;
    double inner_right = (rights.size() > 1) ? rights[1] : 0.0;
    double page_right  = (rights.size() > 2) ? rights[2] : 0.01;

    if (npads == 0) throw std::runtime_error(
        "fatal in tuPadDimSet must request at least one pad");
    if (npads == 1) {
        double pad_width = 1.-first_left-page_left-last_right-page_right;
        if (pad_width<=0) throw std::runtime_error(
                "fatal in tuPadDimSet margins have consumed more than 100\% of TCanvas");
        pads[0] = make_pad( page_left, first_left, pad_width, last_right );
        return pads;
    } 

    double pad_width { (1.-(first_left+page_left+last_right+page_right+
            (inner_left+inner_right)*(npads-1)))/npads };
    if (pad_width<=0) throw std::runtime_error(
            "fatal in tuPadDimSet margins have consumed more than 100\% of TCanvas");

    int index = flip_direction ? npads-1 : 0;
    pads[index] = make_pad(page_left, first_left, pad_width, inner_right);
    double left = pads[index].up;

    for (int i=1;i<npads-1;++i) {
        int index = flip_direction ? npads-i-1 : i;
        pads[index] = make_pad(left, inner_left, pad_width, inner_right);
        left = pads[index].up;
    }
    pads[flip_direction ? 0 : npads-1] = make_pad(left, inner_left, pad_width, last_right);
    return pads;
};

tuPads::tuPads ( int nYpads, vector<double> dimensions, int nXpads ) {
    // build the tuPadDimSet out of dimensions
    tuPadDimSet xPads{ {0.2, 0.0001, 0.01}, {0.0001,0.0,0.01 }};
    tuPadDimSet yPads{ {0.2, 0.0001, 0.01}, {0.0001,0.0,0.01 }};

    int which = 0;
    int cnt   = 0;
    canvas_width  = -1;
    canvas_height = -1;
    for (auto val : dimensions) {
        if (val > 6) { // set first dimensions
            if (canvas_width == -1) canvas_width  = val;
            else                    canvas_height = val;
            continue;
        } else if (val > 1) { 
            which = (int) val;
            cnt = 0;
            continue;
        } else if (which == 0) {
            which = kLeft;
        }
        switch (which) {
            case 6:  // kTop
                yPads.rights[cnt] = val;
                break;
            case 5: // kBottom
                yPads.lefts[cnt] = val;
                break;
            case kLeft:
                xPads.lefts[cnt] = val;
                break;
            case kRight:
                xPads.rights[cnt] = val;
                break;
            default:
                throw std::runtime_error(Form("fatal error: tuPads::tuPads: Error in selection of pad dimensions: was %i but must be (2,3,5,6:kLeft,Right,Bottom,Top)",
                            which));
        }
        ++cnt;
    }
    yPads.nPads = -nYpads;
    xPads.nPads =  nXpads;
    if (canvas_width  == -1) canvas_width = 1200;
    if (canvas_height == -1) canvas_height =  800;
    nCol = TMath::Abs(nXpads);
    nRow = TMath::Abs(nYpads);
    for (auto x_pad : xPads.calc_pads())
        for (auto y_pad : yPads.calc_pads())
            pad_dimensions.push_back( {y_pad, x_pad} );
};

                
tuPads::tuPads ( vector<pair<tuPadDim, tuPadDim>> _pad_dimensions, int
        _canvas_width, int _canvas_height) :
    pad_dimensions{ _pad_dimensions }
{
    if (_canvas_width)  canvas_width  = _canvas_width;
    if (_canvas_height) canvas_height = _canvas_height;
};
/* tuPads::tuPads ( int nYpads, int nXpads, int c_wide, int c_height, */ 
/*             tuPadDimSet Ypads, tuPadDimSet Xpads) { */

/*     if (c_wide)   canvas_width = c_wide; */
/*     if (c_height) canvas_height = c_height; */

/*     Ypads.nPads = nYpads; */
/*     Xpads.nPads = nXpads; */

/*     nCol = TMath::Abs(nXpads); */
/*     nRow = TMath::Abs(nYpads); */
    
/*     for (auto x_pad : Xpads.calc_pads()) */
/*         for (auto y_pad : Ypads.calc_pads()) */
/*             pad_dimensions.push_back( {y_pad, x_pad} ); */
/* }; */
tuPads::tuPads( vector<tuPadDim> y_dim, vector<tuPadDim> x_dim, int c_wide, int c_height) {
    if (x_dim.size()==1) {
        tuPadDim temp_pad;
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
        vector<tuPadDim> temp;
        for (int i{(int)y_dim.size()-1};i>=0;--i) temp.push_back(y_dim[i]);
        y_dim = temp;
    }
    for (auto& x : x_dim) 
        for (auto& y : y_dim)
            pad_dimensions.push_back({y,x});

    nCol = x_dim.size();
    nRow = y_dim.size();
    if (c_wide)   canvas_width = c_wide;
    if (c_height) canvas_height = c_height;
};
TPad* tuPads::operator()(int row, int col) {
    if (pads.size() == 0) init();
    if (row < 0) {
        row = -row;
        col = row % nCol;
        row = row / nCol;
    }
    int i_pad = row+col*nRow;
    if (i_pad >= (int)pads.size()) {
        i_pad = i_pad % (int) pads.size();
    }
    pads[i_pad]->cd();
    return pads[i_pad];
};
void tuPads::init() {
    // make and stylize the TCanvas and pads currently in the list
    const char* t_name = Form("canv_%s",tuUniqueName());
    canvas = new TCanvas(t_name, "",canvas_width, canvas_height);
    tu_fmt(canvas);
    canvas->Draw();
    canvas->cd();

    // add all pads
    add_pad(pad_dimensions);

    canvas->cd();
    /* canvas_pad = new TPad("canvas_pad","",0.,0.,1.,1.); */
    /* tu_fmt(canvas_pad); */
    /* canvas_pad->Draw(); */
};


void tuPads::add_pad(pair<tuPadDim,tuPadDim>& coord){
    canvas->cd();

    if (pads.size()==0) {
        canvas_pad = new TPad(tuUniqueName(),"",0.,0.,1.,1.);
        tu_fmt(canvas_pad);
        canvas_pad->Draw();
        canvas->cd();
    }

    const tuPadDim x { coord.second };
    const tuPadDim y { coord.first  };
    int i{0};
    /* while (gDirectory->FindObjectAny(Form("loc_pad_%i",i))) { ++i; } */
    TPad* p = new TPad(tuUniqueName(),"",x.low,y.low,x.up,y.up);

    // set the boundaries left(l), right(r), top(t), bottom(b)
    p->SetLeftMargin(x.low_margin());
    p->SetRightMargin(x.up_margin());
    p->SetBottomMargin(y.low_margin());
    p->SetTopMargin(y.up_margin());

    tu_fmt(p);
    p->Draw();
    pads.push_back(p);
};

void tuPads::add_pad(vector<pair<tuPadDim,tuPadDim>> input) {
    for (auto& inp : input) add_pad(inp);
};

// Implementation of tuIntList
string tuIntList::make(const char* in_file, bool print) {
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
void tuPads::stamp(const char* msg, tuOptMap options, tuOptMap dict) {
    dict += options;
    canvas_pad->cd();
    /* cout << " x: " << dict["x-loc"] << "  " << dict["y-loc"] << endl; */
    tuDrawTLatex(msg,dict("x-loc"), dict("y-loc"), dict);
};
void tuPads::save(const char* name, const char* tag) {
    TString check {name};
    if (!check.Contains(".")) check.Append(Form("%s.pdf",tag));
    else if (check.EndsWith(".cc")) check.ReplaceAll(".cc",Form("%s.pdf",tag));
    else if (check.EndsWith(".C" )) check.ReplaceAll(".C", Form("%s.pdf",tag));
    canvas->Print(check.Data());
};

/* tuIntList::tuIntList(const char* in_file, ofstream& log, bool print) { */
/*     log << make(in_file, print); */
/* }; */
/* tuIntList::tuIntList(const char* in_file, bool print) { make(in_file, print); }; */

/* bool tuIntList::operator()(int val) { */
/*     return std::binary_search(list.begin(), list.end(), val); */
/* }; */
/* int tuIntList::operator[](int val) { */
/*     return (int)(std::lower_bound(list.begin(), list.end(), val) - list.begin()); */
/* }; */
/* bool tuIntList::has(int val) { return this->operator()(val); }; */
/* bool tuIntList::has_not(int val) { return !(this->operator()(val)); }; */


/* bool tuRunListId::has_run(int id) { */
/*     return static_cast<bool>(map_id.count(id)); */
/* }; */

/* int tuRunListId::size() { return map_id.size(); }; */

/* double tuRunListId::operator()(int run_id) { */
/*     set_id(run_id); */
/*     return id; */
/* }; */

/* double tuRunListId::set_id(int run_id) { */
/*     if (has_run(run_id)) { */
/*         id = static_cast<double>(map_id[run_id]); */
/*     } else { */
/*         id = -1.; */
/*     } */
/*     return id; */
/* }; */

/* // Implementation of tuIntMap */
/* tuIntMap::tuIntMap(const char* file, */
/*         int index_column, */ 
/*         int data_column, */
/*         bool echo_print, */
/*         vector<int> skip_vals */
/*         ) { */
/*     tuIntMap_constructor(file, index_column, data_column, */ 
/*             echo_print, skip_vals); */
/* }; */
/* tuIntMap::tuIntMap(const char* file, */
/*         int index_column, */ 
/*         int data_column, */
/*         bool echo_print, */
/*         ofstream& log, */
/*         vector<int> skip_vals */
/*         ) { */
/*     log << tuIntMap_constructor(file, index_column, data_column, */ 
/*             echo_print, skip_vals); */
/* }; */
/* string tuIntMap::tuIntMap_constructor ( */
/*         const char* in_file, */
/*         int index_column, */ 
/*         int data_column, */
/*         bool echo_print, */
/*         vector<int> skip_vals */
/*         ) { */
/*     set<int> skip_val_set; */
/*     for (auto& v : skip_vals) skip_val_set.insert(v); */
/*     /1* string tuIntMap::make(const char* in_file, bool print) { *1/ */
/*     ostringstream msg; */
/*     ifstream file; */
/*     file.open(in_file); */
/*     if (!file.is_open()) { */
/*         msg << "Could not open int map file \"" << in_file */ 
/*             << "\". No entries entered." << endl; */
/*         cout << msg.str() << endl; */
/*         return msg.str(); */
/*     } else { */
/*         if (echo_print) msg << " Reading file " << in_file << " for map of col(" */
/*             << index_column << ") to col(" << data_column << ")" << endl; */
/*     } */
/*     /1* int n_req { index_column > data_column ? index_column : data_column }; *1/ */
/*     string line; */
/*     while (getline(file,line)) { */
/*         line.append(" "); */
/*         stringstream words(line); */
/*         TString word; */
/*         bool has_index { false }; */
/*         bool has_data  { false }; */
/*         int  val_index{-1}; */
/*         int  val_data{-1}; */
/*         int i {0}; */
/*         bool comment_flag {false}; */
/*         int n_words{0}; */
/*         while (words >> word) { */
/*             /1* cout << " Line: " << line << "  ->(word): " << word << "  IsAlnum: " << word.IsAlnum()<< endl; *1/ */
/*             /1* cout << " word: " << word << endl; *1/ */
/*             if (word.BeginsWith("//") || word.BeginsWith("#")) { */
/*                 comment_flag = true; */
/*                 break; */
/*             } */
/*             // ignore columns that start with a non-number word */
/*             if (n_words == 0 && !word.IsDigit()) { comment_flag = true; break; } */
/*             ++n_words; */
/*             if (i == index_column) { */
/*                 has_index = true; */
/*                 val_index = word.Atoi(); */
/*                 if (skip_val_set.count(val_index)) continue; */
/*             } */ 
/*             if ( i == data_column) { */
/*                 has_data = true; */
/*                 val_data = word.Atoi(); */
/*             } */
/*             if (has_index && has_data) break; */
/*             ++i; */
/*         } */
/*         if (comment_flag) continue; */
/*         if (!has_index || !has_data) { */
/*             /1* ostringstream loc_msg; *1/ */
/*             msg << "fatal error in reading tuIntMap line from file: " */
/*                 << "  -> " << in_file << endl; */
/*             if (!has_index) */ 
/*                 msg << " Couldn't read index column("<<index_column */ 
/*                     <<")" <<endl; */
/*             if (!has_data) */ 
/*                 msg << " Couldn't read data column("<<data_column */ 
/*                     <<")" <<endl; */
/*             msg << "  from line:" << endl */
/*                 << "  -> " << line << endl; */
/*             msg << "  Therefore skipping entries on this line." << endl; */
/*         } else { */
/*             if (echo_print) msg << "  " << val_index */ 
/*                 << " -> " << val_data << endl; */
/*             data_map[val_index] = val_data; */
/*         } */
/*     } */
/*     file.close(); */
/*     msg << " Done reading columns from file \"" << in_file <<"\"" << endl; */
/*     cout << msg.str(); */
/*     return msg.str(); */
/* }; */

/* bool tuIntMap::has(int key) { */
/*     return (bool) data_map.count(key); */
/* }; */

/* int& tuIntMap::operator[](int key) { return data_map[key]; }; */

/* vector<int> tuIntMap::keys() { */
/*     vector<int> vec; */
/*     for (auto m : data_map) vec.push_back(m.first); */
/*     sort(vec.begin(), vec.end()); */
/*     return vec; */
/* }; */
/* int tuIntMap::size() { return data_map.size(); }; */


/* // tuHgStats */
/* tuHgStats::tuHgStats(vector<double> ax_vals, vector<double> _vals, vector<double> _errs, bool cut_zeros) { */
/*     if (ax_vals.size() != _vals.size() || _vals.size() != _errs.size()) */
/*         throw std::runtime_error( */
/*                 "tuHgStats(vec, vec, vec) required vectors of same length"); */
/*     if (ax_vals.size() < 2) */ 
/*         throw std::runtime_error( */
/*                 "tuHgStats(vec, vec, vec) required vectors size > 1"); */
/*     vals = _vals; */
/*     errs = _errs; */
/*     nbins = ax_vals.size(); */

/*     double *p = new double[nbins+1]; */
/*     p[0] = ax_vals[0]-0.5*(ax_vals[1]-ax_vals[0]); */
/*     for (int i{1}; i<nbins; ++i) p[i] = 0.5*(ax_vals[i-1]+ax_vals[i]); */
/*     p[nbins] = ax_vals[nbins-1]+0.5*(ax_vals[nbins-1]-ax_vals[nbins-2]); */
/*     axis = new TAxis(nbins, p); */

/*     if (cut_zeros) { */
/*         for (auto i{0}; i<nbins;++i) { */
/*             weight.push_back( vals[i] == 0 ? 0. : 1.); */
/*         } */
/*     } else { */
/*         for (auto i{0}; i<nbins; ++i) { */
/*             weight.push_back(1.); */
/*         } */
/*     } */
/*     calc_stats(); */
/* }; */
/* tuHgStats::tuHgStats(TH1D* hg, bool cut_zeros) { */
/*     vals = tu_vecBinContent(hg); */
/*     errs = tu_vecBinError  (hg); */
/*     axis=hg->GetXaxis(); */
/*     nbins = axis->GetNbins(); */
/*     if (cut_zeros) { */
/*         for (auto i{0}; i<nbins;++i) { */
/*             weight.push_back( vals[i] == 0 ? 0. : 1.); */
/*         } */
/*     } else { */
/*         for (auto i{0}; i<nbins; ++i) { */
/*             weight.push_back(1.); */
/*         } */
/*     } */
/*     calc_stats(); */
/* }; */
/* tuHgStats::tuHgStats(TProfile* hg, bool cut_zeros, bool weight_by_entries) { */
/*     vals = tu_vecBinContent(hg); */
/*     errs = tu_vecBinError  (hg); */
/*     axis=hg->GetXaxis(); */
/*     nbins = axis->GetNbins(); */
/*     if (weight_by_entries) { */
/*         weight = tu_vecBinEntries(hg); */
/*     } else { */
/*         if (cut_zeros) { */
/*             for (auto i{0}; i<nbins;++i) { */
/*                 weight.push_back( vals[i] == 0 ? 0. : 1.); */
/*             } */
/*         } else { */
/*             for (auto i{0}; i<nbins; ++i) { */
/*                 weight.push_back(1.); */
/*             } */
/*         } */
/*     } */
/*     calc_stats(); */
/* }; */
/* void tuHgStats::calc_stats() { */
/*     mean   = TMath::Mean(vals.begin(), vals.end(), weight.begin()); */
/*     stddev = TMath::StdDev(vals.begin(), vals.end(), weight.begin()); */
/* }; */
/* double tuHgStats::mean_Xsigma(double X) { */
/*     /1* double sumW{0}; *1/ */
/*     /1* double sumV{0}; *1/ */
/*     /1* double sumE{0}; *1/ */
/*     /1* for (auto w : weight) sumW += w; *1/ */
/*     /1* for (auto w : vals) sumV += w; *1/ */
/*     /1* for (auto w : errs) sumE += w; *1/ */
/*     /1* cout << " mean: " << mean << "  stddev " << stddev << "  X " << X << "  val: " << *1/ */
/*     /1* mean+X*stddev << " >> npts" << nbins << Form("stats:%f,%f,%f",sumV,sumE,sumW) <<endl; *1/ */
/*     return mean+ X * stddev; */
/* }; */
/* TLine* tuHgStats::get_horizontal_TLine(double cut, bool cut_times_sigma) { */
/*     double x0 = axis->GetBinLowEdge(1); */
/*     double x1 = axis->GetBinUpEdge(axis->GetNbins()); */
/*     double y = cut_times_sigma ? mean_Xsigma(cut) : cut; */
/*     return new TLine(x0,y,x1,y); */
/* }; */
/* TGraph* tuHgStats::points_above(double cut, bool cut_times_sigma) { */
/*     // count how many points are above */
/*     double y = cut_times_sigma ? mean_Xsigma(cut) : cut; */
/*     int n_above {0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && vals[i]>y) ++n_above; */
/*     } */
/*     double* xpts = new double[n_above]; */
/*     double* ypts = new double[n_above]; */
/*     int n{0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && vals[i]>y) { */
/*             xpts[n] = axis->GetBinCenter(i+1); */
/*             ypts[n] = vals[i]; */
/*             ++n; */
/*         } */
/*     } */
/*     return new TGraph(n_above,xpts,ypts); */
/* }; */
/* TGraph* tuHgStats::points_below(double cut, bool cut_times_sigma) { */
/*     // count how many points are above */
/*     double y = cut_times_sigma ? mean_Xsigma(cut) : cut; */
/*     int n_pts {0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && vals[i]<y) ++n_pts; */
/*     } */
/*     double* xpts = new double[n_pts]; */
/*     double* ypts = new double[n_pts]; */
/*     int n{0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && vals[i]<y) { */
/*             xpts[n] = axis->GetBinCenter(i+1); */
/*             ypts[n] = vals[i]; */
/*             ++n; */
/*         } */
/*     } */
/*     return new TGraph(n_pts,xpts,ypts); */
/* }; */
/* TGraph* tuHgStats::points_between(double locut, double hicut, bool cut_times_sigma) { */
/*     // count how many points are above */
/*     double y_lo = cut_times_sigma ? mean_Xsigma(locut) : locut; */
/*     double y_hi = cut_times_sigma ? mean_Xsigma(hicut) : hicut; */
/*     int n_pts {0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && vals[i]>y_lo && vals[i]<y_hi) ++n_pts; */
/*     } */
/*     double* xpts = new double[n_pts]; */
/*     double* ypts = new double[n_pts]; */
/*     int n{0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && vals[i]<y_hi && vals[i]>y_lo) { */
/*             xpts[n] = axis->GetBinCenter(i+1); */
/*             ypts[n] = vals[i]; */
/*             ++n; */
/*         } */
/*     } */
/*     return new TGraph(n_pts,xpts,ypts); */
/* }; */
/* // cuts */
/* tuHgStats& tuHgStats::cut_above(double cut, bool cut_times_sigma) { */
/*     // count how many points are above */
/*     double y = cut_times_sigma ? mean_Xsigma(cut) : cut; */
/*     int n_above {0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && vals[i]>y) weight[i] = 0.; */
/*     } */
/*     calc_stats(); */
/*     return *this; */
/* }; */
/* tuHgStats& tuHgStats::cut_below(double cut, bool cut_times_sigma) { */
/*     // count how many points are above */
/*     double y = cut_times_sigma ? mean_Xsigma(cut) : cut; */
/*     int n_pts {0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && vals[i]<y) weight[i] = 0.; */
/*     } */
/*     calc_stats(); */
/*     return *this; */
/* }; */
/* vector<double> tuHgStats::unmasked_vals() { */
/*     vector<double> r_vec; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0.) r_vec.push_back(vals[i]); */
/*     } */
/*     return r_vec; */
/* }; */
/* tuHgStats& tuHgStats::mask(const vector<bool>  mask) { */
/*     if ((int)mask.size() != nbins) { */
/*         cout << " error in tuHgStats::mask " << endl */
/*             << " There are " << nbins << " bins, but only " << mask.size() */
/*             << " points in in put mask " << endl << endl */
/*             << " Therefore mask is not applied." << endl; */
/*         return *this; */
/*     } */
/*     for (int i{0}; i<nbins; ++i) { */
/*         if (!mask[i]) weight[i] = 0.; */
/*         calc_stats(); */
/*     } */
/*     return *this; */
/* }; */
/* /1* tuHgStats& tuHgStats::mask(vector<int>& mask, bool mask_keep_true) { *1/ */
/* /1*     if ((int)mask.size() != nbins) { *1/ */
/* /1*         cout << " error in tuHgStats::mask " << endl *1/ */
/* /1*              << " There are " << nbins << " bins, but only " << mask.size() *1/ */
/* /1*              << " points in in put mask " << endl << endl *1/ */
/* /1*              << " Therefore mask is not applied." << endl; *1/ */
/* /1*         return *this; *1/ */
/* /1*     } *1/ */
/* /1*     for (int i{0}; i<nbins; ++i) { *1/ */
/* /1*         if ((mask[i] == 0)) weight[i] = 0.; *1/ */
/* /1*         calc_stats(); *1/ */
/* /1*     } *1/ */
/* /1*     return *this; *1/ */
/* /1* }; *1/ */
/* tuHgStats& tuHgStats::cut_to_range(double locut, double hicut, bool cut_times_sigma) { */
/*     // count how many points are above */
/*     double y_lo = cut_times_sigma ? mean_Xsigma(locut) : locut; */
/*     double y_hi = cut_times_sigma ? mean_Xsigma(hicut) : hicut; */
/*     int n_pts {0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0. && (vals[i]<y_lo || vals[i]>y_hi)) weight[i] = 0.; */
/*     } */
/*     calc_stats(); */
/*     return *this; */
/* }; */
/* TGraph* tuHgStats::points() { */
/*     int n_pts {0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0.) ++n_pts; */
/*     } */
/*     double* xpts = new double[n_pts]; */
/*     double* ypts = new double[n_pts]; */
/*     int n{0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0.) { */
/*             xpts[n] = axis->GetBinCenter(i+1); */
/*             ypts[n] = vals[i]; */
/*             ++n; */
/*         } */
/*     } */
/*     return new TGraph(n_pts,xpts,ypts); */
/* }; */
/* int tuHgStats::count_points() { */
/*     int n_pts {0}; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0.) ++n_pts; */
/*     } */
/*     return n_pts; */
/* }; */
/* vector<int> tuHgStats::bin_indices() { */
/*     vector<int> vec; */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (weight[i]!=0.) vec.push_back(i+1); */
/*     } */
/*     return vec; */
/* }; */
/* tuHgStats& tuHgStats::restore_points() { */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         weight[i] = 1.; */
/*     } */
/*     calc_stats(); */
/*     return *this; */
/* }; */
/* tuHgStats& tuHgStats::cut_zeros() { */
/*     for (auto i{0}; i<nbins; ++i) { */
/*         if (vals[i]==0) weight[i] =0.; */
/*     } */
/*     calc_stats(); */
/*     return *this; */
/* }; */

/* // ------------------------------------------------------ */
/* // | Implementation of tuMsgTree                        | */
/* // ------------------------------------------------------ */
/* tuMsgTree::tuMsgTree(bool set_echo) : */ 
/*     b_msg{""}, */ 
/*     tree{"Messages", "Tree of messages"} */
/* { */
/*     tree.Branch("messages", &b_msg); */
/*     dash(); */
/*     msg(" Start of msg tree "); */
/*     dash(); */
/*     echo_to_cout = set_echo; */
/* }; */
/* void tuMsgTree::msg(string msg) { */
/*     b_msg = msg; */
/*     if (echo_to_cout) cout << b_msg << endl; */
/*     tree.Fill(); */
/* }; */
/* void tuMsgTree::msg(vector<string> messages) { */
/*     for (auto& msg : messages) { */
/*         b_msg = msg; */
/*         if (echo_to_cout) cout << b_msg << endl; */
/*         tree.Fill(); */
/*     } */
/* }; */
/* void tuMsgTree::dash() { */
/*     b_msg = "---------------"; */
/*     if (echo_to_cout) cout << b_msg << endl; */
/*     tree.Fill(); */
/* }; */
/* void tuMsgTree::write(){ */
/*     tree.Write(); */
/* }; */
/* void tuMsgTree::read_messages(const char* f_name){ */
/*     cout << " Reading file: " << f_name << endl; */
/*     /1* TTree *tree; *1/ */
/*     TFile* fin  = new TFile(f_name, "read"); */
/*     if (!fin->IsOpen()) { */
/*         cout << " Input file: " << f_name << " is not open." << endl; */
/*         delete fin; */
/*         return; */
/*     } */
/*     TTreeReader myReader("Messages",fin); */
/*     TTreeReaderValue<string> msg(myReader, "messages"); */
/*     cout << "  Contents of TFile(\""<<f_name<<"\") TTree(\"Messages\"):" << endl; */
/*     while (myReader.Next()) cout << " " << *msg << endl; */
/*     fin->Close(); */
/*     delete fin; */
/*     /1* } *1/ */
/*     }; */
/* void tuMsgTree::slurp_file(const char* which_file) { */
/*     // try and read all lines of which_file into the tree */
/*     msg(Form("--Begin contents of file \"%s\"",which_file)); */
/*     ifstream f_in {which_file}; */
/*     if (!f_in.is_open()) { */
/*         msg(Form("  error: couldn't open file \"%s\"",which_file)); */
/*     } else { */
/*         string line; */
/*         while (getline(f_in,line)) msg(line); */
/*     } */
/*     msg(Form("--End contents of file \"%s\"",which_file)); */
/*     return; */
/* }; */

/* vector<int> tuIntVec::vals(string tag) { */
/*     assert_tag(tag,"vals"); */
/*     int col { i_tag(tag) }; */
/*     vector<int> r_vec {}; */
/*     for (auto i{0}; i<size(); ++i) r_vec.push_back(data[i][col]); */
/*     return r_vec; */
/* }; */

/* vector<int> tuIntVec::vals(string tag, vector<bool> mask) { */
/*     assert_tag(tag,"vals"); */
/*     int col { i_tag(tag) }; */

/*     if (mask.size() != data.size()) */
/*         throw std::runtime_error( */
/*                 Form("fatal error in tuIntVec::vals: size of mask and data do not match")); */

/*     vector<int> r_vec {}; */
/*     for (auto i{0}; i<(int)mask.size(); ++i) if (mask[i]) r_vec.push_back(data[i][col]); */
/*     return r_vec; */
/* }; */

/* bool tuIntVec::flip_tag(string& tag) { */
/*     int t_size = tag.size(); */
/*     if (t_size==0) return true; */
/*     if (tag.substr(0,1)=="!") { */
/*         tag = tag.substr(1,t_size-1); */
/*         return false; */
/*     } */
/*     return true; */
/* }; */

/* vector<bool> tuIntVec::mask(string tag, const char* module) { */
/*     bool on_true = flip_tag(tag); */
/*     assert_tag(tag,module); */
/*     int col { i_tag(tag) }; */
/*     vector<bool> r_vec; */
/*     for (auto& vec : data) r_vec.push_back((vec[col]!=0) == on_true); */
/*     return r_vec; */
/* }; */

/* vector<bool> tuIntVec::is_any(vector<string> _tags){ */
/*     int s = _tags.size(); */
/*     if (!s) return {}; */
/*     auto r_vec = mask(_tags[0]); */
/*     for (int i{1}; i<s; ++i) r_vec = r_vec || mask(_tags[i],"is_any"); */
/*     return r_vec; */
/* }; */
/* vector<bool> tuIntVec::is_all(vector<string> _tags){ */
/*     int s = _tags.size(); */
/*     if (!s) return {}; */
/*     auto r_vec = mask(_tags[0]); */
/*     for (int i{1}; i<s; ++i) r_vec = r_vec && mask(_tags[i],"is_all"); */
/*     return r_vec; */
/* }; */


/* tuIntVec::tuIntVec( const char* file_name, bool echo_print, vector<string>_tags) */ 
/* { tuIntVec_constructor(file_name, echo_print, _tags); }; */

/* tuIntVec::tuIntVec( const char* file_name, ofstream& log, */ 
/*         bool echo_print, vector<string>_tags) */ 
/* { */ 
/*     log << tuIntVec_constructor(file_name, echo_print, _tags); */ 
/*     log << *this; */
/* }; */

/* string tuIntVec::tuIntVec_constructor(const char* in_file, */ 
/*         bool echo_print, vector<string>tags_requested) */ 
/* { */
/*     // read until the tag line is read */
/*     // if given input tags, check they are all present, and use them in that order */
/*     // */ 

/*     ostringstream msg; */
/*     ifstream file; */
/*     file.open(in_file); */
/*     if (!file.is_open()) { */
/*         msg << "Could not open int to vector<int> map file \"" << in_file */ 
/*             << "\"." << endl << "  -> No entries entered." << endl; */
/*         cout << msg.str() << endl; */
/*         return msg.str(); */
/*     } else { */
/*         if (echo_print) msg << " Reading file " << in_file << endl; */
/*     } */

/*     // Read until get the column tags (first uncommented line) */

/*     string line; */
/*     /1* vector<string> tags_read; *1/ */
/*     while (getline(file,line)) { */
/*         /1* cout << " reading line: " << line << endl; *1/ */
/*         line.append(" "); */
/*         stringstream words(line); */
/*         TString word; */
/*         bool is_key_tag {true}; */
/*         while (words >> word) { */
/*             if (word.BeginsWith("//") || word.BeginsWith("#")){ */
/*                 break; */
/*             } */
/*             if (is_key_tag) { */
/*                 is_key_tag = false; */
/*             } else { */
/*                 tags.push_back(word.Data()); */
/*             } */
/*         } */
/*         if (tags.size() != 0) break; */
/*     } */
/*     if (tags_requested.size() == 0) { */
/*         tags_requested = tags; */
/*     } */ 
/*     auto read_cols = tag_cols(tags_requested); */
/*     int max_col = 0; */
/*     for (auto val : read_cols) if (val > max_col) max_col = val; */
/*     tags = tags_requested; */

/*     // read all remaining lines into data */
/*     vector<pair<int,vector<int>>> data_in; */
/*     while(getline(file,line)) { */
/*         vector<int> c_vec; */
/*         int key; */
/*         bool has_key {false}; */

/*         line.append(" "); */
/*         stringstream words(line); */
/*         TString word; */
/*         while (words >> word) { */
/*             if (word.BeginsWith("//") || word.BeginsWith("#")){ */
/*                 break; */
/*             } */
/*             if (!has_key) { */
/*                 has_key = true; */
/*                 key = word.Atoi(); */
/*             } else { */
/*                 c_vec.push_back(word.Atoi()); */
/*             } */
/*         } */
/*         if (has_key) { */
/*             data_in.push_back({key,c_vec}); */
/*             if ((int)c_vec.size() < (max_col+1)) */
/*                 throw std::runtime_error( */
/*                         Form("In tuIntVec needs at least %i entries (+id) in line \"%s\"", */
/*                             max_col+1, line.c_str()) */
/*                         ); */
/*         } */
/*     } */
/*     std::sort(data_in.begin(), data_in.end()); */

/*     for (auto entry : data_in) { */
/*         keys.push_back(entry.first); */
/*         vector<int> vec; */
/*         for (auto i : read_cols) vec.push_back(entry.second[i]); */
/*         data.push_back(vec); */
/*     }; */

/*     file.close(); */
/*     msg << " Done reading columns from file \"" << in_file <<"\"" << endl; */
/*     if (echo_print) { */
/*         cout << msg.str(); */
/*         cout << *this << endl; */
/*     } */
/*     return msg.str(); */
/* }; */

/* int tuIntVec::i_tag(string tag) { */
/*     auto it = std::find(tags.begin(), tags.end(), tag); */
/*     if (it == tags.end()) return -1; */
/*     else return (int)(it-tags.begin()); */
/* }; */

/* void tuIntVec::assert_tag(string tag, const char* module) { */
/*     if (!has_tag(tag)) */ 
/*         throw std::runtime_error(Form("fatal in tuIntVec::%s, couldn't find tag \"%s\"",module,tag.c_str())); */
/* }; */
/* vector<int> tuIntVec::tag_cols(vector<string> _tags, const char* name) { */
/*     vector<int> cols; */
/*     if (_tags.size()==0) { */
/*         for (int i{0}; i<(int)_tags.size(); ++i) cols.push_back(i); */
/*     } else { */
/*         for (auto& tag : _tags) { */
/*             assert_tag(tag,"tag_cols"); */
/*             cols.push_back(i_tag(tag)); */
/*         } */
/*     } */
/*     return cols; */
/* }; */
/* bool tuIntVec::has_tag(string tag) { */
/*     return ( i_tag(tag) != -1 ); */
/* }; */
/* tuIntVec& tuIntVec::swap_tags(string tag0, string tag1) { */
/*     // find the two tags among tags */
/*     int i0 = i_tag(tag0); */
/*     if (i0 == -1) { */ 
/*         cout << "fatal error in tuIntVec::swap_tags " */
/*             << "could not find tag \"" << tag0 << "\"" << endl; */
/*         return *this; */
/*     } */

/*     int i1 = i_tag(tag1); */
/*     if (i1 == -1) { */
/*         cout << "fatal error in tuIntVec::swap_tags " */
/*             << "could not find tag \"" << tag1 << "\"" << endl; */
/*         return *this; */
/*     } */

/*     for (auto& vec : data) { */
/*         auto temp = vec[i0]; */
/*         vec[i0] = vec[i1]; */
/*         vec[i1] = temp; */
/*     } */
/*     tags[i0] = tag1; */
/*     tags[i1] = tag0; */
/*     return *this; */
/* }; */
/* tuIntVec& tuIntVec::subset(vector<string> sub_tags) { */
/*     // cut down the table to only the subset of tags */
/*     auto keep_cols = tag_cols(sub_tags,"subset"); */
/*     for (auto& vec : data) { */
/*         vector<int> copy; */
/*         for (auto& i : keep_cols) { */
/*             copy.push_back(vec[i]); */
/*         } */
/*         vec = copy; */
/*     } */
/*     tags = sub_tags; */
/*     return *this; */
/* }; */
/* tuIntVec& tuIntVec::rm_tags(vector<string> tags) { */
/*     auto rm_cols = tag_cols(tags,"rm_tags"); */
/*     sort(rm_cols.begin(), rm_cols.end()); */
/*     vector<string> keep_tags {}; */
/*     for (int i{0}; i<(int)tags.size(); ++i) { */
/*         if (find(rm_cols.begin(), rm_cols.end(), i) == rm_cols.end()) */
/*             keep_tags.push_back(tags[i]); */
/*     } */
/*     return subset(tags); */
/* }; */
/* tuIntVec& tuIntVec::add_tag(string tag, int def_val) { */
/*     int i_col = i_tag(tag); */
/*     if (i_col != -1) { */
/*         cout << " error in tuIntVec::add_tag; added tag " << tag << " already exists! " << endl; */
/*         return *this; */
/*     } */
/*     tags.push_back(tag); */
/*     cout << " Adding tag: " << tag << endl; */
/*     for (auto& vec : data) vec.push_back(def_val); */
/*     /1* for (auto key : keys()) data_map[key].push_back(def_val); *1/ */
/*     return *this; */
/* }; */
/* void tuIntVec::rename_tag(string tag0, string tag1) { */
/*     assert_tag(tag0,"rename_tag"); */
/*     tags[i_tag(tag0)] = tag1; */
/* }; */
/* bool tuIntVec::has_key(int key) { */ 
/*     return binary_search(keys.begin(), keys.end(), key); */
/* }; */
/* vector<int>& tuIntVec::operator[](int key) { */ 
/*     if (!has_key(key)) { */
/*         throw std::runtime_error(Form("fatal: no key \"%i\" in tuIntVec",key)); */
/*     } */
/*     int i = (int)(std::lower_bound(keys.begin(), keys.end(), key) - keys.begin()); */
/*     return data[i]; */
/* }; */
/* /1* vector<int>    tuIntVec::keys() const { *1/ */
/* /1*     vector<int> vec; *1/ */
/* /1*     for (auto m : data_map) vec.push_back(m.first); *1/ */
/* /1*     sort(vec.begin(), vec.end()); *1/ */
/* /1*     return vec; *1/ */
/* /1* }; *1/ */
/* vector<int> tuIntVec::get_keys(vector<bool> mask, bool keep_on_true) { */
/*     if (mask.size() != data.size()) */
/*         throw std::runtime_error( */
/*                 Form("fatal error in tuIntVec::keys: size of mask and data do not match")); */

/*     vector<int> r_vec {}; */
/*     for (auto i{0}; i<(int)mask.size(); ++i) { */
/*         if (mask[i] == keep_on_true) r_vec.push_back(keys[i]); */
/*     } */
/*     return r_vec; */
/* }; */
/* /1* vector<bool> tuIntVec::is_any(vector<pair<string,bool>> mask_keep_true) { *1/ */
/* /1*     vector<bool>   v_keep_true; *1/ */
/* /1*     vector<string> v_tags; *1/ */
/* /1*     for (auto mpair : mask_keep_true) { *1/ */
/* /1*         v_keep_true.push_back(mpair.second); *1/ */
/* /1*         v_tags.push_back(mpair.first); *1/ */
/* /1*     } *1/ */
/* /1*     auto cols = tag_cols(v_tags, "is_any"); *1/ */
/* /1*     int  n {(int) cols.size()}; *1/ */
/* /1*     vector<bool> r_vec {}; *1/ */
/* /1*     for (auto& vec : data) { *1/ */
/* /1*         bool keep {false}; *1/ */
/* /1*         for (int i{0}; i<n; ++i) { *1/ */
/* /1*             if ( (vec[cols[i]] != 0) == v_keep_true[i]) { *1/ */
/* /1*                 keep = true; *1/ */
/* /1*                 break; *1/ */
/* /1*             } *1/ */
/* /1*         } *1/ */
/* /1*         r_vec.push_back(keep); *1/ */
/* /1*     } *1/ */
/* /1*     return r_vec; *1/ */
/* /1* }; *1/ */
/* int tuIntVec::size(){ return (int)keys.size(); }; */
/* /1* vector<bool> tuIntVec::is_all(vector<pair<string,bool>> mask_keep_true) { *1/ */
/* /1*     vector<bool>   v_keep_true; *1/ */
/* /1*     vector<string> v_tags; *1/ */
/* /1*     for (auto mpair : mask_keep_true) { *1/ */
/* /1*         v_keep_true.push_back(mpair.second); *1/ */
/* /1*         v_tags.push_back(mpair.first); *1/ */
/* /1*     } *1/ */
/* /1*     auto cols = tag_cols(v_tags, "is_all"); *1/ */
/* /1*     int  n {(int) cols.size()}; *1/ */
/* /1*     vector<bool> r_vec {}; *1/ */
/* /1*     for (auto& vec : data) { *1/ */
/* /1*         bool keep {true}; *1/ */
/* /1*         for (int i{0}; i<n; ++i) { *1/ */
/* /1*             if ( (vec[cols[i]] != 0) != v_keep_true[i]) { *1/ */
/* /1*                 keep = false; *1/ */
/* /1*                 break; *1/ */
/* /1*             } *1/ */
/* /1*         } *1/ */
/* /1*         r_vec.push_back(keep); *1/ */
/* /1*     } *1/ */
/* /1*     return r_vec; *1/ */
/* /1* }; *1/ */
/* ostream& operator<<(ostream& os, tuIntVec& tu) { */
/*     // generate the max charactures needed in each column */

/*     // max length id line */
/*     int max_id = 2; */
/*     for (auto key : tu.keys) { */
/*         int n = tu_count_digits(key); */
/*         if (n > max_id) max_id = n; */
/*     } */
/*     os << Form(Form(" %%%is",max_id),"id"); */

/*     vector<int> max_chars; */
/*     for (auto tag : tu.tags) { */
/*         max_chars.push_back(tag.size()); */
/*     } */
/*     for (auto& vec : tu.data) { */
/*         int i{0}; */
/*         for (auto val : vec) { */
/*             int n = tu_count_digits(val); */
/*             if (n > max_chars[i]) max_chars[i] = n; */
/*             ++i; */
/*         } */
/*     } */

/*     int i{0}; */
/*     for (auto tag : tu.tags) { */
/*         os << Form(Form(" %%%is",max_chars[i]),tag.c_str()); */
/*         ++i; */
/*     } */
/*     os << endl; */

/*     const char* fmt_key = Form(" %%%ii",max_id); */

/*     vector<const char*> this_fmt; */
/*     for (auto n : max_chars) { */
/*         this_fmt.push_back(Form(" %%%ii",n)); */
/*     } */

/*     int z { 0}; */
/*     const int n_tags = tu.tags.size(); */
/*     for (int i{0}; i<(int)tu.keys.size(); ++i) { */
/*         os << " " << std::right << std::setw(max_id) << tu.keys[i]; */
/*         for (int k{0}; k<n_tags; ++k) { */
/*             os << " " << std::setw(max_chars[k]) << tu.data[i][k]; */
/*         } */
/*         os << endl; */
/*     } */
/*     return os; */
/* }; */

/* void tuIntVec::write_to_file(const char* which_file, vector<string>comments) { */
/*     ofstream f_out { which_file }; */
/*     if (!f_out.is_open()) { */
/*         cout <<  "  fatal error in tuIntVec::write_to_file : " << endl */
/*             <<  "    Couldn't open file \""<<which_file<<"\""<< endl */
/*             <<  "    -> not file being written to " << endl; */
/*     } */
/*     for (auto comment : comments) f_out << "// " << comment << endl; */
/*     f_out << *this; */
/*     cout << " Write table to file: " << which_file << endl; */
/*     f_out.close(); */
/* }; */
tuIntSet::tuIntSet(vector<int>in_data) {
    list = in_data;
    sort(list.begin(), list.end());
};

tuIntSet& tuIntSet::operator+=(const tuIntSet& rhs) {
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
tuIntSet& tuIntSet::operator-=(const tuIntSet& rhs) {
    vector<int> new_list;
    for (auto val : list) {
        if (!binary_search(rhs.list.begin(),rhs.list.end(),val)) {
            new_list.push_back(val);
        }
    }
    /* for (auto v : new_vals) list.push_back(v); */
    list = new_list;
    return *this;
};
tuIntSet& tuIntSet::operator*=(const tuIntSet& sec) {
    vector<int> new_list;
    for (auto val : sec.list)
        if (binary_search(list.begin(),list.end(),val)) 
            new_list.push_back(val);
    list = new_list;
    return *this;
};

void tuIntSet::write_to_file(const char* file_name, vector<string> comments) {
    ofstream fout;
    fout.open(file_name);
    for (auto& comment : comments) fout << "// " << comment << endl;
    for (auto v : list) fout << v << endl;
    fout.close();
};


int tuIntSet::size() { return list.size(); };
void tuIntSet::clear() { list.clear(); };
ostringstream tuIntSet::read_file(const char* in_file, int col, bool print, bool strip_commas) {
    ostringstream msg;
    if (!strcmp(in_file,"")) return msg;
    if (print) {
        msg << " Reading following values from col " << col << " from " << in_file << endl;
    }
    try {
        auto new_data = tuReadIntVec(in_file, col, true, strip_commas);
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
        cerr << " fatal error in tuIntSet::read_file " << endl;
        cerr << err.what() << endl;
        /* cout << err << endl; */
    }
    return msg;
};
/* tuIntSet::tuIntSet() {}; */
ostream& operator<<(ostream& os, tuIntSet& dt) { 
    for (auto& v : dt.list) cout << v << endl;
    return os;
};
tuIntSet::tuIntSet(const char* in_file, ofstream& log, int col, bool print, bool strip_commas) {
    log << read_file(in_file, col, print, strip_commas).str() << endl;
};
tuIntSet::tuIntSet(const char* in_file, int col, bool print, bool strip_commas) {
    read_file(in_file, col, print, strip_commas);
};
tuIntSet::tuIntSet(const char* file, const char* tag) {
    for (auto val : tuReadValVec(file,tag,{{"sort",true}})) {
        list.push_back((int)val);
    }
};
bool tuIntSet::operator()(int val) { return std::binary_search(list.begin(), list.end(), val); };
bool tuIntSet::has(int i) const { return binary_search(list.begin(),list.end(),i); };
int tuIntSet::operator[](int val) {
    return (int)(std::lower_bound(list.begin(), list.end(), val) - list.begin());
};

/* tuIntBinCnt::tuIntBinCnt(const char* name, vector<int> x_dim, const char* title) : */
/*     tuIntBinCnt{name, x_dim, {}, title} {}; */
/* tuIntBinCnt::tuIntBinCnt(const char* name, vector<int> x_dim, */ 
/*         vector<int> y_dim, const char* title) { */
/*     if (y_dim.size()>0) { */
/*         is2D = true; */
/*         const char* use_title = (strcmp(title,"")) ? name : title; */
/*         hg2 = new TH2D(name, use_title, x_dim.size(), ax_doubleptr(x_dim), */
/*                 y_dim.size(), ax_doubleptr(y_dim)); */
/*     } else { */
/*         const char* use_title = (strcmp(title,"")) ? name : title; */
/*         hg1 = new TH1D(name,use_title, x_dim.size(), ax_doubleptr(x_dim)); */
/*     } */
/* }; */
/* double tuIntBinCnt::getcnt(TH1D* hg, int i) { */
/*     int i_bin = hg->FindBin((double)i); */
/*     return hg->GetBinContent(i_bin); */
/* }; */
/* double tuIntBinCnt::getcnt(TH2D* hg, int i, int j) { */
/*     int i_bin = hg->FindBin((double)i,(double)j); */
/*     return hg->GetBinContent(i_bin); */
/* }; */
/* void tuIntBinCnt::fill(int i, double weight){ */
/*     hg1->Fill( (double)i,  weight); */
/* }; */
/* void tuIntBinCnt::fill(int i, int j, double weight){ */
/*     if (!is2D) { */
/*         cout << " Error: only 1D counter generated " << endl; */
/*         return; */
/*     } */
/*     hg2->Fill( (double)i, (double)j, weight); */
/* }; */
/* void tuIntBinCnt::write() { */
/*     if (is2D) hg2->Write(); */
/*     else hg1->Write(); */
/* }; */

/* tuFnCaller::tuFnCaller(const char* file_data, double(&_fn)(double*,double*)) : */
/*     fn{_fn}, x{new double[2]} */
/* { */
/*     auto p_vals = tuReadValVec(file_data); */
/*     p = new double[p_vals.size()]; */
/*     int i{0}; */
/*     for (auto val : p_vals) p[i++] = val; */
/* }; */
/* double tuFnCaller::operator()(double x0, double x1){ */
/*     x[0] = x0; */
/*     x[1] = x1; */
/*     return fn(x,p); */
/* }; */

/* //------------------------------------------------- */
/* tuXsec::tuXsec( */
/*      const char* tag_file, */
/*      const char* Xsection_tag, */
/*      const char* pthatbin_tag, */
/*      const char* nEvents_tag */
/*  ) : */
/*     Xsection  { tuReadValVec(tag_file, Xsection_tag) }, */
/*     pthatbins { tuReadValVec(tag_file, pthatbin_tag) }, */
/*     nbins_pthat { static_cast<int>( Xsection.size()) }, */
/*     Nevents {}, */
/*     Ncollected ( nbins_pthat, 0 ) */
/* { */
/*     for (auto v : tuReadValVec(tag_file, nEvents_tag)) { */
/*         Nevents.push_back(static_cast<int>(v)); */
/*     } */
/*     tuBinVec pthat_edges { tag_file, {{"tag",pthatbin_tag}} }; */
/*     hg_collected = new TH1D("pthat_bins_collected", */
/*         "Number of events collected in each pthat bin;#hat{p}_{T};n-collected", */
/*         pthat_edges, pthat_edges); */
/* }; */

/* double tuXsec::pthatbin_center(int bin) { */
/*     return 0.5*(pthatbins[bin]+pthatbins[bin+1]); */
/* }; */

/* int tuXsec::pthatbin(pair<double,double> bounds) { */
/*     double first { bounds.first }; */
/*     auto iter = std::lower_bound(pthatbins.begin(), pthatbins.end(), bounds.first); */
/*     if (iter==pthatbins.end() || (*(iter)!=bounds.first)) { */ 
/*         throw std::runtime_error(" fatal in tuXsec::pthatbin: pthatbin not found " ); */
/*         return -1; */
/*     } */
/*     if (*(iter+1) != bounds.second) { */
/*         cout << " warning in tuXsec::pthatbin: pthatbin upperbound sought " << endl */
/*              << " doesnt match pthatbin lower found." << endl; */
/*     } */
/*     return (iter - pthatbins.begin()); */
/* }; */

/* double tuXsec::Xsec(pair<double,double> bounds, int numEvents) { */
/*     return Xsec(pthatbin(bounds), numEvents); */
/* }; */
/* double tuXsec::Xsec(int pthatbin, int numEvents) { */
/*     check_pthatbin(pthatbin); */
/*     if (numEvents != 0) return Xsection[pthatbin] / numEvents; */

/*     if (n_collected_total != 0) { */
/*         if (Ncollected[pthatbin]==0) { */
/*             cout << " fatal error: no events for bin " << pthatbin */ 
/*                 << " have been collected." << endl */
/*                 << " Cannot generate weighted cross section." << endl; */
/*             throw std::runtime_error("Asking for ptbin with no events in tuXsec::Xsec"); */
/*         } else { */
/*             return Xsection[pthatbin] / Ncollected[pthatbin] ; */
/*         } */
/*     } */
/*     return Xsection[pthatbin] / Nevents[pthatbin]; */
/* }; */

/* void tuXsec::collect (int pthatbin) { */
/*     ++n_collected_total; */
/*     ++Ncollected[pthatbin]; */
/* }; */
/* void tuXsec::collect(pair<double,double>bounds) { collect(pthatbin(bounds)); }; */

/* void tuXsec::check_pthatbin(int bin) { */
/*     if (bin >= nbins_pthat) { */
/*         throw std::runtime_error( Form( */
/*         "fatal in tuXsec : asked for pthatbin %i but there are only %i bins", */
/*         bin, nbins_pthat)); */
/*     } */
/* }; */

/* tuIntStrFunctor::tuIntStrFunctor ( const char* file, tuOptMap options, tuOptMap dict) */
/* { */
/*     dict += options; */
/*     data = tuReadIntStrMap(file, dict); */
/* }; */
/* const char* tuIntStrFunctor::operator()(int index) { */
/*     try { */
/*         return data[index].c_str(); */
/*     } */
/*     catch (...) { */
/*         cout << " Key failure in tuIntStrFunctor with index " << index << endl; */
/*         throw; */
/*     }; */
/* }; */


/* //StrStr */
/* tuStrStrFunctor::tuStrStrFunctor ( const char* file, tuOptMap options, tuOptMap dict) */
/* { */
/*     dict += options; */
/*     data = tuReadMapStrStr(file, dict); */
/* }; */
/* const char* tuStrStrFunctor::operator()(const char* key) { */
/*     try { */
/*         return data[key].c_str(); */
/*     } */
/*     catch (...) { */
/*         cout << " Key failure in tuStrStrFunctor with key " << key << endl; */
/*         throw; */
/*     }; */
/* }; */

/* tuFirst::operator bool() { */ 
/*     if (is_first) { */
/*         is_first = false; */
/*         return true; */
/*     } */
/*     return false; */
/* }; */
/* bool tuFirst::operator()() { */
/*     if (is_first) { */
/*         is_first = false; */
/*         return true; */
/*     } */
/*     return false; */
/* }; */

tuCycleTrue::tuCycleTrue(int period_in) :
    period { period_in }, cnt{1}
{};
bool tuCycleTrue::operator()() {
    ++ cnt;
    if (cnt == period) {
        cnt = 0;
        return true;
    } else {
        return false;
    }
};
tuCycleTrue::operator bool() { return this->operator()(); };
void tuCycleTrue::reset() { cnt = 0; };

tuCycleSpacer::tuCycleSpacer(int period, int _n_width, const char* _spacer, const char* _newline_spacer) :
    cycle{period}, spacer{_spacer}, n_width{_n_width}, newline_spacer{_newline_spacer}
{
    cycle.cnt=0;
};

ostream& operator<<(ostream& os, tuCycleSpacer& cs) {
    os << cs.spacer;
    if (cs.cycle) os << endl << cs.newline_spacer;
    if (cs.n_width) os << setw(cs.n_width);
    return os;
};

void tuCycleSpacer::reset() { cycle.cnt=0; };

/* tuXYbounder::tuXYbounder(vector<double> x, vector <double> y, tuOptMap opt) : */
/*     X{x}, Y{y}, */
/*     size { (int) X.size() }, */
/*     lodef{ opt.has("default-lo") ? opt["default-lo"]() : */ 
/*            size > 0 ? Y[0] : 0. */
/*     }, */
/*     hidef{ opt.has("default-hi") ? opt["default-hi"]() : */ 
/*            size > 0 ? Y[size-1] : 0. */
/*     } */
/* {}; */

/* tuXYbounder::tuXYbounder() : X {}, Y{}, size{0}, lodef{0.}, hidef{0.} */
/* {}; */

/* tuXYbounder::tuXYbounder( */
/*     const char* file, const char* tagX, */ 
/*     const char* tagY, tuOptMap opt */
/* ) : */
/*     X { tuReadValVec(file, {{"tag",tagX,"sort",true}}) }, */
/*     Y { tuReadValVec(file, {{"tag",tagY,"sort",false}}) }, */
/*     size { (int) X.size() }, */
/*     lodef{ opt.has("default-lo") ? opt["default-lo"]() : */ 
/*            size > 0 ? Y[0] : 0. */
/*     }, */
/*     hidef{ opt.has("default-hi") ? opt["default-hi"]() : */ 
/*            size > 0 ? Y[size-1] : 0. */
/*     } */
/* {}; */

/* bool tuXYbounder::operator()(double x, double y) { */
/*     if (size == 0) return false; */
/*     int bin = (int)(std::lower_bound(X.begin(), X.end(), x) - X.begin()); */
/*     /1* cout << " bin: " << bin << endl; *1/ */
/*     if (bin == 0) return y>lodef; */
/*     else if (bin == size) return y>hidef; */
/*     else return y>Y[bin]; */
/* }; */

/* double tuXYbounder::operator()(double x) { */
/*     if (size == 0) return -1; */
/*     int bin = (int)(std::lower_bound(X.begin(), X.end(), x) - X.begin()); */
/*     if (bin == 0) return lodef; */
/*     else if (bin == size) return hidef; */
/*     else return Y[bin]; */
/* }; */




tuPtrDbl::tuPtrDbl(TAxis* ax, double bin_loc, bool get_widths) {
    int n_bins = ax->GetNbins();
    if (get_widths) {
        for (int i{1}; i<=n_bins; ++i) vec.push_back(ax->GetBinWidth(i)*bin_loc);
    } else if (bin_loc == 0.5) {
        for (int i{1}; i<=n_bins; ++i) vec.push_back(ax->GetBinCenter(i));
    } else if (bin_loc == 0.) {
        for (int i{1}; i<=n_bins; ++i) vec.push_back(ax->GetBinLowEdge(i));
    } else if (bin_loc == 1.) {
        for (int i{1}; i<=n_bins; ++i) vec.push_back(ax->GetBinUpEdge(i));
    } else {
        for (int i{1}; i<=n_bins; ++i) {
            double W = ax->GetBinWidth(i);
            double L  = ax->GetBinLowEdge(i);
            vec.push_back(L+W*bin_loc);
        }
    }
    build_ptr();
}; 

tuPtrDbl::tuPtrDbl(vector<double> V){
    for (auto v : V) vec.push_back(v);
    build_ptr();
};
tuPtrDbl::tuPtrDbl(TH1* hg, bool get_errors) {
    if (get_errors) {
        for (auto i{1}; i<= hg->GetXaxis()->GetNbins(); ++i) {
            vec.push_back(hg->GetBinError(i));
        }
    } else {
        for (auto i{1}; i<= hg->GetXaxis()->GetNbins(); ++i) {
            vec.push_back(hg->GetBinContent(i));
        }
    }
    build_ptr();
};

tuPtrDbl::tuPtrDbl(const tuPtrDbl& ihs) {
    for (auto v : ihs.vec) vec.push_back(v);
    build_ptr();
};

tuPtrDbl& tuPtrDbl::operator=(const tuPtrDbl& rhs) {
    vec.clear();
    if (ptr) delete ptr;
    for (auto v : rhs.vec) vec.push_back(v);
    build_ptr();
    return *this;
};

/* tuPtrDbl::tuPtrDbl(tuPtrDbl ihs) { */
/*     for (auto v : ihs.vec) vec.push_back(v); */
/*     build_ptr(); */
/* }; */

void tuPtrDbl::build_ptr() {
    size = vec.size();
    ptr = new double[size];
    for (int i{0}; i<size; ++i) ptr[i] = vec[i];
};

tuPtrDbl& tuPtrDbl::update() {
    for (int i{0}; i<size; ++i) ptr[i] = vec[i];
    return *this;
};

tuPtrDbl::operator int () { return size; };
tuPtrDbl::operator double* () { return ptr; };
tuPtrDbl::operator vector<double> () { return vec; };

tuPtrDbl::tuPtrDbl(const char* file, const char* tag) {
    auto read_vec = tuReadValVec(file,tag);
    for (auto v: read_vec) vec.push_back(v);
    build_ptr();
};
tuPtrDbl::tuPtrDbl(int n) {
    for (auto i{0}; i<n; ++i) vec.push_back(0);
    build_ptr();
};

tuPtrDbl::~tuPtrDbl() {
    delete[] ptr;
};
vector<double>::iterator tuPtrDbl::begin() { return vec.begin(); };
vector<double>::iterator tuPtrDbl::end()   { return vec.end(); };

double& tuPtrDbl::operator[](int i) { return vec[i]; };
ostream& operator<<(ostream& os, tuPtrDbl& tu) {
    for (auto v : tu) cout << " " << v;
    return os;
};
int tuPtrDbl::throw_error(const char* msg) {
        throw std::runtime_error(Form(" fatal in tuPtrDbl::%s, sizes don't match",msg));
        return -1;
};

tuPtrDbl& tuPtrDbl::operator+=(const tuPtrDbl& rhs) {
    if (size != rhs.size) throw_error("operator+=");
    for (auto i{0}; i<size; ++i) vec[i] += rhs.vec[i];
    update();
    return *this;
};
tuPtrDbl& tuPtrDbl::operator-=(const tuPtrDbl& rhs) {
    if (size != rhs.size) throw_error("operator-=");
    for (auto i{0}; i<size; ++i) vec[i] -= rhs.vec[i];
    update();
    return *this;
};
tuPtrDbl& tuPtrDbl::operator*=(const tuPtrDbl& rhs) {
    if (size != rhs.size) throw_error("operator-=");
    for (auto i{0}; i<size; ++i) vec[i] *= rhs.vec[i];
    return update();
};
tuPtrDbl& tuPtrDbl::abs() {
    for (auto& v : vec) if (v<0.) v = -v;
    return update();
};
tuPtrDbl& tuPtrDbl::sqrt() {
    for (auto& v : vec) v = TMath::Sqrt(v);
    return update();
};
tuPtrDbl& tuPtrDbl::zero_negatives() {
    for (auto& v : vec) if (v<0.) v = 0;
    return update();
};
tuPtrDbl& tuPtrDbl::abs_diff(const tuPtrDbl& rhs) {
    if (size != rhs.size) throw_error("abs_diff");
    *this -= rhs;
    this->abs();
    return update();
};
tuPtrDbl& tuPtrDbl::square_diff(const tuPtrDbl& rhs) {
    if (size != rhs.size) throw_error("abs_diff");
    *this -= rhs;
    *this *= *this;
    return update();
};
tuPtrDbl& tuPtrDbl::make_min(const tuPtrDbl& rhs) {
    if (size != rhs.size) throw_error("make_min");
    for (auto i{0}; i<size; ++i) {
        if (rhs.vec[i] < vec[i]) vec[i] = rhs.vec[i];
    }
    return update();
};
tuPtrDbl& tuPtrDbl::make_max(const tuPtrDbl& rhs) {
    if (size != rhs.size) throw_error("make_max");
    for (auto i{0}; i<size; ++i) {
        if (rhs.vec[i] > vec[i]) vec[i] = rhs.vec[i];
    }
    return update();
};
tuPtrDbl& tuPtrDbl::operator/=(const tuPtrDbl& rhs) {
    /* cout << " z3 " << size << " " << rhs.size << endl; */
    if (size != rhs.size) throw_error("operator/=");
    for (auto i{0}; i<size; ++i) vec[i] /= rhs.vec[i];
    return update();
};
tuPtrDbl& tuPtrDbl::operator/=(const double rhs) {
    for (auto& v : vec) v /= rhs;
    return update();
};
tuPtrDbl& tuPtrDbl::operator*=(const double rhs) {
    for (auto& v : vec) v *= rhs;
    return update();
};
tuPtrDbl& tuPtrDbl::operator+=(const double rhs) {
    for (auto& v : vec) v += rhs;
    return update();
};
tuPtrDbl& tuPtrDbl::operator-=(const double rhs) {
    for (auto& v : vec) v -= rhs;
    return update();
};

tuPtrDbl  operator+(const tuPtrDbl& lhs, const tuPtrDbl& rhs) {
    if (lhs.size != rhs.size) {
        throw std::runtime_error(" fatal in tuPtrDbl+tuPtrDbl, sizes don't match");
    }
    tuPtrDbl r_val{lhs};
    r_val += rhs;
    return r_val.update();
};
tuPtrDbl  operator-(const tuPtrDbl& lhs, const tuPtrDbl& rhs) {
    if (lhs.size != rhs.size) {
        throw std::runtime_error(" fatal in tuPtrDbl-tuPtrDbl, sizes don't match");
    }
    tuPtrDbl r_val{lhs};
    r_val-=rhs;
    return r_val.update();
};
tuPtrDbl tu_calc_quadrature(vector<tuPtrDbl> data) {
        if (data.size() == 0) return tuPtrDbl{vector<double>{}};
        if (data.size() == 1) return data[0];
        data[0] *= data[0];
        for (unsigned int i=1; i<data.size(); ++i) {
            data[0] += (data[i] *= data[i]);
        }
        return data[0].sqrt();
};
tuPtrDbl tu_calc_quadrature(vector<tuPtrDbl> data, tuPtrDbl mean) {
        if (data.size() == 0) return mean;
        for (auto& dat : data) dat -= mean;
        return tu_calc_quadrature(data);
};
tuPtrDbl tu_calc_mean(vector<tuPtrDbl> data) {
        if (data.size() == 0) return tuPtrDbl{vector<double>{}};
        if (data.size() == 1) return data[0];
        for (unsigned int i=1; i<data.size(); ++i) {
            data[0] += data[i];
        }
        data[0] /= data.size();
        return data[0];
};
tuPtrDbl tu_calc_max_bound(vector<tuPtrDbl> data) {
        if (data.size() == 0) return tuPtrDbl{vector<double>{}};
        if (data.size() == 1) return data[0];
        for (unsigned int i=1; i<data.size(); ++i) {
            data[0].make_max(data[i]);
        }
        return data[0];
};
tuPtrDbl tu_calc_max_berr(vector<tuPtrDbl> data, tuPtrDbl mean) {
        auto max_bound = tu_calc_max_bound(data);
        return max_bound.abs_diff(mean);
};
tuPtrDbl tu_calc_min_bound(vector<tuPtrDbl> data) {
        if (data.size() == 0) return tuPtrDbl{vector<double>{}};
        if (data.size() == 1) return data[0];
        for (unsigned int i=1; i<data.size(); ++i) {
            data[0].make_min(data[i]);
        }
        return data[0];
};
tuPtrDbl tu_calc_min_berr(vector<tuPtrDbl> data, tuPtrDbl mean) {
        return mean.abs_diff(tu_calc_min_bound(data));
};
pair<tuPtrDbl,tuPtrDbl> tu_calc_bounds(vector<tuPtrDbl> data) {
    return { tu_calc_min_bound(data), tu_calc_max_bound(data) };
};

TGraphAsymmErrors* tu_draw_error_boxes(TH1D* mean, tuPtrDbl err, tuOptMap opts, array<double,4>x_set,
        double range_lo, double range_hi) {
    return tu_draw_error_boxes(mean, {err,err}, opts, x_set, range_lo, range_hi);
};
TGraphAsymmErrors* tu_draw_error_boxes(TH1D* mean, array<tuPtrDbl,2> err, 
        tuOptMap dict, array<double,4>x_set,
        double range_lo, double range_hi) 
{
    tuSysErrors pts { mean, x_set, err };
    /* cout << " range_lo: " << range_lo << " " << range_hi << endl; */
    if (range_lo!=0. || range_hi!=0) { pts.set_x_range(range_lo,range_hi); };
    auto tgase = pts.tgase;
    tu_fmt(tgase,dict);
    /* tgase->SetFillColor(kBlue); */
    if ( dict["FillColor"] ) {
        tgase->Draw("E2");
    };
    if ( dict["LineColor"] ) tuDrawBoxErrors(tgase, dict);
    return tgase;
};

tuSysErrors::tuSysErrors(TGraphAsymmErrors* _tgase, array<double,4> x_rat) : tgase{_tgase}
{
    size = tgase->GetN(); 
    if (x_rat[0]>=0) set_rat_xbins(x_rat);
};

tuSysErrors::tuSysErrors(TH1* hg, array<double,4> x_rat, array<tuPtrDbl,2> _err) {
    tuPtrDbl x     {hg->GetXaxis()};
    tuPtrDbl err_x {hg->GetXaxis(), 0.5, true};

    tuPtrDbl y     {hg};
    tuPtrDbl err_y {hg,true};

    tuPtrDbl err_y_lo = _err[0].size==0 ? err_y : _err[0];
    tuPtrDbl err_y_hi = _err[1].size==0 ? err_y : _err[1];
    
    tgase = new TGraphAsymmErrors (x.size,x,y,err_x,err_x,err_y_lo,err_y_hi);
    tgase->SetMarkerColor(hg->GetMarkerColor());
    tgase->SetMarkerSize(hg->GetMarkerSize());
    tgase->SetMarkerStyle(hg->GetMarkerStyle());
    tgase->SetLineColor(hg->GetLineColor());
    tgase->SetLineStyle(hg->GetLineStyle());
    size = x.size;

    if (x_rat[0]>=0) set_rat_xbins(x_rat);
};

/* tuSysErrors& tuSysErrors::add_data(tuPtrDbl data) { vec_data.push_back(data); return *this;}; */
/* tuSysErrors& tuSysErrors::add_data(vector<tuPtrDbl> data) { */ 
    /* for (auto& dat : data) vec_data.push_back(dat); */
    /* return *this; */
/* }; */
tuSysErrors& tuSysErrors::set_x_range (double x_lo, double x_hi) {
    int i_left  = 0;
    int i_right = size-1;
    double* x_old = tgase->GetX();
    double* y_old = tgase->GetY();
    for (int i{0}; i<size; ++i) {
        if (x_old[i] < x_lo) i_left = x_old[i];
        else break;
    }
    for (int i{size-1}; i>=0; --i) {
        if (x_old[i_right] > x_hi) i_right = i;
        else break;
    }
    if (i_left != 0 || i_right != (size-1)) {
        int new_size = i_right-i_left+1;
        double *x = new double[new_size];
        double *y = new double[new_size];
        std::copy(&(x_old[i_left]), &(x_old[i_right+1]), x);
        std::copy(&(y_old[i_left]), &(y_old[i_right+1]), y);
        auto tgase_new = new TGraphAsymmErrors(new_size,x,y);
        size = new_size;
        for (int i{0}; i<new_size; ++i) {
            tgase_new->SetPointError( i,
                 tgase->GetErrorXlow(i_left+i), tgase->GetErrorXhigh(i_left+i),
                 tgase->GetErrorYlow(i_left+i), tgase->GetErrorYhigh(i_left+i));
        }
        delete tgase;
        tgase = tgase_new;
    }
    return *this;
};

tuSysErrors& tuSysErrors::swap_xy () {
    double* x = new double[size];
    double* y = new double[size];

    double* x_old = tgase->GetX();
    double* y_old = tgase->GetY();

    std::copy(x_old, &(x_old[size]), y);
    std::copy(y_old, &(y_old[size]), x);
    auto tgase_new = new TGraphAsymmErrors(size,x,y);

    for (int i{0}; i<size; ++i) {
        tgase_new->SetPointError( i,
                tgase->GetErrorYlow(i), tgase->GetErrorYhigh(i),
                tgase->GetErrorXlow(i), tgase->GetErrorXhigh(i));
    }
    delete tgase;
    tgase = tgase_new;
    return *this;
};
tuSysErrors& tuSysErrors::setYhigh(tuPtrDbl& data) {
    for (int i{0}; i<data.size; ++i) {
        tgase->SetPointEYhigh(i,data[i]);
    }
    return *this;
};
tuSysErrors& tuSysErrors::setYlow(tuPtrDbl& data) {
    for (int i{0}; i<data.size; ++i) {
        tgase->SetPointEYlow(i,data[i]);
    }
    return *this;
};
tuSysErrors& tuSysErrors::setYhilo(tuPtrDbl& data) {
    for (int i{0}; i<data.size; ++i) {
        tgase->SetPointEYlow(i,data[i]);
        tgase->SetPointEYhigh(i,data[i]);
    }
    return *this;
};

tuSysErrors::tuSysErrors(const tuSysErrors& rhs) {
    tgase = (TGraphAsymmErrors*) rhs.tgase->Clone();
    size = tgase->GetN();
};

tuPtrDbl tuSysErrors::getY() {
    tuPtrDbl y {size};
    double* y_dat = tgase->GetY();
    for (int i{0}; i<size; ++i) y[i] = y_dat[i];
    y.update();
    return y;
};
tuPtrDbl tuSysErrors::getYlow() {
    tuPtrDbl y {size};
    for (int i{0}; i<size; ++i) y[i] = tgase->GetErrorYlow(i);
    y.update();
    return y;
};
tuPtrDbl tuSysErrors::getYhigh() {
    tuPtrDbl y {size};
    for (int i{0}; i<size; ++i) y[i] = tgase->GetErrorYhigh(i);
    y.update();
    return y;
};

tuSysErrors& tuSysErrors::set_rat_xbins(array<double,4> rat_rel) {
    double r_left   { rat_rel[0] };
    double r_right  { rat_rel[1] };
    double r_center { rat_rel[2] }; // relative position between left and right
    double offset_c { rat_rel[3] };
    
    double* x = tgase->GetX();
    for (int i{0}; i<size; ++i) {
        double deltaX = tgase->GetErrorXlow(i) + tgase->GetErrorXhigh(i);
        double anchor  = x[i]-tgase->GetErrorXlow(i);
        double p_left  = anchor + r_left  * deltaX;
        double p_right = anchor + r_right * deltaX;
        double p_center = anchor + r_center * deltaX + offset_c;
        tgase->SetPointEXlow (i,  p_center-p_left);
        tgase->SetPointEXhigh(i,  p_right-p_center);
        tgase->SetPoint(i,p_center,tgase->GetY()[i]);
    }
    return *this;
};

TGraphAsymmErrors* tuSysErrors::operator-> () { return tgase; };
tuSysErrors::operator TGraphAsymmErrors* () { return tgase; };

/* tu_TF1fitter::tu_TF1fitter(const char* fnc_str, const char* _name, TH1D* hg, double _lo, double _hi) : */
/*     name { strcmp(_name,"")==0 ? tuUniqueName() : _name } , */
/*     fn_str { fnc_str } */
/* { */
/*     cout << " this name: " << name << endl; */
/*     if (_lo || _hi) { */
/*         lo = _lo; */
/*         hi = _hi; */
/*     } */
    
/*     if (lo || hi) fn = new TF1(name.c_str(), fnc_str,lo,hi); */
/*     else          fn = new TF1(name.c_str(), fnc_str); */
/*     n_par = fn->GetNpar(); */
/*     if (hg)  fit(hg); */
/* }; */
/* void tu_TF1fitter::set_presets(vector<double> pre_set) { */
/*     for (unsigned int i{0}; i<pre_set.size(); ++i) fn->SetParameter(i,pre_set[i]); */
/* }; */
/* vector<double>& tu_TF1fitter::fit(TH1D* hg, double _lo, double _hi) { */
/*     if (_lo || _hi) { */
/*         lo = _lo; */
/*         hi = _hi; */
/*     } */

/*     if (lo || hi) { cout << " LO-HI " << endl; hg->Fit(fn,"","",lo,hi); }//hg->GetXaxis()->SetRangeUser(lo,hi); */
/*     else          hg->Fit(fn); */

/*     /1* gccint n_par = fn->GetNpar(); *1/ */
/*     for (unsigned int i{0}; i<n_par; ++i) { */
/*         if (i >= fval.size()) fval.push_back( fn->GetParameter(i) ); */
/*         else                  fval[i] =       fn->GetParameter(i)  ; */

/*         if (i >= ferr.size()) ferr.push_back( fn->GetParError(i) ); */
/*         else                  ferr[i] =       fn->GetParError(i)  ; */
/*     }; */
/*     return fval; */
/* }; */
/* vector<double>& tu_TF1fitter::operator()(TH1D* hg, double lo, double hi) { */
/*     return fit                          (hg, lo, hi); */
/* }; */
/* vector<double>& tu_TF1fitter::operator()() { return fval; } ; */
/* pair<double,double> tu_TF1fitter::operator()(int i){ */
/*     if (i>=(int)n_par) { */
/*         throw std::runtime_error( */
/*                 Form("Error: requesting paremeter %i of tu_TF1fitter (%s) which has only %i pars", */
/*                     i, name.c_str(), n_par) ); */
/*     } */
/*     return {fval[i], ferr[i]}; */
/* }; */
/* double tu_TF1fitter::operator[](int i){ */
/*     if (i>=(int)n_par) { */
/*         throw std::runtime_error( */
/*                 Form("Error: requesting paremeter %i of tu_TF1fitter (%s) which has only %i pars", */
/*                     i, name.c_str(), n_par) ); */
/*     } */
/*     return fval[i]; */
/* }; */
/* void tu_TF1fitter::fix_match_params(tu_TF1fitter& to_match, std::set<int> which) { */
/*     for (int i=0; i<(int) n_par; ++i) { */
/*         if (which.count(i) != 0) { */
/*             fn->FixParameter(i, to_match[i]); */
/*         } else { */
/*             fn->SetParLimits(i,0,0); */
/*         } */
/*     }; */
/* }; */

/* ostream& operator<<(ostream& os, tu_TF1fitter& fit) { */
/*     os << " Fit for TF1 \"" << fit.name << "\" fn: " << fit.fn_str << endl; */    
/*     double lo, hi; */
/*     for (unsigned int i{0}; i<fit.n_par; ++i) { */
/*         fit.fn->GetParLimits(i,lo,hi); */
/*         /1* if (!lo && !hi) os << " [" << i << "] " << fit.fval[i] << "  err:  " << fit.ferr[i] << endl; *1/ */
/*         /1* else            os << "*[" << i << "] " << fit.fval[i] << "  err:  " << fit.ferr[i] << endl; *1/ */
/*         const char* is_fix = (!lo && !hi) ? " " : "*" ; */
/*         os << Form( "%s [%2i]  %10.5g  err: %10.5g  (err/val: %10.5g)", is_fix, i, fit.fval[i], fit.ferr[i], fit.ferr[i]/fit.fval[i]) << endl; */
/*         /1* else            os << "*[" << i << "] " << fit.fval[i] << "  err:  " << fit.ferr[i] << endl; *1/ */
/*     } */
/*     return os; */
/* }; */

/* tuParticleThrower::tuParticleThrower(TH1D* hg, double multiple, unsigned int _seed) : */
/*     r3{_seed} */ 
/* { */
/*     TAxis* ax { hg->GetXaxis() }; */
/*     for (int i{1}; i<=hg->GetXaxis()->GetNbins(); ++i) { */
/*         auto val { hg->GetBinContent(i) }; */
/*         if (val<=0) continue; */
/*         val *= multiple; */

/*         l_bound.push_back(ax->GetBinLowEdge(i)); */
/*         u_bound.push_back(ax->GetBinUpEdge(i)); */
/*         whole.push_back(static_cast<int>(TMath::Floor(val))); */
/*         remainder.push_back(val-TMath::Floor(val)); */
/*     } */
/*     i_size = l_bound.size(); */
/* }; */
/* ostream& operator<<(ostream& os, tuParticleThrower& tu) { */
/*     for (int i{0}; i<static_cast<int>(tu.whole.size()); ++i) { */
/*         cout << Form("bin[%2i]  [%5.1f-%5.1f] : %i   ", */
/*                 i,tu.l_bound[i],tu.u_bound[i],tu.whole[i]) */ 
/*             << tu.remainder[i] << endl; */
/*     } */
/*     return os; */
/* }; */
/* bool tuParticleThrower::throw_particle() { */
/*     if (i_bin == i_size) { */
/*         i_bin = 0; */
/*         i_whole = 0; */
/*         return false; */
/*     } */
/*     unsigned int n_whole = whole[i_bin]; */
/*     if (i_whole < n_whole) { */
/*         set_rand(i_bin); */
/*         ++i_whole; */
/*         is_thrown = true; */
/*     } else { */
/*         i_whole = 0; */
/*         if (r3.Uniform()<remainder[i_bin]) { */
/*             set_rand(i_bin); */
/*             is_thrown=true; */
/*         } else { */
/*             is_thrown=false; */
/*         } */
/*         ++i_bin; */
/*     } */
/*     return true; */
/* }; */
/* void tuParticleThrower::set_rand(unsigned int i) { */
/*     phi = r3.Uniform(tu_twopi); */
/*     eta = r3.Uniform(-1.,1.); */
/*     pt =  r3.Uniform(l_bound[i],u_bound[i]); */
/* }; */

/* tuPoissonParticleThrower::tuPoissonParticleThrower( */
/*         TRandom3& _r3, TH1D& _dist, double mult) */ 
/*     : r3{_r3}, dist{_dist} */
/* { */
/*     mean = mult * dist.Integral(); */
/*     i_to_throw = r3.Poisson(mean); */
/*     // remove any negative bins */
/*     for (int i=0; i<dist.GetXaxis()->GetNbins(); ++i) { */
/*         if (dist.GetBinContent(i)<0) dist.SetBinContent(i,0.); */
/*     } */
/* }; */

/* bool tuPoissonParticleThrower::throw_particle() { */
/*     if (i_to_throw == i_thrown) { */
/*         i_thrown = 0; */
/*         i_to_throw = r3.Poisson(mean); */
/*         return false; */
/*     } else { */
/*         phi = r3.Uniform(tu_twopi); */
/*         eta = r3.Uniform(-1.,1.); */
/*         pt =  dist.GetRandom(); */
/*         ++i_thrown; */
/*         return true; */
/*     } */
/* }; */

/* tuHopper1D::tuHopper1D(vector<double> _edges) : */
/*     edges { _edges } */
/* { */
/*     nbins = static_cast<int>( edges.size()-1 ); */ 
/*     for (int i{0}; i<nbins; ++i) { */
/*         double center { 0.5*(edges[i]+edges[i+1]) }; */
/*         hopper.push_back({center,0.}); */
/*     } */
/* }; */
/* int tuHopper1D::fill(double val, double weight) { */
/*     auto lb = upper_bound(edges.begin(), edges.end(), val); */
/*     if (lb == edges.begin()) return -1; // don't fill anything */
/*     if (lb == edges.end())   return -1; // don't fill anything */
/*     int bin = static_cast<int>(lb-edges.begin()-1); */
/*     hopper[bin].second += weight; */
/*     return bin; */
/* }; */
/* void tuHopper1D::reset() { */
/*     for (auto& v : hopper) v.second = 0.; */
/* }; */
/* vector<pair<double,double>>::iterator tuHopper1D::begin() { */
/*     return hopper.begin(); */
/* } */
/* vector<pair<double,double>>::iterator tuHopper1D::end() { */
/*     return hopper.end(); */
/* } */
