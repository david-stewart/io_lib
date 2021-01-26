#include "ioClass.h"
#include "io_fmt.h"
#include "io_fnc.h"

#include <sstream>
#include <algorithm>
#include "TString.h" // not used, but present for some header 
#include "TMath.h"

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

ioPads::ioPads ( vector<pair<ioPadDim, ioPadDim>> _pad_dimensions, int
        _canvas_width, int _canvas_height) :
    pad_dimensions{ _pad_dimensions }
{
    if (_canvas_width)  canvas_width  = _canvas_width;
    if (_canvas_height) canvas_height = _canvas_height;
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
            pad_dimensions.push_back({x,y});

    nCol = x_dim.size();
    if (c_wide) canvas_width = c_wide;
    if (c_height) canvas_height = c_height;
};
ioPads::ioPads(int nPads, int c_wide, int c_high){
    if (nPads==2) {
        pad_dimensions.push_back( {{0.,0.12,0.9,0.99},{0.55,0.55,0.95,0.99}} );
        pad_dimensions.push_back( {{0.,0.12,0.9,0.99},{0.00,0.15,0.55     }} );
        canvas_width = 800;
        canvas_height = 800;
    } else if (nPads ==1) {
        pad_dimensions.push_back( {{0.,0.12,0.9,0.99},{0.00,0.15,0.95,0.99}} );
        canvas_width = 800;
        canvas_height = 800;
    } else {
        throw std::runtime_error(" fatal: Called ioPads(int nPads...) with nPads > 2");
    }
};
TPad* ioPads::operator()(int row, int col) {
    if (pads.size() == 0) init();
    int i_pad = row+col*nCol;
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
    int i{0};
    canvas = new TCanvas(ioUniqueName(),"",canvas_width, canvas_height);
    canvas->Draw();

    cout << " name " << canvas->GetName() << endl;
    io_fmt(this->canvas);

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
        while (words >> word) {
            if (word.BeginsWith("//") || word.BeginsWith("#")) break;
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
bool ioIntList::has(int val) { return this->operator()(val); };
bool ioIntList::has_not(int val) { return !(this->operator()(val)); };

// run list id

ioRunListId::ioRunListId(const char* file_name, bool skip_127053_138064) {
    // skip_after_mid_comment for two 
    ifstream file;
    file.open(file_name);
    if (!file.is_open()){
        cout << "fatal error, couldn't open input file \""<< file_name <<"\"" << endl;
        exit(2);
    }
    string line;
    while (getline(file,line)){
        if (line.substr(0,1)=="#") continue;
        int id, runid, sec;
        stringstream words(line);
        words >> id >> runid >> sec;
        if (skip_127053_138064) {
            if (runid == 16127053 || runid == 16138064) continue;
        }
        map_id [runid] = static_cast<double>(id);
        run_sec[runid] = sec;
    }
};

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


// ioHgStats
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
