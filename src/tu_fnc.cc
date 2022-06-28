#include "tu_fnc.h"
#include <numeric>
#include <iostream>
#include <fstream>
#include "tuOptMap.h"
/* #include "tuClass.h" */
/* #include "RooUnfoldBayes.h" */
/* #include "RooUnfoldResponse.h" */
#include "TMath.h"
#include "TProfile2D.h"
#include <numeric>
#include <iostream>
#include <fstream>
#include "tuMiniClass.h"

int tuOpenShape (int i, int cycle) {
    array<int, 10> shapes { kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenStar, kOpenDiamond, kOpenCross, kOpenDoubleDiamond };
    return shapes[i%cycle];
};
int tuFullShape (int i, int cycle) {
    array<int, 10> shapes { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullDiamond, kFullCross, kFullDoubleDiamond };
    return shapes[i%cycle];
};

//*
TH1D*       tuHgBlank     (TH1* _, tuOptMap opt, tuOptMap dict) {
    opt("xAxisTitle") = _->GetXaxis()->GetTitle();
    opt("yAxisTitle") = _->GetYaxis()->GetTitle();
    opt += dict;
    TAxis* ax = _->GetXaxis();
    int nbins = ax->GetNbins();
    double *edges = new double[nbins+1];
    for (int i=1;i<=nbins;++i) edges[i-1] = ax->GetBinLowEdge(i);
    edges[nbins] = ax->GetBinUpEdge(nbins);
    return tuHgBlank(nbins, edges, opt);
};

TH1D*       tuHgBlank     (int nbins, double* edges, tuOptMap opt, tuOptMap dict) {
    TH1D* hg = new TH1D(tuUniqueName(), ";;", nbins, edges);
    opt += dict;
    tu_fmt(hg,opt);
    cout << "dict: " << opt << endl;
    return hg;
};

void tu_fmt_ax (TH1* tu, TH1* from) {
    tu->GetXaxis()->SetTitle(from->GetXaxis()->GetTitle());
    tu->GetYaxis()->SetTitle(from->GetYaxis()->GetTitle());
};

void tuPause(int i)  {
    TCanvas* c = new TCanvas( tuUniqueName(i+100), tuUniqueName(i+100), 100,100);
    c->SetFillColor(kCyan);
    c->SetFillStyle(1001);
    c->Draw();
    c->WaitPrimitive();
};
//*
const char* tuUniqueName(int i) {
    while (gDirectory->FindObjectAny(Form("unique_name__%i",i))!=nullptr) ++i;
    return Form("unique_name__%i",i);
};
//*
vector<int> tuColorVec(int n_colors, int palette, bool print) {
    TCanvas *junk_canvas = new TCanvas(tuUniqueName(),"junk_canvas",1200,800);
    gStyle->SetPalette(palette);
    vector<TProfile*> hg;
    vector<int> colors;
    for (int i{0};i<n_colors;++i){
        if (print) cout << " i: " << i << endl;
        hg.push_back(new TProfile(tuUniqueName(i),";;",n_colors,0.,n_colors) );
        hg[i]->SetStats(0);
        hg[i]->SetMarkerStyle(kFullCircle);
        hg[i]->SetMarkerSize(10.);
        hg[i]->Fill(i,i);
        if (i==0) { 
            hg[i]->GetYaxis()->SetRangeUser(0.,n_colors+2);
            hg[i]->Draw("PLC PMC");
        } else hg[i]->Draw("same PE PLC PMC");
    }
    if (true) junk_canvas->SaveAs(Form("colors_temp.pdf"));
    for (int i{0};i<n_colors; ++i) {
        colors.push_back(hg[i]->GetMarkerColor()) ;
        if (print) cout << " name " << hg[i]->GetName() << "  color " << hg[i]->GetMarkerColor() <<
            " " << hg[i]->GetLineColor() << endl;
    }
    /* if (print) for (auto i : colors)  cout << " color: " << i << endl; */
    cout << " // colors generated: " << endl;
    cout << " vector<int> colors { " << colors[0];
    for (int i{1}; i<n_colors;++i) cout << ", " << colors[i];
    cout << " };" << endl;
    delete junk_canvas;
    return colors;
};

array<double,3> tuPercentAround(TH1* hg, double val) {
    TAxis* ax = hg->GetXaxis();
    int i = ax->FindBin(val);
    double lo = ax->GetBinLowEdge(i);
    double hi = ax->GetBinUpEdge(i);
    double total = hg->Integral();
    array<double,3> r_vals { hg->Integral(1,i-1)/total, hg->GetBinContent(i)/total, hg->Integral(i+1,ax->GetNbins())/total };
    cout << Form(" val: %f in histogram %s:  in bin (%i); location in bin: (%f:%f:%f), perc: (%f:%f:%f) or w/bin (%f,%f)",
            val, hg->GetName(), i, lo,val,hi, r_vals[0], r_vals[1], r_vals[2], r_vals[0]+r_vals[1],r_vals[1]+r_vals[2]) << endl;
    return r_vals;
};
//*
TH1* tuDivide(TH1* num, TH1* den, tuOptMap opt, tuOptMap dict) {
    dict += opt;

    double norm_num {1};
    double norm_den {1};

    if (dict["norm"]) {
        norm_num = num->Integral();
        norm_den = den->Integral();
    };
    if (dict["print"]) cout << " norms: " << norm_num << " and " << norm_den << endl;

    TH1D* ret;
    if (dict["style-den"]) {
        ret = (TH1D*) den->Clone(tuUniqueName());
    } else {
        ret = (TH1D*) num->Clone(tuUniqueName());
    }

    for (int j=1;j<den->GetXaxis()->GetNbins()+1;++j){
        double n = num->GetBinContent(j) / norm_num;
        double d = den->GetBinContent(j) / norm_den;
        if (n == 0 || d == 0) {
            ret->SetBinContent(j,0);
            ret->SetBinError(j,0);
            /* if (i==0) cout << "Jin: " << j << " val: " << den->GetBinContent(j) << endl; */
            continue;
        }
        double n_err = num->GetBinError(j) / norm_num;
        double d_err = den->GetBinError(j) / norm_den;
        double val = n / d;
        double err = val * pow( pow(n_err/n,2)+pow(d_err/d,2),0.5);
        ret->SetBinContent(j,val);
        ret->SetBinError(j,err);
    }
    return ret;
};
TH1* tuMultiply(TH1* num, TH1* den, tuOptMap opt, tuOptMap dict) {
    dict += opt;

    double norm_num {1};
    double norm_den {1};

    if (dict["norm"]) {
        norm_num = num->Integral();
        norm_den = den->Integral();
    };
    cout << " norms: " << norm_num << " and " << norm_den << endl;

    TH1D* ret;
    if (dict["style-den"]) {
        ret = (TH1D*) den->Clone(tuUniqueName());
    } else {
        ret = (TH1D*) num->Clone(tuUniqueName());
    }

    for (int j=1;j<den->GetXaxis()->GetNbins()+1;++j){
        double n = num->GetBinContent(j) / norm_num;
        double d = den->GetBinContent(j) / norm_den;
        if (n == 0 || d == 0) {
            ret->SetBinContent(j,0);
            ret->SetBinError(j,0);
            /* if (i==0) cout << "Jin: " << j << " val: " << den->GetBinContent(j) << endl; */
            continue;
        }
        double n_err = num->GetBinError(j) / norm_num;
        double d_err = den->GetBinError(j) / norm_den;
        double val = n * d;
        double err = val * pow( pow(n_err/n,2)+pow(d_err/d,2),0.5);
        ret->SetBinContent(j,val);
        ret->SetBinError(j,err);
    }
    return ret;
};
//*
void tuDrawTLatex(const char* msg, double x, double y, tuOptMap options, tuOptMap dict
) {
    dict += options;
    TLatex tlatex;
    tlatex.SetTextColorAlpha(dict("TextColor"),dict("TextColorAlpha"));
    tlatex.SetTextAlign(dict("TextAlign"));
    tlatex.SetTextSize (dict("TextSize"));
    tlatex.SetTextFont (dict("TextFont"));
    tlatex.SetTextAngle (dict("TextAngle"));
    tlatex.DrawLatex(x, y, msg);
};
//*
void tuDrawTLine(double x0, double y0, double x1, double y1, 
        tuOptMap options, tuOptMap dict)
{
    dict += options;
    // ndc is needed when getting vertical lines in logaritymic plots
    if (dict["ndc"]) {
        TLine* line = new TLine();
        tu_fmt(line, options);
        line->DrawLineNDC(x0,y0,x1,y1);
    } else {
        TLine* line = new TLine(x0,y0,x1,y1);
        tu_fmt(line, options);
        line->Draw();
    }
};
//*
void tuDrawTLineBox(double x0, double y0, double x1, double y1, tuOptMap options) {
    tuDrawTLine( x0, y0, x1, y0, options);
    tuDrawTLine( x1, y0, x1, y1, options);
    tuDrawTLine( x1, y1, x0, y1, options);
    tuDrawTLine( x0, y1, x0, y0, options);
};
//*
void tuDrawBoxErrors(TGraphAsymmErrors* tgas, tuOptMap options, tuOptMap dict)
{
    dict += options;
    double* x = tgas->GetX();
    double* y = tgas->GetY();
    for (int i=0; i<tgas->GetN(); ++i) {
        double x0 = x[i] - tgas->GetErrorXlow(i);
        double x1 = x[i] + tgas->GetErrorXhigh(i);
        double y0 = y[i] - tgas->GetErrorYlow(i);
        double y1 = y[i] + tgas->GetErrorYhigh(i);
        tuDrawTLineBox(x0,y0,x1,y1,options);
    }
};
//*
double tuPadxRat(double y_in){
    gPad->Update();
    double y0 = gPad->GetUxmin();
    double y1 = gPad->GetUxmax();
    return y0+y_in*(y1-y0);
};
//*
double tuPadyRat(double y_in){
    gPad->Update();
    double y0 = gPad->GetUymin();
    double y1 = gPad->GetUymax();
    return y0+y_in*(y1-y0);
};
//*
void tuDrawTLineHorizontal(double y, tuOptMap options) {
    gPad->Update();
    double x0 { tuPadxRat(0.) };
    double x1 { tuPadxRat(1.) };
    tuDrawTLine(x0,y,x1,y,options);
};
//*
void tuDrawTLineVertical(double x, tuOptMap options) {
    gPad->Update();
    double y0 { tuPadyRat(0.) };
    double y1 { tuPadyRat(1.) };
    tuDrawTLine(x,y0,x,y1,options);
}
//*
void tu_scaleByBinWidth(TH1* hg, double scale_factor) {
    for (int i{1}; i<=(int) hg->GetNbinsX(); ++i) {
        if (hg->GetBinContent(i)) {
            double factor { scale_factor / hg->GetBinWidth(i) };
            hg->SetBinContent(i, hg->GetBinContent(i) * factor);
            hg->SetBinError(i, hg->GetBinError(i) * factor);
        }
    }
};
//*
void tu_scaleByBinWidth(TH2* hg, double scale_factor, bool byXwidth, bool byYwidth) {
    if (!byXwidth && !byYwidth) {
        hg->Scale(scale_factor);
        return;
    }
    for (int x{1}; x<=(int) hg->GetNbinsX(); ++x) {
        double x_wide { byXwidth ? hg->GetXaxis()->GetBinWidth(x) : 1. };
        for (int y{1}; y<=(int) hg->GetNbinsY(); ++y)  {
            double y_wide { byYwidth ? hg->GetYaxis()->GetBinWidth(y) : 1. };
            double factor { scale_factor / y_wide / x_wide };
            hg->SetBinContent(x, y, hg->GetBinContent(x,y) * factor);
            hg->SetBinError  (x, y, hg->GetBinError(x,y)   * factor);
        }
    }
};
//*
vector<string> tu_split_string(string str){
    // split strings like: "first name|| second name|| third name|| fourth name"
    vector<string> vec;
    TString s1 { str };
    TPRegexp r1("\\s*(.*?)\\s*\\|\\|(.*)");
    TPRegexp strip("\\s*([^\\s]+.*?[^\\s]*)\\s*");
    while (s1(r1) != "") {
        TObjArray* subStrL = r1.MatchS(s1);
        const Int_t nrSubStr = subStrL->GetLast()+1;
        /* TString key, value; */
        if (nrSubStr > 1) {
            vec.push_back  (((TObjString *)subStrL->At(1))->GetString().Data());
        } else {
            break;
        }
        if (nrSubStr > 2) {
            s1 = ((TObjString*)subStrL->At(2))->GetString();
        }
    }
    // get a final entry
    if (s1(strip) != "") {
        /* cout << " s1 " << s1 << endl; */
        TObjArray* subStrL = strip.MatchS(s1);
        const Int_t nrSubStr = subStrL->GetLast()+1;
        /* TString key, value; */
        /* cout << " nsub " << nrSubStr << endl; */
        if (nrSubStr > 0) {
            /* cout << " add " << s1 << endl; */
            /* cout << "0:"<<  ((TObjString *)subStrL->At(0))->GetString().Data() << endl; */
            /* cout << "1:"<<  ((TObjString *)subStrL->At(1))->GetString().Data() << endl; */
            vec.push_back  (((TObjString *)subStrL->At(1))->GetString().Data());
            /* vec.push_back  (((TObjString *)subStrL->At(0))->GetString().Data()); */
            /* vec.push_back  (((TObjString *)subStrL->At(2))->GetString().Data()); */
        } 
    }
    return vec;
};
//*
void tu_print(TH1* hg, const char* tag) {
    cout << "  " << tag << endl;
    cout << " Printing histogram: " << hg->GetName() 
         << "  xTitle " << hg->GetXaxis()->GetTitle()
         << "  yTitle " << hg->GetYaxis()->GetTitle() << endl;
    int nbins = hg->GetXaxis()->GetNbins();
    for (int i=1; i<=nbins; ++i) {
        /* cout << i << endl; */
        cout << Form(" bin (%2i)  content(%10.5g)  error(%10.5g)",
                i, hg->GetBinContent(i), hg->GetBinError(i)) << endl;
    }
    cout << " end printing histogram" << endl;
};
//*
TLegend* tuNewTLegend() {
    TLegend *leg = new TLegend(0.6455533,0.6332006,0.8911167,0.8938613,NULL,"brNDC");
    tu_fmt(leg);
    return leg;
};
//*
int tu_geant05(int geantid) {
    switch (geantid) {
        case 8: return 0;
        case 9: return 1;
        case 11: return 2;
        case 12: return 3;
        case 14: return 4;
        case 15: return 5;
    }
    return -1;
};
//*
const char* tu_geant05_ascii(int geantid) {
    switch (geantid) {
        case 8: return "pi";
        case 9: return "antipi";
        case 11: return "K";
        case 12: return "antiK";
        case 14: return "p";
        case 15: return "pbar";

        case 0: return "pi";
        case 1: return "antipi";
        case 2: return "K";
        case 3: return "antiK";
        case 4: return "p";
        case 5: return "pbar";
    }
    return "none";
};
//*
const char* tu_geant05_greek(int geantid) {
    switch (geantid) {
        case 8: return "#pi";
        case 9: return "#pi^{-}";
        case 11: return "K";
        case 12: return "K^{-}";
        case 14: return "p";
        case 15: return "#bar{p}";

        case 0: return "#pi";
        case 1: return "#pi^{-}";
        case 2: return "K";
        case 3: return "K^{-}";
        case 4: return "p";
        case 5: return "#bar{p}";
    }
    return "none";
};
//*
//*
//*
/* TH1D* tu_BayesUnfold(TH1D* data, TH1D* T, TH2D* R, int iRepUnfold, TH1D* M) { */
/*     // return data unfolded with response matrix R */
/*     // note that truth T is included for misses, and if present */
/*     // M will include fakes */
/*     if (M==nullptr) M = R->ProjectionX(tuUniqueName()); */
/*     RooUnfoldResponse* rooUnfRes = new RooUnfoldResponse (M,T,R,tuUniqueName()); */
/*     RooUnfoldBayes*    bayes     = new RooUnfoldBayes(rooUnfRes, data, iRepUnfold); */
/*     TH1D* unfolded = (TH1D*) bayes->Hunfold(); */
/*     unfolded->SetName(tuUniqueName()); */
/*     delete rooUnfRes; */
/*     delete bayes; */
/*     return unfolded; */
/* }; */
//*
TH1D* tuNorm(TH1D* hg, const char which) {
    switch (which) {
        case 'n':
            hg->Scale(1./hg->Integral());
            return hg;
        case 'c':
            hg->Scale(1./hg->Integral()/hg->GetBinWidth(1));
            return hg;
        case 'v':
            {
                double sum=0.;
                TAxis *axis = hg->GetXaxis();
                for (int i{1};i<axis->GetNbins(); ++i) 
                    sum += hg->GetBinContent(i) * axis->GetBinWidth(i);
                hg->Scale(1./sum);
                return hg;
            }
        default:
            cout << "Error in tuNorm, selected option '"<<which<<"' is not valid."<<endl;
            cout << "Only options:" << endl;
            cout << " n : 'none' -- don't scale by any width"<<endl;
            cout << " c : 'constant' -- scale by 1/hg->GetBinWidth(1) width"<<endl;
            cout << " v : 'variable' -- scale by all bin widths individually"<<endl;
            return hg;
    }
};
//*
double tu_get_box_integral(TH2D* hg, pair<double,double>p, pair<double,double>q){
    int x0 { (p.first==0  && q.first==0)  ? 1 : hg->GetXaxis()->FindBin(p.first)};
    int x1 { (p.first==0  && q.first==0)  ? hg->GetXaxis()->GetNbins() : hg->GetXaxis()->FindBin(q.first)};
    int y0 { (p.second==0 && q.second==0) ? 1 : hg->GetYaxis()->FindBin(p.second)};
    int y1 { (p.second==0 && q.second==0) ? hg->GetYaxis()->GetNbins() : hg->GetYaxis()->FindBin(q.second)};
    return hg->Integral(x0,x1,y0,y1);
};
//*
double tu_get_box_integral(TProfile2D* pr, pair<double,double>p, pair<double,double>q){
    int x0 { (p.first==0  && q.first==0)  ? 1 : pr->GetXaxis()->FindBin(p.first)};
    int x1 { (p.first==0  && q.first==0)  ? pr->GetXaxis()->GetNbins() : pr->GetXaxis()->FindBin(q.first)};
    int y0 { (p.second==0 && q.second==0) ? 1 : pr->GetYaxis()->FindBin(p.second)};
    int y1 { (p.second==0 && q.second==0) ? pr->GetYaxis()->GetNbins() : pr->GetYaxis()->FindBin(q.second)};
    return pr->Integral(x0,x1,y0,y1);
};
//*
double tu_get_box_mean(TProfile2D* pr, pair<double,double>p, pair<double,double>q){
    if (p.first  != 0 || q.first != 0)  pr->GetXaxis()->SetRangeUser(p.first, q.first );
    if (p.second != 0 || q.second != 0) pr->GetYaxis()->SetRangeUser(p.second,q.second);
    return pr->GetMean(3);
};
//*
vector<double> tu_vec_BinContent(TH1* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinContent(i));
    return vec;
};
//*
vector<double> tu_vec_BinError(TH1* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinError(i));
    return vec;
};
vector<double> tu_vec_BinCenter(TH1* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    TAxis* ax = hg->GetXaxis();
    for (int i=i0;i<i1;++i) {
        vec.push_back(ax->GetBinCenter(i));
    }
    return vec;
};
/* /1* template<class T> *1/ */
pair<int, double*> tu_vecBinEdge(TH1* hg) {
    int n_bins = hg->GetXaxis()->GetNbins();
    double *edges = new double[n_bins+1];
    for (int i{1};i<=n_bins+1;++i){
        edges[i-1] = hg->GetXaxis()->GetBinLowEdge(i);
    }
    return {n_bins,edges};
};
vector<double> tu_vecBinContent(TH1* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinContent(i));
    return vec;
};
vector<double> tu_vecBinError  (TH1* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinError(i));
    return vec;
};
float tu_dphi(float phi0, float phi1) {
    float rval = phi1-phi0;
    while (rval < -M_PI) rval += 2*M_PI;
    while (rval >  M_PI) rval -= 2*M_PI;
    return rval;
};
float tu_absDphi(float phi0, float phi1) { return TMath::Abs(tu_dphi(phi0,phi1)); };

/* bool tu_isAbsTransPhi(float phi0, float phi1, float lo_bound, float hi_bound){ */
/*     float dphi = TMath::Abs(tu_dphi(phi0,phi1)); */
/*     return (dphi>=lo_bound && dphi<=hi_bound); */
/* }; */
/* bool tu_isAbsTrans358pi(float phi0, float phi1, float lo_bound, float hi_bound){ */
/*     float dphi = TMath::Abs(tu_dphi(phi0,phi1)); */
/*     return (dphi>=lo_bound && dphi<=hi_bound); */
/* }; */

float tu_02pi(float &phi){
    while (phi<0)        phi += 2*M_PI;
    while (phi>2*M_PI) phi -= 2*M_PI;
    return phi;
};
float tu_02pi(float  phi){
    while (phi<0)        phi += 2*M_PI;
    while (phi>2*M_PI) phi -= 2*M_PI;
    return phi;
};

vector<double> tuQuantiles(TH1D* hg, vector<double> percents){
    cout << hg->GetName() << endl;
    int n {(int)percents.size()};
    double *x = new double[n];
    double *q = new double[n];
    for (int i{0}; i<n; ++i) x[i] = percents[i];
    hg->GetQuantiles(n,q,x);
    vector<double> vec;
    for (int i{0}; i<n; ++i) vec.push_back(q[i]);
    delete[] x;
    delete[] q;
    return vec;
};
string  tuStringVec(vector<double> vec, const char* name, const char* formatter){
    if (vec.size()==0) return "";
    ostringstream os;
    const char* fmt = Form("%%%s",formatter);
    os << "vector<double> " << name << " {" << Form(fmt,vec[0]);
    for (int i{1}; i<(int)vec.size();++i) os << ", " << Form(fmt,vec[i]);
    os << "};";
    return os.str();
}

int tu_count_digits(int n, int min_val) {
    int cnt = 0;
    while (n != 0) { 
        n /= 10;
        ++cnt;
    }
    if (cnt < min_val) cnt = min_val;
    return cnt;
};
int tuSum(const vector<int> vec)  { return std::accumulate(vec.begin(), vec.end(), 0); };
int tuSum(const vector<bool> vec) { return std::accumulate(vec.begin(), vec.end(), 0); };

map<int,string>    tuReadIntStrMap(const char* file, tuOptMap options, tuOptMap dict){
    dict += options;
    /* cout << " map =1 : " << options["tag"].str() << endl; */
    auto data = tuReadStrVec(file, options);
    if (data.size() % 2) 
        throw std::runtime_error(
                "Error in tu_ReadIntStrMap : must have even number of entries");
    map<int,string> M;
    for (unsigned int i=0; i<(data.size()/2); ++i) {
        int index = atoi(data[i*2].c_str());
        M[index] = data[i*2+1];
    }
    return M;
};

vector<int> tuReadIntVec(const char* in_file, int col, bool sort, bool strip_commas) {
    vector<int> vec;
    ifstream file;
    file.open(in_file);
    ostringstream msg;
    if (!file.is_open()) {
        msg << "fatal error: could not open int map file \"" 
            << in_file << " in tuReadIntVec" << endl;
        throw std::runtime_error(msg.str());
    }
    /* int n_req { index_column > data_column ? index_column : data_column }; */
    string line;
    bool   read_all_cols = (col==-1);
    while (getline(file,line)) {
        line.append(" ");
        if (strip_commas && (line.find(",") != string::npos)) {
            TString trans = line;
            trans.ReplaceAll(","," ");
            line = trans;
        }
        stringstream words(line);
        TString word;
        int i {-1};
        bool found_val  {false};
        bool is_comment {false};
        while (words >> word) {
            ++i;
            if (!word.IsDigit()) {
                is_comment = true;
                break;
            }
            if (read_all_cols) vec.push_back(word.Atoi());
            if (i == col) {
                found_val = true;
                vec.push_back(word.Atoi());
                break;
            }
        }
        if (is_comment || found_val || read_all_cols) continue;
        msg << "fatal error: could not read col " << col << endl
            << "  in input line \""<<line<<"\""  << endl
            << "  in input file \"" << in_file <<"\" in tuReadIntVec" << endl;
        throw std::runtime_error(msg.str());
    }
    file.close();
    if (sort) std::sort(vec.begin(), vec.end());
    return vec;
};

//----------
map<string,string> tu_VecStrToMapStrStr (vector<string> data) 
{
    if (data.size() % 2) 
        throw std::runtime_error(
                "Error in tu_VecStrToMapStrStrg : must have even number of entries");
    map<string,string> M;
    for (unsigned int i=0; i<(data.size()/2); ++i) {
        M[data[i*2]] = data[i*2+1];
    }
    return M;
};
map<string,string> tuReadMapStrStr(const char* in_file, const char* tag, 
        tuOptMap options, tuOptMap dict) 
{
    dict += options;
    dict("tag") = tag;
    return tu_VecStrToMapStrStr(tuReadStrVec(in_file, dict));
};
map<string,string> tuReadMapStrStr(const char* in_file, tuOptMap options, tuOptMap dict) {
    dict += options;
    return tu_VecStrToMapStrStr(tuReadStrVec(in_file, dict));
};


vector<string> tuReadStrVec(const char* in_file, const char* tag, 
        tuOptMap options, tuOptMap dict) {
    dict += options;
    dict("tag") = tag;
    return tuReadStrVec(in_file, dict);
};
vector<string> tuReadStrVec(const char* in_file, tuOptMap options, tuOptMap dict) {
    dict += options;

    bool use_column = (dict("column").str!="all");
    int  which_column {0};
    if  (use_column) { which_column = dict("column"); }
    /* cout << "a0 : " << dict("tag").str() << endl; */
    bool use_tag {dict("tag").str != "none"};
    string tag = use_tag ? dict("tag") : "";
    bool in_tag { !use_tag };
    bool strip_commas { dict["strip_commas"] };

    /* cout << " Tag: " << tag << endl; */

    // logic:
    // per line:
    //      if strip_commas replace all commas with spaces
    //      if use_tag && !in_tag: ignore everything until the flag is found
    //      if use_tag && in_tag: terminate in </tag> if found
    //      if use_column: only read column entry

    vector<string> vec;
    ifstream file;
    file.open(in_file);
    ostringstream msg;
    if (!file.is_open()) {
        msg << "fatal error: could not open file \"" 
            << in_file << " in tuReadStrVec" << endl;
        throw std::runtime_error(msg.str());
    }
    string line;
    /* bool   read_all_cols = (col==-1); */
    bool found_end_tag {false};
    bool found_start_tag {false};
    while (getline(file,line)) {
        line.append(" ");
        if (strip_commas && (line.find(",") != string::npos)) {
            TString trans = line;
            trans.ReplaceAll(","," ");
            line = trans;
        }
        stringstream words(line);
        TString word;
        int col {-1};
        while (words >> word) {
            if (use_tag) {
                if ( in_tag ) {
                    if (tuWordIsEndTag(word, tag)) { 
                        found_end_tag  = true;
                        break;
                    }
                } else {
                    if (tuWordIsTag(word,tag)) {
                        found_start_tag = true;
                        in_tag = true;
                    }
                    continue;
                }
            }
            if (tuIsAnyTag(word)) continue;
            /* if (word.IsFloat()) { */
            ++col;
            if (!use_column || col==which_column) {
                if (word.Contains("--") && !options.has("keep--")) {
                    word = word.ReplaceAll("--"," ");
                }
                vec.push_back(word.Data());
            }
        }
        if (found_end_tag ) break;
    }
    file.close();
    if (use_tag && !found_start_tag) {
        throw std::runtime_error(
                "fatal: didn't find starting tag \"<"
             + tag + ">\" in file \"" + in_file + "\"");
    };
    if (use_tag && !found_end_tag) {
        throw std::runtime_error(
                "fatal: didn't find end tag \"</"
             + tag + ">\" in file \"" + in_file +"\"");
    };
             
    if ( dict["sort"] ) std::sort(vec.begin(), vec.end());
    return vec;
};

vector<double> tuReadValVec(const char* in_file, const char* tag, 
        tuOptMap options, tuOptMap dict) {
    dict += options;
    dict("tag") = tag;
    return tuReadValVec(in_file, dict);
};
vector<double> tuReadValVec(const char* in_file, tuOptMap options, tuOptMap dict) {
    dict += options;

    bool use_column = (dict("column").str!="all");
    int  which_column {0};
    if  (use_column) { which_column = dict("column"); }
    bool use_tag {dict("tag").str != "none"};
    string tag = use_tag ? dict("tag") : "";
    bool in_tag { !use_tag };
    bool strip_commas { dict["strip_commas"] };

    // logic:
    // per line:
    //      if strip_commas replace all commas with spaces
    //      if use_tag && !in_tag: ignore everything until the flag is found
    //      if use_tag && in_tag: terminate in </tag> if found
    //      if use_column: only read column entry

    vector<double> vec;
    ifstream file;
    file.open(in_file);
    ostringstream msg;
    if (!file.is_open()) {
        msg << "fatal error: could not open int map file \"" 
            << in_file << " in tuReadIntVec" << endl;
        throw std::runtime_error(msg.str());
    }
    /* int n_req { index_column > data_column ? index_column : data_column }; */
    string line;
    /* bool   read_all_cols = (col==-1); */
    bool found_end_tag {false};
    bool found_start_tag {false};
    while (getline(file,line)) {
        line.append(" ");
        if (strip_commas && (line.find(",") != string::npos)) {
            TString trans = line;
            trans.ReplaceAll(","," ");
            line = trans;
        }
        stringstream words(line);
        TString word;
        int col {-1};
        bool found_val  {false};
        while (words >> word) {
            if (word.IsFloat()) {
                ++col;
                if (in_tag) {
                    if (!use_column || col==which_column) vec.push_back(word.Atof());
                    if (use_column && col == which_column) break;
                }
                continue;
            } 
            bool is_a_tag { tuIsAnyTag(word) };
            if (!is_a_tag) break; // neither number or tag, so is a comment
            // here: is_a_tag == true
            if (!use_tag) continue;
            // use_tag == true && this is a tag
            if ( in_tag ) {
                if (tuWordIsEndTag(word, tag)) { 
                    found_end_tag  = true;
                    break;
                }
            } else {
                if (tuWordIsTag(word,tag)) {
                    found_start_tag = true;
                    in_tag = true;
                }
            }
        }
        if (found_end_tag ) break;
    }
    file.close();
    if (use_tag && !found_start_tag) {
        throw std::runtime_error(
                "fatal: didn't find starting tag \"<"
             + tag + ">\" in file \"" + in_file + "\"");
    };
    if (use_tag && !found_end_tag) {
        throw std::runtime_error(
                "fatal: didn't find end tag \"</"
             + tag + ">\" in file \"" + in_file +"\"");
    };
             
    if (  dict("sort").val == 1) std::sort(vec.begin(), vec.end());
    return vec;
};

pair<int,double*> tuReadValsPtr(const char* file, tuOptMap options, tuOptMap dict) {
    dict += options;
    auto vec = tuReadValVec(file, dict);
    int i0 = dict("begin_index");
    int i1 = dict("end_index");
    int size_data {(int) vec.size()};
    int size_req = (i1 != -1) 
        ? i1-i0
        : size_data-i0;
    if (i0!=0 || i1!=-1) {
        if (size_req<=0 || size_req > size_data) {
            ostringstream msg;
            msg << "tuReadValsPtr in file " << file << " requested start and end indices " 
                << i0 << " and " << i1
                << " for total size of " << size_req 
                << " out of total available " << size_data << endl;
            throw std::runtime_error(msg.str());
        }
    }
    if (size_data == 0) {
        throw std::runtime_error(
        (string)"tuReadValsPtr  in file " + file + " found no data. ");
    }
    double* vals = new double[size_req];
    for (int i{0}; i<size_req; ++i) {
        vals[i] = vec[i0+i];
    }
    return {size_req, vals};
};

TGraph* tuMakeTGraph(TH1D* hg, bool invert, bool skip_zeros, bool normalize) {
    if (normalize) hg->Scale(1./hg->Integral());
    vector<double> x, y;
    TAxis* axis = hg->GetXaxis();
    for (int i{1}; i<=hg->GetXaxis()->GetNbins(); ++i) {
        if (hg->GetBinContent(i) == 0 && hg->GetBinError(i) == 0 && skip_zeros) continue;
        x.push_back(axis->GetBinCenter(i));
        y.push_back(hg->GetBinContent(i));
    }
    TGraph* gr;
    double lo = hg->GetXaxis()->GetBinLowEdge(1);
    double hi = hg->GetXaxis()->GetBinUpEdge( hg->GetXaxis()->GetNbins());
    if (invert) {
        gr = tuMakeTGraph(y,x);
        // set the limits
        gr->SetMinimum(lo);
        gr->SetMaximum(hi);
    } else {
        gr = tuMakeTGraph(x,y);
        gr->GetXaxis()->SetLimits(lo,hi);
    }
    return gr;
};

TGraph* tuMakeTGraph(TProfile* pr, bool invert, bool skip_zeros, bool normalize) {
    TH1D* hg = (TH1D*) pr->ProjectionX(tuUniqueName());
    if (normalize) hg->Scale(1./hg->Integral());
    return tuMakeTGraph(hg,invert,skip_zeros);
};

TGraphErrors* tuMakeTGraphErrors(TH1D* hg, bool invert, bool skip_zeros, bool normalize ) {
    vector<double> x, x_err, y, y_err;
    TAxis* axis = hg->GetXaxis();
    if (normalize) hg->Scale(1./hg->Integral());
    for (int i{1}; i<=hg->GetXaxis()->GetNbins(); ++i) {
        if (hg->GetBinContent(i) == 0 && hg->GetBinError(i) == 0 && skip_zeros) continue;
        x.push_back(axis->GetBinCenter(i));
        x_err.push_back(0.);
        y.push_back(hg->GetBinContent(i));
        y_err.push_back(hg->GetBinError(i));
        
    }
    TGraphErrors* gr;
    double lo = hg->GetXaxis()->GetBinLowEdge(1);
    double hi = hg->GetXaxis()->GetBinUpEdge( hg->GetXaxis()->GetNbins());
    if (invert) {
        gr = tuMakeTGraphErrors(y,x,x_err,y_err);
        gr->SetMinimum(lo);
        gr->SetMaximum(hi);
    } else {
        gr = tuMakeTGraphErrors(x,y,y_err,x_err);
        gr->GetXaxis()->SetLimits(lo,hi);
    }
    return gr;
};

/* TGraph* tuMakeTGraph(vector<double>& x, vector<double>& y) { */
/*     if (x.size() != y.size()) */ 
/*         throw std::runtime_error("tuMakeTGraph(vec, vec) required vectors of same length"); */
/*     const unsigned int n = x.size(); */
/*     double* xpts = new double[n]; */
/*     double* ypts = new double[n]; */
/*     for (unsigned int i{0}; i<n; ++i) { */
/*         xpts[i] = x[i]; */
/*         ypts[i] = y[i]; */
/*     } */
/*     return new TGraph(n,xpts,ypts); */
/* }; */
TGraphErrors* tuMakeTGraphErrors(vector<double> x, vector<double> y, vector<double> y_err, vector<double> x_err, vector<bool> use) {
    tuMinMax nsize;
    nsize.fill(x.size());
    nsize.fill(y.size());
    nsize.fill(y_err.size());
    if (x_err.size() != 0) nsize.fill(x_err.size());
    if (use.size() != 0)   nsize.fill(use.size());

    if (x.size() != y.size() || x.size() != y_err.size()) 
        cout << " Warning: tuMakeTGraphErrors(vec, vec, vec={}) don't have equal input sizes" << endl
             << " Therefor using shortest vector of size " <<  nsize.max << endl;
    if (x_err.size() != 0 && x_err.size() != x.size())
        cout << " Warning: tuMakeTGraphErrors(vec, vec, vec={}) don't have equal input sizes" << endl
             << " Therefor using shortest vector of size " <<  nsize.max << endl;
    if (use.size() != 0 && use.size() != x.size())
        cout << " Warning: tuMakeTGraphErrors(vec, vec, vec={}) don't have equal input sizes" << endl
             << " Therefor using shortest vector of size " <<  nsize.max << endl;

    if (use.size() != 0) {
        vector<double> x0, y0, y0_err, x0_err;
        for (int i{0}; i<(int)use.size();++i) {
            if (use[i]) {
                x0.push_back(x[i]);
                y0.push_back(y[i]);
                x0_err.push_back( x_err.size() == 0 ? 0. : x_err[i] );
                y0_err.push_back( y_err[i] );
            }
        }
        x=x0;
        y=y0;
        x_err=x0_err;
        y_err=y0_err;
    }

    const unsigned int n = nsize.max;
    /* cout << " size: " << n << " " << x.size() << " " << y.size() << " " << x_err.size() << " " << y_err.size() << endl; */
    double* xpts = new double[n];
    double* ypts = new double[n];
    double* yerr = new double[n];
    double* xerr = new double[n];
    bool use_x_err { x_err.size() != 0 };

    for (unsigned int i{0}; i<n; ++i) {
        xpts[i] = x[i];
        ypts[i] = y[i];
        yerr[i] = y_err[i];
        xerr[i] = ( use_x_err ? x_err[i] : 0. );
    }
    return new TGraphErrors(n,xpts,ypts,xerr,yerr);
};

TGraph* tuMakeTGraph(vector<double> x, vector<double> y) {
    if (x.size() != y.size()) 
        throw std::runtime_error("tuMakeTGraph(vec, vec) required vectors of same length");
    const unsigned int n = x.size();
    double* xpts = new double[n];
    double* ypts = new double[n];
    for (unsigned int i{0}; i<n; ++i) {
        xpts[i] = x[i];
        ypts[i] = y[i];
    }
    return new TGraph(n,xpts,ypts);
};

// find the bin in a vector
int tuwhichbin0(double val, vector<double>& vec) {
    return (int)(std::upper_bound(vec.begin(), vec.end(), val) - vec.begin()-1);
} ;
int tuwhichbin0(double val, TH1D* hg) {
    return (int)(hg->GetXaxis()->FindBin(val)-1);
} ;
int tuwhichbin1(double val, vector<double>& vec) {
    return (int)(std::upper_bound(vec.begin(), vec.end(), val) - vec.begin());
} ;
int tuwhichbin1(double val, TH1D* hg) {
    return (int)(hg->GetXaxis()->FindBin(val));
} ;

double* ax_doubleptr(vector<int> vals) {
    sort(vals.begin(), vals.end());
    int size { static_cast<int>(vals.size()) };
    double* x = new double[size + 1];
    if (size == 0) return x;
    if (size == 1) {
        x[0] = vals[0]-0.5;
        x[1] = vals[0]+0.5;
        return x;
    }
    x[0] = vals[0] - (vals[1] - vals[0])/2.;
    for (int i{1}; i<size; ++i) {
        x[i] = (vals[i]+vals[i-1])/2.;
    }
    x[size] = vals[size-1]+(vals[size-1]-vals[size-2])/2.;
    return x;
};
double tu_D(double x0,double y0,double x1,double y1) 
{ return TMath::Sqrt( TMath::Sq(x1-x0)+TMath::Sq(y1-y0)); };
double tu_D2(double x0,double y0,double x1,double y1) 
{ return TMath::Sq(x1-x0)+TMath::Sq(y1-y0); };
double tu_R(double x0,double y0,double x1,double y1) 
{ return TMath::Sqrt( TMath::Sq(x1-x0)+TMath::Sq(tu_dphi(y1,y0))); };
double tu_R2(double x0,double y0,double x1,double y1) 
{ return TMath::Sq(x1-x0)+TMath::Sq(tu_dphi(y1,y0)); };



/* RooUnfoldResponse tuMakeRooUnfoldResponse( */
/*      int nbins, double lo_bin, double hi_bin, */ 
/*      const char* tag, const char* title */
/* ) { */
/*     const char* this_tag { */
/*         (strcmp(tag,"")==0) */
/*             ? tuUniqueName() */
/*             : tag */ 
/*     }; */
/*     TH1D truth { */
/*             Form("%s_truth",this_tag), */
/*             Form("%s;truth;N",title), */
/*             nbins, lo_bin, hi_bin }; */
/*     TH1D measured { */
/*             Form("%s_measured",this_tag), */
/*             Form("%s;measured;N",title), */
/*             nbins, lo_bin, hi_bin }; */
/*     return {&measured, &truth, Form("%s_RooUnfR",this_tag), title}; */
/* }; */

/* RooUnfoldResponse tuMakeRooUnfoldResponse( */
/*      int nbins, double* edges, */ 
/*      const char* tag, const char* title */
/* ) { */
/*     const char* this_tag { */
/*         (strcmp(tag,"")==0) */
/*             ? tuUniqueName() */
/*             : tag */ 
/*     }; */
/*     TH1D truth { */
/*             Form("%s_truth",this_tag), */
/*             Form("%s;truth;N",title), */
/*             nbins, edges }; */
/*     TH1D measured { */
/*             Form("%s_measured",this_tag), */
/*             Form("%s;measured;N",title), */
/*             nbins,  edges }; */
/*     return {&measured, &truth, Form("%s_RooUnfR",this_tag), title}; */
/* }; */


/* RooUnfoldResponse tuMakeRooUnfoldResponse( */
/*      int nb_measured, double lo_measured, double hi_measured, */
/*      int nb_truth, double lo_truth, double hi_truth, */
/*      const char* tag, const char* title */
/* ) { */
/*     const char* this_tag { */
/*         (strcmp(tag,"")==0) */
/*             ? tuUniqueName() */
/*             : tag */ 
/*     }; */
/*     TH1D truth{ */
/*             Form("%s_truth",this_tag), */
/*             Form("%s;truth;N",title), */
/*             nb_truth, lo_truth, hi_truth }; */
/*     TH1D measured { */
/*             Form("%s_measured",this_tag), */
/*             Form("%s;measured;N",title), */
/*             nb_measured, lo_measured, hi_measured }; */
/*     return {&measured, &truth, Form("%s_RooUnfR",this_tag), title}; */
/* }; */
/* RooUnfoldResponse tuMakeRooUnfoldResponse( */
/*      int nb_measured, double* edges_measured, */
/*      int nb_truth, double* edges_truth, */
/*      const char* tag, const char* title */
/* ) { */
/*     const char* this_tag { */
/*         (strcmp(tag,"")==0) */
/*             ? tuUniqueName() */
/*             : tag */ 
/*     }; */
/*     TH1D truth { */
/*             Form("%s_truth",this_tag), */
/*             Form("%s;truth;N",title), */
/*             nb_truth, edges_truth}; */
/*     TH1D measured { */
/*             Form("%s_measured",this_tag), */
/*             Form("%s;measured;N",title), */
/*             nb_measured, edges_measured}; */
/*     return {&measured, &truth, Form("%s_RooUnfR",this_tag), title}; */
/* }; */
// return which bin (starting from 0) the data is in: lower bound <= val < upper bound
/* int tuwhichbin(int, double*); */
/* int tuwhichbin(vector<double>, vector<int> remap); // return which bin (starting from 0) the data is in: lower bound <= val < upper bound */
/* int tuwhichbin(int, double*,   vector<int> remap); */

double tuRatCircleOverLine (double R, double d) {
    R = TMath::Abs(R);
    d = TMath::Abs(d);
    if (R <= d) return 0.;
    return TMath::ACos(d/R)/M_PI;
};
double tuRatCircleInTwoParallelLines (const double d0,const double d1,double C,double R) {
    if (C<d0) C = d0;
    if (C>d1) C = d1;
    return 1 - tuRatCircleOverLine(R,C-d0)
             - tuRatCircleOverLine(R,d1-C);
};
double tuPolyP6_a0_a1x_a2xx_a3y_a4yy_a5xy(double* x, double *p){
    return p[0]      + p[1]*x[0] + p[2]*x[0]*x[0]
                     + p[3]*x[1] + p[4]*x[1]*x[1] + p[5]*x[0]*x[1];
};

bool tuWordIsTag   (string  word, string tag)  { return (word == ("<" +tag+">")); };
bool tuWordIsEndTag(string  word, string tag)  { return (word == ("</" +tag+">")); };
bool tuWordIsTag   (TString word, string tag) { return (word == ("<" +tag+">")); };
bool tuWordIsEndTag(TString word, string tag) { return (word == ("</" +tag+">")); };
bool tuIsAnyTag    (string  word) { return tuIsAnyTag((TString)word); };
bool tuIsAnyTag    (TString word) {
    if (!word.BeginsWith("<")) return false;
    return word.EndsWith(">");
};

void tu_normByRow(TH2D* hg, double factor, bool use_max_val) {
    int nCols = hg->GetXaxis()->GetNbins();
    int nRows   = hg->GetYaxis()->GetNbins(); 

    for (int row{1}; row<=nRows; ++row) {
        double mult;
        if (use_max_val) {
            double vmax { hg->GetBinContent(1,row) };
            for (int col{1}; col <= nCols; ++col) {
                if (hg->GetBinContent(col,row) > vmax) {
                    vmax = hg->GetBinContent(col,row);
                }
            }
            mult = 1./vmax;
        } else {
            mult = factor / hg->Integral(1,nCols,row,row);
        }


        for (int col {1}; col <= nCols; ++col) {
            hg->SetBinContent(col, row, hg->GetBinContent(col, row) * mult);
            hg->SetBinError  (col, row, hg->GetBinError  (col, row) * mult);
        }
    };
};

vector<double> tu_print_first_blank(TH2D* hg) {
    vector<double> vec;
    int nY = hg->GetYaxis()->GetNbins();
    int nX = hg->GetXaxis()->GetNbins();
    for (int y = 1; y<=nY; ++y) {
        bool found = false;
        for (int x = 1; x<=nX; ++x) {
            int ibin = hg->GetBin(x,y);
            if (hg->GetBinContent(ibin) == 0) {
                vec.push_back(hg->GetXaxis()->GetBinCenter(x));
                found = true;
                break;
            }
        }
        if (!found) vec.push_back(hg->GetXaxis()->GetBinCenter(nX));
    }
    int i=0;
    for (auto v : vec) cout << Form("  %3i  -> %5.2f",i++,v) << endl;
    return vec;
};


/* void tuTrimSmallBins(TH2D* hg, int Nmin, bool cut_underover_flow=true) { */
/*     TAxis* Xaxis = hg->GetXaxis(); */
/*     TAxis* Yaxis = hg->GetYaxis(); */
/*     int nXbins = Xaxis->GetNbins(); */
/*     int nYbins = Xaxis->GetNbins(); */
/*     if (cut_underover_flow) { */
        
/* tuS_keepcut_stats tu_keepcutbin(TH1* h, int ibin, double minval, double zero_val, double zero_err) { */
/*     tuS_keepcut_stats stats{}; */
/*     double val = h->GetBinContent(ibin); */
/*     if (val == 0 && h->GetBinError(ibin)==0) return {}; */
/*     if (val < minval) { */
/*         ++stats.n_cut; */
/*         stats.sum_cut = val-zero_val; */
/*         stats.was_cut = true; */
/*         h->SetBinContent(ibin, zero_val); */

/*         if (zero_val != 0 && zero_err == 0.)  h->SetBinError(ibin, sqrt(zero_val)); */ 
/*         else  h->SetBinError(ibin, 0.); */ 

/*         int x, y, z; */
/*         h->GetBinXYZ(ibin,x,y,z); */
/*         stats.x = h->GetXaxis()->GetBinCenter(x); */
/*         stats.y = h->GetYaxis()->GetBinCenter(y); */
/*     } else { */
/*         ++stats.n_keep; */
/*         stats.sum_keep = val; */
/*     } */
/*     /1* if (stats.was_cut) cout << " - " << stats.x << " " << stats.y << endl; *1/ */
/* return stats; */
/* }; */

/* tuS_keepcut_stats tu_cullsmallbins( */
/*         TH1* h, */ 
/*         double min_val, */ 
/*         tu area, */ 
/*         bool remove_overunderflow) */ 
/* { */
/*     if (remove_overunderflow) */ 
/*         for (auto i : tu_binvec(h, tu::overunderflow)) tu_setbinzero(h,i); */

/*     tuS_keepcut_stats stats; */
/*     for (auto i : tu_binvec(h, area)) { */
/*         stats += tu_keepcutbin(h, i, min_val); */
/*         /1* if (stats.was_cut) cout << stats.x << " " << stats.y << endl; *1/ */
/*         /1* cout << " stats: " << stats.was_cut << endl; *1/ */
/*     } */
/*     return stats; */
/* }; */
/* bool tu_cullsmallbins_( */
/*         TH1* h, */ 
/*         double min_val) */
/* { */
/*     for (auto i : tu_binvec(h, tu::in)) { */
/*         tu_keepcutbin(h, i, min_val); */
/*     } */
/*     return true; */
/* }; */

/* vector<int> tu_binvec(TH1* h, tu loc) { */
/*     vector<int> vec{}; */
/*     if (h->GetYaxis()->GetNbins() == 1) { */
/*         switch (loc) { */
/*             case tu::overunderflow: */
/*                 vec = { 0, h->GetXaxis()->GetNbins()+1 }; */
/*                 break; */
/*             case tu::all: */
/*                 vec = { 0, h->GetXaxis()->GetNbins()+1 }; */
/*                 // no break -- continue to next block */
/*             case tu::in: */
/*                 for (int i{1}; i<=h->GetXaxis()->GetNbins(); ++i) */ 
/*                     vec.push_back(i); */
/*                 break; */
/*             default: break; */
/*         } */
/*     } else { */
/*         int nX = h->GetXaxis()->GetNbins(); */
/*         int nY = h->GetYaxis()->GetNbins(); */
/*         switch (loc) { */
/*             case tu::overunderflow: */ 
/*                 for (auto x{0}; x<=nX+1; ++x) vec.push_back(h->GetBin(x,0)); */
/*                 for (auto x{0}; x<=nX+1; ++x) vec.push_back(h->GetBin(x,nY+1)); */
/*                 for (auto y{1}; y<=nY;   ++y) vec.push_back(h->GetBin(0,y)); */
/*                 for (auto y{1}; y<=nY;   ++y) vec.push_back(h->GetBin(nX+1,y)); */
/*                 break; */
/*             case tu::all: */
/*                 vec = tu_binvec(h, tu::overunderflow); */
/*                 // no break -- continue to next block */
/*             case tu::in: */
/*                 for (auto y{1}; y<nY+1; ++y) */
/*                 for (auto x{1}; x<nX+1; ++x) */
/*                         vec.push_back(h->GetBin(x,y)); */
/*                 break; */
/*             default: break; */
/*         } */
/*     } */
/*     return vec; */
/* }; */

TH2D* tu_cut_high_sigmaX(TH2D* hg, double n_sigma, double offset, bool print) {
    // cut all bin content above n_sigma from mean

    if (n_sigma==0) return hg;
    TAxis* y_axis = hg->GetYaxis();
    TAxis* x_axis = hg->GetXaxis();
    double pre_sum = hg->Integral();
    for (int y=1;y<y_axis->GetNbins();++y) {
        TH1D* proj = hg->ProjectionX(tuUniqueName(),y,y);
        double mu = proj->GetMean();
        double sigma = proj->GetStdDev();
        int i_bin = x_axis->FindBin(mu+n_sigma*sigma+offset);
        for (int i=i_bin+1;i<=x_axis->GetNbins();++i) {
            if (hg->GetBinContent(i,y)!=0) {
                hg->SetBinContent(i,y,0.);
                hg->SetBinError(i,y,0.);
            }
        }

    }
    double post_sum = hg->Integral();
    if (print) cout << " percent cut: " << 100.*(pre_sum-post_sum)/pre_sum << endl;
    return hg;
};
TH1D* tu_cut_high_sigmaX(TH1D* hg, double n_sigma, double offset, bool print) {
    // cut all bin content above n_sigma from mean
    if (n_sigma==0) return hg;
    TAxis* x_axis = hg->GetXaxis();
    double pre_sum = hg->Integral();
    double mu = hg->GetMean();
    double sigma = hg->GetStdDev();
    int i_bin = x_axis->FindBin(mu+n_sigma*sigma+offset);
    for (int i=i_bin+1;i<=x_axis->GetNbins();++i) {
      if (hg->GetBinContent(i)!=0) {
          hg->SetBinContent(i,0.);
          hg->SetBinError(i,0.);
      }
    }
    double post_sum = hg->Integral();
    if (print) cout << " percent cut: " << 100.*(pre_sum-post_sum)/pre_sum << endl;
    return hg;
};

double tu_setbinzero(TH1* hg, int bin, double val, double err)
{ 
    double p_val = hg->GetBinContent(bin);
    hg->SetBinContent(bin, val);
    if (val != 0. && err == 0.) err = sqrt(val);
    hg->SetBinError(bin, err);
    return p_val;
};


const char* tu_cutdiff(int a, int b, const char* fmt) {
    int diff {b-a};
    return Form(fmt,b,diff,(double)diff/a);
};

TH1* tuSetCntErrors(TH1* hg) {
    for (int i{0}; i<hg->GetNcells(); ++i) {
        if (hg->GetBinContent(i) != 0) {
            hg->SetBinError(i, TMath::Sqrt(hg->GetBinContent(i)));
        }
    }
    return hg;
};
TH1* tuAddBinCnt(TH1* hg_to, TH1* hg_fr, bool set_cnt_errors, bool rm_underoverflow) {
    // Add all content from bin centers in hg_from to values in corresponding bins in hg_to
    if (rm_underoverflow) hg_fr->ClearUnderflowAndOverflow();
    int x0, y0, z0;
    TAxis* x_ax = hg_fr->GetXaxis();
    TAxis* y_ax = hg_fr->GetYaxis();
    TAxis* z_ax = hg_fr->GetZaxis();

    if (true) for (int i{0}; i<hg_fr->GetNcells(); ++i) {
        double val { hg_fr->GetBinContent(i) }; 
        if (val==0) continue;
        hg_fr->GetBinXYZ(i,x0,y0,z0);
        int i_to { hg_to->FindBin(
                x_ax->GetBinCenter(x0),
                y_ax->GetBinCenter(y0),
                z_ax->GetBinCenter(z0)) };
        hg_to->SetBinContent( i_to, hg_to->GetBinContent(i_to) + val );
    }
    if (rm_underoverflow) hg_to->ClearUnderflowAndOverflow();
    if (set_cnt_errors) tuSetCntErrors(hg_to);
    return hg_to;
};
TH1* tu_build_CDF(TH1D* _in, int first_bin, int last_bin, double alpha) {
    // build the bins 
    //
    TAxis* ax_in = _in->GetXaxis();
    vector<double> bins_edge;
    if (last_bin==0) last_bin = ax_in->GetNbins();
    for (int i=first_bin; i<=last_bin; ++i) bins_edge.push_back(ax_in->GetBinLowEdge(i));
    bins_edge.push_back(ax_in->GetBinUpEdge(last_bin));
    
    int     nbins = bins_edge.size()-1;
    double* p_edges = new double[nbins+1];
    for (auto i=0;i<=nbins;++i) p_edges[i] = bins_edge[i];

    /* tuBinVec bins { bins_edge} ; */
    /* cout << " bins: " << bins << endl; */
    TH1D* hg_cdf = new TH1D(tuUniqueName(), Form("%s;%s;%s",
                _in->GetTitle(), _in->GetXaxis()->GetTitle(), _in->GetYaxis()->GetTitle()),
            nbins, p_edges );
    
    delete[] p_edges;
    hg_cdf->SetMarkerStyle(_in->GetMarkerStyle());
    hg_cdf->SetMarkerColorAlpha(_in->GetMarkerColor(),alpha);
    hg_cdf->SetMarkerSize(_in->GetMarkerSize());
    hg_cdf->SetLineColorAlpha(_in->GetLineColor(),alpha);

    // warning: this only does 
    /* TH1D* hg_cdf = (TH1D*) _in->Clone(tuUniqueName()); */ 
    /* hg_cdf->Reset(); */
    TAxis* ax = hg_cdf->GetXaxis();
    /* if (last_bin==0) last_bin = ax->GetNbins(); */

    double sum_all {0};
    double err_sq_all {0};
    for (int i{first_bin}; i<=last_bin; ++i) {
        double width { ax->GetBinWidth(i) };
        sum_all += _in->GetBinContent(i) * width;
        err_sq_all += TMath::Sq(_in->GetBinError(i)*width);
    }
    double rel_err_all_sq = err_sq_all / TMath::Sq(sum_all);

    double sum_cdf {0};
    double err_cdf_sq {0};
    int i_to = 1;
    for (int i{first_bin}; i<=last_bin; ++i) {
        double width { ax->GetBinWidth(i) };
        sum_cdf += _in->GetBinContent(i) * width;
        err_cdf_sq += TMath::Sq(_in->GetBinError(i)*width);

        double in_val = err_cdf_sq / TMath::Sq(sum_cdf) + rel_err_all_sq - 2*err_cdf_sq/sum_cdf/sum_all;
        /* if (in_val<0) */ 

        double bin_err = (in_val<0) ? 0. : sum_cdf/sum_all * TMath::Sqrt(in_val);
                /* err_cdf_sq / TMath::Sq(sum_cdf) */ 
              /* + rel_err_all_sq */ 
              /* - 2*err_cdf_sq/sum_cdf/sum_all); */
        /* if (in_val < 0) cout << " FOUND IT: " << in_val << "  " << bin_err << endl; */
        hg_cdf->SetBinContent(i_to,sum_cdf/sum_all);
        hg_cdf->SetBinError(i_to, bin_err);
        ++i_to;
    }
    return hg_cdf;
};

tuOptMap tuCalcRowStats(TH1* hg, double q0, double q1, double nsig, bool b_cut) {
    // return a tuOptMap like:
    // q0 = q0
    // q1 = q1
    // x0 = quantile location q0
    // x1 = quantile location q1
    // i0 = bin x0
    // i1 = bin i1
    // mean = mean from i0 to i1
    // var  = var  from i0 to i1
    // cut  = mean + nsig*var
    tuOptMap dict {{ "q0", q0, "q1", q1, "nsig", nsig }};
    if (hg->Integral()==0) {
        dict("nothing_here") = 1.;
        return dict;
    };
    double q[2]; q[0] = q0; q[1]=q1;
    double x[2];

    hg->GetQuantiles(2,x,q);
    TAxis *ax = hg->GetXaxis();
    int nbins = ax->GetNbins();
    double i0 = ax->FindBin(x[0]);
    double i1 = ax->FindBin(x[1]);
    ax->SetRange(i0,i1);
    double mean = hg->GetMean();
    double var  = hg->GetStdDev();
    double x_cut = mean + nsig*var;
    int    i_cut = ax->FindBin(x_cut)+1;
    if (i_cut>(nbins+1)) i_cut = nbins+1;
    if (b_cut) dict("n_cut") = tuScrubBlock(hg, i_cut, nbins);
    tuInflate(hg);
    
    dict += {{ "x0", x[0], "x1", x[1], "i0", i0, "i1", i1, "mean", mean, "var", var, "x_cut", x_cut,
        "n_cut", dict("n_cut",0), "var", var, "i_cut", i_cut }};
    return dict;
};
double tuScrubBlock(TH2* hg, int x0, int x1, int y0, int y1) {
    if (x0 < 1) x0 = 1;
    if (x1 > hg->GetXaxis()->GetNbins()) x1 = hg->GetXaxis()->GetNbins();
    if (y0 < 1) y0 = 1;
    if (y1 > hg->GetYaxis()->GetNbins()) y1 = hg->GetYaxis()->GetNbins();
    if ( (x1<x0) || (y1<y0) ) return 0.;
    double start = hg->Integral(x0,x1,y0,y1);
    if (start ==0) return 0.;
    for (int x=x0;x<=x1;++x) {
        for (int y=y0;y<=y1;++y) {
            if (hg->GetBinContent(x,y) != 0) {
                hg->SetBinContent(x,y,0.);
                hg->SetBinError(x,y,0.);
            }
        }
    }
    return start - hg->Integral();
};
double tuScrubBlock(TH1* hg, int x0, int x1){
    if (x0 < 1) x0 = 1;
    if (x1 > hg->GetXaxis()->GetNbins()) x1 = hg->GetXaxis()->GetNbins();
    if (x1<x0) return 0.;
    double start = hg->Integral(x0,x1);
    if (start ==0) return 0.;
    for (int x=x0;x<=x1;++x) {
        if (hg->GetBinContent(x) != 0) {
            hg->SetBinContent(x,0.);
            hg->SetBinError(x,0.);
        }
    }
    return start-hg->Integral(x0,x1);
};

double tuScrubbBins(TH1* hg, int min_entries) {
    double start = hg->Integral();
    if (start == 0) return start;
    for (int i = 0; i<hg->GetNcells(); ++i) {
        if (hg->GetBinContent(i)>0. && hg->GetBinContent(i)<min_entries) {
            hg->SetBinContent(i,0.);
            hg->SetBinError  (i,0.);
        }
    }
    return start-hg->Integral();
};

double tuScrubIslands(TH2* hg,  bool isX, int nblank_max, double max_scrub_rat) {
    double total = hg->Integral();
    if (total==0.) return 0.;
    /* int nbins_Y = hg->GetYaxis()->GetNbins(); */
    /* int nbins_X = hg->GetXaxis()->GetNbins(); */
    int nbins_i = isX ? hg->GetNbinsY() : hg->GetNbinsX();
    int nbins_j = isX ? hg->GetNbinsX() : hg->GetNbinsY();
    int i; // inner loop
    int j; // outer loop
    int one {1};
    int& x = isX ? j : i;
    int& y = isX ? i : j;
    int& x_min = isX ? one : x;
    int& y_min = isX ? y : one;
    int& x_max = isX ? nbins_j : x;
    int& y_max = isX ? y : nbins_j;

    for (i=1;i<nbins_i; ++i) {
        // see if this row is blank
        double sum_strip = hg->Integral(x_min,x_max,y_min,y_max);
        if (sum_strip == 0.) continue;

        int jbin_max = 0;
        double content_max = 0;

        for (j=1;j<=nbins_j;++j) {
            if (hg->GetBinContent(x,y) > content_max) {
                content_max = hg->GetBinContent(x,y);
                jbin_max = j;
            }
        }
        if (jbin_max == 0) continue;
        // remove all islands below jbin_max
        int n_blank=0;
        for (j=jbin_max-1; j>=1; --j) {
            if (hg->GetBinContent(x,y) == 0.) {
                ++n_blank;
                if (n_blank == nblank_max) {
                    // check percentage above the blanks
                    if ( hg->Integral(x_min, x, y_min, y)/sum_strip <= max_scrub_rat) {
                        tuScrubBlock(hg, x_min, x, y_min, y);
                        break;
                    }
                }
            } else {
                n_blank = 0;
            }
        }
        // remove all islands above jbin_max
        j = jbin_max;
        n_blank=0;
        for (j=jbin_max+1; j<=nbins_j; ++j) {
            if (hg->GetBinContent(x,y) == 0.) {
                ++n_blank;
                if (n_blank == nblank_max) {
                    // check percentage above the blanks
                    if ( hg->Integral(x, x_max, y, y_max)/sum_strip <= max_scrub_rat) {
                        tuScrubBlock(hg, x, x_max, y, y_max);
                        break;
                    }
                }
            } else {
                n_blank = 0;
            }
        }
    }
    return total - hg->Integral();
};

void tuInflate(TH1* h) {
    h->GetXaxis()->SetRange(1,h->GetXaxis()->GetNbins());
    h->GetYaxis()->SetRange(1,h->GetYaxis()->GetNbins());
};
