#include "io_fnc.h"
#include "ioClass.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "TProfile2D.h"

#include <numeric>
#include <iostream>
#include <fstream>



const char* ioUniqueName(int i) {
    while (gDirectory->FindObjectAny(Form("unique_name__%i",i))!=nullptr) ++i;
    return Form("unique_name__%i",i);
};

void ioWaitPrimitive(int i)  {
    TCanvas* c = new TCanvas( ioUniqueName(i+100), ioUniqueName(i+100), 500,500);
    c->WaitPrimitive();
};

vector<int> ioColorVec(int n_colors, int palette, bool print) {
    TCanvas *junk_canvas = new TCanvas(ioUniqueName(),"junk_canvas",1200,800);
    gStyle->SetPalette(palette);
    vector<TProfile*> hg;
    vector<int> colors;
    for (int i{0};i<n_colors;++i){
        if (print) cout << " i: " << i << endl;
        hg.push_back(new TProfile(ioUniqueName(i),";;",n_colors,0.,n_colors) );
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
    if (print) for (auto i : colors)  cout << " color: " << i << endl;
    delete junk_canvas;
    return colors;
};


void ioDrawTLatex(const char* msg, double x, double y, ioOptMap options, ioOptMap dict
) {
    dict += options;
    TLatex tlatex;
    tlatex.SetTextColorAlpha(dict["TextColor"],dict["TextColorAlpha"]);
    tlatex.SetTextAlign(dict["TextAlign"]);
    tlatex.SetTextSize (dict["TextSize"]);
    tlatex.SetTextFont (dict["TextFont"]);
    tlatex.SetTextAngle (dict["TextAngle"].val());
    tlatex.DrawLatex(x, y, msg);
};

void ioDrawTPaveText(double x0, double y0, double x1, double y1,
        vector<const char*> msg, ioOptMap options, ioOptMap dict)
{
    dict += options;
    TPaveText *pav = new TPaveText(x0,y0,x1,y1, dict["BorderOpt"]);
    for (auto m : msg)  pav->AddText(m);

    pav->SetTextColorAlpha(dict["TextColor"],dict["TextColorAlpha"].val());
    pav->SetFillColorAlpha(dict["FillColor"],dict["FillColorAlpha"].val());
    pav->SetTextAlign(dict["TextAlign"]);
    pav->SetTextSize (dict["TextSize"].val());
    pav->SetTextFont (dict["TextFont"]);
    pav->SetTextAngle(dict["TextAngle"].val());
    pav->SetLineWidth(dict["LineWidth"]);
    pav->SetFillStyle(dict["FillStyle"]);

    pav->Draw();
};

void ioDrawTLine(double x0, double y0, double x1, double y1, 
        ioOptMap options, ioOptMap dict)
{
    dict += options;
    TLine* line = new TLine(x0,y0,x1,y1);
    io_fmt(line, options);
    /* cout << options << endl; */
    line->Draw("same");
};

void ioDrawTLineBox(double x0, double y0, double x1, double y1, 
        ioOptMap options) {
    ioDrawTLine( x0, y0, x1, y0, options);
    ioDrawTLine( x1, y0, x1, y1, options);
    ioDrawTLine( x1, y1, x0, y1, options);
    ioDrawTLine( x0, y1, x0, y0, options);
};

double ioPadxRat(double x_in){
    gPad->Update();
    /* cout << " x_in " << x_in << endl; */
    double x0 = gPad->GetUxmin();
    double x1 = gPad->GetUxmax();
    /* cout << " x0 " << x0 << endl; */
    /* cout << " x1 " << x1 << endl; */
    return x0+x_in*(x1-x0);
};
double ioPadyRat(double y_in){
    gPad->Update();
    double y0 = gPad->GetUymin();
    double y1 = gPad->GetUymax();
    return y0+y_in*(y1-y0);
};

void ioDrawTLineHorizontal(double y, ioOptMap options) {
    gPad->Update();
    double x0 { ioPadxRat(0.) };
    double x1 { ioPadxRat(1.) };
    ioDrawTLine(x0,y,x1,y,options);
};

void ioDrawTLineVertical(double x, ioOptMap options) {
    gPad->Update();
    double y0 { ioPadyRat(0.) };
    double y1 { ioPadyRat(1.) };
    ioDrawTLine(x,y0,x,y1,options);
};

double io_get_box_integral(TH2D* hg, pair<double,double>p, pair<double,double>q){
    int x0 { (p.first==0  && q.first==0)  ? 1 : hg->GetXaxis()->FindBin(p.first)};
    int x1 { (p.first==0  && q.first==0)  ? hg->GetXaxis()->GetNbins() : hg->GetXaxis()->FindBin(q.first)};
    int y0 { (p.second==0 && q.second==0) ? 1 : hg->GetYaxis()->FindBin(p.second)};
    int y1 { (p.second==0 && q.second==0) ? hg->GetYaxis()->GetNbins() : hg->GetYaxis()->FindBin(q.second)};
    return hg->Integral(x0,x1,y0,y1);
};
double io_get_box_integral(TProfile2D* pr, pair<double,double>p, pair<double,double>q){
    int x0 { (p.first==0  && q.first==0)  ? 1 : pr->GetXaxis()->FindBin(p.first)};
    int x1 { (p.first==0  && q.first==0)  ? pr->GetXaxis()->GetNbins() : pr->GetXaxis()->FindBin(q.first)};
    int y0 { (p.second==0 && q.second==0) ? 1 : pr->GetYaxis()->FindBin(p.second)};
    int y1 { (p.second==0 && q.second==0) ? pr->GetYaxis()->GetNbins() : pr->GetYaxis()->FindBin(q.second)};
    return pr->Integral(x0,x1,y0,y1);
};
double io_get_box_mean(TProfile2D* pr, pair<double,double>p, pair<double,double>q){
    if (p.first  != 0 || q.first != 0)  pr->GetXaxis()->SetRangeUser(p.first, q.first );
    if (p.second != 0 || q.second != 0) pr->GetYaxis()->SetRangeUser(p.second,q.second);
    return pr->GetMean(3);
};

    
vector<double> io_vec_BinContent(TH1* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinContent(i));
    return vec;
};

vector<double> io_vec_BinError(TH1* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinError(i));
    return vec;
};

/* template<class T> */
vector<double> io_vec_BinCenter(TH1* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    TAxis* ax = hg->GetXaxis();
    for (int i=i0;i<i1;++i) {
        vec.push_back(ax->GetBinCenter(i));
    }
    return vec;
};

/* template<class T> */
pair<int, double*> io_vec_BinEdge(TH1* hg) {
    int n_bins = hg->GetXaxis()->GetNbins();
    double *edges = new double[n_bins+1];
    for (int i{1};i<=n_bins+1;++i){
        edges[i-1] = hg->GetXaxis()->GetBinLowEdge(i);
    }
    return {n_bins,edges};
};

vector<string> io_split_string(string str){
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

TH1D* ioDivideTH1(TH1* num, TH1* den, ioOptMap opt, ioOptMap dict) {
    dict += opt;

    double norm_num {1};
    double norm_den {1};

    if (dict["norm"]()) {
        norm_num = num->Integral();
        norm_den = den->Integral();
    };

    TH1D* ret;
    if (dict["style-den"]==1) {
        ret = (TH1D*) den->Clone(ioUniqueName());
    } else {
        ret = (TH1D*) num->Clone(ioUniqueName());
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

TH2* ioDivideTH2byTH1(TH2* num, TH1* den, bool scale_by_cols) {
    int nbins = den->GetXaxis()->GetNbins();
    if (scale_by_cols) {
        if (num->GetXaxis()->GetNbins() != nbins) {
             cout << " error in ioDivideTH2byTH1: TH2->xAxis->nBins="
             << num->GetXaxis()->GetNbins() << " != TH1->xAxis->nBins="
             << nbins <<"!" <<endl
             << "Division not taking place."<<endl;
             return num;
        }
    } else {
        if (num->GetYaxis()->GetNbins() != nbins) {
            cout << " error in ioDivideTH2byTH1: TH2->xYxis->nBins="
             << num->GetYaxis()->GetNbins() << " != TH1->xAxis->nBins="
             << nbins <<"!" <<endl
             << "Division not taking place."<<endl;
            return num;
        }
    }

    int n_per_loop = scale_by_cols ? num->GetYaxis()->GetNbins() : num->GetXaxis()->GetNbins();
    for (int n=0;n<=nbins+1;++n) {
        double d = den->GetBinContent(n);
        /* cout << " n: "<<n<<"  -> " << den->GetBinContent(n) << endl; */
        for (int j=0;j<=n_per_loop+1;++j) {
            int i;
            if (n_per_loop) {
                i=n;
            } else {
                i=j;
                j=n;
            }
            double n = num->GetBinContent(i,j);
            if (n == 0 || d == 0) {
                num->SetBinContent(i,j,0);
                num->SetBinError  (i,j,0);
                continue;
            } 
            double n_err = num->GetBinError(i,j);
            double d_err = den->GetBinError(n)  ;
            double val = n / d;
            double err = val * pow( pow(n_err/n,2)+pow(d_err/d,2),0.5);
            num->SetBinContent(i,j,val);
            num->SetBinError(i,j,err);
        }
    }
    return num;
};

void io_scaleByBinWidth(TH1D* hg, double scale_factor) {
    for (int i{1}; i<=(int) hg->GetNbinsX(); ++i) {
        if (hg->GetBinContent(i)) {
            double factor { scale_factor / hg->GetBinWidth(i) };
            hg->SetBinContent(i, hg->GetBinContent(i) * factor);
            hg->SetBinError(i, hg->GetBinError(i) * factor);
        }
    }
};
void io_scaleByBinWidth(TH2D* hg, double scale_factor, bool byXwidth, bool byYwidth) {
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


vector<double> io_vecBinContent(TH1D* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinContent(i));
    return vec;
};
vector<double> io_vecBinError  (TH1D* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinError(i));
    return vec;
};
vector<double> io_vecBinContent(TProfile* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinContent(i));
    return vec;
};
vector<double> io_vecBinError  (TProfile* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinError(i));
    return vec;
};
vector<double> io_vecBinEntries(TProfile* hg, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? hg->GetNcells() : hg->GetNcells()-1;
    for (int i=i0;i<i1;++i) vec.push_back(hg->GetBinEntries(i));
    return vec;
};
vector<double> io_vecAxisBinCenter (TAxis* axis, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? axis->GetNbins()+1 : axis->GetNbins();
    for (int i=i0;i<=i1;++i) {
        vec.push_back(axis->GetBinCenter(i));
    }
    return vec;
};
vector<double> io_vecAxisBinEdges  (TAxis* axis, bool under_over_flow){
    vector<double> vec;
    int i0 = under_over_flow ? 0 : 1;
    int i1 = under_over_flow ? axis->GetNbins()+2: axis->GetNbins()+1;
    for (int i=i0;i<=i1;++i) vec.push_back(axis->GetBinLowEdge(i));
    return vec;
};

int io_geant05(int geantid) {
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
const char* io_geant05_ascii(int geantid) {
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
const char* io_geant05_greek(int geantid) {
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


TF1*  io_TsallisFit(double m0, double A, double T, double n, double x_min, double x_max) {
    TF1* fun = new TF1( ioUniqueName(), "[0] / TMath::Power((1.+"
            "(TMath::Sqrt(x*x+[1]*[1])-[1])/([2]*[3])), [3])", x_min, x_max);
    fun->SetParameter(0, A);
    fun->SetParameter(1, m0);
    fun->SetParameter(2, T);
    fun->SetParameter(3, n);
    return fun;
};
TF1* io_dAu_200GeV_TsallisFit(const char* name, double x_min, double x_max) {
    if ("pi") return io_TsallisFit( 0.135, 13.62511532585, 
            0.1505198389696, 10.0949606259, x_min, x_max);
    else if (!strcmp(name,"antipi")) return io_TsallisFit( 0.135, 13.62511532585, 
                0.1505198389696, 10.0949606259, x_min, x_max);
    else if (!strcmp(name,"p")) return io_TsallisFit( 0.938272, 0.2873036847272, 
            0.2114254602693, 11.44109908352, x_min, x_max);
    else if (!strcmp(name,"pbar")) return io_TsallisFit( 0.938272, 0.2257531137243, 
            0.2195209346999, 13.02603055528, x_min, x_max);
    else if (!strcmp(name,"K")) return io_TsallisFit( 0.493677, 0.461550522549, 
            0.2215729988205, 11.61375838522, x_min, x_max);
    else if (!strcmp(name,"antiK")) return io_TsallisFit( 0.493677, 0.4675391506256, 
            0.2192846302987, 11.44758460939, x_min, x_max);
    else {
        cout << "fatal error: particle name \"" << name << "\" not valid." << endl;
        throw std::runtime_error("Particle name must be pi, antipi, p, pbar, K, antiK");
        return nullptr;
    }
};
TF1* io_pp_200GeV_TsallisFit(const char* name, double x_min, double x_max) {
    if (!strcmp(name,"pi")) return io_TsallisFit( 0.135, 5.503706201866,
        0.1277746601471, 9.759524221237, x_min, x_max);
    else if (!strcmp(name,"antipi")) return io_TsallisFit( 0.135, 5.58650216135,
        0.1270067076772, 9.73258646095,  x_min, x_max);
    else if (!strcmp(name,"p")) return io_TsallisFit( 0.938272, 0.07499819855665,
        0.1758814595453, 10.54524945929, x_min, x_max);
    else if (!strcmp(name,"pbar")) return io_TsallisFit( 0.938272, 0.06427254263927,
        0.1691919676508, 10.06359691238, x_min, x_max);
    else if (!strcmp(name,"K")) return io_TsallisFit( 0.493677, 0.06934945866744,
        0.1929677543032, 11.82478291302, x_min, x_max);
    else if (!strcmp(name,"antiK")) return io_TsallisFit( 0.493677, 0.06934945866744,
        0.1929677543032, 11.82478291302, x_min, x_max);
    else {
        cout << "fatal error: particle name \"" << name << "\" not valid." << endl;
        throw std::runtime_error("Particle name must be pi, antipi, p, pbar, K, antiK");
        return nullptr;
    }
};
void io_apply_prior(TF1* fn, TH1D* T) {
    // respone: x-axis=Meas. y-axis=True
    // Note: x-axis underflow may contain all the misses
    /* TH1D* T = response->ProjectionY(); */
    TAxis* Y_ax = T->GetXaxis();
    /* TAxis* X_ax = MT->GetXaxis(); */
    for (int y=0; y<Y_ax->GetNbins()+2; ++y){
        if (T->GetBinContent(y) == 0) continue;
        double w = fn->Integral(Y_ax->GetBinLowEdge(y), Y_ax->GetBinUpEdge(y)) /
                T->GetBinContent(y);
            double val = T->GetBinContent(y)*w;
            double err = T->GetBinError(y)*w;
            T->SetBinContent(y,val);
            T->SetBinError(y,err);
    };
};
void io_apply_prior(TF1* fn, TH2D* MT, TH1D* T, bool both) {
    // Note: this only applies the prior to MT, NOT to T
    // respone: x-axis=Meas. y-axis=True
    // NB: assumed that the x-axis underflow bin contains all the misses
    /* TH1D* T = response->ProjectionY(); */
    TAxis* Y_ax = MT->GetYaxis();
    TAxis* X_ax = MT->GetXaxis();

    for (int y=0; y<Y_ax->GetNbins()+2; ++y){
        if (T->GetBinContent(y) == 0) continue;
        double w = fn->Integral(Y_ax->GetBinLowEdge(y), Y_ax->GetBinUpEdge(y)) /
                T->GetBinContent(y);
        for (int x = 0; x<X_ax->GetNbins()+2; ++x) {
            double val = MT->GetBinContent(x,y) * w;
            double err = MT->GetBinError(x,y) * w;
            MT->SetBinContent(x,y,val);
            MT->SetBinError(x,y,err);
        }
        if (both) {
            double val = T->GetBinContent(y)*w;
            double err = T->GetBinError(y)*w;
            T->SetBinContent(y,val);
            T->SetBinError(y,err);
        }
    };
};

TH1D* io_BayesUnfold(TH1D* data, RooUnfoldResponse* response, int iRepUnfold) {
    RooUnfoldBayes*    bayes     = new RooUnfoldBayes(response, data, iRepUnfold);
    TH1D* unfolded = (TH1D*) bayes->Hreco();
    unfolded->SetName(ioUniqueName());
    return unfolded;
};

TH1D* io_BayesUnfold(TH1D* data, TH1D* T, TH2D* R, int iRepUnfold, TH1D* M) {
    // return data unfolded with response matrix R
    // note that truth T is included for misses, and if present
    // M will include fakes
    if (M==nullptr) M = R->ProjectionX(ioUniqueName());
    RooUnfoldResponse* rooUnfRes = new RooUnfoldResponse (M,T,R,ioUniqueName());
    RooUnfoldBayes*    bayes     = new RooUnfoldBayes(rooUnfRes, data, iRepUnfold);
    TH1D* unfolded = (TH1D*) bayes->Hreco();
    unfolded->SetName(ioUniqueName());
    delete rooUnfRes;
    delete bayes;
    return unfolded;
};

TLegend* ioNewTLegend() {
    TLegend *leg = new TLegend(0.6455533,0.6332006,0.8911167,0.8938613,NULL,"brNDC");
    io_fmt(leg);
    return leg;
};

float io_dphi(float phi0, float phi1) {
    float rval = phi1-phi0;
    while (rval < -IO_pi) rval += IO_twopi;
    while (rval >  IO_pi) rval -= IO_twopi;
    return rval;
};

bool io_isAbsTransPhi(float phi0, float phi1, float lo_bound, float hi_bound){
    float dphi = TMath::Abs(io_dphi(phi0,phi1));
    return (dphi>=lo_bound && dphi<=hi_bound);
};

float io_02pi(float &phi){
    while (phi<0)       phi += IO_twopi;
    while (phi>IO_twopi) phi -= IO_twopi;
    return phi;
};
float io_02pi(float  phi){
    while (phi<0)        phi += IO_twopi;
    while (phi>IO_twopi) phi -= IO_twopi;
    return phi;
};

vector<double> ioQuantiles(TH1D* hg, vector<double> percents){
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
string  ioStringVec(vector<double> vec, const char* name, const char* formatter){
    if (vec.size()==0) return "";
    ostringstream os;
    const char* fmt = Form("%%%s",formatter);
    os << "vector<double> " << name << " {" << Form(fmt,vec[0]);
    for (int i{1}; i<(int)vec.size();++i) os << ", " << Form(fmt,vec[i]);
    os << "};";
    return os.str();
}

int io_count_digits(int n, int min_val) {
    int cnt = 0;
    while (n != 0) { 
        n /= 10;
        ++cnt;
    }
    if (cnt < min_val) cnt = min_val;
    return cnt;
};
int ioSum(const vector<int> vec)  { return std::accumulate(vec.begin(), vec.end(), 0); };
int ioSum(const vector<bool> vec) { return std::accumulate(vec.begin(), vec.end(), 0); };

map<int,string>    ioReadIntStrMap(const char* file, ioOptMap options, ioOptMap dict){
    dict += options;
    /* cout << " map =1 : " << options["tag"].str() << endl; */
    auto data = ioReadStrVec(file, options);
    if (data.size() % 2) 
        throw std::runtime_error(
                "Error in io_ReadIntStrMap : must have even number of entries");
    map<int,string> M;
    for (unsigned int i=0; i<(data.size()/2); ++i) {
        int index = atoi(data[i*2].c_str());
        M[index] = data[i*2+1];
    }
    return M;
};

vector<int> ioReadIntVec(const char* in_file, int col, bool sort, bool strip_commas) {
    vector<int> vec;
    ifstream file;
    file.open(in_file);
    ostringstream msg;
    if (!file.is_open()) {
        msg << "fatal error: could not open int map file \"" 
            << in_file << " in ioReadIntVec" << endl;
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
            << "  in input file \"" << in_file <<"\" in ioReadIntVec" << endl;
        throw std::runtime_error(msg.str());
    }
    file.close();
    if (sort) std::sort(vec.begin(), vec.end());
    return vec;
};

//----------
map<string,string> io_VecStrToMapStrStr (vector<string> data) 
{
    if (data.size() % 2) 
        throw std::runtime_error(
                "Error in io_VecStrToMapStrStrg : must have even number of entries");
    map<string,string> M;
    for (unsigned int i=0; i<(data.size()/2); ++i) {
        M[data[i*2]] = data[i*2+1];
    }
    return M;
};
map<string,string> ioReadMapStrStr(const char* in_file, const char* tag, 
        ioOptMap options, ioOptMap dict) 
{
    dict += options;
    dict["tag"] = tag;
    return io_VecStrToMapStrStr(ioReadStrVec(in_file, dict));
};
map<string,string> ioReadMapStrStr(const char* in_file, ioOptMap options, ioOptMap dict) {
    dict += options;
    return io_VecStrToMapStrStr(ioReadStrVec(in_file, dict));
};


vector<string> ioReadStrVec(const char* in_file, const char* tag, 
        ioOptMap options, ioOptMap dict) {
    dict += options;
    dict["tag"] = tag;
    return ioReadStrVec(in_file, dict);
};
vector<string> ioReadStrVec(const char* in_file, ioOptMap options, ioOptMap dict) {
    dict += options;

    bool use_column = (dict["column"].str()!="all");
    int  which_column {0};
    if  (use_column) { which_column = dict["column"]; }
    /* cout << "a0 : " << dict["tag"].str() << endl; */
    bool use_tag {dict["tag"].str() != "none"};
    string tag = use_tag ? dict["tag"] : "";
    bool in_tag { !use_tag };
    bool strip_commas { (bool) dict["strip_commas"]() };

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
            << in_file << " in ioReadStrVec" << endl;
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
                    if (ioWordIsEndTag(word, tag)) { 
                        found_end_tag  = true;
                        break;
                    }
                } else {
                    if (ioWordIsTag(word,tag)) {
                        found_start_tag = true;
                        in_tag = true;
                    }
                    continue;
                }
            }
            if (ioIsAnyTag(word)) continue;
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
             
    if ( (bool) dict["sort"]() ) std::sort(vec.begin(), vec.end());
    return vec;
};




vector<double> ioReadValVec(const char* in_file, const char* tag, 
        ioOptMap options, ioOptMap dict) {
    dict += options;
    dict["tag"] = tag;
    return ioReadValVec(in_file, dict);
};
vector<double> ioReadValVec(const char* in_file, ioOptMap options, ioOptMap dict) {
    dict += options;

    bool use_column = (dict["column"].str()!="all");
    int  which_column {0};
    if  (use_column) { which_column = dict["column"]; }
    bool use_tag {dict["tag"].str() != "none"};
    string tag = use_tag ? dict["tag"] : "";
    bool in_tag { !use_tag };
    bool strip_commas { (bool) dict["strip_commas"]() };

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
            << in_file << " in ioReadIntVec" << endl;
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
            bool is_a_tag { ioIsAnyTag(word) };
            if (!is_a_tag) break; // neither number or tag, so is a comment
            // here: is_a_tag == true
            if (!use_tag) continue;
            // use_tag == true && this is a tag
            if ( in_tag ) {
                if (ioWordIsEndTag(word, tag)) { 
                    found_end_tag  = true;
                    break;
                }
            } else {
                if (ioWordIsTag(word,tag)) {
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
             
    if ( (bool) dict["sort"]() ) std::sort(vec.begin(), vec.end());
    return vec;
};

pair<int,double*> ioReadValsPtr(const char* file, ioOptMap options, ioOptMap dict) {
    dict += options;
    auto vec = ioReadValVec(file, dict);
    int i0 = dict["begin_index"];
    int i1 = dict["end_index"];
    int size_data {(int) vec.size()};
    int size_req = (i1 != -1) 
        ? i1-i0
        : size_data-i0;
    if (i0!=0 || i1!=-1) {
        if (size_req<=0 || size_req > size_data) {
            ostringstream msg;
            msg << "ioReadValsPtr in file " << file << " requested start and end indices " 
                << i0 << " and " << i1
                << " for total size of " << size_req 
                << " out of total available " << size_data << endl;
            throw std::runtime_error(msg.str());
        }
    }
    if (size_data == 0) {
        throw std::runtime_error(
        (string)"ioReadValsPtr  in file " + file + " found no data. ");
    }
    double* vals = new double[size_req];
    for (int i{0}; i<size_req; ++i) {
        vals[i] = vec[i0+i];
    }
    return {size_req, vals};
};

/* TGraph* ioMakeTGraph(vector<double>& x, vector<double>& y) { */
/*     if (x.size() != y.size()) */ 
/*         throw std::runtime_error("ioMakeTGraph(vec, vec) required vectors of same length"); */
/*     const unsigned int n = x.size(); */
/*     double* xpts = new double[n]; */
/*     double* ypts = new double[n]; */
/*     for (unsigned int i{0}; i<n; ++i) { */
/*         xpts[i] = x[i]; */
/*         ypts[i] = y[i]; */
/*     } */
/*     return new TGraph(n,xpts,ypts); */
/* }; */
TGraphErrors* ioMakeTGraphErrors(vector<double> x, vector<double> y, vector<double> y_err, vector<double> x_err, vector<bool> use) {
    ioMinMax nsize;
    nsize.fill(x.size());
    nsize.fill(y.size());
    nsize.fill(y_err.size());
    if (x_err.size() != 0) nsize.fill(x_err.size());
    if (use.size() != 0)   nsize.fill(use.size());

    if (x.size() != y.size() || x.size() != y_err.size()) 
        cout << " Warning: ioMakeTGraphErrors(vec, vec, vec={}) don't have equal input sizes" << endl
             << " Therefor using shortest vector of size " <<  nsize.max << endl;
    if (x_err.size() != 0 && x_err.size() != x.size())
        cout << " Warning: ioMakeTGraphErrors(vec, vec, vec={}) don't have equal input sizes" << endl
             << " Therefor using shortest vector of size " <<  nsize.max << endl;
    if (use.size() != 0 && use.size() != x.size())
        cout << " Warning: ioMakeTGraphErrors(vec, vec, vec={}) don't have equal input sizes" << endl
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

TGraph* ioMakeTGraph(vector<double> x, vector<double> y) {
    if (x.size() != y.size()) 
        throw std::runtime_error("ioMakeTGraph(vec, vec) required vectors of same length");
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
int iowhichbin0(double val, vector<double>& vec) {
    return (int)(std::lower_bound(vec.begin(), vec.end(), val) - vec.begin());
} ;
int iowhichbin0(double val, TH1D* hg) {
    return (int)(hg->GetXaxis()->FindBin(val)-1);
} ;
int iowhichbin1(double val, vector<double>& vec) {
    return (int)(std::lower_bound(vec.begin(), vec.end(), val) - vec.begin())+1;
} ;
int iowhichbin1(double val, TH1D* hg) {
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
double io_D(double x0,double y0,double x1,double y1) 
{ return TMath::Sqrt( TMath::Sq(x1-x0)+TMath::Sq(y1-y0)); };
double io_D2(double x0,double y0,double x1,double y1) 
{ return TMath::Sq(x1-x0)+TMath::Sq(y1-y0); };
double io_R(double x0,double y0,double x1,double y1) 
{ return TMath::Sqrt( TMath::Sq(x1-x0)+TMath::Sq(io_dphi(y1,y0))); };
double io_R2(double x0,double y0,double x1,double y1) 
{ return TMath::Sq(x1-x0)+TMath::Sq(io_dphi(y1,y0)); };



RooUnfoldResponse ioMakeRooUnfoldResponse(
     int nbins, double lo_bin, double hi_bin, 
     const char* tag, const char* title
) {
    const char* this_tag {
        (strcmp(tag,"")==0)
            ? ioUniqueName()
            : tag 
    };
    TH1D truth {
            Form("%s_truth",this_tag),
            Form("%s;truth;N",title),
            nbins, lo_bin, hi_bin };
    TH1D measured {
            Form("%s_measured",this_tag),
            Form("%s;measured;N",title),
            nbins, lo_bin, hi_bin };
    return {&measured, &truth, Form("%s_RooUnfR",this_tag), title};
};

RooUnfoldResponse ioMakeRooUnfoldResponse(
     int nbins, double* edges, 
     const char* tag, const char* title
) {
    const char* this_tag {
        (strcmp(tag,"")==0)
            ? ioUniqueName()
            : tag 
    };
    TH1D truth {
            Form("%s_truth",this_tag),
            Form("%s;truth;N",title),
            nbins, edges };
    TH1D measured {
            Form("%s_measured",this_tag),
            Form("%s;measured;N",title),
            nbins,  edges };
    return {&measured, &truth, Form("%s_RooUnfR",this_tag), title};
};


RooUnfoldResponse ioMakeRooUnfoldResponse(
     int nb_measured, double lo_measured, double hi_measured,
     int nb_truth, double lo_truth, double hi_truth,
     const char* tag, const char* title
) {
    const char* this_tag {
        (strcmp(tag,"")==0)
            ? ioUniqueName()
            : tag 
    };
    TH1D truth{
            Form("%s_truth",this_tag),
            Form("%s;truth;N",title),
            nb_truth, lo_truth, hi_truth };
    TH1D measured {
            Form("%s_measured",this_tag),
            Form("%s;measured;N",title),
            nb_measured, lo_measured, hi_measured };
    return {&measured, &truth, Form("%s_RooUnfR",this_tag), title};
};
RooUnfoldResponse ioMakeRooUnfoldResponse(
     int nb_measured, double* edges_measured,
     int nb_truth, double* edges_truth,
     const char* tag, const char* title
) {
    const char* this_tag {
        (strcmp(tag,"")==0)
            ? ioUniqueName()
            : tag 
    };
    TH1D truth {
            Form("%s_truth",this_tag),
            Form("%s;truth;N",title),
            nb_truth, edges_truth};
    TH1D measured {
            Form("%s_measured",this_tag),
            Form("%s;measured;N",title),
            nb_measured, edges_measured};
    return {&measured, &truth, Form("%s_RooUnfR",this_tag), title};
};
// return which bin (starting from 0) the data is in: lower bound <= val < upper bound
/* int iowhichbin(int, double*); */
/* int iowhichbin(vector<double>, vector<int> remap); // return which bin (starting from 0) the data is in: lower bound <= val < upper bound */
/* int iowhichbin(int, double*,   vector<int> remap); */

double ioRatCircleOverLine (double R, double d) {
    R = TMath::Abs(R);
    d = TMath::Abs(d);
    if (R <= d) return 0.;
    return TMath::ACos(d/R)/IO_pi;
};
double ioRatCircleInTwoParallelLines (const double d0,const double d1,double C,double R) {
    if (C<d0) C = d0;
    if (C>d1) C = d1;
    return 1 - ioRatCircleOverLine(R,C-d0)
             - ioRatCircleOverLine(R,d1-C);
};
double ioPolyP6_a0_a1x_a2xx_a3y_a4yy_a5xy(double* x, double *p){
    return p[0]      + p[1]*x[0] + p[2]*x[0]*x[0]
                     + p[3]*x[1] + p[4]*x[1]*x[1] + p[5]*x[0]*x[1];
};

bool ioWordIsTag   (string  word, string tag)  { return (word == ("<" +tag+">")); };
bool ioWordIsEndTag(string  word, string tag)  { return (word == ("</" +tag+">")); };
bool ioWordIsTag   (TString word, string tag) { return (word == ("<" +tag+">")); };
bool ioWordIsEndTag(TString word, string tag) { return (word == ("</" +tag+">")); };
bool ioIsAnyTag    (string  word) { return ioIsAnyTag((TString)word); };
bool ioIsAnyTag    (TString word) {
    if (!word.BeginsWith("<")) return false;
    return word.EndsWith(">");
};

void io_normByRow(TH2D* hg, double factor, bool use_max_val) {
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

vector<double> io_print_first_blank(TH2D* hg) {
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


pair<TF1*, ioOptMap> ioFitJESJER(TH1D* hg, double pt_jet, 
        double quant_lo, double quant_hi, 
        const char* tag, double sigma_rat_guess) 
{
    string my_tag = (strcmp(tag,"")==0) ? ioUniqueName() : tag;
    double center = hg->GetXaxis()->GetBinCenter(hg->GetMaximumBin());
    auto quant = ioQuantiles(hg,{quant_lo,quant_hi});
    if (quant[0] >= center || quant[1] <= center) {
        cout << "Can't fit in ioFitJESJER : the max value bin at "
             << Form("%f doesn't lie between quantiles (%f,%f) at (%f,%f)",
                     center, quant_lo, quant_hi, quant[0], quant[1]) << endl;
        return {nullptr, {{
            "a0",0.,
            "a1",0.,
            "a2",0.,
            "pt_jet", pt_jet,
            "JES", 0.,
            "JES_err",0.,
            "JER", 0.,
            "JER_err", 0.,
            "bound_lo", 0.,
            "bound_hi", 0.,
            "quant_lo", 0.,
            "quant_hi", 0.}}
        };
    }
    TF1* f1 = new TF1(my_tag.c_str(),"gaus(0)",quant[0],quant[1]);
    f1->SetParameter(0, hg->GetMaximum());
    f1->SetParameter(1, center);
    f1->SetParameter(2, center*0.1);
    hg->Fit(my_tag.c_str(),"","",quant[0],quant[1]);
    return {f1, {{
            "a0",f1->GetParameter(0),
            "a1",f1->GetParameter(1),
            "a2",f1->GetParameter(2),
            "pt_jet", pt_jet,
            "JES", (f1->GetParameter(1) - pt_jet),
            "JES_err", (f1->GetParError(1)),
            "JER", (f1->GetParameter(2) ),
            "JER_err", (f1->GetParError(2)),
            "bound_lo", quant[0],
            "bound_hi", quant[1],
            "quant_lo", quant_lo,
            "quant_hi", quant_hi}}};
};

/* void ioTrimSmallBins(TH2D* hg, int Nmin, bool cut_underover_flow=true) { */
/*     TAxis* Xaxis = hg->GetXaxis(); */
/*     TAxis* Yaxis = hg->GetYaxis(); */
/*     int nXbins = Xaxis->GetNbins(); */
/*     int nYbins = Xaxis->GetNbins(); */
/*     if (cut_underover_flow) { */
        
IOS_keepcut_stats io_keepcutbin(TH1* h, int ibin, double minval, double zero_val, double zero_err) {
    IOS_keepcut_stats stats{};
    double val = h->GetBinContent(ibin);
    if (val == 0 && h->GetBinError(ibin)==0) return {};
    if (val < minval) {
        ++stats.n_cut;
        stats.sum_cut = val-zero_val;
        stats.was_cut = true;
        h->SetBinContent(ibin, zero_val);

        if (zero_val != 0 && zero_err == 0.)  h->SetBinError(ibin, sqrt(zero_val)); 
        else  h->SetBinError(ibin, 0.); 

        int x, y, z;
        h->GetBinXYZ(ibin,x,y,z);
        stats.x = h->GetXaxis()->GetBinCenter(x);
        stats.y = h->GetYaxis()->GetBinCenter(y);
    } else {
        ++stats.n_keep;
        stats.sum_keep = val;
    }
    /* if (stats.was_cut) cout << " - " << stats.x << " " << stats.y << endl; */
return stats;
};

IOS_keepcut_stats io_cullsmallbins(
        TH1* h, 
        double min_val, 
        IO area, 
        bool remove_overunderflow) 
{
    if (remove_overunderflow) 
        for (auto i : io_binvec(h, IO::overunderflow)) io_setbinzero(h,i);

    IOS_keepcut_stats stats;
    for (auto i : io_binvec(h, area)) {
        stats += io_keepcutbin(h, i, min_val);
        /* if (stats.was_cut) cout << stats.x << " " << stats.y << endl; */
        /* cout << " stats: " << stats.was_cut << endl; */
    }
    return stats;
};

vector<int> io_binvec(TH1* h, IO loc) {
    vector<int> vec{};
    if (h->GetYaxis()->GetNbins() == 1) {
        switch (loc) {
            case IO::overunderflow:
                vec = { 0, h->GetXaxis()->GetNbins()+1 };
                break;
            case IO::all:
                vec = { 0, h->GetXaxis()->GetNbins()+1 };
                // no break -- continue to next block
            case IO::in:
                for (int i{1}; i<=h->GetXaxis()->GetNbins(); ++i) 
                    vec.push_back(i);
                break;
            default: break;
        }
    } else {
        int nX = h->GetXaxis()->GetNbins();
        int nY = h->GetYaxis()->GetNbins();
        switch (loc) {
            case IO::overunderflow: 
                for (auto x{0}; x<=nX+1; ++x) vec.push_back(h->GetBin(x,0));
                for (auto x{0}; x<=nX+1; ++x) vec.push_back(h->GetBin(x,nY+1));
                for (auto y{1}; y<=nY;   ++y) vec.push_back(h->GetBin(0,y));
                for (auto y{1}; y<=nY;   ++y) vec.push_back(h->GetBin(nX+1,y));
                break;
            case IO::all:
                vec = io_binvec(h, IO::overunderflow);
                // no break -- continue to next block
            case IO::in:
                for (auto y{1}; y<nY+1; ++y)
                for (auto x{1}; x<nX+1; ++x)
                        vec.push_back(h->GetBin(x,y));
                break;
            default: break;
        }
    }
    return vec;
};

double io_setbinzero(TH1* hg, int bin, double val, double err)
{ 
    double p_val = hg->GetBinContent(bin);
    hg->SetBinContent(bin, val);
    if (val != 0. && err == 0.) err = sqrt(val);
    hg->SetBinError(bin, err);
    return p_val;
};


const char* io_cutdiff(int a, int b, const char* fmt) {
    int diff {b-a};
    return Form(fmt,b,diff,(double)diff/a);
};

TH1* ioSetCntErrors(TH1* hg) {
    for (int i{0}; i<hg->GetNcells(); ++i) {
        if (hg->GetBinContent(i) != 0) {
            hg->SetBinError(i, TMath::Sqrt(hg->GetBinContent(i)));
        }
    }
    return hg;
};
TH1* ioAddBinCnt(TH1* hg_to, TH1* hg_fr, bool set_cnt_errors, bool rm_underoverflow) {
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
    if (set_cnt_errors) ioSetCntErrors(hg_to);
    return hg_to;
};

