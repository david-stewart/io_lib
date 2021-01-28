#include "io_fnc.h"

const char* ioUniqueName(int i) {
    while (gDirectory->FindObjectAny(Form("unique_name__%i",i))!=nullptr) ++i;
    return Form("unique_name__%i",i);
};

void ioWaitPrimitive(int i)  {
    TCanvas* c = new TCanvas( ioUniqueName(i), "", 500,500);
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

double io_get_box_integral(TH2D* hg, vector<double>p, vector<double>q){
    int x0 { hg->GetXaxis()->FindBin(p[0])};
    int x1 { hg->GetXaxis()->FindBin(q[0])};
    int y0 { hg->GetXaxis()->FindBin(p[1])};
    int y1 { hg->GetXaxis()->FindBin(q[1])};
    return hg->Integral(x0,x1,y0,y1);
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

TH1D* ioDivideTH1(TH1* num, TH1* den, bool norm) {
    double norm_num {1};
    double norm_den {1};

    if (norm) {
        norm_num = num->Integral();
        norm_den = den->Integral();
    };

    TH1D* ret = (TH1D*) num->Clone(ioUniqueName());

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
        case 15: return 6;
    }
    return -1;
};
string io_geant05_ascii(int geantid) {
    switch (geantid) {
        case 8: return "pi";
        case 9: return "antipi";
        case 11: return "K";
        case 12: return "antiK";
        case 14: return "P";
        case 15: return "pbar";

        case 0: return "pi";
        case 1: return "antipi";
        case 2: return "K";
        case 3: return "antiK";
        case 4: return "P";
        case 5: return "pbar";
    }
    return "none";
};
string io_geant05_greek(int geantid) {
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
