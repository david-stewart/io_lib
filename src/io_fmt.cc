#include "io_fmt.h"
#include "io_fnc.h"

TLine* io_fmt(TLine* line, ioOptMap options, ioOptMap dict)
{
    dict += options;
    line->SetLineColorAlpha(dict["LineColor"],dict["LineColorAlpha"].val());
    line->SetLineStyle(dict["LineStyle"]);
    line->SetLineWidth(dict["LineWidth"]);
    return line;
};



ioOptMap io_fmt__hg_dict() {
    return {{
        "MarkerStyle", kFullCircle,
        "MarkerColor", kBlack,
        "MarkerAlpha", 1.,
        "MarkerSize",  1.,

        "LineWidth",   1 ,
        "LineStyle",   1 ,
        /* {"LineColor",   kBlack }, */
        /* {"LineAlpha",   1. }, */

        "xAxisTitleFont", 43,
        "xAxisTitleSize", 22,
        "xAxisTitleOffset", 1.6 ,

        "xAxisLabelFont", 43,
        "xAxisLabelSize", 18,
        "xAxisLabelOffset", 0.02 ,

        "yAxisTitleFont", 43,
        "yAxisTitleSize", 22,
        "yAxisTitleOffset", 1.85 ,

        "yAxisLabelFont", 43,
        "yAxisLabelSize", 18,
        "yAxisLabelOffset", 0.02 ,

        "zAxisTitleFont", 43,
        "zAxisTitleSize", 22,
        "zAxisTitleOffset", 2.4 ,

        "zAxisLabelFont", 43,
        "zAxisLabelSize", 18,
        "zAxisLabelOffset", 0.02 ,

        "SetStats", 0
    }};
};


void io_fmt_ranges(vector<TH1D*> hgs, ioOptMap dict) {
    // Used to set common ranges on all hgs (assumed to have common x_range)
    // 1. Find yMinVal:min(all hgs minimima) and yMaxVal:max(all hgs maxima)
    // 2. Options:
    //      "yRangeLowEdge":double Abs value
    //      "yRangeLowRatio":double   = yMinVal-(yMaxVal-yMinVal)*yRangeUpRatio
    //      "yRangeUpEdge":double  Abs value
    //      "yRangeUpRatio":double    = yMinVal+(yMaxVal-yMinVal)*yRangeUpRatio
    //      "yRange":no-option =  Default to 0. to 1.2 * (ymax - min)
    // 3. For xRange:
    //      "xRangeLowEdge":double Abs value
    //      "xRangeLowRatio":double = bin of xMinVal-(xMaxVal-xMinVal)*xRangeUpRatio
    //      "xRangeUpEdge":double  Abs value
    //      "xRangeUpRatio":double    = yMinVal+(yMaxVal-yMinVal)*yRangeUpRatio
    if   (     dict("yRangeLowEdge") 
            || dict("yRangeLowRatio")
            || dict("yRangeUpEdge")
            || dict("yRangeUpRatio")
            || dict("yRange") )
    {
        double yMinVal { hgs[0]->GetBinContent(hgs[0]->GetMinimumBin()) };
        double yMaxVal { hgs[0]->GetBinContent(hgs[0]->GetMaximumBin()) };
        for (unsigned int i{1}; i<hgs.size(); ++i) {
            double minTemp { hgs[i]->GetBinContent(hgs[i]->GetMinimumBin()) };
            double maxTemp { hgs[i]->GetBinContent(hgs[i]->GetMaximumBin()) };
            if (minTemp < yMinVal) yMinVal = minTemp;
            if (maxTemp > yMaxVal) yMaxVal = maxTemp;
        }

        double yMinRange;
        if      (dict("yRangeLowEdge" ) ) yMinRange = dict["yRangeLowEdge"].val();
        else if (dict("yRangeLowRatio") ) yMinRange = 
            yMinVal - (yMaxVal-yMinVal) * dict["yRangeLowRatio"].val();
        else if (yMinVal > 0) yMinRange = 0;
        else                  yMinRange = yMinVal;

        double yMaxRange;
        if      (dict("yRangeUpEdge" ) ) yMaxRange = dict["yRangeUpEdge"].val();
        else if (dict("yRangeUpRatio") ) yMaxRange = 
            yMaxVal + (yMaxVal-yMaxVal) * dict["yRangeUpRatio"].val();
        else yMaxRange = yMaxVal + (yMaxVal-yMinVal) * 0.2;

        /* cout << " beta yrange " << yMinVal << ","<< yMinRange << " " << yMaxVal<<","<<yMaxRange << endl; */
        for (auto hg :hgs) hg->GetYaxis()->SetRangeUser(yMinRange,yMaxRange);
    }

    if   (     dict("xRangeLowEdge") 
            || dict("xRangeLowRatio")
            || dict("xRangeUpEdge")
            || dict("xRangeUpRatio")
            || dict("xRangeCommon")
            || dict("xRange") )
    {
        // For x_Range, find the smallest bin that has non-zero entries
        int nbins { hgs[0]->GetXaxis()->GetNbins() };
        int i0 { nbins };
        int i1 { 1 };

        int i0_common{1};
        int i1_common{nbins};

        for (auto hg : hgs) {
            int i0_s { nbins };
            int i1_s { 1 };
            for (int k{1}; k<=nbins; ++k) {
                if (hg->GetBinContent(k) != 0) { 
                    i0_s = k;
                    break;
                }
            }
            if (i0_s < i0) i0 = i0_s;
            if (i0_s > i0_common) i0_common = i0_s;

            for (int k{nbins}; k>=1; --k) {
                if (hg->GetBinContent(k) != 0) { 
                    i1_s = k;
                    break;
                }
            }
            if (i1_s > i1) i1 = i1_s;
            if (i1_s < i1_common) i1_common = i1_s;
        }

        // if limiting to common range, use that and ignore all else
        if (dict["xRangeCommon"]) {
            for (auto hg : hgs) {
                hg->GetXaxis()->SetRange(i0_common, i1_common);
                /* cout << " cutting: " << hg->GetXaxis()->GetBinLowEdge(i0_common) */ 
                /* << " to " <<        hg->GetXaxis()->GetBinUpEdge(i1_common) << endl; */
            }
            return;
        };

        TAxis* axis = hgs[0]->GetXaxis();

        double xMinVal { axis->GetBinLowEdge(i0) };
        double xMaxVal { axis->GetBinUpEdge( i1) };

        double xMinRange;
        if      (dict("xRangeLowEdge" ) ) xMinRange = dict["xRangeLowEdge"].val();
        else if (dict("xRangeLowRatio") ) xMinRange = 
            xMinVal - (xMaxVal-xMinVal) * dict["xRangeLowRatio"].val();
        else xMinRange = axis->GetBinLowEdge(1);

        int iMinRange = axis->FindBin(xMinRange);
        xMinRange = axis->GetBinLowEdge(iMinRange);

        double xMaxRange;
        if      (dict("xRangeUpEdge" ) ) xMaxRange = dict["xRangeUpEdge"].val();
        else if (dict("xRangeUpRatio") ) xMaxRange = 
            xMaxVal + (xMaxVal-xMaxVal) * dict["xRangeUpRatio"].val();
        else xMaxRange = axis->GetBinUpEdge(nbins);

        int iMaxRange = axis->FindBin(xMaxRange);
        double b_edge = axis->GetBinLowEdge(iMaxRange);
        double width  = axis->GetBinWidth  (iMaxRange);
        double ratio  = (xMaxRange - b_edge) / width;
        if (ratio > 0) xMaxRange = axis->GetBinUpEdge(iMaxRange);
        else           xMaxRange = axis->GetBinLowEdge(iMaxRange);

        for (auto hg :hgs) hg->GetXaxis()->SetRangeUser(xMinRange,xMaxRange);

    }
};
void io_fmt_ranges(vector<TH1D*> hgs, vector<TH1D*> hgs_2, ioOptMap dict) {
    vector<TH1D*> tmp;
    for (auto& hg : hgs)   tmp.push_back(hg);
    for (auto& hg : hgs_2) tmp.push_back(hg);
    io_fmt_ranges(tmp, dict);
}

TCanvas* io_fmt (TCanvas* canv, ioOptMap dict) {
    int FillStyle = dict("TCanvasFillStyle") ?
        dict["TCanvasFillStyle"].get_int() : 4000;

    int FrameFillStyle = dict("TCanvasFrameFillStyle")  ?
        dict["TCanvasFrameFillStyle"].get_int() : 4000;

    canv->SetFillStyle(FillStyle);
    canv->SetFrameFillStyle(FrameFillStyle);
    return canv;
};
TPad* io_fmt (TPad* pad, ioOptMap dict) {
    int FillStyle = dict("TCanvasFillStyle") ?
        dict["TCanvasFillStyle"].get_int() : 4000;

    int FrameFillStyle = dict("TCanvasFrameFillStyle")  ?
        dict["TCanvasFrameFillStyle"].get_int() : 4000;

    pad->SetFillStyle(FillStyle);
    pad->SetFrameFillStyle(FrameFillStyle);
    return pad;
};
TLegend* io_fmt (TLegend* leg, ioOptMap dict) {
    int LineColor =  dict("LineColor") ? dict["LineColor"]   : kWhite;
    int FillColor =  dict("FillColor") ? dict["FillColor"]   : kWhite;
    double alpha  =  dict("alpha")     ? dict["alpha"].val() : 0.;
    leg->SetLineColorAlpha(LineColor, alpha);
    leg->SetFillColorAlpha(FillColor, alpha);
    return leg;
};
TProfile* io_fmt (TProfile* hg, ioOptMap _override, ioOptMap dict) {
    dict += _override;

    /* cout << " format " << endl; */
    /* cout << dict << endl; */

    if (dict("normalize")) hg->Scale(1./hg->Integral());
    if (dict("Rebin"))        hg->Rebin(dict["Rebin"]);

    if (!dict("MarkerAlpha")) dict["MarkerAlpha"] = 1.;
    if (dict("noTitle"))      hg->SetTitle("");
    if (dict("SetStats"))     hg->SetStats(dict["SetStats"]);


    if (dict("MarkerStyle"))
        hg->SetMarkerStyle(dict["MarkerStyle"]);
    if (dict("MarkerColor"))
        hg->SetMarkerColorAlpha(dict["MarkerColor"], dict["MarkerAlpha"].val());
    if (dict("MarkerSize"))
        hg->SetMarkerSize(dict["MarkerSize"].val());

    if (dict("LineWidth"))
        hg->SetLineWidth(dict["LineWidth"]);
    if (dict("LineStyle"))
        hg->SetLineStyle(dict["LineStyle"]);
    if (dict("MarkerColor"))
        hg->SetLineColorAlpha(dict["MarkerColor"], dict["MarkerAlpha"].val());
    
    if (dict("LineColor")) {
        double LineAlpha = dict("LineAlpha") ? dict["LineAlpha"].val() : dict["MarkerAlpha"].val();
        hg->SetLineColorAlpha(dict["LineColor"], LineAlpha);
    }

    if (dict("FillColor")) {
        double FillAlpha = dict("FillAlpha") ? dict["FillAlpha"].val() : dict["MarkerAlpha"].val();
        hg->SetFillColorAlpha(dict["FillColor"], FillAlpha);
    }

    // Set titles
    if (dict("Title")) hg->SetTitle(dict["Title"]);
    if (dict("xAxisTitle")) hg->GetXaxis()->SetTitle(dict["xAxisTitle"]);
    if (dict("yAxisTitle")) hg->GetYaxis()->SetTitle(dict["yAxisTitle"]);

    // Set axes styles

    // x Axis
    if (dict("xAxisTitleFont"))
        hg->GetXaxis()->SetTitleFont(dict["xAxisTitleFont"]);
    if (dict("xAxisTitleSize"))
        hg->GetXaxis()->SetTitleSize(dict["xAxisTitleSize"].val());
    /* cout << " beta " << endl; */
    /* if (dict("xAxisTitleOffset")) cout << " Setting xaxistitleoffset " <<  dict["xAxisTitleOffset"].val() << endl; */
    if (dict("xAxisTitleOffset"))
        hg->GetXaxis()->SetTitleOffset(dict["xAxisTitleOffset"].val());
    if (dict("xAxisLabelFont"))
        hg->GetXaxis()->SetLabelFont(dict["xAxisLabelFont"]);
    if (dict("xAxisLabelSize"))
        hg->GetXaxis()->SetLabelSize(dict["xAxisLabelSize"].val());
    if (dict("xAxisLabelOffset"))
        hg->GetXaxis()->SetLabelOffset(dict["xAxisLabelOffset"].val());

    // y Axis
    if (dict("yAxisTitleFont"))
        hg->GetYaxis()->SetTitleFont(dict["yAxisTitleFont"]);
    if (dict("yAxisTitleSize"))
        hg->GetYaxis()->SetTitleSize(dict["yAxisTitleSize"].val());
    if (dict("yAxisTitleOffset"))
        hg->GetYaxis()->SetTitleOffset(dict["yAxisTitleOffset"].val());
    if (dict("yAxisLabelFont"))
        hg->GetYaxis()->SetLabelFont(dict["yAxisLabelFont"]);
    if (dict("yAxisLabelSize"))
        hg->GetYaxis()->SetLabelSize(dict["yAxisLabelSize"].val());
    if (dict("yAxisLabelOffset"))
        hg->GetYaxis()->SetLabelOffset(dict["yAxisLabelOffset"].val());


    // z Axis
    if (dict("zAxisTitleFont"))
        hg->GetZaxis()->SetTitleFont(dict["zAxisTitleFont"]);
    if (dict("zAxisTitleSize"))
        hg->GetZaxis()->SetTitleSize(dict["zAxisTitleSize"].val());
    if (dict("zAxisTitleOffset"))
        hg->GetZaxis()->SetTitleOffset(dict["zAxisTitleOffset"].val());
    if (dict("zAxisLabelFont"))
        hg->GetZaxis()->SetLabelFont(dict["zAxisLabelFont"]);
    if (dict("zAxisLabelSize"))
        hg->GetZaxis()->SetLabelSize(dict["zAxisLabelSize"].val());
    if (dict("zAxisLabelOffset"))
        hg->GetZaxis()->SetLabelOffset(dict["zAxisLabelOffset"].val());
    return hg;
};

TH1D* io_fmt (TH1D* hg, TPad* pad, ioOptMap _override, ioOptMap dict) {
    // cancel the x and y-axis if there is no margin for them
    if (!pad->GetBottomMargin()) {
        io_fmt(hg, {"no-xAxis:0;;"},{});
    }
    if (!pad->GetLeftMargin()) {
        io_fmt(hg, {"no-yAxis:0;;"},{});
    }
    io_fmt(hg, _override, dict);
    return hg;
};

TH1D* io_fmt (TH1D* hg, ioOptMap _override, ioOptMap dict) {
    /* cout << " alpha " << endl; */
    dict += _override;

    /* cout << " format " << endl; */
    /* cout << dict << endl; */

    if (dict("normalize")) hg->Scale(1./hg->Integral());
    if (dict("Rebin"))        hg->Rebin(dict["Rebin"]);

    if (!dict("MarkerAlpha")) dict["MarkerAlpha"] = 1.;
    if (dict("noTitle"))      hg->SetTitle("");
    if (dict("SetStats"))     hg->SetStats(dict["SetStats"]);


    if (dict("MarkerStyle"))
        hg->SetMarkerStyle(dict["MarkerStyle"]);
    if (dict("MarkerColor"))
        hg->SetMarkerColorAlpha(dict["MarkerColor"], dict["MarkerAlpha"].val());
    if (dict("MarkerSize"))
        hg->SetMarkerSize(dict["MarkerSize"].val());

    if (dict("LineWidth"))
        hg->SetLineWidth(dict["LineWidth"]);
    if (dict("LineStyle"))
        hg->SetLineStyle(dict["LineStyle"]);
    if (dict("MarkerColor"))
        hg->SetLineColorAlpha(dict["MarkerColor"], dict["MarkerAlpha"].val());
    
    if (dict("LineColor")) {
        double LineAlpha = dict("LineAlpha") ? dict["LineAlpha"].val() : dict["MarkerAlpha"].val();
        hg->SetLineColorAlpha(dict["LineColor"], LineAlpha);
    }

    if (dict("FillColor")) {
        double FillAlpha = dict("FillAlpha") ? dict["FillAlpha"].val() : dict["MarkerAlpha"].val();
        hg->SetFillColorAlpha(dict["FillColor"], FillAlpha);
    }

    // Set titles
    if (dict("Title")) hg->SetTitle(dict["Title"]);
    if (dict("xAxisTitle")) hg->GetXaxis()->SetTitle(dict["xAxisTitle"]);
    if (dict("yAxisTitle")) hg->GetYaxis()->SetTitle(dict["yAxisTitle"]);

    if (dict("no-xAxis")) { 
        hg->GetXaxis()->SetLabelColor(kWhite,0.);
        hg->GetXaxis()->SetTitle("");
    }
    if (dict("no-yAxis")) { 
        hg->GetYaxis()->SetLabelColor(kWhite,0.);
        hg->GetYaxis()->SetTitle("");
    }

    // x Axis
    if (dict("xAxisTitleFont"))
        hg->GetXaxis()->SetTitleFont(dict["xAxisTitleFont"]);
    if (dict("xAxisTitleSize"))
        hg->GetXaxis()->SetTitleSize(dict["xAxisTitleSize"].val());
    if (dict("xAxisTitleOffset"))
        hg->GetXaxis()->SetTitleOffset(dict["xAxisTitleOffset"].val());
    if (dict("xAxisLabelFont"))
        hg->GetXaxis()->SetLabelFont(dict["xAxisLabelFont"]);
    if (dict("xAxisLabelSize"))
        hg->GetXaxis()->SetLabelSize(dict["xAxisLabelSize"].val());
    if (dict("xAxisLabelOffset"))
        hg->GetXaxis()->SetLabelOffset(dict["xAxisLabelOffset"].val());

    // y Axis
    if (dict("yAxisTitleFont"))
        hg->GetYaxis()->SetTitleFont(dict["yAxisTitleFont"]);
    if (dict("yAxisTitleSize"))
        hg->GetYaxis()->SetTitleSize(dict["yAxisTitleSize"].val());
    if (dict("yAxisTitleOffset"))
        hg->GetYaxis()->SetTitleOffset(dict["yAxisTitleOffset"].val());
    if (dict("yAxisLabelFont"))
        hg->GetYaxis()->SetLabelFont(dict["yAxisLabelFont"]);
    if (dict("yAxisLabelSize"))
        hg->GetYaxis()->SetLabelSize(dict["yAxisLabelSize"].val());
    if (dict("yAxisLabelOffset"))
        hg->GetYaxis()->SetLabelOffset(dict["yAxisLabelOffset"].val());
    
    // Set Axis ranges with {x||y||z}AxisRange{Lo||Hi}
    if (dict("xAxisRangeLo") || dict("xAxisRangeHi")) {
        if (!dict("xAxisRangeLo") || !dict("xAxisRangeHi")) {
            cout << " Warning in io_fmt: has xAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting xAxisRange." << endl;
        } else {
            hg->GetXaxis()->SetRangeUser(dict["xAxisRangeLo"].val(), dict["xAxisRangeHi"].val());
        }
    }
    if (dict("yAxisRangeLo") || dict("yAxisRangeHi")) {
        if (!dict("yAxisRangeLo") || !dict("yAxisRangeHi")) {
            cout << " Warning in io_fmt: has yAxisRange{Lo||Hi} but not both. Needs both."<<endl;
            cout << " -> Not setting yAxisRange." << endl;
        } else {
            hg->GetYaxis()->SetRangeUser(dict["yAxisRangeLo"].val(), dict["yAxisRangeHi"].val());
        }
    }
    if (dict("zAxisRangeLo") || dict("zAxisRangeHi")) {
        if (!dict("zAxisRangeLo") || !dict("zAxisRangeHi")) {
            cout << " Warning in io_fmt: has zAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting zAxisRange." << endl;
        } else {
            hg->GetZaxis()->SetRangeUser(dict["zAxisRangeLo"].val(), dict["zAxisRangeHi"].val());
        }
    }
    // Set Ndivisions
    if (dict("xAxisNdivisions")) hg->GetXaxis()->SetNdivisions(dict["xAxisNdivisions"]);
    if (dict("yAxisNdivisions")) hg->GetYaxis()->SetNdivisions(dict["yAxisNdivisions"]);
    if (dict("zAxisNdivisions")) hg->GetZaxis()->SetNdivisions(dict["zAxisNdivisions"]);

    // z Axis
    if (dict("zAxisTitleFont"))
        hg->GetZaxis()->SetTitleFont(dict["zAxisTitleFont"]);
    if (dict("zAxisTitleSize"))
        hg->GetZaxis()->SetTitleSize(dict["zAxisTitleSize"].val());
    if (dict("zAxisTitleOffset"))
        hg->GetZaxis()->SetTitleOffset(dict["zAxisTitleOffset"].val());
    if (dict("zAxisLabelFont"))
        hg->GetZaxis()->SetLabelFont(dict["zAxisLabelFont"]);
    if (dict("zAxisLabelSize"))
        hg->GetZaxis()->SetLabelSize(dict["zAxisLabelSize"].val());
    if (dict("zAxisLabelOffset"))
        hg->GetZaxis()->SetLabelOffset(dict["zAxisLabelOffset"].val());

    if (dict("Normalize")) hg->Scale(1./hg->Integral());
    if (dict("ScaleByBinWidth")) io_scaleByBinWidth(hg);
    /* cout << " MarkerColor: " << hg->GetMarkerColor() << endl; */
    return hg;
};

TH2D* io_fmt (TH2D* hg, ioOptMap _override, ioOptMap dict) {
    /* cout << " alpha " << endl; */
    dict += _override;

    /* cout << " format " << endl; */
    /* cout << dict << endl; */

    if (dict("normalize")) hg->Scale(1./hg->Integral());
    /* if (dict("Rebin"))        hg->Rebin(dict["Rebin"]); */

    /* if (!dict("MarkerAlpha")) dict["MarkerAlpha"] = 1.; */
    if (dict("noTitle"))      hg->SetTitle("");
    if (dict("SetStats"))     hg->SetStats(dict["SetStats"]);


    /* if (dict("MarkerStyle")) */
        /* hg->SetMarkerStyle(dict["MarkerStyle"]); */
    /* if (dict("MarkerColor")) */
        /* hg->SetMarkerColorAlpha(dict["MarkerColor"], dict["MarkerAlpha"].val()); */
    /* if (dict("MarkerSize")) */
        /* hg->SetMarkerSize(dict["MarkerSize"].val()); */

    /* if (dict("LineWidth")) */
        /* hg->SetLineWidth(dict["LineWidth"]); */
    /* if (dict("LineStyle")) */
        /* hg->SetLineStyle(dict["LineStyle"]); */
    /* if (dict("MarkerColor")) */
        /* hg->SetLineColorAlpha(dict["MarkerColor"], dict["MarkerAlpha"].val()); */
    
    /* if (dict("LineColor")) { */
        /* double LineAlpha = dict("LineAlpha") ? dict["LineAlpha"].val() : dict["MarkerAlpha"].val(); */
        /* hg->SetLineColorAlpha(dict["LineColor"], LineAlpha); */
    /* } */

    // Set titles
    if (dict("Title")) hg->SetTitle(dict["Title"]);
    if (dict("xAxisTitle")) hg->GetXaxis()->SetTitle(dict["xAxisTitle"]);
    if (dict("yAxisTitle")) hg->GetYaxis()->SetTitle(dict["yAxisTitle"]);

    if (dict("no-xAxis")) { 
        hg->GetXaxis()->SetLabelColor(kWhite,0.);
        hg->GetXaxis()->SetTitle("");
    }
    if (dict("no-yAxis")) { 
        hg->GetYaxis()->SetLabelColor(kWhite,0.);
        hg->GetYaxis()->SetTitle("");
    }

    // x Axis
    if (dict("xAxisTitleFont"))
        hg->GetXaxis()->SetTitleFont(dict["xAxisTitleFont"]);
    if (dict("xAxisTitleSize"))
        hg->GetXaxis()->SetTitleSize(dict["xAxisTitleSize"].val());
    if (dict("xAxisTitleOffset"))
        hg->GetXaxis()->SetTitleOffset(dict["xAxisTitleOffset"].val());
    if (dict("xAxisLabelFont"))
        hg->GetXaxis()->SetLabelFont(dict["xAxisLabelFont"]);
    if (dict("xAxisLabelSize"))
        hg->GetXaxis()->SetLabelSize(dict["xAxisLabelSize"].val());
    if (dict("xAxisLabelOffset"))
        hg->GetXaxis()->SetLabelOffset(dict["xAxisLabelOffset"].val());

    // y Axis
    if (dict("yAxisTitleFont"))
        hg->GetYaxis()->SetTitleFont(dict["yAxisTitleFont"]);
    if (dict("yAxisTitleSize"))
        hg->GetYaxis()->SetTitleSize(dict["yAxisTitleSize"].val());
    if (dict("yAxisTitleOffset"))
        hg->GetYaxis()->SetTitleOffset(dict["yAxisTitleOffset"].val());
    if (dict("yAxisLabelFont"))
        hg->GetYaxis()->SetLabelFont(dict["yAxisLabelFont"]);
    if (dict("yAxisLabelSize"))
        hg->GetYaxis()->SetLabelSize(dict["yAxisLabelSize"].val());
    if (dict("yAxisLabelOffset"))
        hg->GetYaxis()->SetLabelOffset(dict["yAxisLabelOffset"].val());
    
    // Set Axis ranges with {x||y||z}AxisRange{Lo||Hi}
    if (dict("xAxisRangeLo") || dict("xAxisRangeHi")) {
        if (!dict("xAxisRangeLo") || !dict("xAxisRangeHi")) {
            cout << " Warning in io_fmt: has xAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting xAxisRange." << endl;
        } else {
            hg->GetXaxis()->SetRangeUser(dict["xAxisRangeLo"].val(), dict["xAxisRangeHi"].val());
        }
    }
    if (dict("yAxisRangeLo") || dict("yAxisRangeHi")) {
        if (!dict("yAxisRangeLo") || !dict("yAxisRangeHi")) {
            cout << " Warning in io_fmt: has yAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting yAxisRange." << endl;
        } else {
            hg->GetYaxis()->SetRangeUser(dict["yAxisRangeLo"].val(), dict["yAxisRangeHi"].val());
        }
    }
    if (dict("zAxisRangeLo") || dict("zAxisRangeHi")) {
        if (!dict("zAxisRangeLo") || !dict("zAxisRangeHi")) {
            cout << " Warning in io_fmt: has zAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting zAxisRange." << endl;
        } else {
            hg->GetZaxis()->SetRangeUser(dict["zAxisRangeLo"].val(), dict["zAxisRangeHi"].val());
        }
    }
    // Set Ndivisions
    if (dict("xAxisNdivisions")) hg->GetXaxis()->SetNdivisions(dict["xAxisNdivisions"]);
    if (dict("yAxisNdivisions")) hg->GetYaxis()->SetNdivisions(dict["yAxisNdivisions"]);
    if (dict("zAxisNdivisions")) hg->GetZaxis()->SetNdivisions(dict["zAxisNdivisions"]);

    // z Axis
    if (dict("zAxisTitleFont"))
        hg->GetZaxis()->SetTitleFont(dict["zAxisTitleFont"]);
    if (dict("zAxisTitleSize"))
        hg->GetZaxis()->SetTitleSize(dict["zAxisTitleSize"].val());
    if (dict("zAxisTitleOffset"))
        hg->GetZaxis()->SetTitleOffset(dict["zAxisTitleOffset"].val());
    if (dict("zAxisLabelFont"))
        hg->GetZaxis()->SetLabelFont(dict["zAxisLabelFont"]);
    if (dict("zAxisLabelSize"))
        hg->GetZaxis()->SetLabelSize(dict["zAxisLabelSize"].val());
    if (dict("zAxisLabelOffset"))
        hg->GetZaxis()->SetLabelOffset(dict["zAxisLabelOffset"].val());
    /* cout << " MarkerColor: " << hg->GetMarkerColor() << endl; */
    return hg;
};
// repeat of the above TGraph*
TGraph* io_fmt (TGraph* hg, TPad* pad, ioOptMap _override, ioOptMap dict) {
    // cancel the x and y-axis if there is no margin for them
    if (!pad->GetBottomMargin()) {
        io_fmt(hg, {"no-xAxis:0;;"},{});
    }
    if (!pad->GetLeftMargin()) {
        io_fmt(hg, {"no-yAxis:0;;"},{});
    }
    io_fmt(hg, _override, dict);
    return hg;
};

TGraph* io_fmt (TGraph* hg, ioOptMap _override, ioOptMap dict) {
    /* cout << " alpha " << endl; */
    dict += _override;

    /* cout << " format " << endl; */
    /* cout << dict << endl; */

    /* if (dict("normalize")) hg->Scale(1./hg->Integral()); */
    /* if (dict("Rebin"))        hg->Rebin(dict["Rebin"]); */

    if (!dict("MarkerAlpha")) dict["MarkerAlpha"] = 1.;
    if (dict("noTitle"))      hg->SetTitle("");
    /* if (dict("SetStats"))     hg->SetStats(dict["SetStats"]); */


    if (dict("MarkerStyle"))
        hg->SetMarkerStyle(dict["MarkerStyle"]);
    if (dict("MarkerColor"))
        hg->SetMarkerColorAlpha(dict["MarkerColor"], dict["MarkerAlpha"].val());
    if (dict("MarkerSize"))
        hg->SetMarkerSize(dict["MarkerSize"].val());

    if (dict("LineWidth"))
        hg->SetLineWidth(dict["LineWidth"]);
    if (dict("LineStyle"))
        hg->SetLineStyle(dict["LineStyle"]);
    if (dict("MarkerColor"))
        hg->SetLineColorAlpha(dict["MarkerColor"], dict["MarkerAlpha"].val());
    
    if (dict("LineColor")) {
        double LineAlpha = dict("LineAlpha") ? dict["LineAlpha"].val() : dict["MarkerAlpha"].val();
        hg->SetLineColorAlpha(dict["LineColor"], LineAlpha);
    }

    // Set titles
    if (dict("Title")) hg->SetTitle(dict["SetTitle"]);
    if (dict("xAxisTitle")) hg->GetXaxis()->SetTitle(dict["xAxisTitle"]);
    if (dict("yAxisTitle")) hg->GetYaxis()->SetTitle(dict["yAxisTitle"]);

    if (dict("no-xAxis")) { 
        hg->GetXaxis()->SetLabelColor(kWhite,0.);
        hg->GetXaxis()->SetTitle("");
    }
    if (dict("no-yAxis")) { 
        hg->GetYaxis()->SetLabelColor(kWhite,0.);
        hg->GetYaxis()->SetTitle("");
    }

    

    // Set axes styles

    // x Axis
    if (dict("xAxisTitleFont"))
        hg->GetXaxis()->SetTitleFont(dict["xAxisTitleFont"]);
    if (dict("xAxisTitleSize"))
        hg->GetXaxis()->SetTitleSize(dict["xAxisTitleSize"].val());
    /* cout << " beta " << endl; */
    /* if (dict("xAxisTitleOffset")) cout << " Setting xaxistitleoffset " <<  dict["xAxisTitleOffset"].val() << endl; */
    if (dict("xAxisTitleOffset"))
        hg->GetXaxis()->SetTitleOffset(dict["xAxisTitleOffset"].val());
    if (dict("xAxisLabelFont"))
        hg->GetXaxis()->SetLabelFont(dict["xAxisLabelFont"]);
    if (dict("xAxisLabelSize"))
        hg->GetXaxis()->SetLabelSize(dict["xAxisLabelSize"].val());
    if (dict("xAxisLabelOffset"))
        hg->GetXaxis()->SetLabelOffset(dict["xAxisLabelOffset"].val());

    // y Axis
    if (dict("yAxisTitleFont"))
        hg->GetYaxis()->SetTitleFont(dict["yAxisTitleFont"]);
    if (dict("yAxisTitleSize"))
        hg->GetYaxis()->SetTitleSize(dict["yAxisTitleSize"].val());
    if (dict("yAxisTitleOffset"))
        hg->GetYaxis()->SetTitleOffset(dict["yAxisTitleOffset"].val());
    if (dict("yAxisLabelFont"))
        hg->GetYaxis()->SetLabelFont(dict["yAxisLabelFont"]);
    if (dict("yAxisLabelSize"))
        hg->GetYaxis()->SetLabelSize(dict["yAxisLabelSize"].val());
    if (dict("yAxisLabelOffset"))
        hg->GetYaxis()->SetLabelOffset(dict["yAxisLabelOffset"].val());
    return hg;
};

