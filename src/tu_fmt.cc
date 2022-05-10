#include "tu_fmt.h"
/* #include "tu_fnc.h" */

TLine* tu_fmt(TLine* line, tuOptMap options, tuOptMap dict)
{
    dict += options;
    line->SetLineColorAlpha(dict("LineColor"),dict("LineColorAlpha"));
    line->SetLineStyle(dict("LineStyle"));
    line->SetLineWidth(dict("LineWidth"));
    return line;
};

tuOptMap tu_fmt__hg_dict() {
    return {{
        "MarkerStyle", kFullCircle,
        "MarkerColor", kBlack,
        "MarkerAlpha", 1.,
        "MarkerSize",  1.,

        "LineWidth",   1 ,
        "LineStyle",   1 ,

        "FillStyle", 4000,
        /* {"LineColor",   kBlack }, */
        /* {"LineAlpha",   1. }, */

        "xAxisTitleFont",     43,
        "xAxisTitleSize",     22,
        "xAxisTitleOffset", 1.6 ,

        "xAxisLabelFont",      43,
        "xAxisLabelSize",      18,
        "xAxisLabelOffset",  0.02,

        "yAxisTitleFont",     43,
        "yAxisTitleSize",     22,
        "yAxisTitleOffset", 1.85,

        "yAxisLabelFont",     43,
        "yAxisLabelSize",     18,
        "yAxisLabelOffset",  0.02,

        "zAxisTitleFont", 43,
        "zAxisTitleSize", 22,
        "zAxisTitleOffset", 2.4 ,

        "zAxisLabelFont", 43,
        "zAxisLabelSize", 18,
        "zAxisLabelOffset", 0.02 ,

        "SetStats", 0
    }};
};

TH1* tu_fmt (TH1* hg, tuOptMap _override, tuOptMap dict) {
    // if don't want defaults, add {} to override value
    dict += _override;

    if (dict["normalize"])    hg->Scale(1./hg->Integral());
    if (dict["Rebin"])        hg->Rebin(dict("Rebin"));

    if (dict["noTitle"])      hg->SetTitle("");
    if (dict["SetStats"])     hg->SetStats((int)dict("SetStats"));


    if (dict["MarkerStyle"])
        hg->SetMarkerStyle(dict("MarkerStyle"));
    if (dict["MarkerColor"])
        hg->SetMarkerColorAlpha(dict("MarkerColor"), dict("MarkerAlpha",1.));
    if (dict["MarkerSize"])
        hg->SetMarkerSize(dict("MarkerSize"));

    if (dict["LineWidth"])
        hg->SetLineWidth(dict("LineWidth"));
    if (dict["LineStyle"])
        hg->SetLineStyle(dict("LineStyle"));
    if (dict["LineColor"])
        hg->SetLineColorAlpha(dict("LineColor"), dict("LineAlpha",dict("MarkerAlpha",1.)));
    else if (dict["MarkerColor"])
        hg->SetLineColorAlpha(dict("MarkerColor"), dict("LineAlpha",dict("MarkerAlpha",1.)));
    
    if (dict["FillColor"]) 
        hg->SetFillColorAlpha(dict("FillColor"), dict("FillAlpha",1.));

    // Set titles
    if (dict["Title"]) hg->SetTitle(dict("Title"));
    if (dict["xAxisTitle"]) hg->GetXaxis()->SetTitle(dict("xAxisTitle"));
    if (dict["yAxisTitle"]) hg->GetYaxis()->SetTitle(dict("yAxisTitle"));

    // Set axes styles

    if (dict["TitleSize"]) hg->SetTitleSize(dict("TitleSize"));
    if (dict["TitleFont"]) hg->SetTitleFont(dict("TitleFont"));
    // x Axis
    if (dict["xAxisTitleFont"])
        hg->GetXaxis()->SetTitleFont(dict("xAxisTitleFont"));
    if (dict["xAxisTitleSize"])
        hg->GetXaxis()->SetTitleSize(dict("xAxisTitleSize"));
    /* cout << " beta " << endl; */
    /* if (dict["xAxisTitleOffset"]) cout << " Setting xaxistitleoffset " <<  dict("xAxisTitleOffset") << endl; */
    if (dict["xAxisTitleOffset"])
        hg->GetXaxis()->SetTitleOffset(dict("xAxisTitleOffset"));
    if (dict["xAxisLabelFont"])
        hg->GetXaxis()->SetLabelFont(dict("xAxisLabelFont"));
    if (dict["xAxisLabelSize"])
        hg->GetXaxis()->SetLabelSize(dict("xAxisLabelSize"));
    if (dict["xAxisLabelOffset"])
        hg->GetXaxis()->SetLabelOffset(dict("xAxisLabelOffset"));

    // y Axis
    if (dict["yAxisTitleFont"])
        hg->GetYaxis()->SetTitleFont(dict("yAxisTitleFont"));
    if (dict["yAxisTitleSize"])
        hg->GetYaxis()->SetTitleSize(dict("yAxisTitleSize"));
    if (dict["yAxisTitleOffset"])
        hg->GetYaxis()->SetTitleOffset(dict("yAxisTitleOffset"));
    if (dict["yAxisLabelFont"])
        hg->GetYaxis()->SetLabelFont(dict("yAxisLabelFont"));
    if (dict["yAxisLabelSize"])
        hg->GetYaxis()->SetLabelSize(dict("yAxisLabelSize"));
    if (dict["yAxisLabelOffset"])
        hg->GetYaxis()->SetLabelOffset(dict("yAxisLabelOffset"));


    // z Axis
    if (dict["zAxisTitleFont"])
        hg->GetZaxis()->SetTitleFont(dict("zAxisTitleFont"));
    if (dict["zAxisTitleSize"])
        hg->GetZaxis()->SetTitleSize(dict("zAxisTitleSize"));
    if (dict["zAxisTitleOffset"])
        hg->GetZaxis()->SetTitleOffset(dict("zAxisTitleOffset"));
    if (dict["zAxisLabelFont"])
        hg->GetZaxis()->SetLabelFont(dict("zAxisLabelFont"));
    if (dict["zAxisLabelSize"])
        hg->GetZaxis()->SetLabelSize(dict("zAxisLabelSize"));
    if (dict["zAxisLabelOffset"])
        hg->GetZaxis()->SetLabelOffset(dict("zAxisLabelOffset"));
    
    // Set Axis ranges with {x||y||z}AxisRange{Lo||Hi}
    if (dict["xAxisRangeLo"] || dict["xAxisRangeHi"]) {
        if (!dict["xAxisRangeLo"] || !dict["xAxisRangeHi"]) {
            cout << " Warning in tu_fmt: has xAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting xAxisRange." << endl;
        } else {
            hg->GetXaxis()->SetRangeUser(dict("xAxisRangeLo"), dict("xAxisRangeHi"));
        }
    }
    if (dict["yAxisRangeLo"] || dict["yAxisRangeHi"]) {
        if (!dict["yAxisRangeLo"] || !dict["yAxisRangeHi"]) {
            cout << " Warning in tu_fmt: has yAxisRange{Lo||Hi} but not both. Needs both."<<endl;
            cout << " -> Not setting yAxisRange." << endl;
        } else {
            hg->GetYaxis()->SetRangeUser(dict("yAxisRangeLo"), dict("yAxisRangeHi"));
        }
    }
    if (dict["zAxisRangeLo"] || dict["zAxisRangeHi"]) {
        if (!dict["zAxisRangeLo"] || !dict["zAxisRangeHi"]) {
            cout << " Warning in tu_fmt: has zAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting zAxisRange." << endl;
        } else {
            hg->GetZaxis()->SetRangeUser(dict("zAxisRangeLo"), dict("zAxisRangeHi"));
        }
    }
    // Set Ndivisions
    if (dict["xAxisNdivisions"]) hg->GetXaxis()->SetNdivisions(dict("xAxisNdivisions"));
    if (dict["yAxisNdivisions"]) hg->GetYaxis()->SetNdivisions(dict("yAxisNdivisions"));
    if (dict["zAxisNdivisions"]) hg->GetZaxis()->SetNdivisions(dict("zAxisNdivisions"));
    return hg;
};

/* TH1D* tu_fmt (TH1D* hg, TPad* pad, tuOptMap _override, tuOptMap dict) { */
/*     // cancel the x and y-axis if there is no margin for them */
/*     if (!pad->GetBottomMargin()) { */
/*         tu_fmt(hg, {"no-xAxis:0;;"},{}); */
/*     } */
/*     if (!pad->GetLeftMargin()) { */
/*         tu_fmt(hg, {"no-yAxis:0;;"},{}); */
/*     } */
/*     tu_fmt(hg, _override, dict); */
/*     return hg; */
/* }; */

// repeat of the above TGraph*
/* TGraph* tu_fmt (TGraph* hg, TPad* pad, tuOptMap _override, tuOptMap dict) { */
/*     // cancel the x and y-axis if there is no margin for them */
/*     if (!pad->GetBottomMargin()) { */
/*         tu_fmt(hg, {"no-xAxis:0;;"},{}); */
/*     } */
/*     if (!pad->GetLeftMargin()) { */
/*         tu_fmt(hg, {"no-yAxis:0;;"},{}); */
/*     } */
/*     tu_fmt(hg, _override, dict); */
/*     return hg; */
/* }; */


TGraph* tu_fmt (TGraph* hg, tuOptMap _override, tuOptMap dict) {
    /* cout << " alpha " << endl; */
    dict += _override;

    if (dict["Title"])      hg->SetTitle(dict("Title"));

    if (dict["MarkerStyle"])
        hg->SetMarkerStyle(dict("MarkerStyle"));
    if (dict["MarkerColor"])
        hg->SetMarkerColorAlpha(dict("MarkerColor"), dict("MarkerAlpha",1.));
    if (dict["MarkerSize"])
        hg->SetMarkerSize(dict("MarkerSize"));

    if (dict["LineWidth"])
        hg->SetLineWidth(dict("LineWidth"));
    if (dict["LineStyle"])
        hg->SetLineStyle(dict("LineStyle"));
    if (dict["LineColor"])
        hg->SetLineColorAlpha(dict("LineColor"), dict("LineAlpha",dict("MarkerAlpha",1.)));
    else if (dict["MarkerColor"])
        hg->SetLineColorAlpha(dict("MarkerColor"), dict("LineAlpha",dict("MarkerAlpha",1.)));

    // Set titles
    if (dict["Title"]) hg->SetTitle(dict("SetTitle"));
    if (dict["xAxisTitle"]) hg->GetXaxis()->SetTitle(dict("xAxisTitle"));
    if (dict["yAxisTitle"]) hg->GetYaxis()->SetTitle(dict("yAxisTitle"));

    if (dict["no-xAxis"]) { 
        hg->GetXaxis()->SetLabelColor(kWhite,0.);
        hg->GetXaxis()->SetTitle("");
    }
    if (dict["no-yAxis"]) { 
        hg->GetYaxis()->SetLabelColor(kWhite,0.);
        hg->GetYaxis()->SetTitle("");
    }

    // Set axes styles
    // x Axis
    if (dict["xAxisTitleFont"])
        hg->GetXaxis()->SetTitleFont(dict("xAxisTitleFont"));
    if (dict["xAxisTitleSize"])
        hg->GetXaxis()->SetTitleSize(dict("xAxisTitleSize"));
    /* cout << " beta " << endl; */
    /* if (dict["xAxisTitleOffset"]) cout << " Setting xaxistitleoffset " <<  dict("xAxisTitleOffset") << endl; */
    if (dict["xAxisTitleOffset"])
        hg->GetXaxis()->SetTitleOffset(dict("xAxisTitleOffset"));
    if (dict["xAxisLabelFont"])
        hg->GetXaxis()->SetLabelFont(dict("xAxisLabelFont"));
    if (dict["xAxisLabelSize"])
        hg->GetXaxis()->SetLabelSize(dict("xAxisLabelSize"));
    if (dict["xAxisLabelOffset"])
        hg->GetXaxis()->SetLabelOffset(dict("xAxisLabelOffset"));

    // y Axis
    if (dict["yAxisTitleFont"])
        hg->GetYaxis()->SetTitleFont(dict("yAxisTitleFont"));
    if (dict["yAxisTitleSize"])
        hg->GetYaxis()->SetTitleSize(dict("yAxisTitleSize"));
    if (dict["yAxisTitleOffset"])
        hg->GetYaxis()->SetTitleOffset(dict("yAxisTitleOffset"));
    if (dict["yAxisLabelFont"])
        hg->GetYaxis()->SetLabelFont(dict("yAxisLabelFont"));
    if (dict["yAxisLabelSize"])
        hg->GetYaxis()->SetLabelSize(dict("yAxisLabelSize"));
    if (dict["yAxisLabelOffset"])
        hg->GetYaxis()->SetLabelOffset(dict("yAxisLabelOffset"));
    if (dict["xAxisRangeLo"] || dict["xAxisRangeHi"]) {
        if (!dict["xAxisRangeLo"] || !dict["xAxisRangeHi"]) {
            cout << " Warning in tu_fmt: has xAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting xAxisRange." << endl;
        } else {
            hg->GetXaxis()->SetLimits(dict("xAxisRangeLo"), dict("xAxisRangeHi"));
        }
    }
    if (dict["xAxisNdivisions"]) hg->GetXaxis()->SetNdivisions(dict("xAxisNdivisions"));
    if (dict["yAxisNdivisions"]) hg->GetYaxis()->SetNdivisions(dict("yAxisNdivisions"));

    if (dict["yAxisRangeLo"]) hg->SetMinimum(dict("yAxisRangeLo"));
    if (dict["yAxisRangeHi"]) hg->SetMaximum(dict("yAxisRangeHi"));
    return hg;
};

/* TGraphErrors* tu_fmt (TGraphErrors* hg, tuOptMap _override, tuOptMap dict) { */
/*     /1* cout << " alpha " << endl; *1/ */
/*     dict += _override; */

/*     /1* cout << " format " << endl; *1/ */
/*     /1* cout << dict << endl; *1/ */

/*     /1* if (dict["normalize"]) hg->Scale(1./hg->Integral()); *1/ */
/*     /1* if (dict["Rebin"])        hg->Rebin(dict("Rebin")); *1/ */

/*     if (!dict["MarkerAlpha"]) dict("MarkerAlpha") = 1.; */
/*     if (dict["noTitle"])      hg->SetTitle(""); */
/*     /1* if (dict["SetStats"])     hg->SetStats(dict("SetStats")); *1/ */


/*     if (dict["MarkerStyle"]) */
/*         hg->SetMarkerStyle(dict("MarkerStyle")); */
/*     if (dict["MarkerColor"]) */
/*         hg->SetMarkerColorAlpha(dict("MarkerColor"), dict("MarkerAlpha")); */
/*     if (dict["MarkerSize"]) */
/*         hg->SetMarkerSize(dict("MarkerSize")); */

/*     if (dict["LineWidth"]) */
/*         hg->SetLineWidth(dict("LineWidth")); */
/*     if (dict["LineStyle"]) */
/*         hg->SetLineStyle(dict("LineStyle")); */
/*     if (dict["LineColor"]) */
/*         hg->SetLineColorAlpha(dict("LineColor"), dict("LineAlpha",dict("MarkerAlpha",1.))); */
/*     else if (dict["MarkerColor"]) */
/*         hg->SetLineColorAlpha(dict("MarkerColor"), dict("LineAlpha",dict("MarkerAlpha",1.))); */

/*     // Set titles */
/*     if (dict["Title"]) hg->SetTitle(dict("SetTitle")); */
/*     if (dict["xAxisTitle"]) hg->GetXaxis()->SetTitle(dict("xAxisTitle")); */
/*     if (dict["yAxisTitle"]) hg->GetYaxis()->SetTitle(dict("yAxisTitle")); */

/*     if (dict["no-xAxis"]) { */ 
/*         hg->GetXaxis()->SetLabelColor(kWhite,0.); */
/*         hg->GetXaxis()->SetTitle(""); */
/*     } */
/*     if (dict["no-yAxis"]) { */ 
/*         hg->GetYaxis()->SetLabelColor(kWhite,0.); */
/*         hg->GetYaxis()->SetTitle(""); */
/*     } */

/*     // Set axes styles */
/*     // x Axis */
/*     if (dict["xAxisTitleFont"]) */
/*         hg->GetXaxis()->SetTitleFont(dict("xAxisTitleFont")); */
/*     if (dict["xAxisTitleSize"]) */
/*         hg->GetXaxis()->SetTitleSize(dict("xAxisTitleSize")); */
/*     /1* cout << " beta " << endl; *1/ */
/*     /1* if (dict["xAxisTitleOffset"]) cout << " Setting xaxistitleoffset " <<  dict("xAxisTitleOffset") << endl; *1/ */
/*     if (dict["xAxisTitleOffset"]) */
/*         hg->GetXaxis()->SetTitleOffset(dict("xAxisTitleOffset")); */
/*     if (dict["xAxisLabelFont"]) */
/*         hg->GetXaxis()->SetLabelFont(dict("xAxisLabelFont")); */
/*     if (dict["xAxisLabelSize"]) */
/*         hg->GetXaxis()->SetLabelSize(dict("xAxisLabelSize")); */
/*     if (dict["xAxisLabelOffset"]) */
/*         hg->GetXaxis()->SetLabelOffset(dict("xAxisLabelOffset")); */

/*     // y Axis */
/*     if (dict["yAxisTitleFont"]) */
/*         hg->GetYaxis()->SetTitleFont(dict("yAxisTitleFont")); */
/*     if (dict["yAxisTitleSize"]) */
/*         hg->GetYaxis()->SetTitleSize(dict("yAxisTitleSize")); */
/*     if (dict["yAxisTitleOffset"]) */
/*         hg->GetYaxis()->SetTitleOffset(dict("yAxisTitleOffset")); */
/*     if (dict["yAxisLabelFont"]) */
/*         hg->GetYaxis()->SetLabelFont(dict("yAxisLabelFont")); */
/*     if (dict["yAxisLabelSize"]) */
/*         hg->GetYaxis()->SetLabelSize(dict("yAxisLabelSize")); */
/*     if (dict["yAxisLabelOffset"]) */
/*         hg->GetYaxis()->SetLabelOffset(dict("yAxisLabelOffset")); */

/*     if (dict["xAxisRangeLo"] || dict["xAxisRangeHi"]) { */
/*         if (!dict["xAxisRangeLo"] || !dict["xAxisRangeHi"]) { */
/*             cout << " Warning in tu_fmt: has xAxisRange{lo||hi} but not both. Needs both."<<endl; */
/*             cout << " -> Not setting xAxisRange." << endl; */
/*         } else { */
/*             hg->GetXaxis()->SetLimits(dict("xAxisRangeLo"), dict("xAxisRangeHi")); */
/*         } */
/*     } */
/*     if (dict["yAxisRangeLo"]) hg->SetMinimum(dict("yAxisRangeLo")); */
/*     if (dict["yAxisRangeHi"]) hg->SetMaximum(dict("yAxisRangeHi")); */

/*     if (dict["xAxisNdivisions"]) hg->GetXaxis()->SetNdivisions(dict("xAxisNdivisions")); */
/*     if (dict["yAxisNdivisions"]) hg->GetYaxis()->SetNdivisions(dict("yAxisNdivisions")); */
/*     return hg; */
/* }; */
/* TMultiGraph* tu_fmt (TMultiGraph* hg, tuOptMap _override, tuOptMap dict) { */
/*     /1* cout << " alpha " << endl; *1/ */
/*     dict += _override; */

/*     if (!dict["MarkerAlpha"]) dict("MarkerAlpha") = 1.; */
/*     if (dict["noTitle"])      hg->SetTitle(""); */
/*     /1* if (dict["SetStats"])     hg->SetStats(dict("SetStats")); *1/ */

/*     // Set titles */
/*     if (dict["Title"]) hg->SetTitle(dict("SetTitle")); */
/*     if (dict["xAxisTitle"]) hg->GetXaxis()->SetTitle(dict("xAxisTitle")); */
/*     if (dict["yAxisTitle"]) hg->GetYaxis()->SetTitle(dict("yAxisTitle")); */

/*     if (dict["no-xAxis"]) { */ 
/*         hg->GetXaxis()->SetLabelColor(kWhite,0.); */
/*         hg->GetXaxis()->SetTitle(""); */
/*     } */
/*     if (dict["no-yAxis"]) { */ 
/*         hg->GetYaxis()->SetLabelColor(kWhite,0.); */
/*         hg->GetYaxis()->SetTitle(""); */
/*     } */
    
/*     // x Axis */
/*     if (dict["xAxisTitleFont"]) */
/*         hg->GetXaxis()->SetTitleFont(dict("xAxisTitleFont")); */
/*     if (dict["xAxisTitleSize"]) */
/*         hg->GetXaxis()->SetTitleSize(dict("xAxisTitleSize")); */
/*     if (dict["xAxisTitleOffset"]) */
/*         hg->GetXaxis()->SetTitleOffset(dict("xAxisTitleOffset")); */
/*     if (dict["xAxisLabelFont"]) */
/*         hg->GetXaxis()->SetLabelFont(dict("xAxisLabelFont")); */
/*     if (dict["xAxisLabelSize"]) */
/*         hg->GetXaxis()->SetLabelSize(dict("xAxisLabelSize")); */
/*     if (dict["xAxisLabelOffset"]) */
/*         hg->GetXaxis()->SetLabelOffset(dict("xAxisLabelOffset")); */

/*     // y Axis */
/*     if (dict["yAxisTitleFont"]) */
/*         hg->GetYaxis()->SetTitleFont(dict("yAxisTitleFont")); */
/*     if (dict["yAxisTitleSize"]) */
/*         hg->GetYaxis()->SetTitleSize(dict("yAxisTitleSize")); */
/*     if (dict["yAxisTitleOffset"]) */
/*         hg->GetYaxis()->SetTitleOffset(dict("yAxisTitleOffset")); */
/*     if (dict["yAxisLabelFont"]) */
/*         hg->GetYaxis()->SetLabelFont(dict("yAxisLabelFont")); */
/*     if (dict["yAxisLabelSize"]) */
/*         hg->GetYaxis()->SetLabelSize(dict("yAxisLabelSize")); */
/*     if (dict["yAxisLabelOffset"]) */
/*         hg->GetYaxis()->SetLabelOffset(dict("yAxisLabelOffset")); */

/*     if (dict["xAxisRangeLo"] || dict["xAxisRangeHi"]) { */
/*         if (!dict["xAxisRangeLo"] || !dict["xAxisRangeHi"]) { */
/*             cout << " Warning in tu_fmt: has xAxisRange{lo||hi} but not both. Needs both."<<endl; */
/*             cout << " -> Not setting xAxisRange." << endl; */
/*         } else { */
/*             hg->GetXaxis()->SetLimits(dict("xAxisRangeLo"), dict("xAxisRangeHi")); */
/*         } */
/*     } */
/*     if (dict["yAxisRangeLo"]) hg->SetMinimum(dict("yAxisRangeLo")); */
/*     if (dict["yAxisRangeHi"]) hg->SetMaximum(dict("yAxisRangeHi")); */
/*     return hg; */
/* }; */
/* TGraphAsymmErrors* tu_fmt (TGraphAsymmErrors* hg, tuOptMap _override, tuOptMap dict) { */
/*     dict += _override; */

/*     if (!dict["MarkerAlpha"]) dict("MarkerAlpha") = 1.; */
/*     if (dict["noTitle"])      hg->SetTitle(""); */
/*     if (!dict["FillStyle"])   dict("FillStyle") = 1001; */
/*     if (!dict["FillAlpha"])   dict("FillAlpha") = 1; */
/*     /1* if (dict["SetStats"])     hg->SetStats(dict("SetStats")); *1/ */


/*     if (dict["MarkerStyle"]) */
/*         hg->SetMarkerStyle(dict("MarkerStyle")); */
/*     if (dict["MarkerColor"]) */
/*         hg->SetMarkerColorAlpha(dict("MarkerColor"), dict("MarkerAlpha")); */
/*     if (dict["MarkerSize"]) */
/*         hg->SetMarkerSize(dict("MarkerSize")); */
/*     if (dict["FillColor"]) */
/*         hg->SetFillColorAlpha(dict("FillColor"), dict("FillAlpha")); */

/*     hg->SetFillStyle(dict("FillStyle")); */

/*     if (dict["LineWidth"]) */
/*         hg->SetLineWidth(dict("LineWidth")); */
/*     if (dict["LineStyle"]) */
/*         hg->SetLineStyle(dict("LineStyle")); */
/*     if (dict["LineColor"]) */
/*         hg->SetLineColorAlpha(dict("LineColor"), dict("LineAlpha",dict("MarkerAlpha",1.))); */
/*     else if (dict["MarkerColor"]) */
/*         hg->SetLineColorAlpha(dict("MarkerColor"), dict("LineAlpha",dict("MarkerAlpha",1.))); */


/*     // Set titles */
/*     if (dict["Title"]) hg->SetTitle(dict("SetTitle")); */
/*     if (dict["xAxisTitle"]) hg->GetXaxis()->SetTitle(dict("xAxisTitle")); */
/*     if (dict["yAxisTitle"]) hg->GetYaxis()->SetTitle(dict("yAxisTitle")); */

/*     if (dict["no-xAxis"]) { */ 
/*         hg->GetXaxis()->SetLabelColor(kWhite,0.); */
/*         hg->GetXaxis()->SetTitle(""); */
/*     } */
/*     if (dict["no-yAxis"]) { */ 
/*         hg->GetYaxis()->SetLabelColor(kWhite,0.); */
/*         hg->GetYaxis()->SetTitle(""); */
/*     } */

/*     // Set axes styles */
/*     // x Axis */
/*     if (dict["xAxisTitleFont"]) */
/*         hg->GetXaxis()->SetTitleFont(dict("xAxisTitleFont")); */
/*     if (dict["xAxisTitleSize"]) */
/*         hg->GetXaxis()->SetTitleSize(dict("xAxisTitleSize")); */
/*     /1* cout << " beta " << endl; *1/ */
/*     /1* if (dict["xAxisTitleOffset"]) cout << " Setting xaxistitleoffset " <<  dict("xAxisTitleOffset") << endl; *1/ */
/*     if (dict["xAxisTitleOffset"]) */
/*         hg->GetXaxis()->SetTitleOffset(dict("xAxisTitleOffset")); */
/*     if (dict["xAxisLabelFont"]) */
/*         hg->GetXaxis()->SetLabelFont(dict("xAxisLabelFont")); */
/*     if (dict["xAxisLabelSize"]) */
/*         hg->GetXaxis()->SetLabelSize(dict("xAxisLabelSize")); */
/*     if (dict["xAxisLabelOffset"]) */
/*         hg->GetXaxis()->SetLabelOffset(dict("xAxisLabelOffset")); */

/*     // y Axis */
/*     if (dict["yAxisTitleFont"]) */
/*         hg->GetYaxis()->SetTitleFont(dict("yAxisTitleFont")); */
/*     if (dict["yAxisTitleSize"]) */
/*         hg->GetYaxis()->SetTitleSize(dict("yAxisTitleSize")); */
/*     if (dict["yAxisTitleOffset"]) */
/*         hg->GetYaxis()->SetTitleOffset(dict("yAxisTitleOffset")); */
/*     if (dict["yAxisLabelFont"]) */
/*         hg->GetYaxis()->SetLabelFont(dict("yAxisLabelFont")); */
/*     if (dict["yAxisLabelSize"]) */
/*         hg->GetYaxis()->SetLabelSize(dict("yAxisLabelSize")); */
/*     if (dict["yAxisLabelOffset"]) */
/*         hg->GetYaxis()->SetLabelOffset(dict("yAxisLabelOffset")); */

/*     if (dict["xAxisRangeLo"] || dict["xAxisRangeHi"]) { */
/*         if (!dict["xAxisRangeLo"] || !dict["xAxisRangeHi"]) { */
/*             cout << " Warning in tu_fmt: has xAxisRange{lo||hi} but not both. Needs both."<<endl; */
/*             cout << " -> Not setting xAxisRange." << endl; */
/*         } else { */
/*             hg->GetXaxis()->SetLimits(dict("xAxisRangeLo"), dict("xAxisRangeHi")); */
/*         } */
/*     } */
/*     if (dict["yAxisRangeLo"]) hg->SetMinimum(dict("yAxisRangeLo")); */
/*     if (dict["yAxisRangeHi"]) hg->SetMaximum(dict("yAxisRangeHi")); */

/*     if (dict["xAxisNdivisions"]) hg->GetXaxis()->SetNdivisions(dict("xAxisNdivisions")); */
/*     if (dict["yAxisNdivisions"]) hg->GetYaxis()->SetNdivisions(dict("yAxisNdivisions")); */
/*     return hg; */
/* }; */

/* void tu_fmt_ranges(vector<TH1D*> hgs, tuOptMap dict) { */
/*     // Used to set common ranges on all hgs (assumed to have common x_range) */
/*     // 1. Find yMinVal:min(all hgs minimima) and yMaxVal:max(all hgs maxima) */
/*     // 2. Options: */
/*     //      "yRangeLowEdge":double Abs value */
/*     //      "yRangeLowRatio":double   = yMinVal-(yMaxVal-yMinVal)*yRangeUpRatio */
/*     //      "yRangeUpEdge":double  Abs value */
/*     //      "yRangeUpRatio":double    = yMinVal+(yMaxVal-yMinVal)*yRangeUpRatio */
/*     //      "yRange":no-option =  Default to 0. to 1.2 * (ymax - min) */
/*     // 3. For xRange: */
/*     //      "xRangeLowEdge":double Abs value */
/*     //      "xRangeLowRatio":double = bin of xMinVal-(xMaxVal-xMinVal)*xRangeUpRatio */
/*     //      "xRangeUpEdge":double  Abs value */
/*     //      "xRangeUpRatio":double    = yMinVal+(yMaxVal-yMinVal)*yRangeUpRatio */
/*     if   (     dict("yRangeLowEdge") */ 
/*             || dict("yRangeLowRatio") */
/*             || dict("yRangeUpEdge") */
/*             || dict("yRangeUpRatio") */
/*             || dict("yRange") ) */
/*     { */
/*         double yMinVal { hgs[0]->GetBinContent(hgs[0]->GetMinimumBin()) }; */
/*         double yMaxVal { hgs[0]->GetBinContent(hgs[0]->GetMaximumBin()) }; */
/*         for (unsigned int i{1}; i<hgs.size(); ++i) { */
/*             double minTemp { hgs[i]->GetBinContent(hgs[i]->GetMinimumBin()) }; */
/*             double maxTemp { hgs[i]->GetBinContent(hgs[i]->GetMaximumBin()) }; */
/*             if (minTemp < yMinVal) yMinVal = minTemp; */
/*             if (maxTemp > yMaxVal) yMaxVal = maxTemp; */
/*         } */

/*         double yMinRange; */
/*         if      (dict("yRangeLowEdge" ) ) yMinRange = dict("yRangeLowEdge"); */
/*         else if (dict["yRangeLowRatio"] ) yMinRange = */ 
/*             yMinVal - (yMaxVal-yMinVal) * dict("yRangeLowRatio"); */
/*         else if (yMinVal > 0) yMinRange = 0; */
/*         else                  yMinRange = yMinVal; */

/*         double yMaxRange; */
/*         if      (dict("yRangeUpEdge" ) ) yMaxRange = dict("yRangeUpEdge"); */
/*         else if (dict["yRangeUpRatio"] ) yMaxRange = */ 
/*             yMaxVal + (yMaxVal-yMaxVal) * dict("yRangeUpRatio"); */
/*         else yMaxRange = yMaxVal + (yMaxVal-yMinVal) * 0.2; */

/*         /1* cout << " beta yrange " << yMinVal << ","<< yMinRange << " " << yMaxVal<<","<<yMaxRange << endl; *1/ */
/*         for (auto hg :hgs) hg->GetYaxis()->SetRangeUser(yMinRange,yMaxRange); */
/*     } */

/*     if   (     dict("xRangeLowEdge") */ 
/*             || dict("xRangeLowRatio") */
/*             || dict("xRangeUpEdge") */
/*             || dict("xRangeUpRatio") */
/*             || dict("xRangeCommon") */
/*             || dict("xRange") ) */
/*     { */
/*         // For x_Range, find the smallest bin that has non-zero entries */
/*         int nbins { hgs[0]->GetXaxis()->GetNbins() }; */
/*         int i0 { nbins }; */
/*         int i1 { 1 }; */

/*         int i0_common{1}; */
/*         int i1_common{nbins}; */

/*         for (auto hg : hgs) { */
/*             int i0_s { nbins }; */
/*             int i1_s { 1 }; */
/*             for (int k{1}; k<=nbins; ++k) { */
/*                 if (hg->GetBinContent(k) != 0) { */ 
/*                     i0_s = k; */
/*                     break; */
/*                 } */
/*             } */
/*             if (i0_s < i0) i0 = i0_s; */
/*             if (i0_s > i0_common) i0_common = i0_s; */

/*             for (int k{nbins}; k>=1; --k) { */
/*                 if (hg->GetBinContent(k) != 0) { */ 
/*                     i1_s = k; */
/*                     break; */
/*                 } */
/*             } */
/*             if (i1_s > i1) i1 = i1_s; */
/*             if (i1_s < i1_common) i1_common = i1_s; */
/*         } */

/*         // if limiting to common range, use that and ignore all else */
/*         if (dict["xRangeCommon"]) { */
/*             for (auto hg : hgs) { */
/*                 hg->GetXaxis()->SetRange(i0_common, i1_common); */
/*                 /1* cout << " cutting: " << hg->GetXaxis()->GetBinLowEdge(i0_common) *1/ */ 
/*                 /1* << " to " <<        hg->GetXaxis()->GetBinUpEdge(i1_common) << endl; *1/ */
/*             } */
/*             return; */
/*         }; */

/*         TAxis* axis = hgs[0]->GetXaxis(); */

/*         double xMinVal { axis->GetBinLowEdge(i0) }; */
/*         double xMaxVal { axis->GetBinUpEdge( i1) }; */

/*         double xMinRange; */
/*         if      (dict("xRangeLowEdge" ) ) xMinRange = dict("xRangeLowEdge"); */
/*         else if (dict["xRangeLowRatio"] ) xMinRange = */ 
/*             xMinVal - (xMaxVal-xMinVal) * dict("xRangeLowRatio"); */
/*         else xMinRange = axis->GetBinLowEdge(1); */

/*         int iMinRange = axis->FindBin(xMinRange); */
/*         xMinRange = axis->GetBinLowEdge(iMinRange); */

/*         double xMaxRange; */
/*         if      (dict("xRangeUpEdge" ) ) xMaxRange = dict("xRangeUpEdge"); */
/*         else if (dict["xRangeUpRatio"] ) xMaxRange = */ 
/*             xMaxVal + (xMaxVal-xMaxVal) * dict("xRangeUpRatio"); */
/*         else xMaxRange = axis->GetBinUpEdge(nbins); */

/*         int iMaxRange = axis->FindBin(xMaxRange); */
/*         double b_edge = axis->GetBinLowEdge(iMaxRange); */
/*         double width  = axis->GetBinWidth  (iMaxRange); */
/*         double ratio  = (xMaxRange - b_edge) / width; */
/*         if (ratio > 0) xMaxRange = axis->GetBinUpEdge(iMaxRange); */
/*         else           xMaxRange = axis->GetBinLowEdge(iMaxRange); */

/*         for (auto hg :hgs) hg->GetXaxis()->SetRangeUser(xMinRange,xMaxRange); */

/*     } */
/* }; */
/* void tu_fmt_ranges(vector<TH1D*> hgs, vector<TH1D*> hgs_2, tuOptMap dict) { */
/*     vector<TH1D*> tmp; */
/*     for (auto& hg : hgs)   tmp.push_back(hg); */
/*     for (auto& hg : hgs_2) tmp.push_back(hg); */
/*     tu_fmt_ranges(tmp, dict); */
/* } */


TCanvas* tu_fmt (TCanvas* canv, tuOptMap dict) {
    canv->SetFillStyle(dict("FillStyle",4000));
    canv->SetFrameFillStyle(dict("FrameFillStyle",4000));
    return canv;
};
TPad* tu_fmt (TPad* canv, tuOptMap dict) {
    canv->SetFillStyle(dict("FillStyle",4000));
    canv->SetFrameFillStyle(dict("FrameFillStyle",4000));
    return canv;
};
/* TPad* tu_fmt (TPad* pad, tuOptMap dict) { */
/*     int FillStyle = dict("TCanvasFillStyle") ? */
/*         dict("TCanvasFillStyle").get_int() : 4000; */

/*     int FrameFillStyle = dict("TCanvasFrameFillStyle")  ? */
/*         dict("TCanvasFrameFillStyle").get_int() : 4000; */

/*     pad->SetFillStyle(FillStyle); */
/*     pad->SetFrameFillStyle(FrameFillStyle); */
/*     return pad; */
/* }; */

tuOptMap tu_fmt__leg_dict() {
    return {{
    "LineColor", kWhite,
    "LineColorAlpha",0.,
    "FillColor",kWhite,
    "FillColorAlpha",0.,
    "FillStyle",4000
    }};
};

TLegend* tu_fmt (TLegend* leg, tuOptMap _override, tuOptMap dict) {
    dict += _override;
    if (dict["LineColor"]) 
        leg->SetLineColorAlpha(dict("LineColor"),dict("LineColorAlpha",0.));
    if (dict["FillColor"]) 
        leg->SetLineColorAlpha(dict("FillColor"),dict("FillColorAlpha",0.));
    if (dict["FillStyle"])
        leg->SetFillStyle(dict("FillStyle"));
    if (dict["TextSize"]) leg->SetTextSize(dict("TextSize"));
    return leg;
};

/* TBox* tu_fmt(TBox* box, tuOptMap dict) { */
/*     int line_color = dict("Color",kBlack); */
/*     int fill_color = dict("Color",kWhite); */
/*     if (dict["LineColor"] && dict("FillColor")) { */
/*         line_color = dict("LineColor"); */
/*         fill_color = dict("FillColor"); */
/*     } else if (dict["LineColor"]) { */
/*         line_color = dict("LineColor"); */
/*         fill_color = dict("LineColor"); */
/*     } else if (dict["FillColor"]) { */
/*         line_color = dict("FillColor"); */
/*         fill_color = dict("FillColor"); */
/*     } */
/*     double fill_alpha = dict("FillAlpha") ? dict("FillAlpha") :  0.5; */
/*     double line_alpha = dict("LineAlpha") ? dict("FillAlpha") :  1.0; */
/*     cout << "Fill : " << fill_color << " " << fill_alpha << endl; */
/*     cout << "Fill : " << line_color << " " << line_alpha << endl; */
/*     box->SetLineColorAlpha(line_color, line_alpha); */
/*     box->SetFillColorAlpha(fill_color, fill_alpha); */
/*     return box; */
/* }; */
