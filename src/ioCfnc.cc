#include "ioCfnc.h"
std::tuple<ioPads, TH2D*, TH1D*, TH1D*, TH1D*> ioMakeClosure(
        const char* file, 
        const char* rooResponseName,
        ioOptMap opt,
        ioOptMap dict
    ){

    dict += opt;
    ioPads pads { 
        {{0.4,0.5,0.99,0.99},{.25,.4},{0,0.1,0.25,0.25}},
        {{0.,0.15,0.80,0.99}} 
    };
    pads(0);
    ioGetter got;
    RooUnfoldResponse* response = (RooUnfoldResponse*) got(file, rooResponseName);

    string name_A = dict["name_A"].str()=="" 
        ? Form("%s_A",rooResponseName) : dict["name_A"].str();
    RooUnfoldResponse* resp_A = (RooUnfoldResponse*) got(file, name_A);

    string name_B = dict["name_B"].str()=="" 
        ? Form("%s_B",rooResponseName) : dict["name_B"].str();
    RooUnfoldResponse* resp_B = (RooUnfoldResponse*) got(file, name_B);
    
    pads(0)->SetLogz();
    TH2D* hg_response = (TH2D*) resp_A->Hresponse();
    io_fmt(hg_response,{{
       "xAxisTitle",dict["xTitleResponse"],
       "yAxisTitle",dict["yTitleResponse"],
       "zAxisRangeLo",dict["zAxisRangeLo"],
       "zAxisRangeHi",dict["zAxisRangeHo"],
       "Stats",dict["StatsResponse"],
       "Title",dict["TitleResponse"]}});
    if (dict("zAxisRangeLo") || dict("zAxisRangeHi")) {
        if (!dict("zAxisRangeLo") || !dict("zAxisRangeHi")) {
            cout << " Warning in ioCfnc.cc: has zAxisRange{lo||hi} but not both. Needs both."<<endl;
            cout << " -> Not setting zAxisRange." << endl;
        } else {
            hg_response->GetZaxis()->SetRangeUser(dict["zAxisRangeLo"].val(), dict["zAxisRangeHi"].val());
        }
    }

    hg_response->Draw("colz");

    pads(1)->SetLogy();

    TH1D* truth = (TH1D*) resp_A->Htruth();
    string truth_yAxisTitle = dict["yTitleTruth"].str()=="" 
        ? truth->GetYaxis()->GetTitle() : dict["yTitleTruth"].str();
    io_fmt(truth,{{
            "MarkerStyle",kOpenCircle,
            "Title",dict["TitleTruth"],
            "MarkerColor",kBlack,
            "yAxisTitle",dict["yTitleTruth"],
            /* "Normalize","-", */
            /* "ScaleByBinWidth","-" */
            }}
    );
    truth->Draw("PE");

    TH1D* measured = (TH1D*) resp_A->Hmeasured();
    TH1D* unf = io_BayesUnfold(measured, resp_B, dict["iter_unfold"]);
    io_fmt(unf,{{
            "MarkerStyle",kOpenSquare,
            "MarkerColor",kRed,
            /* "Normalize","-", */
            /* "ScaleByBinWidth","-" */
            }}
    );
    unf->Draw("PEsame");

    pads(2);
    auto rat = ioDivideTH1(unf,truth);
    io_fmt(rat,{{
            "yAxisTitle",dict["yTitleRatio"],
            "xAxisTitle",dict["xTitleRatio"],
            "Title","",
            "MarkerStyle",kOpenSquare,
            "MarkerColor",kRed}}); 
    rat->Draw("PE");

    pads.pads.push_back(pads.canvas_pad);
    return std::tie(pads, hg_response, truth, unf, rat);
};
RooUnfoldResponse ioMakeRooUnfoldResponse(
        const char* name, const char* file, const char* tag_M, const char* tag_T ) {
    ioBinVec bin_M { file, tag_M };
    ioBinVec bin_T { file, tag_T };
    TH1D truth {
            Form("%s_truth",ioUniqueName()),
            Form("%s;truth;N",ioUniqueName()),
            bin_T, bin_T};
    TH1D measured {
            Form("%s_measured",ioUniqueName()),
            Form("%s;measured;N",ioUniqueName()),
            bin_M, bin_M};
    return {&measured, &truth, Form("%s_RUR",name), "RooUnfoldResponse"};
};

TH1D* ioBlankTH1D(ioBinVec bins, ioOptMap dict, bool draw_vert_lines) {
    TH1D* hg = new TH1D(ioUniqueName(), ";;", bins, bins);
    hg->SetMarkerStyle(kDot);
    hg->SetMarkerColorAlpha(kWhite,0.);
    hg->SetLineColorAlpha(kWhite,0.);
    hg->SetTitle("");
    if (dict.dict.size()>0) {
        io_fmt(hg,dict);
        hg->Draw();
    }
    if (draw_vert_lines) {
        double y_lo = (dict("yAxisRangeLo")) ? dict["yAxisRangeLo"]() : 0.;
        double y_hi = (dict("yAxisRangeHi")) ? dict["yAxisRangeHi"]() : 1.;
        cout << "ylo " << y_lo << " y_hi " << y_hi << endl;
        for (int i=1; i<bins.size-1; ++i) {
            int line_width = dict("LineWidth",1);
            int line_style = dict("LineStyle",2);
            float alpha = (dict("LineColorAlpha")) ? dict["LineColorAlpha"].val() : 0.85;
            int line_color = dict("LineColor",kGray);
            ioDrawTLine( bins[i], y_lo, bins[i], y_hi, 
                    {{"LineWidth", line_width,
                      "LineStyle", line_style,
                      "LineColorAlpha", alpha,
                      "LineColor", line_color}} );
        };
    };
    return hg;
};
