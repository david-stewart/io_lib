root -l <<EOF
    .x io_loadlibs.C
    .L loc_fnc.cc

    // Get the EAbbc spectra for triggers 8-12 GeV and 12-30 GeV

    /* auto data_trig = (THnSparseD*) got("hadd_zdcXall_vzall_trig8_30.root","data_trig"); */
    ioGetter got{};

    auto data_jet = (THnSparseD*) got("hadd_zdcXall_vzall_trig8_30.root","data_jet");
    auto data_trig = (THnSparseD*) got("hadd_zdcXall_vzall_trig8_30.root","data_trigs");
    ioJetSpectraSparse sdata ( data_jet, data_trig );



/*     auto hg3 = (TH3D*) got("hadd_zdcXall_vzall_trig8_30.root","jet_8pi_EAbbc"); */

/*     hg3->GetZaxis()->SetRange(3,3); // high EA */
/*     auto hg2_hi = (TH2D*) hg3->Project3D("yx"); */ 
/*     hg2_hi->SetName(ioUniqueName()); */

/*     hg3->GetZaxis()->SetRange(1,1); // low EA */
/*     auto hg2_lo = (TH2D*) hg3->Project3D("yx"); */ 
/*     hg2_lo->SetName(ioUniqueName()); */

/*     auto trigs = (TH2D*) got("hadd_zdcXall_vzall_trig8_30.root", "nTrig_EAbbc"); */

    array<array<TH1D*,8>,2> jets_hi, jets_lo, jets_rat; // index 1 AbsDeltaPhi, index 2 -- EtTrigger
    ioOptMap fmt {{ "yAxisTitle", "S = #frac{1}{N_{trig}} #frac{dN_{jet}}{d#it{p}_{T}^{jet}}",
        "xAxisTitle","#it{p}_{T}^{jet} [GeV/#it{c}]",
        "yAxisTitleSize", 20, "yAxisTitleOffset",1.9,
        "xAxisTitleSize", 20, "xAxisTitleOffset",3.3,
        "yAxisRangeLo", 1.01e-6, "yAxisRangeHi",0.99,
        "Title",""
    }};

    ioPads pads{ {{0.55,0.99},{0.15,0.25,0.55}},{{0.01,0.15,0.99,0.999}},800,600};

    /* auto bins = jetbins_pt2; */
    auto bins = jetbins_pt2;
    /* auto ntrigs = trigs->Integral(); */

    /* auto ntrigs_hi = trigs->Integral(3,4,3,3); // (3,4 == TrigEt 8-30 GeV, 3,3 = high EA */
    /* auto ntrigs_lo = trigs->Integral(3,4,1,1); // (3,4 == TrigEt 8-30 GeV, 1,1 = low EA */

   TLegend *leg = new TLegend(0.8017748,0.7346907,0.9680068,0.9438131,NULL,"brNDC");
   io_fmt(leg);

    pads(0)->SetLogy();
    double r_range{50.};
    double l_range{8.};

    for (auto k : vector<int> {0,1}) {
    for (auto i : vector<int> {0,7}) {
        if (k==0) {
            sdata.range_TrigEt(9,12); // lower EtTrig value 8-12
        } else {
            sdata.range_TrigEt(13,30); // lower EtTrig value 8-12
        }

        sdata.range_EAbbc(3,3); // high BBC value
        fmt["MarkerSize"]  = 0.6*sizes8[i];
        fmt["MarkerStyle"] = k == 0 ? fullshapes8[i] : i == 0 ? kFullStar : kFullTriangleUp ;
        fmt["MarkerColor"] = k == 0 ? colors8[i]     : i == 0 ? kCyan : kMagenta;
        fmt["MarkerAlpha"] = 0.8;
        auto hg_hi = sdata.hg_JetPt8(i+1);
        jets_hi[i][k] = (TH1D*) hg_hi->Rebin(bins, ioUniqueName(), bins);
        io_scaleByBinWidth(jets_hi[i][k]);
        /* jets_hi[i]->Scale(1./ntrigs_hi); */
        jets_hi[i][k]->GetXaxis()->SetRangeUser(l_range,r_range);
        io_fmt(jets_hi[i][k],fmt);

        sdata.range_EAbbc(1,1);
        fmt["MarkerStyle"] =  k == 0 ? openshapes8[i] : i == 0 ? kOpenStar : kOpenTriangleUp ;
        auto hg_lo = sdata.hg_JetPt8(i+1);
        jets_lo[i][k] = (TH1D*) hg_lo->Rebin(bins, ioUniqueName(), bins);
        io_scaleByBinWidth(jets_lo[i][k]);
        /* jets_lo[i][k]->Scale(1./ntrigs_lo); */
        io_fmt(jets_lo[i][k],fmt);

        pads(0);
        jets_hi[i][k]->Draw( i==0 ? "PE" : "PE same" );
        ioWaitPrimitive();
        jets_lo[i][k]->Draw( "PE same" );
        leg->AddEntry(jets_hi[i][k],Deltaphi8[i]);

        jets_rat[i][k] = ioDivideTH1(jets_hi[i][k],jets_lo[i][k]);
        jets_rat[i][k]->GetYaxis()->SetTitle("#frac{S_{EA.BBC(0-30\%)}}{S_{EA.BBC(70-100\%)}}");
        pads(1);
        jets_rat[i][k]->SetMarkerStyle( k == 0 ? openshapes8[i] : i==0 ? kOpenStar : kOpenTraingleUp);
        jets_rat[i][k]->GetXaxis()->SetRangeUser(l_range,r_range);
        jets_rat[i][k]->GetYaxis()->SetNdivisions(7);
        jets_rat[i][k]->GetYaxis()->SetRangeUser(0.,1.29);
        jets_rat[i][k]->Draw(i==0 ? "PE" : "PE same" );
    }
    data.range_EAbbc(1,1); 
    TH1D* _hg = data.hg_JetPt8(8);
    hg = (TH1D*) _hg->Rebin(bins, ioUniqueName(), bins);
    io_scaleByBinWidth(hg);
    io_fmt(hg,{{"MarkerStyle",kOpenCircle,"MarkerSize",3.}});
    pads(0);
    hg->Draw("PE same");

    data.range_EAbbc(3,3);
    TH1D* _zg = data.hg_JetPt8(1);
    zg = (TH1D*) _zg->Rebin(bins, ioUniqueName(), bins);
    io_scaleByBinWidth(zg);
    io_fmt(zg,{{"MarkerStyle",kOpenSquare,"MarkerSize",3.}});
    pads(0);
    zg->Draw("PE same");
    /* cout << " CHECK " << data.data_trig->Projection(1)->Integral() << " vs " << ntrigs_hi << endl; */

    pads(0);
    leg->Draw();
    pads.stamp("`date` `pwd`/$0");
    pads(1);
    ioDrawTLine(l_range,1.,r_range,1., {{"LineColor",kGray+2,"LineStyle",7}});

    pads.canvas_pad->cd();
    ioDrawTLatex("Open Markers: 70-90\% EA", 0.20, 0.7,{{"TextSize",18}});
    ioDrawTLatex("Full Markers:  0-30\% EA", 0.20, 0.66,{{"TextSize",18}});
    /* pads.canvas->SaveAs("TrigTrans_EAbbc.pdf"); */
    /* pads.canvas->SaveAs("all_EA_all_dphi.pdf"); */
    ioWaitPrimitive();

    
EOF
