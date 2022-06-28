root -l <<EOF
    .x tu_loadlibs.C
    tuGetter got{};

    auto h1 = (TH1D*) got("in_test.root", "h1");
    auto g2 = (TH1D*) h1->Clone("g1");
    auto stats = tuCalcRowStats(h1, 0.6, 0.8, 200.3,true);
    tu_fmt(h1,{{"MarkerColor",kRed,"MarkerStyle",kOpenCircle}});
    tu_fmt(g1,{{"MarkerColor",kBlue,"MarkerStyle",kOpenSquare}});
    g1->Draw("PE");
    h1->Draw("PEsame");
    cout << stats << endl;
    tuPause();

EOF
