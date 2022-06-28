root -l <<EOF
    .x tu_loadlibs.C
    tuGetter got{};

    auto h1 = (TH1D*) got("in_test.root", "h1");
    auto h2 = (TH2D*) got("in_test.root", "h2");

    auto g1 = (TH1D*) h1->Clone("g1");
    auto g2 = (TH2D*) h2->Clone("g2");

    tuPads pads{2};
    pads(0)->SetLogz();
    h2->Draw("colz");

    pads(2)->SetLogz();
    tuScrubBlock(g2,10,30,10,30);
    g2->Draw("colz");

    tuPause();

EOF
