root -l /home/dsjohnny/pc2021/pico-to-trees/tower_qa_study/hadd_500004.root<<EOF
    .X io_lib/io_loadlibs.C
    .ls
    /* hit8->Draw("PE"); */

    hit8->Draw("PE");
    ioHgStats hgs{hit8,true};
    double sigLo {-1.0};
    double sigHi { 5.0};
    io_fmt(
            hgs.points_above(sigHi,true),
            {{"MarkerStyle",kOpenCircle,"MarkerColor",kRed}}
            )->Draw("PEsame");
    /* io_fmt( */
    /*         hgs.points_below(sigLo,true), */
    /*         {{"MarkerStyle",kOpenCircle,"MarkerColor",kRed}} */
    /*         )->Draw("PEsame"); */
    io_fmt( hgs.get_horizontal_TLine(sigHi,true),{{"LineStyle",2,"LineColor",kBlue}})->Draw();
    /* io_fmt( hgs.get_horizontal_TLine(sigLo,true),{{"LineStyle",2,"LineColor",kBlue}})->Draw(); */
    ioDrawTLatex("#mu_{non-masked towers}+5#sigma", 500., hgs.mean_Xsigma(sigHi+0.7));

    ioWaitPrimitive();
EOF
