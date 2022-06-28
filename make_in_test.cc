root -l <<EOF
    .x tu_loadlibs.C
    TRandom3 _rand;
    TFile fout{"in_test.root","recreate"};

    TH1D h1 { "h1",";;",100,0.,300.};
    TH2D h2 { "h2",";;",100,0.,300., 100, 0., 300.};
    
    for (int i=0; i<1000000; ++i) {
        h1.Fill(_rand.Gaus(30., 15.));
        h1.Fill(_rand.Gaus(70., 15.));

        h2.Fill(_rand.Gaus(30., 15.), _rand.Gaus(40., 69.));
        h2.Fill(_rand.Gaus(70., 15.), _rand.Gaus(20., 23.));
    }

    h1.Write();
    h2.Write();

    fout.Save();
    fout.Close();

EOF

