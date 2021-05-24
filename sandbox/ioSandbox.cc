root -l << EOF
    gSystem->Load("../lib/liboiJetMaker.so")
    oiJetMaker areajets { {{"calc_areas",1.}} };
    oiJetMaker jets  { };

EOF
