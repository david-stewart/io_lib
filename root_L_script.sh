root -l <<EOF
gSystem->Load("${ROOUNFOLD}/libRooUnfold.so");
cout << ".L src/ioTHnSparse.cc+" << endl;
.L src/ioTHnSparse.cc+
cout << ".L src/io_fmt.cc+" << endl;
.L src/io_fmt.cc+
cout << ".L src/io_fnc.cc+" << endl;
.L src/io_fnc.cc+
cout << ".L src/ioOptMap.cc+" << endl;
.L src/ioOptMap.cc+
cout << ".L src/io_IOS.cc+" << endl;
.L src/io_IOS.cc+
cout << ".L src/io_operators.cc+" << endl;
.L src/io_operators.cc+
cout << ".L src/ioClass.cc+" << endl;
.L src/ioClass.cc+
cout << ".L src/ioCfnc.cc+" << endl;
.L src/ioCfnc.cc+
cout << ".L src/ioXsec_pAu2015.cc+" << endl;
.L src/ioXsec_pAu2015.cc+
cout << ".L src/ioTowerLoc.cc+" << endl;
.L src/ioTowerLoc.cc+
cout << ".L src/ioJetMatcher.cc+" << endl;
.L src/ioJetMatcher.cc+
cout << ".L src/ioJetMatcherGoodBins.cc+" << endl;
.L src/ioJetMatcherGoodBins.cc+
cout << ".L src/io_test.cc+" << endl;
.L src/io_test.cc+
EOF
