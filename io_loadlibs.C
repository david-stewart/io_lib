#include "TSystem.h"

extern TSystem* gSystem;

void io_loadlibs() {
    for (auto iolib : vector<string> {"Class","Cfnc","OptMap","_fnc","_fmt","TowerLoc","JetMatcher","JetMatcher100","_operators","Xsec_pAu2015", "_test", "_IOS","THnSparse","_pAu2015"}) {
        gSystem->Load(Form("io_lib/lib/libio%s.so",iolib.c_str())); 
    }
    /* gSystem->Load("io_lib/lib/liboiJetMaker.so"); */
};

