#include "TSystem.h"

extern TSystem* gSystem;

void io_loadlibs() {
    for (auto iolib : vector<string> {"Class","OptMap","_fnc","_fmt","TowerLoc","JetMatcher"}) {
        gSystem->Load(Form("io_lib/lib/libio%s.so",iolib.c_str())); 
    }
};
