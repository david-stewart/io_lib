#include "TSystem.h"
extern TSystem* gSystem;
void tu_loadlibs() {
    gSystem->Load("$ROOUNFOLD/libRooUnfold.so");
    string line;
    ifstream f_in;
    f_in.open(Form("%s/tu_lib_list",gSystem->Getenv("IO_LIB")));
    while (getline(f_in,line)) {
        istringstream iss(line);
        string word;
        while (iss >> word) {
            if (word[0] == '#') break;
            gSystem->Load(Form("$IO_LIB/lib/libtu%s.so",word.c_str())); 
        }
    };
};

//test end of text to puth to github
