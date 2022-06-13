#include "TSystem.h"
extern TSystem* gSystem;
void tu_loadlibs() {
    string line;
    ifstream f_in;
    f_in.open("/Users/hl7947/software/io_lib/tu_lib_list");
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
