#ifndef RunListId__h
#define RunListId__h

#include <map>
#include "TProfile.h"
using namespace std;

struct RunListId {
    RunListId(const char* file_name, bool skip_127053_138064=true);
    map<int,double> map_id;  // runid -> map_id
    map<int,int>    run_sec; // runid -> duration
    bool has_run(int);
    double set_id(int run_id);
    double id; // last set id
    int  size();
    double operator()(int);  // return map_id[runid], unless not present, then -1.
};

#endif
