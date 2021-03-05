#include "RunListId.h"
#include <iostream>
#include <sstream>
#include <fstream>

RunListId::RunListId(const char* file_name, bool skip_127053_138064) {
    // skip_after_mid_comment for two 
    ifstream file;
    file.open(file_name);
    if (!file.is_open()){
        cout << "fatal error, couldn't open input file \""<< file_name <<"\"" << endl;
        exit(2);
    }
    string line;
    while (getline(file,line)){
        if (line.substr(0,1)=="#") continue;
        int id, runid, sec;
        stringstream words(line);
        words >> id >> runid >> sec;
        if (skip_127053_138064) {
            if (runid == 16127053 || runid == 16138064) continue;
        }
        map_id [runid] = static_cast<double>(id);
        run_sec[runid] = sec;
    }
};

bool RunListId::has_run(int id) {
    return static_cast<bool>(map_id.count(id));
};

int RunListId::size() { return map_id.size(); };

double RunListId::operator()(int run_id) {
    set_id(run_id);
    return id;
};

double RunListId::set_id(int run_id) {
    if (has_run(run_id)) {
        id = static_cast<double>(map_id[run_id]);
    } else {
        id = -1.;
    }
    return id;
};


