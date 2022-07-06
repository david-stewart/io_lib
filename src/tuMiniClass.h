#ifndef tuMiniClass__h
#define tuMiniClass__h

#include <iostream>
#include "TFile.h"
#include <vector>
#include "TPad.h"
#include <fstream>
#include "TCanvas.h"
#include <map>
#include <string>
#include "TLine.h"
#include "TGraph.h"
#include "TTree.h"
#include <string>
#include "tuOptMap.h"

struct tuMinMaxPtr {
    tuMinMaxPtr(){};
    double min_val;
    double max_val;
    bool has_data{false};
    void* min_ptr{nullptr};
    void* max_ptr{nullptr};
    void operator()(double, void* ptr=nullptr);
};

struct tuMinMax {
    tuMinMax(string _name="");

    tuMinMax(double*,       string _name="");
    tuMinMax(int*,          string _name="");
    tuMinMax(unsigned int*, string _name="");
    tuMinMax(short*,        string _name="");
    tuMinMax(char*,         string _name="");

    tuMinMax(double*,       int& _index, string _name="");
    tuMinMax(int*,          int& _index, string _name="");
    tuMinMax(unsigned int*, int& _index, string _name="");
    tuMinMax(short*,        int& _index, string _name="");
    tuMinMax(char*,        int& _index, string _name="");

    tuMinMax(double&,       string _name="");
    tuMinMax(int&,          string _name="");
    tuMinMax(unsigned int&, string _name="");
    tuMinMax(short&,        string _name="");
    tuMinMax(char&,        string _name="");

    int fill_option{-1};
    string name;
    long long int n_entries{0};
    double min{0.};
    double max{0.};
    int    nbins(); // return (int)(max-min)+1
    long long int operator()();
    long long int operator()(int); // to fill an array then
    long long int fill(double val);

    double* ptr_double{nullptr}; // option 0
    int*    ptr_int{nullptr};    // option 1
    unsigned int* ptr_uint{nullptr};   // option 2
    short*  ptr_short{nullptr};  // option 3
    char*   ptr_char{nullptr};  // option 3
    int*    index{nullptr};

    friend ostream& operator<<(ostream& os, tuMinMax& self);
};
 

#endif
