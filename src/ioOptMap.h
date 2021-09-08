#ifndef ioOptMap__h
#define ioOptMap__h

#include <sstream>
#include <iostream>
#include <map>

#include "TString.h"
#include "TPRegexp.h"
#include <iostream>
#include "TClonesArray.h"
#include "TObjString.h"

using std::cout;
using std::endl;
using std::string;
using std::ostream;
using std::map;
using std::vector;


// Takes entries as either strings, values, or pairs of values,
// and can return values accordingly
struct ioOpt {
    //statuses
    bool is_str {false};

    // data
    string my_str      {};
    double my_val      {1.};
    vector<double> vec {};

    // use:
    // interface:

    //type-casts
    operator const char*();// { return my_str.c_str(); };
    operator string();//
    operator int   ();//
    friend ostream& operator<< (ostream& os, const ioOpt& opt);// {
    // accessors:
    string str()             ;
    double val(int i=-1)     ; // my_val or vec[i]
    double operator()(int i=-1);
    int    get_int(int i=-1) ;
    int    operator[](int i) ;
    const char*  c_str()     ;

    // constructors
    ioOpt(const char* _);
    ioOpt(string      _);
    ioOpt(int         _);
    ioOpt(double      _);
    ioOpt(             );
    // to make a vector
    void push_back(double val); // fills my_vec
};


struct ioOptMap {
    /* ioOpt null_opt{0.}; */
    map<string,ioOpt> dict;

    void add(ioOpt); // { cout << " adding ioOpt " << endl; };
    void add(vector<ioOpt>); // { cout << " adding vector<ioOpt>" << endl; };

    ioOptMap(ioOpt _);
    ioOptMap(vector<ioOpt> _);
    ioOptMap();
    ioOpt& operator[](string _);

    bool has(string _);
    bool operator()(string _);
    ioOpt operator()(string name, ioOpt defVal);

    friend ostream& operator<< (ostream& os, const ioOptMap& _);
    ioOptMap operator+=(ioOptMap _);
    friend ioOptMap operator+ (ioOptMap lhs, const ioOptMap& rhs);

    ioOptMap& update(ioOptMap _);

};

#endif
