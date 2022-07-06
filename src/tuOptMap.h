#ifndef tuOptMap__h
#define tuOptMap__h

#include <sstream>
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
struct tuOpt {
    // values
    bool is_str {false};
    string str     {""};
    double val     {0.};

    // constructors
    tuOpt(const char* _);
    tuOpt(string      _);
    tuOpt(int         _);
    tuOpt(double      _);
    tuOpt(             );
    /* tuOpt(short        ); */
    /* tuOpt(float        ); */

    // friend os
    friend ostream& operator<< (ostream& os, const tuOpt& opt);// {
                                                               //
    //type-casts
    operator const char*();// { return my_str.c_str(); };
    operator string     ();//
    operator int        ();//
    operator double     ();
    operator short      ();
    operator float      ();

    const    char* c_str()     ;

    // to make a vector
    /* void push_back(double val); // fills my_vec */
};


struct tuOptMap {
    tuOpt blank_option {};
    map<string,tuOpt> dict;
    void add(vector<tuOpt>); // { cout << " adding vector<tuOpt>" << endl; };

    /* void add(tuOpt); // { cout << " adding tuOpt " << endl; }; */

    /* tuOptMap(tuOpt _); */
    tuOptMap(vector<tuOpt> _);
    tuOptMap();

    // check for existence of key
    bool operator[](string _);
    bool has(string _);

    // return values with ()
    tuOpt& operator()(string _);
    tuOpt  operator()(string name, tuOpt defVal);

    friend ostream& operator<< (ostream& os, const tuOptMap& _);
    tuOptMap operator+=(tuOptMap _);
    friend tuOptMap operator+ (tuOptMap lhs, const tuOptMap& rhs);

    tuOptMap& update(tuOptMap _);

};

#endif
