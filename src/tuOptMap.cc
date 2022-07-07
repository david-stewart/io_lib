#include "tuOptMap.h"

using std::cout;
using std::endl;


// implement tuOpt --------------------------------------------
//    constructors
tuOpt::operator const char*() { return str.c_str(); }; 
tuOpt::operator string     () { return str; };
tuOpt::operator int        () { return val; };
tuOpt::operator double     () { return val; };
tuOpt::operator short      () { return (short)val; };
tuOpt::operator float      () { return (float)val; };

// friend os
ostream& operator<< (ostream& os, const tuOpt& opt) {
    if (opt.is_str) os << opt.str;
    else os << opt.val;
    /* os << (opt.is_str ? opt.str : opt.val) << endl; */
    return os;
};

/* string tuOpt::str()             { return str; }; */
const char*  tuOpt::c_str()     { return str.c_str(); };

//type-casts
tuOpt::tuOpt(const char* _) : is_str{true}, str{_} {};
tuOpt::tuOpt(string      _) : is_str{true}, str{_} {};
tuOpt::tuOpt(int         _) :               val{(double)_} {};
tuOpt::tuOpt(double      _) :               val{_} {};
tuOpt::tuOpt(             ) {};

// implement tuOptMap
/* tuOptMap::tuOptMap(tuOpt _)         { add(_); }; */
tuOptMap::tuOptMap(vector<tuOpt> _) { add(_); };
tuOptMap::tuOptMap()                {};


// check for existence of key
bool tuOptMap::operator[](string _) { return has(_); }; //if (dict.count(_)!=0) return dict[_]; else return blank; };
bool tuOptMap::has(string _) { return dict.count(_) != 0; };
                                                        //
// return values with ()
tuOpt& tuOptMap::operator()(string _) {  return dict[_]; };
tuOpt  tuOptMap::operator()(string name, tuOpt defVal) { return has(name) ? dict[name] : defVal; };

ostream& operator<< (ostream& os, const tuOptMap& _) {
        for (auto& v : _.dict)  os << Form("  %-12s -> ",v.first.c_str()) << v.second << endl;
        return os;
    };

tuOptMap tuOptMap::operator+=(tuOptMap _) { 
    for (auto& key : _.dict) dict[key.first] = key.second; return *this; 
};
tuOptMap operator+ (tuOptMap lhs, const tuOptMap& rhs) { lhs += rhs; return lhs; };
tuOptMap& tuOptMap::update(tuOptMap _) {
    for (auto& key : _.dict) dict[key.first] = key.second; return *this; 
}



/* void tuOptMap::add(tuOpt _) { */
/*     // only acceptable syntax is string with "key:value;;[key2:value2;; ...]" */
/*     if (!_.is_str) { */ 
/*         cout << " !: warning in tuOptMap.add: bad syntax to add option: " << _ */ 
/*              << endl << "    skipping entry" << endl; */
/*         return; */
/*     } */
    
/*     // parse the string */
/*     TString s1 { _.str() }; */
/*     /1* vector<pair<TString,TString>> rval; *1/ */
/*     TPRegexp r1("\\s*([\\w-]+):\\s*(.*?)\\s*;;(.*)"); */
/*     bool has_first = false; */
/*     while (s1(r1) != "") { */
/*         TObjArray* subStrL = r1.MatchS(s1); */
/*         const Int_t nrSubStr = subStrL->GetLast()+1; */
/*         /1* TString key, value; *1/ */
/*         if (nrSubStr > 2) { */
/*             string  key = ((TObjString *)subStrL->At(1))->GetString().Data(); */
/*             TString value  = ((TObjString *)subStrL->At(2))->GetString().Data(); */

/*             if (value.IsFloat()) dict[key] = value.Atof(); */
/*             else                 dict[key] = value.Data(); */
/*             has_first = true; */
/*             /1* add(key,value); *1/ */
/*             /1* rval.push_back({key,value}); *1/ */
/*         } else { */
/*             break; */
/*         } */
/*         if (nrSubStr > 3) { */
/*             s1 = ((TObjString*)subStrL->At(3))->GetString(); */
/*         } */
/*     } */
/*     if (has_first) { // check for a trailing value of key:value (without key:value;;) */
/*         TPRegexp strip("\\s*([\\w-]+):\\s*(\\S+(\\s*\\S+)*)"); */
/*         /1* TPRegexp strip("\\s*([\\w-]+):\\s*([^\\s]+.*[^\\s]+)\\s*"); *1/ */
/*         /1* TPRegexp strip("\\s*([\\w-]+):\\s*([^\\s]+.*[^\\s]+)\\s*"); *1/ */
/*         if (s1(strip) != "") { */
/*             cout << " in : " << endl; */
/*             TObjArray* subStrL = strip.MatchS(s1); */
/*             const Int_t nrSubStr = subStrL->GetLast()+1; */
/*             if (nrSubStr > 2) { */
/*                 cout << "0:"<<  ((TObjString *)subStrL->At(0))->GetString().Data() << endl; */
/*                 cout << "1:"<<  ((TObjString *)subStrL->At(1))->GetString().Data() << endl; */
/*                 cout << "2:"<<  ((TObjString *)subStrL->At(2))->GetString().Data() << endl; */

/*                 string  key    = ((TObjString *)subStrL->At(1))->GetString().Data(); */
/*                 TString value  = ((TObjString *)subStrL->At(2))->GetString().Data(); */
/*                 if (value.IsFloat()) dict[key] = value.Atof(); */
/*                 else                 dict[key] = value.Data(); */
/*             } */ 
/*         } */
/*     } */

/*     return; */
/* }; */

void tuOptMap::add(vector<tuOpt> args) {
    // input must of the form: string, val, string, val, etc...
    if (args.size() % 2)  
       throw std::runtime_error("Fatal: tuOptMap::add(args) : number of arguments added not even.");
    for (int i=0; i<(int)args.size(); i+= 2) {
        if (!args[i].is_str)
           throw std::runtime_error(Form("Fatal: in tuOptMap::add(args), arg[%i]=(%f) must be a string to be a key, but is not",i,args[i].val));
        dict[args[i].str]=args[i+1];
    }
};

