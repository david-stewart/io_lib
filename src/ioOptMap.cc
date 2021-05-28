#include "ioOptMap.h"

using std::cout;
using std::endl;

ioOpt::operator const char*() { return my_str.c_str(); }; 
ioOpt::operator string     () { return my_str; };
ioOpt::operator int   ()      { return my_val; };

ostream& operator<< (ostream& os, const ioOpt& opt) {
    if (opt.is_str) 
        os << opt.my_str;
    else if (opt.vec.size() > 0) 
        for (auto &v : opt.vec) os << " " << v;
    else
        os << opt.my_val;
    return os;
};


string ioOpt::str()             { return my_str; };
double ioOpt::operator()(int i) { return val(i);};
int    ioOpt::get_int(int i   ) { return val(i); }; 
int    ioOpt::operator[](int i) { return val(i);};
const char*  ioOpt::c_str()     { return my_str.c_str(); };

// constructors
ioOpt::ioOpt(const char* _) : is_str{true}, my_str{_} {};
ioOpt::ioOpt(string      _) : is_str{true}, my_str{_} {};
ioOpt::ioOpt(int         _) :               my_val{(double)_} {};
ioOpt::ioOpt(double      _) :               my_val{_} {};
ioOpt::ioOpt(             ) {};
/* ioOpt(const ioOpt& _) { is_str = _.is_str; my_str = _.my_str; vec = _.vec; my_val = _.my_val; }; */

// to make a vector
void push_back(double val); // fills my_vec

double ioOpt::val(int i) {
    if (i==-1) return my_val;
    else if (i >= (int)vec.size()) {
        cout << " !: warning: requested element " << i << " in ioOpt with vector size " << vec.size() << endl;
        cout << " returning vec[0] instead" << endl;
        return vec[0];
    } else {
        return vec[i];
    }
};

void ioOpt::push_back(double val) {
    my_val = val; // my_val is therefore always the last value in
    vec.push_back(val);
};

ioOptMap::ioOptMap(ioOpt _)         { add(_); };
ioOptMap::ioOptMap(vector<ioOpt> _) { add(_); };
ioOptMap::ioOptMap()                {};

ioOpt& ioOptMap::operator[](string _) { return dict[_]; }; //if (dict.count(_)!=0) return dict[_]; else return null_opt; };
bool ioOptMap::has(string _) { return dict.count(_) != 0; };
bool ioOptMap::operator()(string _) { return has(_); };

ostream& operator<< (ostream& os, const ioOptMap& _) {
        for (auto& v : _.dict)  os << Form("  %-12s -> ",v.first.c_str()) << v.second << endl;
        return os;
    };

ioOptMap ioOptMap::operator+=(ioOptMap _) { 
    for (auto& key : _.dict) dict[key.first] = key.second; return *this; 
};
ioOptMap operator+ (ioOptMap lhs, const ioOptMap& rhs) { lhs += rhs; return lhs; };

void ioOptMap::add(ioOpt _) {
    // only acceptable syntax is string with "key:value;;[key2:value2;; ...]"
    if (!_.is_str) { 
        cout << " !: warning in ioOptMap.add: bad syntax to add option: " << _ 
             << endl << "    skipping entry" << endl;
        return;
    }
    
    // parse the string
    TString s1 { _.str() };
    /* vector<pair<TString,TString>> rval; */
    TPRegexp r1("\\s*([\\w-]+):\\s*(.*?)\\s*;;(.*)");
    bool has_first = false;
    while (s1(r1) != "") {
        TObjArray* subStrL = r1.MatchS(s1);
        const Int_t nrSubStr = subStrL->GetLast()+1;
        /* TString key, value; */
        if (nrSubStr > 2) {
            string  key = ((TObjString *)subStrL->At(1))->GetString().Data();
            TString value  = ((TObjString *)subStrL->At(2))->GetString().Data();

            if (value.IsFloat()) dict[key] = value.Atof();
            else                 dict[key] = value.Data();
            has_first = true;
            /* add(key,value); */
            /* rval.push_back({key,value}); */
        } else {
            break;
        }
        if (nrSubStr > 3) {
            s1 = ((TObjString*)subStrL->At(3))->GetString();
        }
    }
    if (has_first) { // check for a trailing value of key:value (without key:value;;)
        TPRegexp strip("\\s*([\\w-]+):\\s*(\\S+(\\s*\\S+)*)");
        /* TPRegexp strip("\\s*([\\w-]+):\\s*([^\\s]+.*[^\\s]+)\\s*"); */
        /* TPRegexp strip("\\s*([\\w-]+):\\s*([^\\s]+.*[^\\s]+)\\s*"); */
        if (s1(strip) != "") {
            cout << " in : " << endl;
            TObjArray* subStrL = strip.MatchS(s1);
            const Int_t nrSubStr = subStrL->GetLast()+1;
            if (nrSubStr > 2) {
                cout << "0:"<<  ((TObjString *)subStrL->At(0))->GetString().Data() << endl;
                cout << "1:"<<  ((TObjString *)subStrL->At(1))->GetString().Data() << endl;
                cout << "2:"<<  ((TObjString *)subStrL->At(2))->GetString().Data() << endl;

                string  key    = ((TObjString *)subStrL->At(1))->GetString().Data();
                TString value  = ((TObjString *)subStrL->At(2))->GetString().Data();
                if (value.IsFloat()) dict[key] = value.Atof();
                else                 dict[key] = value.Data();
            } 
        }
    }

    return;
};

void ioOptMap::add(vector<ioOpt> args) {
    int n_args { (int) args.size() };
    int i {0};
    int n {n_args-i}; // n = number remaining
    while (i<n_args) {
        if (!args[i].is_str) {
            cout << " !: warning in ioOptMap: input argument in vector " << args[i] << endl
                 << "    cannot be a map key because it is not a string." << endl
                 << "    Terminating population of ioOptMap." << endl;
            return;
        }
        int n {n_args-i-1}; // number of arguments left
        if (args[i].str().find(";;") != string::npos) {
            // This is is a string map of "key:value;; ..."
            add(args[i]);
            ++i;
            continue;
        } //  args[i] is a key and must be followed by:
          //  * string
          //  * double (x1)
          //  * double double ... (>1)
          //  If string, use the next entry at the value
          //  If val (x1) enter it as a single value
          //  If val (>1) enter it as a vector
        /* cout << " a0:  " << args[i] << " " << n << endl; */
        if (n == 0) {
            cout << " !: warning in ioOptMap: input key " << args[i] << endl
                 << "    not followed by an entry. Entering as a blank key" << endl;
            dict[args[i].str()]={};
            ++i;
        } else if (n == 1) {
            dict[args[i].str()]=args[i+1];
            i += 2;
        } else {
            int i_lastval{i+1};
            while (i_lastval < n_args && !args[i_lastval].is_str) ++i_lastval;
            int n_vals { i_lastval - i - 1};
            if (n_vals < 2) {
                dict[args[i].str()] = args[i+1];
                i += 2;
            } else {
                for (int k{i+1}; k<i_lastval; ++k) {
                    dict[args[i]].push_back(args[k].val());
                }
                i += n_vals+1;
            }
        }

    }
};

