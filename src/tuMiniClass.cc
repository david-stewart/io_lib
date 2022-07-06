#include "tuMiniClass.h"

void tuMinMaxPtr::operator()(double _, void* p){
    if (!has_data) {
        has_data = true;
        min_ptr = p;
        max_ptr = p;
        min_val = _;
        max_val = _;
    } else {
        if (_ > max_val) {
            max_val = _;
            max_ptr = p;
        }
        if (_ < min_val) {
            min_val = _;
            min_ptr = p;
        }
    }
};

tuMinMax::tuMinMax(string _name) : name{_name} {};
long long int tuMinMax::fill(double val) {
    if (n_entries == 0) {
        min = val;
        max = val;
    } else {
        if (val < min) min = val;
        if (val > max) max = val;
    }
    ++n_entries;
    return n_entries;
};
tuMinMax::tuMinMax(double& ptr, string _name) : 
    name{_name}, fill_option{5}, ptr_double{&ptr} {};
tuMinMax::tuMinMax(int& ptr, string _name) : 
    name{_name}, fill_option{6}, ptr_int{&ptr} {};
tuMinMax::tuMinMax(unsigned int& ptr, string _name) : 
    name{_name}, fill_option{7}, ptr_uint{&ptr} {};
tuMinMax::tuMinMax(short& ptr, string _name) : 
    name{_name}, fill_option{8}, ptr_short{&ptr} {};
tuMinMax::tuMinMax(char& ptr, string _name) : 
    name{_name}, fill_option{9}, ptr_char{&ptr} {};

tuMinMax::tuMinMax(double* ptr, string _name) : 
    name{_name}, fill_option{10}, ptr_double{ptr} {};
tuMinMax::tuMinMax(int* ptr, string _name) : 
    name{_name}, fill_option{11}, ptr_int{ptr} {};
tuMinMax::tuMinMax(unsigned int* ptr, string _name) : 
    name{_name}, fill_option{12}, ptr_uint{ptr} {};
tuMinMax::tuMinMax(short* ptr, string _name) : 
    name{_name}, fill_option{13}, ptr_short{ptr} {};
tuMinMax::tuMinMax(char* ptr, string _name) : 
    name{_name}, fill_option{14}, ptr_char{ptr} {};

tuMinMax::tuMinMax(double* ptr, int& _index, string _name) : 
    name{_name}, fill_option{0}, ptr_double{ptr}, index{&_index} {};
tuMinMax::tuMinMax(int* ptr, int& _index, string _name) : 
    name{_name}, fill_option{1}, ptr_int{ptr}, index{&_index} {};
tuMinMax::tuMinMax(unsigned int* ptr, int& _index, string _name) : 
    name{_name}, fill_option{2}, ptr_uint{ptr}, index{&_index} {};
tuMinMax::tuMinMax(short* ptr, int& _index, string _name) : 
    name{_name}, fill_option{3}, ptr_short{ptr}, index{&_index} {};
tuMinMax::tuMinMax(char* ptr, int& _index, string _name) : 
    name{_name}, fill_option{4}, ptr_char{ptr}, index{&_index} {};

long long int tuMinMax::operator()() {
    switch (fill_option) {
        case 0:
            return fill(ptr_double[*index]);
        case 1:
            return fill(ptr_int[*index]);
        case 2:
            return fill(ptr_uint[*index]);
        case 3:
            return fill(ptr_short[*index]);
        case 4:
            return fill(ptr_char[*index]);

        case 5:
            return fill(*ptr_double);
        case 6:
            return fill(*ptr_int);
        case 7:
            return fill(*ptr_uint);
        case 8:
            return fill(*ptr_short);
        case 9:
            return fill(*ptr_char);
        default:
            throw std::runtime_error("tuMinMax::operator()() called with no pointer set");
    }
};
long long int tuMinMax::operator()(int index) {
    switch (fill_option) {
        case 10:
            return fill(ptr_double[index]);
        case 11:
            return fill(ptr_int[index]);
        case 12:
            return fill(ptr_uint[index]);
        case 13:
            return fill(ptr_short[index]);
        case 14:
            return fill(ptr_char[index]);
        default:
            throw std::runtime_error("tuMinMax::operator()(int) called with no pointer set");
    }
};
int tuMinMax::nbins() { return (int)(max-min)+1; };
ostream& operator<<(ostream& os, tuMinMax& self) {
    if (self.name != "") cout << self.name << ": ";
    os << self.min << " " << self.max;
    return os;
};
