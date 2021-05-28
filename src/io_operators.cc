#include "io_operators.h"
vector<bool> operator&&(const vector<bool>& A, const vector<bool>& B) {
    if (A.size() != B.size()) 
        throw std::runtime_error("Fatal error: trying to && two vectors of different sizes");
    vector<bool> r_vec;
    for (unsigned int i{0}; i<A.size(); ++i) r_vec.push_back(A[i] && B[i]);
    return r_vec;
};
vector<bool> operator||(const vector<bool>& A, const vector<bool>& B) {
    if (A.size() != B.size()) 
        throw std::runtime_error("Fatal error: trying to || two vectors of different sizes");
    vector<bool> r_vec;
    for (unsigned int i{0}; i<A.size(); ++i) r_vec.push_back(A[i] || B[i]);
    return r_vec;
};
ostream& operator<<(ostream& os, const vector<bool>& vec) { 
    for (auto v : vec) os << " " << v;
    return os;
}
ostream& operator<<(ostream& os, const vector<int>& vec) { 
    for (auto v : vec) os << " " << v;
    return os;
}

vector<bool> operator!(const vector<bool>& vec) {
    vector<bool> r_val;
    for (auto v : vec) r_val.push_back(!v);
    return r_val;
};
/* vector<bool> operator!(const vector<bool>  vec) { */
/*     vector<bool> r_val; */
/*     for (auto v : vec) r_val.push_back(!v); */
/*     return r_val; */
/* }; */
