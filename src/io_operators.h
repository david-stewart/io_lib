#ifndef io_operators__h
#define io_operators__h
#include <vector>
#include <ostream>
using std::vector;
using std::ostream;

vector<bool> operator&&(const vector<bool>&, const vector<bool>&);
vector<bool> operator||(const vector<bool>&, const vector<bool>&);

vector<bool> operator!(const vector<bool>&);
/* vector<bool> operator!(const vector<bool>); */

ostream& operator<<(ostream& os, const vector<bool>&);
ostream& operator<<(ostream& os, const vector<int>&);



#endif
