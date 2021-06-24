#ifndef io_test__h
#define io_test__h

#include "io_enum.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include <iostream>
using namespace std;


// next two functions remove all bins with entries < min_entries



/*
// study
struct io_intIterBase {
    int iend;
    int index{0};
    io_intIterBase(int _iend) : iend{_iend} {};
    io_intIterBase(int _index, int _iend) : iend{_iend}, index{_index} {};
    bool operator!=(const io_intIterBase& rhs) { return index != rhs.iend; };
    void operator++() { ++index; };
};
struct io_iiter : public io_intIterBase {
    io_iiter(int _iend) : io_intIterBase{_iend} {};
    io_iiter(int _index, int _iend) :io_intIterBase{_index, _iend} {};
    io_iiter begin() { return io_iiter{index, iend}; };
    io_iiter end()   { return io_iiter{iend,  iend}; };
    int operator*() { return index; };
};
*/

#endif
