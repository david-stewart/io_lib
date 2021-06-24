#ifndef io_IOS__h
#define io_IOS__h

#include "TGraph.h"
    // these are classes, but their primary purpose is to be return structure data types
struct IOS_keepcut_stats {
    int  n_keep {0};
    long double sum_keep {0.};
    int  n_cut {0};
    long double sum_cut {0.};
    // stats on last point cut
    bool was_cut{false};
    double x{0.};
    double y{0.};

    // values for total
    vector<double> vec_x{}; // collect X and Y values
    vector<double> vec_y{};
    IOS_keepcut_stats& operator+=(const IOS_keepcut_stats& rhs);

    // functions
    TGraph* graph(); // make graph from vec_x and vec_y
    double n_cut_ratio();
    double sum_cut_ratio();

    friend ostream& operator<<(ostream&, IOS_keepcut_stats&);

};

#endif
