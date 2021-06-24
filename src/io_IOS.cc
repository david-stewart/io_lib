#include "io_IOS.h"
#include "io_fnc.h"
#include <iostream>

IOS_keepcut_stats& IOS_keepcut_stats::operator+=(const IOS_keepcut_stats& rhs) {
    n_cut    += rhs.n_cut;
    sum_cut  += rhs.sum_cut;
    n_keep   += rhs.n_keep;
    sum_keep += rhs.sum_keep;

    was_cut   = rhs.was_cut;
    if (was_cut) {
        vec_x.push_back(rhs.x);
        vec_y.push_back(rhs.y);
    }
    return *this;
};
TGraph* IOS_keepcut_stats::graph() { return ioMakeTGraph(vec_x,vec_y); };
double IOS_keepcut_stats::n_cut_ratio()   { return (double)n_cut/(n_cut+n_keep); };
double IOS_keepcut_stats::sum_cut_ratio() { return (double)(sum_cut/(sum_cut+sum_keep)); };

ostream& operator<<(ostream& os, IOS_keepcut_stats& obj) {
    os << Form("nKeep:%6i  sum_keep:%8.4f  nCut:%6i  sum_cut:%8.4f Nrat_cut:%5.2f  sum_rat_cut:%5.2f  hvec size: %4i",
            obj.n_keep, (double)obj.sum_keep, 
            obj.n_cut,  (double)obj.sum_cut, 
            obj.n_cut_ratio(), obj.sum_cut_ratio(),
            (int)obj.vec_x.size())
             << endl;
    return os;
};
