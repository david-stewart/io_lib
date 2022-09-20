#ifndef ioXsec_pAu2015__h
#define ioXsec_pAu2015__h

#include "TH1D.h"

struct ioXsec_pAu2015 {
    static constexpr double f_n10xPyth6[9] {
        3920155 , 2101168 , 1187058 , 1695141 , 
        4967075 , 1797387 , 260676 , 261926 , 262366
    };
    static constexpr double f_XsecPyth6[9] {
        0.107509,   0.0190967,  0.00475202,
        0.00198812, 0.000361282, 9.65463E-06,
        4.71077E-07, 2.68464E-08, 1.38211E-09
    };
    static double w10xPyth6(int i) {
        return f_XsecPyth6[i] / f_n10xPyth6[i];
    };
    // There are 9 pthat bins:
    // 0: 5-7
    // 1: 7-9
    // 2: 9-11
    // 3: 11-15
    // 4: 15-25
    // 5: 25-35
    // 6: 35-45
    // 7: 45-55
    // 8: 55-65

    // return the bin (0-9) of above.
    // fail if not one of the above sets
    unsigned int pthatbin(int pthat_min, int pthat_max);

    // Return the Xsection according to Pythia8/6
    //   divided by number_of_events.
    //   If number_of_events == 0:
    //      if collected events, use collected events
    //      else use default numbers (from my local trees)
    /* double XsecPyth8(int pthat_min, int pthat_max, int number_of_events=0); */
    double XsecPyth8(pair<int,int>, int number_of_events=0);
    double XsecPyth8(int pthat_bin, int number_of_events=0);
    double XsecPyth6(pair<int,int>, int number_of_events=0);
    double XsecPyth6(int pthat_bin, int number_of_events=0);
    void collect(int pthat_bin);
    void collect(int pthat_min, int pthat_max);
    long int n_collected[9];
    long int n_collected_total;
    TH1D* hg_collected;

    ioXsec_pAu2015();
};


#endif
