#ifndef ioXsec_pAu2015__h
#define ioXsec_pAu2015__h

#include "TH1D.h"

struct ioXsec_pAu2015 {
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
    double XsecPyth8(int pthat_bin, int number_of_events=0);
    /* double XsecPyth6(int pthat_min, int pthat_max, int number_of_events=0); */
    double XsecPyth6(int pthat_bin, int number_of_events=0);
    void collect(int pthat_bin);
    void collect(int pthat_min, int pthat_max);
    long int n_collected[9];
    long int n_collected_total;
    TH1D* hg_collected;

    ioXsec_pAu2015();
};


#endif
