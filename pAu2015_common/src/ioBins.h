#ifndef ioBins__h
#define ioBins__h

// for binning edges
double* ioBins_jetpt_prelim_13bins(); // 14 edges for 13 bins
double* ioBins_BBCES_3bin(); // 3 bins, four edges
double* ioBins_BBCES_4bin(); // 3 bins, four edges

/* struct ioBins_BBCES_3bins_binner { */
/*     ioBins_BBCES_3bins_binner (); */
/*     int last_val; */
/*     int operator()(); */
/*     int operator(double ea); // convert EA to which bin */
/* } */

#endif
