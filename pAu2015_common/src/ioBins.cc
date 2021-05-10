#include "ioBins.h"

double* ioBins_prelim_jetpt_13bins(){
    return  new double[14]{0.,1.,2.,3.,4.,5.,6.,8.,10.,15.,20.,25.,36.,50.};
};
double* ioBins_BBCES_3bin(){
    return new double[4]{2750.,8292.,24159.,57433};
};
double* ioBins_BBCES_4bin(){
    return new double[5]{0., 2750.,8292.,24159.,57433};
};
