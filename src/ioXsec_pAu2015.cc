#include "ioXsec_pAu2015.h"
#include "TH1D.h"
#include <iostream>

using std::cout;
using std::endl;


ioXsec_pAu2015::ioXsec_pAu2015() : 
    n_collected{0,0,0,0,0,0,0,0,0},
    n_collected_total {0} 
{};

unsigned int ioXsec_pAu2015::pthatbin(int pthat_min, int pthat_max){
	int pt_bin { pthat_min * 1000 + pthat_max };
	switch ( pthat_min * 1000 + pthat_max ) {
		case 5007: return 0;
		case 7009: return 1;
		case 9011: return 2;
		case 11015: return 3;
		case 15025: return 4;
		case 25035: return 5;
		case 35045: return 6;
		case 45055: return 7;
		case 55065: return 8;
		default:
			throw std::runtime_error(Form("Not a valid pthat range %i-%i",pthat_min,pthat_max));
			return 900;
	}
};

// Return the Xsection according to Pythia8/6
//   divided by number_of_events.
//   If number_of_events == 0:
//      if collected events, use collected events
//      else use default numbers (from my local trees)
/* double ioXsec_pAu2015::XsecPyth8(int pthat_min, int pthat_max, int number_of_events=0){ */
/*     return XsecPyth8(pthatbin(pthat_min, pthat_max), number_of_events); */
/* }; */

double ioXsec_pAu2015::XsecPyth6(pair<int,int> ptrange, int number_of_events){
    return XsecPyth6(pthatbin(ptrange.first, ptrange.second), number_of_events);
};
double ioXsec_pAu2015::XsecPyth8(pair<int,int> ptrange, int number_of_events){
    return XsecPyth8(pthatbin(ptrange.first, ptrange.second), number_of_events);
};

double ioXsec_pAu2015::XsecPyth8(int pthat_bin, int number_of_events){
	const static double Xsec[9] {
		0.1604997783173907, 0.0279900193730690,     0.006924398431,
	    0.0028726057079642, 0.0005197051748372, 0.0000140447879818,
	    0.0000006505378525, 0.0000000345848665, 0.0000000016149182 };
    const static int nEvents[9] { // for my trees
        375735, 217715, 110680,
        168886, 518397, 177438,
         59805,  60212,  60270
    };
    if (number_of_events) return Xsec[pthat_bin] / number_of_events;
    else if (n_collected_total) {
        if (!n_collected[pthat_bin]) {
            cout << " fatal error: no events for bin " << pthat_bin 
                 << " have been collected." << endl
                 << " Cannot generate weighted cross section." << endl;
			throw std::runtime_error("Asking for ptbin with no events");
        }
        return Xsec[pthat_bin] / n_collected[pthat_bin];
    } else {
        return Xsec[pthat_bin] / nEvents[pthat_bin];
    }
};

double ioXsec_pAu2015::XsecPyth6(int pthat_bin, int number_of_events){
	const static double Xsec[9] {
        // these versions were pulled from the actual embedding files
        0.107509,   0.0190967,  0.00475202,
      0.00198812, 0.000361282, 9.65463E-06,
     4.71077E-07, 2.68464E-08, 1.38211E-09 };
    // The following values were generated independently running Pythia6
		/* 0.10849815,      0.019171747,  0.0047063602, */
	    /* 0.0020085739,  0.00035981473, 9.6174432e-06, */
		/* 4.722213e-07,  2.6900098e-08, 1.3887308e-09 }; */
    const static int nEvents[9] { // for my trees
        375735, 217715, 110680,
        168886, 518397, 177438,
         59805,  60212,  60270
    };
    if (number_of_events) return Xsec[pthat_bin] / number_of_events;
    else if (n_collected_total) {
        if (!n_collected[pthat_bin]) {
            cout << " fatal error: no events for bin " << pthat_bin 
                 << " have been collected." << endl
                 << " Cannot generate weighted cross section." << endl;
			throw std::runtime_error("Asking for ptbin with no events");
        }
        return Xsec[pthat_bin] / n_collected[pthat_bin];
    } else {
        return Xsec[pthat_bin] / nEvents[pthat_bin];
    }
};

void   ioXsec_pAu2015::collect(int pthat_min, int pthat_max){
    collect(pthatbin(pthat_min, pthat_max));
};
void   ioXsec_pAu2015::collect(int pthat_bin){
    ++n_collected[pthat_bin];
    ++n_collected_total;
    if (!hg_collected){
        double* edges = new double[10]{5.,7.,9.,11.,15.,25.,35.,45.,55.,65.};
        hg_collected = new TH1D("pthat_bins_collected",
            "Number of events collected in each pthat bin;#hat{p}_{T};n-collected",
            9, edges);
    }
    static const double x_in_bin[9]{ 5.5,7.5,9.5,11.5,15.5,25.5,35.5,45.5,55.5 };
    hg_collected->Fill(x_in_bin[pthat_bin]);
};
