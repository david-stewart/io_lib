#include "THnSparse.h"
#include "TH1D.h"

class ioJetSpectraSparse {
    private:
    int nbins[6];
    double hopper[6];

    public:
    ioJetSpectraSparse(const char* bin_file, const char* tag);
    void write();

    THnSparseD* data_trig;
    THnSparseD* data_jet;
    // axes:
    // 0 : EAbbc  - 3 bins
    // 1 : EAtpc  - 3 bins
    // 2 : TrigEt - 30 bins
    // 3 : ZDCx   - 3 bins
    // 4 : Vz     - 3 bins
    // 5 : JetPt  - 70 bins // for data_jet only -- triggers have one less axis
    
    void fill_trig(double EAbbc, double EAtpc, double TrigEt, double ZDCx, double Vz);
    void fill_jetpt(double pt); // note: will fill with last values in hopper from fill_trig;

    void range_EAbbc (int, int);
    void range_EAtpc (int, int);
    void range_TrigEt(int,int);
    void range_ZDCx  (int,int);
    void range_Vz    (int,int);
    void range_JetPt (int,int);

    TH1D* hg_EAbbc ();
    TH1D* hg_EAtpc ();
    TH1D* hg_TrigEt();
    TH1D* hg_ZDCx  ();
    TH1D* hg_Vz    ();
    TH1D* hg_JetPt ();
    double n_triggers();
};
