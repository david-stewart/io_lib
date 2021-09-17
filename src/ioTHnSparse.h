#include "THnSparse.h"
#include "TH1D.h"
#include "ioClass.h"

class ioJetSpectraSparse {
    private:
    int nbins[7];
    double hopper[7];
    double n_triggers{-1};

    public:
    ioBinVec* bins {nullptr};
    ioJetSpectraSparse(const char* bin_file, const char* tag="");
    ioJetSpectraSparse(THnSparseD* _data_jet, THnSparseD* data_trig);
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
    // 6 : JetAbsDeltaPhi - 64 bins
    
    void fill_trig(double EAbbc, double EAtpc, double TrigEt, double ZDCx, double Vz);
    void fill_jetpt_absDphi(double pt, double absdeltaphi); // note: will fill with last values in hopper from fill_trig;

    void range_axes (int i_axis, int i0, int i1);
    void range_EAbbc (int i0, int i1) { range_axes(0,i0,i1); };
    void range_EAtpc (int i0, int i1) { range_axes(1,i0,i1); };
    void range_TrigEt(int i0, int i1) { range_axes(2,i0,i1); };
    void range_ZDCx  (int i0, int i1) { range_axes(3,i0,i1); };
    void range_Vz    (int i0, int i1) { range_axes(4,i0,i1); };
    void range_JetPt  (int i0, int i1) { range_axes(5,i0,i1); };
    void range_absDphi(int i0, int i1) { range_axes(6,i0,i1); };
    void range8_absDphi(int i) { range_axes (6, 8*(i-1)+1, 8*i); };

    TH1D* hg_axis  (int i_axis, double norm=0.);
    TH1D* hg_EAbbc (double norm=0.) { return hg_axis(0, norm); }; // if norm = 0, then use triggers
    TH1D* hg_EAtpc (double norm=0.) { return hg_axis(1, norm); };
    TH1D* hg_TrigEt(double norm=0.) { return hg_axis(2, norm); };
    TH1D* hg_ZDCx  (double norm=0.) { return hg_axis(3, norm); };
    TH1D* hg_Vz    (double norm=0.) { return hg_axis(4, norm); };
    TH1D* hg_JetPt (double norm=0.) { return hg_axis(5, norm); };
    TH1D* hg_JetPt64 (int i0, int i1, double norm=0.);
    TH1D* hg_JetPt8  (int i0, double norm=0.);
    TH1D* hg_absDphi (double norm=0.) { return hg_axis(6, norm); };
    double get_n_triggers();
};
