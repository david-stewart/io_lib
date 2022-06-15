// some comment
#include "THnSparse.h"
#include "TH1D.h"
#include "tuClass.h"
/* #include "fastjet/PseudoJet.hh" */
#include "TRandom3.h"

/* using fastjet::PseudoJet; */

class tuTrackSparse {
    // axes:
    // EVENT:
    // 0 : bbc     - 3 bins
    // 1 : TrigEt  - 3 bins
    // 2 : TrigEta - 3 bins
    // 3 : ZDCx    - 3 bins -- 4-10, 10-16, 16-22 kHz
    // 4 : Vz      - 10 bins (because, why not?)
    // TRACK:
    // 5 : pt 
    // 6 : abs_dphi - 3 bins -- 0 - trans, 1, same-size, 2, recoil 
    // 7 : eta      - 3 bins -- east, mid, west
    private:
    int nbins[8];
    double n_triggers{-1};

    public:
    bool scaleByBinWidth=true;

    double weight{1.};
    double hopper[8];
    tuBinVec* bins {nullptr};
    tuTrackSparse(const char* tag="", bool debug_print=false);
    tuTrackSparse(THnSparseD* _data_track, THnSparseD* _data_trig, THnSparseD* _data_PU, bool debug_print=false);
    void write();

    THnSparseD* data_trig;
    THnSparseD* data_PU;
    THnSparseD* data_track;

    bool debug_print {false};
    
    void fill_trig(double EAbbc, double TrigEt, double TrigEta, double ZDCx, double Vz, double rhoPU, 
            double rhoPU_east=0., double rhoPU_mid=0., double rhoPU_west=0.);
    void fill_pt_eta_dphi(double pt, double eta, double dphi=0.); 
    // note: will fill with last values in hopper from fill_trig;

    void range_axes     (int i_axis, int i0, int i1);
    void range_bbc      (int i0, int i1) { range_axes(0,i0,i1); };
    void range_TrigEt   (int i0, int i1) { range_axes(1,i0,i1); };
    void range_TrigEta  (int i0, int i1) { range_axes(2,i0,i1); };
    void range_ZDCx     (int i0, int i1) { range_axes(3,i0,i1); };
    void range_vz       (int i0, int i1) { range_axes(4,i0,i1); };
    void range_pt       (int i0, int i1) { range_axes(5,i0,i1); };
    void range_abs_dphi (int i0, int i1) { range_axes(6,i0,i1); };
    void range_eta      (int i0, int i1) { range_axes(7,i0,i1); };

    TH1D* hg_axis     (int i_axis, double norm=1., bool use_data_track=true);
    TH1D* hg_bbc      (double norm=1., bool isjets=true) { return hg_axis(0, norm,isjets); }; // if norm = 0, then use triggers
    TH1D* hg_TrigEt   (double norm=1., bool isjets=true) { return hg_axis(1, norm,isjets); };
    TH1D* hg_TrigEta  (double norm=1., bool isjets=true) { return hg_axis(2, norm,isjets); };
    TH1D* hg_ZDCx     (double norm=1., bool isjets=true) { return hg_axis(3, norm,isjets); };
    TH1D* hg_vz       (double norm=1., bool isjets=true) { return hg_axis(4, norm,isjets); };
    TH1D* hg_pt       (double norm=1.)                   { return hg_axis(5, norm); };
    TH1D* hg_abs_dphi (double norm=1.)                   { return hg_axis(6, norm); };
    TH1D* hg_eta      (double norm=1.)                   { return hg_axis(7, norm); };
    double get_n_triggers();
    double get_sum_PU(int i_bin=4);
};


