#ifndef ioTHnSparse__h
#define ioTHnSparse__h
// some comment
#include "THnSparse.h"
#include "TH1D.h"
#include "ioClass.h"
#include "fastjet/PseudoJet.hh"
#include "TRandom3.h"

using fastjet::PseudoJet;

struct ioSparseBinVecSet {
    ioBinVec EAbbc   {{0., 0., 10., 64000.}};
    ioBinVec EAtpc   {{0., 0., 20., 2.    }};
    ioBinVec TrigEt  {{0., 4., 8., 30.    }};
    ioBinVec vz      {{-10., -10., 20, 10. }};
    ioBinVec ZDCx    {{0., 0., 10., 50000. }};
    ioBinVec absDphi {{0., 0., 12, M_PI}};
    ioBinVec JetPt   {{0., 0., 70, 70 }};

    ioBinVec leadPt  {{0., 0., 70, 70  }};
    ioBinVec matchPt {{0., 0., 70, 70  }};
    ioBinVec AJ      {{0., 0., 200, 1. }};
    ioSparseBinVecSet () {};
    friend ostream& operator<<(ostream& os, ioSparseBinVecSet& val);
};
ostream& operator <<(ostream& os, ioSparseBinVecSet& bins);

class ioTrackSparse {
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
    bool scaleByBinWidth=true;

    public:
    double weight{1.};
    double hopper[8];
    ioBinVec* bins {nullptr};
    ioTrackSparse(const char* tag="", bool debug_print=false);
    ioTrackSparse(THnSparseD* _data_jet, THnSparseD* data_trig, THnSparseD* data_PU, bool debug_print=false);
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

    TH1D* hg_axis  (int i_axis, double norm=0., bool use_jet_data=true);
    TH1D* hg_bbc      (double norm=0., bool isjets=true) { return hg_axis(0, norm,isjets); }; // if norm = 0, then use triggers
    TH1D* hg_TrigEt   (double norm=0., bool isjets=true) { return hg_axis(1, norm,isjets); };
    TH1D* hg_TrigEta  (double norm=0., bool isjets=true) { return hg_axis(2, norm,isjets); };
    TH1D* hg_ZDCx     (double norm=0., bool isjets=true) { return hg_axis(3, norm,isjets); };
    TH1D* hg_vz       (double norm=0., bool isjets=true) { return hg_axis(4, norm,isjets); };
    TH1D* hg_pt       (double norm=0.)                   { return hg_axis(5, norm); };
    TH1D* hg_abs_dphi (double norm=0.)                   { return hg_axis(6, norm); };
    TH1D* hg_eta      (double norm=0.)                   { return hg_axis(7, norm); };
    double get_n_triggers();
};


class ioJetSpectraSparse {
    private:
    int nbins[7];
    double n_triggers{-1};

    public:
    bool scaleByBinWidth=true;
    double weight{1.};
    double hopper[7];
    ioBinVec* bins {nullptr};
    ioJetSpectraSparse(ioSparseBinVecSet& _inp, const char* tag="",bool debug_print=false);
    ioJetSpectraSparse(THnSparseD* _data_jet, THnSparseD* data_trig, 
            bool debug_print=false);
    void write();

    THnSparseD* data_trig;
    THnSparseD* data_jet;
    bool debug_print {false};
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

    void range_axes    (int i_axis, int i0, int i1);
    void range_EAbbc   (int i0, int i1) { range_axes(0,i0,i1); };
    void range_EAtpc   (int i0, int i1) { range_axes(1,i0,i1); };
    void range_TrigEt  (int i0, int i1) { range_axes(2,i0,i1); };
    void range_ZDCx    (int i0, int i1) { range_axes(3,i0,i1); };
    void range_Vz      (int i0, int i1) { range_axes(4,i0,i1); };
    void range_JetPt   (int i0, int i1) { range_axes(5,i0,i1); };
    void range_absDphi (int i0, int i1) { range_axes(6,i0,i1); };
    void range8_absDphi(int i) { range_axes (6, 8*(i-1)+1, 8*i); };

    TH1D* hg_axis  (int i_axis, double norm=0., bool use_jet_data=true);
    TH1D* hg_EAbbc (double norm=0., bool isjets=true) { return hg_axis(0, norm,isjets); }; // if norm = 0, then use triggers
    TH1D* hg_EAtpc (double norm=0., bool isjets=true) { return hg_axis(1, norm,isjets); };
    TH1D* hg_TrigEt(double norm=0., bool isjets=true) { return hg_axis(2, norm,isjets); };
    TH1D* hg_ZDCx  (double norm=0., bool isjets=true) { return hg_axis(3, norm,isjets); };
    TH1D* hg_Vz    (double norm=0., bool isjets=true) { return hg_axis(4, norm,isjets); };
    TH1D* hg_JetPt (double norm=0.) { return hg_axis(5, norm); };
    TH1D* hg_JetPt64 (int i0, int i1, double norm=0.);
    TH1D* hg_JetPt8  (int i0, double norm=0.);
    TH1D* hg_absDphi (double norm=0.) { return hg_axis(6, norm); };
    double get_n_triggers();
};


class ioAjSparse {
    private:
    int nbins[8];
    double weight{1.};
    double hopper[8];
    double dphi_min{0};

    // add a possibility to re-weight the ioAjSparse with a correction to the jets
    double mean_pTCorr { 1. };
    double sig_pTCorr  { 0. };
    bool   flag_pTCorr { false };
    TRandom3* pTrand { nullptr };
    double jet_area { 0.50265482457436691815402294132472 }; // assume R=0.4 jets

    public:
    void set_pTCorr ( double _mean_pTCorr=1., double _sig_pTCorr=0. );

    bool scaleByBinWidth=true;
    ioBinVec* bins {nullptr};
    ioAjSparse(ioSparseBinVecSet& bins, const char* tag="", double _recoil_match=0.4);
    ioAjSparse(THnSparseD* _data);
    double recoil_match {0.4};
    void write();

    THnSparseD* data;
    double Aj();
    
    // axes:
    // 0 : EAbbc  - X bins
    // 1 : EAtpc  - X bins
    // 2 : TrigEt - X bins
    // 3 : ZDCx   - X bins
    // 4 : Vz     - X bins
    // 6 : LeadPt
    // 7 : MatchPt
    // 8 : AJ
    
    bool fill(double EAbbc, double EAtpc, double TrigEt, double ZDCx, 
              double Vz, double leadPt, double matchPt);
    bool fill(double* _in, vector<PseudoJet>& jets);

    void range_axes    (int i_axis, int i0, int i1);
    void range_EAbbc   (int i0, int i1) { range_axes(0,i0,i1); };
    void range_EAtpc   (int i0, int i1) { range_axes(1,i0,i1); };
    void range_TrigEt  (int i0, int i1) { range_axes(2,i0,i1); };
    void range_ZDCx    (int i0, int i1) { range_axes(3,i0,i1); };
    void range_Vz      (int i0, int i1) { range_axes(4,i0,i1); };
    void range_leadPt  (int i0, int i1) { range_axes(5,i0,i1); };
    void range_matchPt (int i0, int i1) { range_axes(6,i0,i1); };
    void range_AJ      (int i0, int i1) { range_axes(7,i0,i1); };

    TH1D* hg_axis   (int i_axis, double norm=0.);
    TH1D* hg_EAbbc  (double norm=0.) { return hg_axis(0, norm); }; // if norm = 0, then use triggers
    TH1D* hg_EAtpc  (double norm=0.) { return hg_axis(1, norm); };
    TH1D* hg_TrigEt (double norm=0.) { return hg_axis(2, norm); };
    TH1D* hg_ZDCx   (double norm=0.) { return hg_axis(3, norm); };
    TH1D* hg_Vz     (double norm=0.) { return hg_axis(4, norm); };
    TH1D* hg_leadPt (double norm=0.) { return hg_axis(5, norm); };
    TH1D* hg_matchPt(double norm=0.) { return hg_axis(6, norm); };
    TH1D* hg_AJ     (double norm=0.) { return hg_axis(7, norm); };
};

/* class ioAjSparse_dPhi { */
/*     private: */
/*     int    nbins[9]; */
/*     double weight{1.}; */
/*     double hopper[9]; */
/*     double dphi_min{0}; */

/*     // add a possibility to re-weight the ioAjSparse_dPhi with a correction to the jets */
/*     double mean_pTCorr { 1. }; */
/*     double sig_pTCorr  { 0. }; */
/*     bool   flag_pTCorr { false }; */
/*     TRandom3* pTrand { nullptr }; */
/*     double jet_area { 0.50265482457436691815402294132472 }; // assume R=0.4 jets */

/*     public: */
/*     void set_pTCorr ( double _mean_pTCorr=1., double _sig_pTCorr=0. ); */

/*     bool scaleByBinWidth=true; */
/*     ioBinVec* bins {nullptr}; */
/*     ioAjSparse_dPhi(const char* bin_file, const char* tag="", double _recoil_match=IO_halfpi); */
/*     ioAjSparse_dPhi(THnSparseD* _data); */
/*     double recoil_match {0.4}; */
/*     void write(); */

/*     THnSparseD* data; */
/*     double Aj(); */
    
/*     // axes: */
/*     // 0 : EAbbc  - X bins */
/*     // 1 : EAtpc  - X bins */
/*     // 2 : TrigEt - X bins */
/*     // 3 : ZDCx   - X bins */
/*     // 4 : Vz     - X bins */
/*     // 5 : LeadPt */
/*     // 6 : MatchPt */
/*     // 7 : AJ */
/*     // 8 : dPhi */
    
/*     bool fill(double EAbbc, double EAtpc, double TrigEt, double ZDCx, */ 
/*               double Vz, double leadPt, double matchPt, double dPhi); */
/*     bool fill(double* _in, vector<PseudoJet>& jets); */

/*     void range_axes    (int i_axis, int i0, int i1); */
/*     void range_EAbbc   (int i0, int i1) { range_axes(0,i0,i1); }; */
/*     void range_EAtpc   (int i0, int i1) { range_axes(1,i0,i1); }; */
/*     void range_TrigEt  (int i0, int i1) { range_axes(2,i0,i1); }; */
/*     void range_ZDCx    (int i0, int i1) { range_axes(3,i0,i1); }; */
/*     void range_Vz      (int i0, int i1) { range_axes(4,i0,i1); }; */
/*     void range_leadPt  (int i0, int i1) { range_axes(5,i0,i1); }; */
/*     void range_matchPt (int i0, int i1) { range_axes(6,i0,i1); }; */
/*     void range_AJ      (int i0, int i1) { range_axes(7,i0,i1); }; */
/*     void range_dPhi    (int i0, int i1) { range_axes(8,i0,i1); }; */

/*     TH1D* hg_axis   (int i_axis, double norm=0.); */
/*     TH1D* hg_EAbbc  (double norm=0.) { return hg_axis(0, norm); }; // if norm = 0, then use triggers */
/*     TH1D* hg_EAtpc  (double norm=0.) { return hg_axis(1, norm); }; */
/*     TH1D* hg_TrigEt (double norm=0.) { return hg_axis(2, norm); }; */
/*     TH1D* hg_ZDCx   (double norm=0.) { return hg_axis(3, norm); }; */
/*     TH1D* hg_Vz     (double norm=0.) { return hg_axis(4, norm); }; */
/*     TH1D* hg_leadPt (double norm=0.) { return hg_axis(5, norm); }; */
/*     TH1D* hg_matchPt(double norm=0.) { return hg_axis(6, norm); }; */
/*     TH1D* hg_AJ     (double norm=0.) { return hg_axis(7, norm); }; */
/*     TH1D* hg_dPhi   (double norm=0.) { return hg_axis(8, norm); }; */

/*     TH2D* hg2_axis  (int xAxis, int yAxis, double norm=0); */
/* }; */
#endif
