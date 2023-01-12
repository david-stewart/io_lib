class tuTrSparsely {
    // axes:
    // EVENT:
    // 0 : bbc     - 3 bins
    // 1 : TrigPt  - 3 bins  -- 4-8, 8-12, 12-30
    // 4 : ZDCx    - 6 bins  -- 4-30 kHz
    // 5 : nGlob   - 6 : 300-1600
    // 5 : runId   - 3 bins  -- 1, 260, 780
    // TRACK:
    // 6 : pt       - 18 : .2, 2.0
    // 7 : abs_dphi - 2  bins -- 0 - trans, 1 not-transe
    // 8 : nhit poss- 45 bins
    // 9 : nhit fit - 45 bins
    // 10 : dca     - 0, 1, 3,
    
    private:
    int nbins[11];
    double n_triggers{-1};
    bool scaleByBinWidth=true;

    public:
    double weight{1.};
    double hopper[11];
    double trig_phi;
    tuBinVec* bins {nullptr};
    tuTrSparsely(const char* tag="", bool debug_print=false);
    tuTrSparsely(THnSparseD* _data_jet, THnSparseD* data_trig, bool debug_print=false);
    void write();

    THnSparseD* data_trig;
    THnSparseD* data_track;

    bool debug_print {false};
    
    void fill_trig (double EAbbc, double TrigEt, double TrigEta, double trigPhi, double ZDCx, double runId);
    void fill_track(double pt, double eta, double phi, double dca); 
    // note: will fill with last values in hopper from fill_trig;

    void range_axes       (int i_axis, int i0, int i1);
    void range_axes_float (int i_axis, double f0, double f1);
    void range_bbc        (int i0, int i1) { range_axes(0,i0,i1); };

    void range_EAbbc    (int i0, int i1) { range_axes(0,i0,i1); };
    void range_TrigEt   (int i0, int i1) { range_axes(1,i0,i1); };
    void range_TrigEta  (int i0, int i1) { range_axes(2,i0,i1); };
    void range_TrigPhi  (int i0, int i1) { range_axes(3,i0,i1); };
    void range_ZDCx     (int i0, int i1) { range_axes(4,i0,i1); };
    void range_runId    (int i0, int i1) { range_axes(5,i0,i1); };

    void range_pt       (int i0, int i1) { range_axes(6, i0,i1); };
    void range_abs_dphi (int i0, int i1) { range_axes(7, i0,i1); };
    void range_phi      (int i0, int i1) { range_axes(8, i0,i1); };
    void range_eta      (int i0, int i1) { range_axes(9, i0,i1); };
    void range_dca      (int i0, int i1) { range_axes(10, i0,i1); };

    TH1D* hg_axis     (int i_axis,  double norm=0., bool is_trig=true);
    TH2D* hg2_axis    (int i0_axis, int i1_axis, double norm=0., bool is_trigger=true);

    TH1D* hg_EAbbc    (double norm=1., bool is_trig=true) { return hg_axis(0, norm, is_trig); };
    TH1D* hg_TrigEt   (double norm=1., bool is_trig=true) { return hg_axis(1, norm, is_trig); };
    TH1D* hg_TrigEta  (double norm=1., bool is_trig=true) { return hg_axis(2, norm, is_trig); };
    TH1D* hg_TrigPhi  (double norm=1., bool is_trig=true) { return hg_axis(3, norm, is_trig); };
    TH1D* hg_ZDCx     (double norm=1., bool is_trig=true) { return hg_axis(4, norm, is_trig); };
    TH1D* hg_runId    (double norm=1., bool is_trig=true) { return hg_axis(5, norm, is_trig); };

    TH1D* hg_pt       (double norm=1.) { return hg_axis(6, norm, false); };
    TH1D* hg_abs_dphi (double norm=1.) { return hg_axis(7, norm, false); };
    TH1D* hg_phi      (double norm=1.) { return hg_axis(8, norm, false); };
    TH1D* hg_eta      (double norm=1.) { return hg_axis(9, norm, false); };
    TH1D* hg_dca      (double norm=1.) { return hg_axis(10, norm, false); };

    double get_n_triggers();
};

