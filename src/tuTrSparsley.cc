tuTrSparsely::tuTrSparsely(
        THnSparseD* _data_track, THnSparseD* _data_trig, bool _dbprint
) :
    data_trig { _data_trig }, data_track { _data_track }, debug_print{_dbprint} 
{};

tuTrSparsely::tuTrSparsely( const char* tag, bool _dbprint) :
    debug_print{_dbprint} 
{
    tuBinVec bin_EAbbc3  {{ 0.0000,  8315.4606, 24292.8207, 64000.0000 }};
    tuBinVec bin_TrigEt  {{ 0.,4.,8.,12.,30.}};
    tuBinVec bin_TrigEta {{ -1.0, 1.0001 }};
    tuBinVec bin_TrigPhi {{ 0., 2*M_PI+0.0001 }};
    tuBinVec bin_eta     {{ -1.0, -0.3, 0.3, 1. }};
    tuBinVec bin_phi     {{ 0., 0., 60, 2*M_PI }};
    tuBinVec bin_ZDCx    {{ 4000., 4000., 18, 22000. }};
    tuBinVec bin_runId   {{ 16125035, 16142058.5, 16149001.5, 16159024 }};
    tuBinVec bin_dca     {{ 0., 1., 2., 3. }};

    const tuBinVec bin_trpt {{ // track pT bins
     0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,
     1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,
     3.0,  3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,  4.0,  4.1,  4.2,  4.3,  4.4,
     4.6,  4.8,  5.0,  5.2,  5.4,  5.6,  5.8,  6.0,  6.2,  6.4,  6.6,  6.8,  7.1,  7.4,  7.7,
     8.0,  8.3,  8.6,  9.0,  9.4,  9.8, 10.3, 10.8, 11.3, 11.9, 12.5, 13.2, 14.0, 15.0
    }};

    tuBinVec bin_abs_dphi {{ 0., M_PI/3., 2*M_PI/3., M_PI+0.00001 }};
    //tuBinVec bin_phi      {{ 0., 0., 60, 2*M_PI }};
    //tuBinVec bin_eta      {{ -1.0, -0.3, 0.3, 1. }};

    // trigger 
    nbins[0] = bin_EAbbc3;
    nbins[1] = bin_TrigEt;
    nbins[2] = bin_TrigEta;
    nbins[3] = bin_TrigPhi;
    nbins[4] = bin_ZDCx;
    nbins[5] = bin_runId;

    // track
    nbins[6] = bin_trpt;
    nbins[7] = bin_abs_dphi;
    nbins[8] = bin_phi;
    nbins[9] = bin_eta;
    nbins[10] = bin_dca;

    if (debug_print) {
        cout << " debug_print, nbins: " << endl;
        for (int i{0}; i<7; ++i) cout << " nbins["<<i<<"] " << nbins[i] << endl;
    }

    string trig_string= "bbc;E_{T};Trig #eta;Trig #phi;ZDCx;runId;";
    string word_trig = "triggers;";
           
    data_trig = new THnSparseD(Form("data_trig2%s",tag),
            (word_trig+trig_string).c_str(),
            6, nbins, NULL, NULL);
    data_trig->SetBinEdges(0,bin_EAbbc3);
    data_trig->SetBinEdges(1,bin_TrigEt);
    data_trig->SetBinEdges(2,bin_TrigEta);
    data_trig->SetBinEdges(3,bin_TrigPhi);
    data_trig->SetBinEdges(4,bin_ZDCx);
    data_trig->SetBinEdges(5,bin_runId);
    data_trig->Sumw2();

    string track_string= "track-pT;track #Delta#phi;track #phi;track #eta;dca;";
    string word_track ="track;";
    data_track = new THnSparseD(Form("data_track2%s",tag),
            (word_track+trig_string+track_string).c_str(),
            11, nbins, NULL, NULL);
    data_track->SetBinEdges(0,bin_EAbbc3);
    data_track->SetBinEdges(1,bin_TrigEt);
    data_track->SetBinEdges(2,bin_TrigEta);
    data_track->SetBinEdges(3,bin_TrigPhi);
    data_track->SetBinEdges(4,bin_ZDCx);
    data_track->SetBinEdges(5,bin_runId);
    data_track->SetBinEdges(6,bin_trpt);
    data_track->SetBinEdges(7,bin_abs_dphi);
    data_track->SetBinEdges(8,bin_phi);
    data_track->SetBinEdges(9,bin_eta);
    data_track->SetBinEdges(10,bin_dca);
    data_track->Sumw2();
};
void tuTrSparsely::write() { 
    if (debug_print) {
        cout << " entries for " << data_trig->GetName() << " trig-entries: " <<
            data_trig->GetEntries() << " jet-entries: " << data_track->GetEntries() << endl;
        auto t_proj = data_trig->Projection(0);
        t_proj->SetName("p_trip");
        cout << " First axes: (mean) trig " << t_proj->GetMean() << endl;
        delete t_proj;

        auto j_proj = data_track->Projection(0);
        j_proj->SetName("j_trip");
        cout << " First axes: (mean) jet  " << j_proj->GetMean() << endl;
        delete j_proj;
    }
    data_trig->Write();
    data_track->Write();
};
void tuTrSparsely::fill_trig (double EAbbc, double TrigEt, double TrigEta, double TrigPhi, double ZDCx, double runId) {
    hopper[0] = EAbbc;
    hopper[1] = TrigEt;
    hopper[2] = TrigEta;
    hopper[3] = TrigPhi;
    hopper[4] = ZDCx;
    hopper[5] = runId;
    data_trig->Fill(hopper,weight);
};
void tuTrSparsely::fill_track(double pt, double eta, double phi, double dca){
    hopper[6] = pt;
    hopper[7] = tu_absDphi(hopper[3],phi);
    hopper[8] = phi;
    hopper[9] = eta;
    hopper[10] = dca;
    data_track->Fill(hopper,weight);
};
void tuTrSparsely::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 9) throw std::runtime_error(
        Form("fatal: error in tuTrSparsely, axis(%i) not valid, but by <7",
        i_axis)
    );

    n_triggers = -1.;
    data_track ->GetAxis(i_axis)->SetRange(i0, i1);
    if (i_axis < 6) { data_trig->GetAxis(i_axis)->SetRange(i0, i1); }
};
void tuTrSparsely::range_axes_float (int i_axis, double f0, double f1) {
    if (i_axis > 9) throw std::runtime_error(
        Form("fatal: error in tuTrSparsely, axis(%i) not valid, but by <7",
        i_axis)
    );
    TAxis *ax = data_track->GetAxis(i_axis);
    int i0 = ax->FindBin(f0);
    int i1 = ax->FindBin(f1);
    range_axes(i_axis,i0,i1);
};

TH1D* tuTrSparsely::hg_axis (int i_axis, double norm, bool is_trig){
    TH1D* hg = (TH1D*) (is_trig ? data_trig : data_track)->Projection(i_axis,"E");
    int i = 0;
    /* while (gDirectory->FindObjectAny(Form("THnSparse_unique_name__%i",i))!=nullptr) ++i; */
    hg->SetName(tuUniqueName(0,"THnSparse2__"));
    if (norm == 0.) {
        hg->Scale(1./get_n_triggers());
    } else if (norm != 1) {
        hg->Scale(1./norm);
    }
    return hg;
};
double tuTrSparsely::get_n_triggers() { 
    if (n_triggers != -1.) return n_triggers;
    auto hg = (TH1D*) data_trig->Projection(0);
    n_triggers = hg->Integral();
    delete hg;
    return n_triggers;
};
