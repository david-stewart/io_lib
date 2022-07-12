#include "tuTHnSparse.h"
#include "tuClass.h"
#include "tu_fnc.h"


tuTrackSparse2::tuTrackSparse2(
        THnSparseD* _data_track, THnSparseD* _data_trig, bool _dbprint
) :
    data_trig { _data_trig }, data_track { _data_track }, debug_print{_dbprint} 
{};

tuTrackSparse2::tuTrackSparse2( const char* tag, bool _dbprint) :
    debug_print{_dbprint} 
{
    tuBinVec bin_EAbbc3 {{ 0.0000,  8315.4606, 24292.8207, 64000.0000 }};
    tuBinVec bin_TrigEt  {{ 0.,4.,8.,12.,30.}};
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
    nbins[2] = bin_eta;
    nbins[3] = bin_phi;
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

    string trig_string= "bbc;E_{T};Trig #eta;Trig #phi;ZDCx;runId";
    string word_trig = "triggers;";
           
    data_trig = new THnSparseD(Form("data_trig2%s",tag),
            (word_trig+trig_string).c_str(),
            6, nbins, NULL, NULL);
    data_trig->SetBinEdges(0,bin_EAbbc3);
    data_trig->SetBinEdges(1,bin_TrigEt);
    data_trig->SetBinEdges(2,bin_eta);
    data_trig->SetBinEdges(3,bin_phi);
    data_trig->SetBinEdges(4,bin_ZDCx);
    data_trig->SetBinEdges(5,bin_runId);
    data_trig->Sumw2();

    string track_string= "track-pT;track #Delta#phi;track #phi;track #eta;dca";
    string word_track ="track;";
    data_track = new THnSparseD(Form("data_track2%s",tag),
            (word_track+trig_string+track_string).c_str(),
            11, nbins, NULL, NULL);
    data_track->SetBinEdges(0,bin_EAbbc3);
    data_track->SetBinEdges(1,bin_TrigEt);
    data_track->SetBinEdges(2,bin_eta);
    data_track->SetBinEdges(3,bin_phi);
    data_track->SetBinEdges(4,bin_ZDCx);
    data_track->SetBinEdges(5,bin_runId);
    data_track->SetBinEdges(6,bin_trpt);
    data_track->SetBinEdges(7,bin_abs_dphi);
    data_track->SetBinEdges(8,bin_phi);
    data_track->SetBinEdges(9,bin_eta);
    data_track->SetBinEdges(10,bin_dca);
    data_track->Sumw2();
};
void tuTrackSparse2::write() { 
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
void tuTrackSparse2::fill_trig (double EAbbc, double TrigEt, double TrigEta, double TrigPhi, double ZDCx, double runId) {
    hopper[0] = EAbbc;
    hopper[1] = TrigEt;
    hopper[2] = TrigEta;
    hopper[3] = TrigPhi;
    hopper[4] = ZDCx;
    hopper[5] = runId;
    data_trig->Fill(hopper,weight);
};
void tuTrackSparse2::fill_track(double pt, double eta, double phi, double dca){
    hopper[6] = pt;
    hopper[7] = tu_absDphi(hopper[3],phi);
    hopper[8] = phi;
    hopper[9] = eta;
    hopper[10] = dca;
    data_track->Fill(hopper,weight);
};
void tuTrackSparse2::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 9) throw std::runtime_error(
        Form("fatal: error in tuTrackSparse2, axis(%i) not valid, but by <7",
        i_axis)
    );

    n_triggers = -1.;
    data_track ->GetAxis(i_axis)->SetRange(i0, i1);
    if (i_axis < 6) { data_trig->GetAxis(i_axis)->SetRange(i0, i1); }
};
void tuTrackSparse2::range_axes_float (int i_axis, double f0, double f1) {
    if (i_axis > 9) throw std::runtime_error(
        Form("fatal: error in tuTrackSparse2, axis(%i) not valid, but by <7",
        i_axis)
    );
    TAxis *ax = data_track->GetAxis(i_axis);
    int i0 = ax->FindBin(f0);
    int i1 = ax->FindBin(f1);
    range_axes(i_axis,i0,i1);
};

TH1D* tuTrackSparse2::hg_axis (int i_axis, double norm, bool is_trig){
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
double tuTrackSparse2::get_n_triggers() { 
    if (n_triggers != -1.) return n_triggers;
    auto hg = (TH1D*) data_trig->Projection(0);
    n_triggers = hg->Integral();
    delete hg;
    return n_triggers;
};
//END tuTrackSprase2 LION
// BEAR0
tuTrackSparse::tuTrackSparse(
        THnSparseD* _data_track, THnSparseD* _data_trig, THnSparseD* _data_rhoPU, bool _dbprint
) :
    data_trig { _data_trig }, data_track { _data_track }, data_PU { _data_rhoPU}, debug_print{_dbprint} 
{};

tuTrackSparse::tuTrackSparse( const char* tag, bool _dbprint) :
    debug_print{_dbprint} 
{
    tuBinVec bin_TrigEt   {{ 0.,4.,8.,12.,30.}};
    tuBinVec bin_ZDCx     {{ 4000., 10000., 16000., 22000. }};
    tuBinVec bin_vz       {{ -10., -10., 20, 10. }};
    tuBinVec bin_abs_dphi {{ 0., M_PI/3., 2*M_PI/3., M_PI+0.00001 }};
    tuBinVec bin_eta      {{ -1.0, -0.3, 0.3, 1. }};
    tuBinVec bin_eta_PU    {{ -1.0, -0.3, 0.3, 1., 2.}};

    const tuBinVec bin_trpt {{ // track pT bins
     0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,
     1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,
     3.0,  3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,  4.0,  4.1,  4.2,  4.3,  4.4,
     4.6,  4.8,  5.0,  5.2,  5.4,  5.6,  5.8,  6.0,  6.2,  6.4,  6.6,  6.8,  7.1,  7.4,  7.7,
     8.0,  8.3,  8.6,  9.0,  9.4,  9.8, 10.3, 10.8, 11.3, 11.9, 12.5, 13.2, 14.0, 15.0
    }};
    const tuBinVec bin_EAbbc10 {{ 0, 2767.15, 5397.98, 8333.35, 11610.5,
        15280.9, 19440.2, 24219.7, 29959, 37534.5, 64000.0 }};

    // get the ZDCx bins:

    // trigger 
    nbins[0] = bin_EAbbc10;
    nbins[1] = bin_TrigEt;
    nbins[2] = bin_eta;
    nbins[3] = bin_ZDCx;
    nbins[4] = bin_vz;

    // track
    nbins[5] = bin_trpt;
    nbins[6] = bin_abs_dphi;
    nbins[7] = bin_eta;

    if (debug_print) {
        cout << " debug_print, nbins: " << endl;
        for (int i{0}; i<7; ++i) cout << " nbins["<<i<<"] " << nbins[i] << endl;
    }

    data_trig = new THnSparseD(Form("data_trig%s",tag),
            "triggers;EAbbc;TrigEt;TrigEta;ZDCx;Vz;",
            5, nbins, NULL, NULL);
    data_trig->SetBinEdges(0,bin_EAbbc10);
    data_trig->SetBinEdges(1,bin_TrigEt);
    data_trig->SetBinEdges(2,bin_eta);
    data_trig->SetBinEdges(3,bin_ZDCx);
    data_trig->SetBinEdges(4,bin_vz);
    data_trig->Sumw2();

    nbins[5] = bin_eta_PU;
    data_PU = new THnSparseD(Form("data_PU%s",tag),
            "triggers;EAbbc;TrigEt;TrigEta;ZDCx;Vz;eta_PU",
            6, nbins, NULL, NULL);
    data_PU->SetBinEdges(0,bin_EAbbc10);
    data_PU->SetBinEdges(1,bin_TrigEt);
    data_PU->SetBinEdges(2,bin_eta);
    data_PU->SetBinEdges(3,bin_ZDCx);
    data_PU->SetBinEdges(4,bin_vz);
    data_PU->SetBinEdges(5,bin_eta_PU); // this is the PU eta, and will get filled in each bins once per event
    data_PU->Sumw2();

    if (debug_print) {
        for (int k = 0; k<5; ++k) {
            TAxis* x = data_PU->GetAxis(k);
            cout << k << " debug print: THnSparse Axis: " << x->GetTitle() 
                 << " axis:  " << x->GetTitle() << " nbins: " << x->GetNbins() << endl;
            for (int j = 1; j<=x->GetNbins()+1; ++j) 
                cout << "bin["<<j<<"] " << x->GetBinLowEdge(j) << " ";
            cout << endl;
            cout << k << " ------------- " << endl;
        }
    }

    nbins[5] = bin_trpt;
    data_track = new THnSparseD(Form("data_track%s",tag),
            "jets;EAbbc;TrigEt;TrigEta;ZDCx;Vz;track pT;|#Delta#phi|;#eta",
            8, nbins, NULL, NULL);
    data_track->SetBinEdges(0,bin_EAbbc10);
    data_track->SetBinEdges(1,bin_TrigEt);
    data_track->SetBinEdges(2,bin_eta);
    data_track->SetBinEdges(3,bin_ZDCx);
    data_track->SetBinEdges(4,bin_vz);
    data_track->SetBinEdges(5,bin_trpt);
    data_track->SetBinEdges(6,bin_abs_dphi);
    data_track->SetBinEdges(7,bin_eta);
    data_track->Sumw2();
};
void tuTrackSparse::write() { 
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
    data_PU->Write();
    data_track->Write();
};
void tuTrackSparse::fill_trig(
        double EAbbc, double TrigEt, double TrigEta, double ZDCx, double Vz, 
        double PU, double east, double mid, double west){
    hopper[0] = EAbbc;
    hopper[1] = TrigEt;
    hopper[2] = TrigEta;
    hopper[3] = ZDCx;
    hopper[4] = Vz;
    data_trig->Fill(hopper,weight);

    hopper[5] = 1.5; // will fill in PU overall bin (bin 3)
    data_PU->Fill(hopper,weight*PU);
    if (east != 0) {
        hopper[5] = -0.9;
        data_PU->Fill(hopper,weight*east);
    }
    if (mid != 0) {
        hopper[5] = 0.;
        data_PU->Fill(hopper,weight*mid);
    }
    if (west != 0) {
        hopper[5] = 0.9;
        data_PU->Fill(hopper,weight*west);
    }
};
void tuTrackSparse::fill_pt_eta_dphi(double pt, double eta, double dphi) {
    hopper[5] = pt;
    hopper[6] = dphi;
    hopper[7] = eta;
    data_track->Fill(hopper,weight);
};
void tuTrackSparse::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 7) throw std::runtime_error(
        Form("fatal: error in tuTrackSparse, axis(%i) not valid, but by <7",
        i_axis)
    );

    n_triggers = -1.;
    data_track ->GetAxis(i_axis)->SetRange(i0, i1);
    if (i_axis < 5) {
        data_trig->GetAxis(i_axis)->SetRange(i0, i1);
        data_PU  ->GetAxis(i_axis)->SetRange(i0, i1);
    }
};

TH1D* tuTrackSparse::hg_axis (int i_axis, double norm, bool use_jet_data){
    TH1D* hg = (TH1D*) (use_jet_data ? data_track : data_trig)->Projection(i_axis,"E");
    int i = 0;
    while (gDirectory->FindObjectAny(Form("THnSparse_unique_name__%i",i))!=nullptr) ++i;
    hg->SetName(Form("THnSparse_unique_name__%i",i));
    if (norm == 0.) {
        hg->Scale(1./get_n_triggers());
    } else {
        hg->Scale(1./norm);
    }
    return hg;
};
double tuTrackSparse::get_n_triggers() { 
    if (n_triggers != -1.) return n_triggers;
    auto hg = (TH1D*) data_trig->Projection(0);
    n_triggers = hg->Integral();
    delete hg;
    return n_triggers;
};
double tuTrackSparse::get_sum_PU(int i_bin) { // bins 1 thorugh 4
    data_PU->GetAxis(5)->SetRange(i_bin,i_bin);
    auto hg = (TH1D*) data_PU->Projection(5);
    double val = hg->Integral();
    delete hg;
    return val;
};
// BEAR1

tuJetSpectraSparse::tuJetSpectraSparse(
        THnSparseD* _data_jet, THnSparseD* _data_trig,
        bool _dbprint
) :
    data_trig { _data_trig }, data_jet { _data_jet }, debug_print{_dbprint} 
{};

tuJetSpectraSparse::tuJetSpectraSparse(
        const char* bin_file, const char* tag, bool _dbprint) :
    debug_print{_dbprint} 
{
    if (debug_print) cout << " bin_file: " << bin_file << endl;
    TString s_tag = tag;
    int i_bins = (s_tag.Contains("_10")) ? 10 : 3;
    tuBinVec bin_EAbbc { bin_file, Form("EAbbc_%ibin",i_bins) };
    tuBinVec bin_EAtpc { bin_file, Form("EAtpc_%ibin",i_bins) }; //"EAtpc_3bin" };
    tuBinVec bin_TrigEt    {{ 0.,0.,30,30.}};

    tuBinVec info_ZDCx { bin_file, "zdcX_mu_sigma_206" };
    tuBinVec bin_ZDCx  {{ info_ZDCx[6], info_ZDCx[2], info_ZDCx[3], info_ZDCx[7] }};

    tuBinVec info_vz   { bin_file, "vz_mu_sigma_206" };
    tuBinVec bin_vz    {{ info_vz  [6], info_vz  [2], info_vz  [3], info_vz  [7] }};

    tuBinVec bin_JetPt {{ 0., 0., 70, 70.}};
    tuBinVec bin_absDphi {{ 0., 0., 64., M_PI }};

    // get the ZDCx bins:

    nbins[0] = bin_EAbbc;
    nbins[1] = bin_EAtpc;
    nbins[2] = bin_TrigEt;
    nbins[3] = bin_ZDCx;
    nbins[4] = bin_vz;
    nbins[5] = bin_JetPt;
    nbins[6] = bin_absDphi;

    if (debug_print) {
        cout << " debug_print, nbins: " << endl;
        for (int i{0}; i<7; ++i) cout << " nbins["<<i<<"] " << nbins[i] << endl;
    }

    data_trig = new THnSparseD(Form("data_trig%s",tag),
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz;",
            5, nbins, NULL, NULL);
    data_trig->SetBinEdges(0,bin_EAbbc);
    data_trig->SetBinEdges(1,bin_EAtpc);
    data_trig->SetBinEdges(2,bin_TrigEt);
    data_trig->SetBinEdges(3,bin_ZDCx);
    data_trig->SetBinEdges(4,bin_vz);

    if (debug_print) {
        for (int k = 0; k<5; ++k) {
            TAxis* x = data_trig->GetAxis(k);
            cout << k << " debug print: THnSparse Axis: " << x->GetTitle() 
                 << " axis:  " << x->GetTitle() << " nbins: " << x->GetNbins() << endl;
            for (int j = 1; j<=x->GetNbins()+1; ++j) 
                cout << "bin["<<j<<"] " << x->GetBinLowEdge(j) << " ";
            cout << endl;
            cout << k << " ------------- " << endl;
        }
    }

    data_jet = new THnSparseD(Form("data_jet%s",tag),
            "jets;EAbbc;EAtpc;TrigEt;ZDCx;Vz;Jet #it{p}_{T};|#Delta#phi|;",
            7, nbins, NULL, NULL);
    data_jet->SetBinEdges(0,bin_EAbbc);
    if (debug_print) cout << "bin_EAbbc: " << bin_EAbbc << endl;
    data_jet->SetBinEdges(1,bin_EAtpc);
    data_jet->SetBinEdges(2,bin_TrigEt);
    data_jet->SetBinEdges(3,bin_ZDCx);
    data_jet->SetBinEdges(4,bin_vz);
    data_jet->SetBinEdges(5,bin_JetPt);
    data_jet->SetBinEdges(6,bin_absDphi);

    data_jet->Sumw2();
};
void tuJetSpectraSparse::write() { 
    if (debug_print) {
        cout << " entries for " << data_trig->GetName() << " trig-entries: " <<
            data_trig->GetEntries() << " jet-entries: " << data_jet->GetEntries() << endl;
        auto t_proj = data_trig->Projection(0);
        t_proj->SetName("p_trip");
        cout << " First axes: (mean) trig " << t_proj->GetMean() << endl;
        delete t_proj;

        auto j_proj = data_jet->Projection(0);
        j_proj->SetName("j_trip");
        cout << " First axes: (mean) jet  " << j_proj->GetMean() << endl;
        delete j_proj;
    }
    data_trig->Write();
    data_jet->Write();
};
/* void tuJetSpectraSparse::fill_trig( */
/*         double EAbbc, double EAtpc, double TrigEt, double ZDCx, double Vz){ */
/*     hopper[0] = EAbbc; */
/*     hopper[1] = EAtpc; */
/*     hopper[2] = TrigEt; */
/*     hopper[3] = ZDCx; */
/*     hopper[4] = Vz; */
/*     data_trig->Fill(hopper,weight); */
/* }; */
/* void tuJetSpectraSparse::fill_jetpt_absDphi(double jetpt, double absDphi) { */
/*     hopper[5] = jetpt; */
/*     hopper[6] = absDphi; */
/*     data_jet->Fill(hopper,weight); */
/* }; */
void tuJetSpectraSparse::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 6) throw std::runtime_error(
        Form("fatal: error in tuJetSpectraSparse, axis(%i) not valid, but by <7",
        i_axis)
    );

    n_triggers = -1.;
    data_jet ->GetAxis(i_axis)->SetRange(i0, i1);
    if (i_axis < 5) data_trig->GetAxis(i_axis)->SetRange(i0, i1);
};

TH1D* tuJetSpectraSparse::hg_axis (int i_axis, double norm, bool use_jet_data){
    TH1D* hg;
    TH1D* _hg = (TH1D*) (use_jet_data ? data_jet : data_trig)->Projection(i_axis,"E");
    if (bins != nullptr) {
        hg = (TH1D*) _hg->Rebin(*bins, tuUniqueName(), *bins);
        delete _hg;
        if (scaleByBinWidth) tu_scaleByBinWidth(hg);
    } else {
        hg = _hg;
        hg->SetName(tuUniqueName());
    }
    if (norm == 0.) {
        hg->Scale(1./get_n_triggers());
    } else {
        hg->Scale(1./norm);
    }
    return hg;
};
TH1D* tuJetSpectraSparse::hg_JetPt64 (int i0, int i1, double norm){
    range_absDphi(i0,i1); 
    return hg_axis(5,norm);
};
TH1D* tuJetSpectraSparse::hg_JetPt8 (int i, double norm) {
    range8_absDphi(i); 
    return hg_axis(5,norm);
};
double tuJetSpectraSparse::get_n_triggers() { 
    if (n_triggers != -1.) return n_triggers;
    auto hg = (TH1D*) data_trig->Projection(0);
    n_triggers = hg->Integral();
    delete hg;
    return n_triggers;
};



tuAjSparse::tuAjSparse(THnSparseD* _data ) : data { _data } {};
tuAjSparse::tuAjSparse(const char* bin_file, const char* tag, double _rec_match) :
    recoil_match {_rec_match} {
    dphi_min = M_PI-recoil_match;

    TString s_tag = tag;
    /* int i_bins = (s_tag.Contains("_10")) ? 10 : 3; */
    tuBinVec bin_EAbbc { bin_file, "EAbbc_10bin" };
    tuBinVec bin_EAtpc { bin_file, "EAtpc_10bin" }; //"EAtpc_3bin" };
    tuBinVec bin_TrigEt    {{ 0.,0.,30,30.}};

    tuBinVec info_ZDCx { bin_file, "zdcX_mu_sigma_206" };
    tuBinVec bin_ZDCx  {{ info_ZDCx[6], info_ZDCx[2], info_ZDCx[3], info_ZDCx[7] }};

    tuBinVec info_vz   { bin_file, "vz_mu_sigma_206" };
    tuBinVec bin_vz    {{ info_vz  [6], info_vz  [2], info_vz  [3], info_vz  [7] }};

    tuBinVec bin_leadPt  {{ 0., 0., 70, 70.  }};
    tuBinVec bin_matchPt {{ 0., 0., 70, 70.  }};
    tuBinVec bin_AJ      {{ 0., 0., 200., 1. }};

    // get the ZDCx bins:
    nbins[0] = bin_EAbbc;
    nbins[1] = bin_EAtpc;
    nbins[2] = bin_TrigEt;
    nbins[3] = bin_ZDCx;
    nbins[4] = bin_vz;
    nbins[5] = bin_leadPt; 
    nbins[6] = bin_matchPt;
    nbins[7] = bin_AJ;
    data = new THnSparseD(Form("data_%s",tag),
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz;leadPt;matchPt;AJ;",
            8, nbins, NULL, NULL);
    data->SetBinEdges(0,bin_EAbbc);
    data->SetBinEdges(1,bin_EAtpc);
    data->SetBinEdges(2,bin_TrigEt);
    data->SetBinEdges(3,bin_ZDCx);
    data->SetBinEdges(4,bin_vz);
    data->SetBinEdges(5,bin_leadPt);
    data->SetBinEdges(6,bin_matchPt);
    data->SetBinEdges(7,bin_AJ);
    data->Sumw2();
};
void tuAjSparse::write() { 
    cout << " entries for " << data->GetName() << ": " << data->GetEntries() << endl;
    data->Write();
};

void tuAjSparse::set_pTCorr ( double _mean_pTCorr, double _sig_pTCorr ) {
    mean_pTCorr = _mean_pTCorr;
    sig_pTCorr = _sig_pTCorr;
    flag_pTCorr = true;
    pTrand = new TRandom3();
};

double tuAjSparse::Aj() {
    return hopper[7];
};
/* bool tuAjSparse::fill( */
/*         double EAbbc, double EAtpc, double TrigEt, double ZDCx, double Vz, double leadPt, double matchPt){ */
/*     hopper[0] = EAbbc; */
/*     hopper[1] = EAtpc; */
/*     hopper[2] = TrigEt; */
/*     hopper[3] = ZDCx; */
/*     hopper[4] = Vz; */
/*     hopper[5] = leadPt; */
/*     hopper[6] = matchPt; */
/*     hopper[7] = (leadPt - matchPt)/(leadPt+matchPt); */
/*     data->Fill(hopper,weight); */
/*     return (leadPt >= 20 && matchPt >= 10); */
/* }; */
/* bool tuAjSparse::fill(double* _in, vector<PseudoJet>& jets) { */
/*     if (jets.size() < 2) return false; */
/*     double lead_phi = jets[0].phi(); */
/*     double matchPt = -1.; */
/*     for (unsigned int i{1}; i<jets.size(); ++i) { */
/*         auto& jet = jets[i]; */
/*         if (tu_absDphi(lead_phi,jet.phi())>dphi_min) { */
/*             for (int k=0; k<5; ++k) hopper[k] = _in[k]; */
/*             hopper[5] = jets[0].perp(); */
/*             hopper[6] = jet.perp(); */
/*             double pT_0 = hopper[5]; */
/*             double pT_1 = hopper[6]; */
/*             if (flag_pTCorr) { */
/*                 pT_0 -= _in[1] * jet_area * (pTrand->Gaus(mean_pTCorr, sig_pTCorr)); */
/*                 pT_1 -= _in[1] * jet_area * (pTrand->Gaus(mean_pTCorr, sig_pTCorr)); */
/*             } */
/*             /1* } else { *1/ */
/*             hopper[7] = (pT_0-pT_1)/(pT_0+pT_1); */
/*             /1* } *1/ */
/*             data->Fill(hopper,weight); */

/*             return (hopper[5] >= 20 && hopper[6] >= 10); */
/*         } */
/*     } */
/*     return false; */
/* }; */

void tuAjSparse::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 7) throw std::runtime_error(
        Form("fatal: error in tuAjSparse, axis(%i) not valid, but by <7",
        i_axis)
    );
    data ->GetAxis(i_axis)->SetRange(i0, i1);
};

TH1D* tuAjSparse::hg_axis (int i_axis, double norm){
    TH1D* hg;
    TH1D* _hg = (TH1D*) data->Projection(i_axis,"E");
    if (bins != nullptr) {
        hg = (TH1D*) _hg->Rebin(*bins, tuUniqueName(), *bins);
        delete _hg;
        if (scaleByBinWidth) tu_scaleByBinWidth(hg);
    } else {
        hg = _hg;
        hg->SetName(tuUniqueName());
    }
    if (norm != 0.) {
        hg->Scale(1./norm);
    }
    return hg;
};

tuAjSparse_dPhi::tuAjSparse_dPhi(THnSparseD* _data ) : data { _data } {};
tuAjSparse_dPhi::tuAjSparse_dPhi(const char* bin_file, const char* tag, double _rec_match) :
    recoil_match {_rec_match} {
    dphi_min = M_PI-recoil_match;

    TString s_tag = tag;
    /* int i_bins = (s_tag.Contains("_10")) ? 10 : 3; */
    tuBinVec bin_EAbbc { bin_file, "EAbbc_10bin" };
    tuBinVec bin_EAtpc { bin_file, "EAtpc_10bin" }; //"EAtpc_3bin" };
    tuBinVec bin_TrigEt    {{ 0.,0.,30, 30.}};

    tuBinVec info_ZDCx { bin_file, "zdcX_mu_sigma_206" };
    tuBinVec bin_ZDCx  {{ info_ZDCx[6], info_ZDCx[2], info_ZDCx[3], info_ZDCx[7] }};

    tuBinVec info_vz   { bin_file, "vz_mu_sigma_206" };
    tuBinVec bin_vz    {{ info_vz  [6], info_vz  [2], info_vz  [3], info_vz  [7] }};

    tuBinVec bin_leadPt  {{ 0., 0., 70, 70.  }};
    tuBinVec bin_matchPt {{ 0., 0., 70, 70.  }};
    tuBinVec bin_AJ      {{ 0., 0., 200., 1. }};
    tuBinVec bin_dPhi    {{ 0., 0., 100., 2*M_PI }};

    // get the ZDCx bins:
    nbins[0] = bin_EAbbc;
    nbins[1] = bin_EAtpc;
    nbins[2] = bin_TrigEt;
    nbins[3] = bin_ZDCx;
    nbins[4] = bin_vz;
    nbins[5] = bin_leadPt; 
    nbins[6] = bin_matchPt;
    nbins[7] = bin_AJ;
    nbins[8] = bin_dPhi;

    data = new THnSparseD(Form("data_%s",tag),
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz;leadPt;matchPt;AJ;#Delta#phi",
            9, nbins, NULL, NULL);
    data->SetBinEdges(0,bin_EAbbc);
    data->SetBinEdges(1,bin_EAtpc);
    data->SetBinEdges(2,bin_TrigEt);
    data->SetBinEdges(3,bin_ZDCx);
    data->SetBinEdges(4,bin_vz);
    data->SetBinEdges(5,bin_leadPt);
    data->SetBinEdges(6,bin_matchPt);
    data->SetBinEdges(7,bin_AJ);
    data->SetBinEdges(8,bin_dPhi);
    data->Sumw2();
};
void tuAjSparse_dPhi::write() { 
    cout << " entries for " << data->GetName() << ": " << data->GetEntries() << endl;
    data->Write();
};

void tuAjSparse_dPhi::set_pTCorr ( double _mean_pTCorr, double _sig_pTCorr ) {
    mean_pTCorr = _mean_pTCorr;
    sig_pTCorr = _sig_pTCorr;
    flag_pTCorr = true;
    pTrand = new TRandom3();
};

double tuAjSparse_dPhi::Aj() {
    return hopper[7];
};
/* bool tuAjSparse_dPhi::fill( */
/*         double EAbbc, double EAtpc, double TrigEt, */ 
/*         double ZDCx, double Vz, double leadPt, double matchPt, */
/*         double dPhi){ */
/*     hopper[0] = EAbbc; */
/*     hopper[1] = EAtpc; */
/*     hopper[2] = TrigEt; */
/*     hopper[3] = ZDCx; */
/*     hopper[4] = Vz; */
/*     hopper[5] = leadPt; */
/*     hopper[6] = matchPt; */
/*     hopper[7] = (leadPt - matchPt)/(leadPt+matchPt); */
/*     hopper[8] = dPhi; */
/*     data->Fill(hopper,weight); */
/*     return (leadPt >= 20 && matchPt >= 10); */
/* }; */
/* bool tuAjSparse_dPhi::fill(double* _in, vector<PseudoJet>& jets) { */
/*     if (jets.size() < 2) return false; */
/*     double lead_phi = jets[0].phi(); */
/*     double matchPt = -1.; */
/*     for (unsigned int i{1}; i<jets.size(); ++i) { */
/*         auto& jet = jets[i]; */
/*         double dPhi {tu_dphi(lead_phi,jet.phi())}; */
/*         if (TMath::Abs(dPhi)>dphi_min) { */
/*             for (int k=0; k<5; ++k) hopper[k] = _in[k]; */
/*             hopper[5] = jets[0].perp(); */
/*             hopper[6] = jet.perp(); */
/*             double pT_0 = hopper[5]; */
/*             double pT_1 = hopper[6]; */
/*             if (flag_pTCorr) { */
/*                 pT_0 -= _in[1] * jet_area * (pTrand->Gaus(mean_pTCorr, sig_pTCorr)); */
/*                 pT_1 -= _in[1] * jet_area * (pTrand->Gaus(mean_pTCorr, sig_pTCorr)); */
/*             } */
/*             hopper[7] = (pT_0-pT_1)/(pT_0+pT_1); */
/*             if (dPhi<0) dPhi += tu_twopi; */
/*             hopper[8] = dPhi; */
/*             data->Fill(hopper,weight); */
/*             return (hopper[5] >= 10 && hopper[6] >= 10); */
/*         } */
/*     } */
/*     return false; */
/* }; */

void tuAjSparse_dPhi::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 7) throw std::runtime_error(
        Form("fatal: error in tuAjSparse_dPhi, axis(%i) not valid, but by <7",
        i_axis)
    );
    data ->GetAxis(i_axis)->SetRange(i0, i1);
};

TH1D* tuAjSparse_dPhi::hg_axis (int i_axis, double norm){
    TH1D* hg;
    TH1D* _hg = (TH1D*) data->Projection(i_axis,"E");
    if (bins != nullptr) {
        hg = (TH1D*) _hg->Rebin(*bins, tuUniqueName(), *bins);
        delete _hg;
        if (scaleByBinWidth) tu_scaleByBinWidth(hg);
    } else {
        hg = _hg;
        hg->SetName(tuUniqueName());
    }
    if (norm != 0.) {
        hg->Scale(1./norm);
    }
    return hg;
};

TH2D* tuAjSparse_dPhi::hg2_axis (int xAxis, int yAxis, double norm){
    TH2D* hg = (TH2D*) data->Projection(yAxis, xAxis);
    if (norm) hg->Scale(norm);
    return hg;
};

// BREAK v2
