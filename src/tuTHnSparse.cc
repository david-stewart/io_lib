#include "tuTHnSparse.h"
#include "tuClass.h"
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
