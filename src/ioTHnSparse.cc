#include "ioTHnSparse.h"
#include "io_fnc.h"
#include "ioClass.h"
// BEAR0
ostream& operator <<(ostream& os, ioSparseBinVecSet& bins) {
    os << "EAbbc:   " << bins.EAbbc   << endl;
    os << "EAtpc:   " << bins.EAtpc   << endl;
    os << "TrigEt:  " << bins.TrigEt  << endl;
    os << "vz:      " << bins.vz      << endl;
    os << "ZDCx:    " << bins.ZDCx    << endl;
    os << "absDphi: " << bins.absDphi << endl;
    os << "JetPt:   " << bins.JetPt   << endl;
    os << "leadPt:  " << bins.leadPt  << endl;
    os << "matchPt: " << bins.matchPt << endl;
    os << "AJ:      " << bins.AJ      << endl;
    return os;
};
ioTrackSparse::ioTrackSparse(
        THnSparseD* _data_track, THnSparseD* _data_trig, THnSparseD* _data_rhoPU, bool _dbprint
) :
    data_trig { _data_trig }, data_track { _data_track }, data_PU { _data_rhoPU}, debug_print{_dbprint} 
{};

ioTrackSparse::ioTrackSparse( const char* tag, bool _dbprint) :
    debug_print{_dbprint} 
{
    ioBinVec bin_TrigEt   {{ 0.,4.,8.,12.,30.}};
    ioBinVec bin_ZDCx     {{ 4000., 10000., 16000., 22000. }};
    ioBinVec bin_vz       {{ -10., -10., 20, 10. }};
    ioBinVec bin_abs_dphi {{ 0., M_PI/3., 2*M_PI/3., M_PI+0.00001 }};
    ioBinVec bin_etabin_      {{ -1.0, -0.3, 0.3, 1. }};
    ioBinVec bin_eta_PU    {{ -1.0, -0.3, 0.3, 1., 2.}};// the eta bin is for overall PU density
    ioBinVec bin_eta       {{ -1.0, -1.0, 10, 1.}};

    const ioBinVec bin_trpt {{ // track pT bins
     0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,
     1.5,  1.6,  1.7,  1.8,  1.9,  2.0,  2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,  2.9,
     3.0,  3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,  3.9,  4.0,  4.1,  4.2,  4.3,  4.4,
     4.6,  4.8,  5.0,  5.2,  5.4,  5.6,  5.8,  6.0,  6.2,  6.4,  6.6,  6.8,  7.1,  7.4,  7.7,
     8.0,  8.3,  8.6,  9.0,  9.4,  9.8, 10.3, 10.8, 11.3, 11.9, 12.5, 13.2, 14.0, 15.0
    }};
    const ioBinVec bins_EAbbc10 {{ 0, 2752.52, 5378.27, 8315.51, 11603.5, 15290.3, 19479.4, 24293.3, 30060.2, 37651.3, 64000 }};

    // get the ZDCx bins:

    // trigger 
    nbins[0] = bins_EAbbc10;
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
    data_trig->SetBinEdges(0,bins_EAbbc10);
    data_trig->SetBinEdges(1,bin_TrigEt);
    data_trig->SetBinEdges(2,bin_eta);
    data_trig->SetBinEdges(3,bin_ZDCx);
    data_trig->SetBinEdges(4,bin_vz);
    data_trig->Sumw2();

    nbins[5] = bin_eta_PU;
    data_PU = new THnSparseD(Form("data_PU%s",tag),
            "triggers;EAbbc;TrigEt;TrigEta;ZDCx;Vz;#eta(PU)",
            6, nbins, NULL, NULL);
    data_PU->SetBinEdges(0,bins_EAbbc10);
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
    data_track->SetBinEdges(0,bins_EAbbc10);
    data_track->SetBinEdges(1,bin_TrigEt);
    data_track->SetBinEdges(2,bin_eta);
    data_track->SetBinEdges(3,bin_ZDCx);
    data_track->SetBinEdges(4,bin_vz);
    data_track->SetBinEdges(5,bin_trpt);
    data_track->SetBinEdges(6,bin_abs_dphi);
    data_track->SetBinEdges(7,bin_eta);
    data_track->Sumw2();
};
void ioTrackSparse::write() { 
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
void ioTrackSparse::fill_trig(
        double EAbbc, double TrigEt, double TrigEta, double ZDCx, double Vz, 
        double PU, double PU_east, double PU_mid, double PU_west){
    hopper[0] = EAbbc;
    hopper[1] = TrigEt;
    hopper[2] = TrigEta;
    hopper[3] = ZDCx;
    hopper[4] = Vz;
    data_trig->Fill(hopper,weight);

    hopper[5] = 1.5; // will fill in PU overall bin (bin 3)
    data_PU->Fill(hopper,weight*PU);
    if (PU_east != 0) {
        hopper[5] = -0.9;
        data_PU->Fill(hopper,weight*PU_east);
    }
    if (PU_mid != 0) {
        hopper[5] = 0.;
        data_PU->Fill(hopper,weight*PU_mid);
    }
    if (PU_west != 0) {
        hopper[5] = 0.9;
        data_PU->Fill(hopper,weight*PU_west);
    }
};
void ioTrackSparse::fill_pt_eta_dphi(double pt, double eta, double dphi) {
    hopper[5] = pt;
    hopper[6] = dphi;
    hopper[7] = eta;
    data_track->Fill(hopper,weight);
};
void ioTrackSparse::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 7) throw std::runtime_error(
        Form("fatal: error in ioTrackSparse, axis(%i) not valid, but by <7",
        i_axis)
    );

    n_triggers = -1.;
    data_track ->GetAxis(i_axis)->SetRange(i0, i1);
    if (i_axis < 4) data_trig->GetAxis(i_axis)->SetRange(i0, i1);
};

TH1D* ioTrackSparse::hg_axis (int i_axis, double norm, bool use_jet_data){
    TH1D* hg;
    TH1D* _hg = (TH1D*) (use_jet_data ? data_track : data_trig)->Projection(i_axis,"E");
    if (bins != nullptr) {
        hg = (TH1D*) _hg->Rebin(*bins, ioUniqueName(), *bins);
        delete _hg;
        if (scaleByBinWidth) io_scaleByBinWidth(hg);
    } else {
        hg = _hg;
        hg->SetName(ioUniqueName());
    }
    if (norm == 0.) {
        hg->Scale(1./get_n_triggers());
    } else {
        hg->Scale(1./norm);
    }
    return hg;
};
double ioTrackSparse::get_n_triggers() { 
    if (n_triggers != -1.) return n_triggers;
    auto hg = (TH1D*) data_trig->Projection(0);
    n_triggers = hg->Integral();
    delete hg;
    return n_triggers;
};
// BEAR1

ioJetSpectraSparse::ioJetSpectraSparse(
        THnSparseD* _data_jet, THnSparseD* _data_trig,
        bool _dbprint
) :
    data_trig { _data_trig }, data_jet { _data_jet }, debug_print{_dbprint} 
{};

ioJetSpectraSparse::ioJetSpectraSparse(ioSparseBinVecSet& bins, const char* tag, bool _dbprint) :
        /* const char* bin_file, const char* tag, bool _dbprint) : */
    debug_print{_dbprint} 
{
    /* if (debug_print) cout << " bin_file: " << bin_file << endl; */
    /* TString s_tag = tag; */
    /* int i_bins = (s_tag.Contains("_10")) ? 10 : 3; */

    /* ioBinVec bin_EAbbc  { bin_file, Form("bins_EAbbc%s",tag)  }; */
    /* ioBinVec bin_EAtpc  { bin_file, Form("bins_EAtpc%s",tag)  }; */
    /* ioBinVec bin_TrigEt { bin_file, Form("bins_TrigEt%s",tag) }; */
    /* ioBinVec bin_ZDCx   { bin_file, Form("bins_ZDCx%s"  ,tag) }; */
    /* ioBinVec bin_vz     { bin_file, Form("bins_vz%s"    ,tag) }; */
    /* ioBinVec bin_JetPt  { bin_file, Form("bins_JetPt"   ,tag) }; */
    /* ioBinVec bin_absDphi{ bin_file, Form("bins_AbsDphi" , tag)}; */

    // get the ZDCx bins:

    nbins[0] = bins.EAbbc;
    nbins[1] = bins.EAtpc;
    nbins[2] = bins.TrigEt;
    nbins[3] = bins.ZDCx;
    nbins[4] = bins.vz;
    nbins[5] = bins.JetPt;
    nbins[6] = bins.absDphi;

    if (debug_print) {
        cout << " debug_print, nbins: " << endl;
        for (int i{0}; i<7; ++i) cout << " nbins["<<i<<"] " << nbins[i] << endl;
    }

    data_trig = new THnSparseD(Form("data_trig%s",tag),
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz;",
            5, nbins, NULL, NULL);
    data_trig->SetBinEdges(0,bins.EAbbc);
    data_trig->SetBinEdges(1,bins.EAtpc);
    data_trig->SetBinEdges(2,bins.TrigEt);
    data_trig->SetBinEdges(3,bins.ZDCx);
    data_trig->SetBinEdges(4,bins.vz);

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
    data_jet->SetBinEdges(0,bins.EAbbc);
    if (debug_print) cout << "bins.EAbbc: " << bins.EAbbc << endl;
    data_jet->SetBinEdges(1,bins.EAtpc);
    data_jet->SetBinEdges(2,bins.TrigEt);
    data_jet->SetBinEdges(3,bins.ZDCx);
    data_jet->SetBinEdges(4,bins.vz);
    data_jet->SetBinEdges(5,bins.JetPt);
    data_jet->SetBinEdges(6,bins.absDphi);

    data_jet->Sumw2();
};
void ioJetSpectraSparse::write() { 
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
void ioJetSpectraSparse::fill_trig(
        double EAbbc, double EAtpc, double TrigEt, double ZDCx, double Vz){
    hopper[0] = EAbbc;
    hopper[1] = EAtpc;
    hopper[2] = TrigEt;
    hopper[3] = ZDCx;
    hopper[4] = Vz;
    data_trig->Fill(hopper,weight);
};
void ioJetSpectraSparse::fill_jetpt_absDphi(double jetpt, double absDphi) {
    hopper[5] = jetpt;
    hopper[6] = absDphi;
    data_jet->Fill(hopper,weight);
};
void ioJetSpectraSparse::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 6) throw std::runtime_error(
        Form("fatal: error in ioJetSpectraSparse, axis(%i) not valid, but by <7",
        i_axis)
    );

    n_triggers = -1.;
    data_jet ->GetAxis(i_axis)->SetRange(i0, i1);
    if (i_axis < 5) data_trig->GetAxis(i_axis)->SetRange(i0, i1);
};

TH1D* ioJetSpectraSparse::hg_axis (int i_axis, double norm, bool use_jet_data){
    TH1D* hg;
    TH1D* _hg = (TH1D*) (use_jet_data ? data_jet : data_trig)->Projection(i_axis,"E");
    if (bins != nullptr) {
        hg = (TH1D*) _hg->Rebin(*bins, ioUniqueName(), *bins);
        delete _hg;
        if (scaleByBinWidth) io_scaleByBinWidth(hg);
    } else {
        hg = _hg;
        hg->SetName(ioUniqueName());
    }
    if (norm == 0.) {
        hg->Scale(1./get_n_triggers());
    } else {
        hg->Scale(1./norm);
    }
    return hg;
};
TH1D* ioJetSpectraSparse::hg_JetPt64 (int i0, int i1, double norm){
    range_absDphi(i0,i1); 
    return hg_axis(5,norm);
};
TH1D* ioJetSpectraSparse::hg_JetPt8 (int i, double norm) {
    range8_absDphi(i); 
    return hg_axis(5,norm);
};
double ioJetSpectraSparse::get_n_triggers() { 
    if (n_triggers != -1.) return n_triggers;
    auto hg = (TH1D*) data_trig->Projection(0);
    n_triggers = hg->Integral();
    delete hg;
    return n_triggers;
};



ioAjSparse::ioAjSparse(THnSparseD* _data ) : data { _data } {};
ioAjSparse::ioAjSparse(ioSparseBinVecSet& bins, const char* tag, double _rec_match) :
    recoil_match {_rec_match} {
    dphi_min = IO_pi-recoil_match;

    TString s_tag = tag;
    /* int i_bins = (s_tag.Contains("_10")) ? 10 : 3; */
    /* ioBinVec bin_EAbbc { bin_file, "EAbbc_10bin" }; */
    /* ioBinVec bin_EAtpc { bin_file, "EAtpc_10bin" }; //"EAtpc_3bin" }; */
    /* ioBinVec bin_TrigEt    {{ 0.,0.,30,30.}}; */

    /* ioBinVec info_ZDCx { bin_file, "zdcX_mu_sigma_206" }; */
    /* ioBinVec bin_ZDCx  {{ info_ZDCx[6], info_ZDCx[2], info_ZDCx[3], info_ZDCx[7] }}; */

    /* ioBinVec info_vz   { bin_file, "vz_mu_sigma_206" }; */
    /* ioBinVec bin_vz    {{ info_vz  [6], info_vz  [2], info_vz  [3], info_vz  [7] }}; */

    /* ioBinVec bin_leadPt  {{ 0., 0., 70, 70.  }}; */
    /* ioBinVec bin_matchPt {{ 0., 0., 70, 70.  }}; */
    /* ioBinVec bin_AJ      {{ 0., 0., 200., 1. }}; */

    // get the ZDCx bins:
    nbins[0] = bins.EAbbc;
    nbins[1] = bins.EAtpc;
    nbins[2] = bins.TrigEt;
    nbins[3] = bins.ZDCx;
    nbins[4] = bins.vz;
    nbins[5] = bins.leadPt; 
    nbins[6] = bins.matchPt;
    nbins[7] = bins.AJ;
    data = new THnSparseD(Form("data_%s",tag),
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz;leadPt;matchPt;AJ;",
            8, nbins, NULL, NULL);
    data->SetBinEdges(0,bins.EAbbc);
    data->SetBinEdges(1,bins.EAtpc);
    data->SetBinEdges(2,bins.TrigEt);
    data->SetBinEdges(3,bins.ZDCx);
    data->SetBinEdges(4,bins.vz);
    data->SetBinEdges(5,bins.leadPt);
    data->SetBinEdges(6,bins.matchPt);
    data->SetBinEdges(7,bins.AJ);
    data->Sumw2();
};
void ioAjSparse::write() { 
    cout << " entries for " << data->GetName() << ": " << data->GetEntries() << endl;
    data->Write();
};

void ioAjSparse::set_pTCorr ( double _mean_pTCorr, double _sig_pTCorr ) {
    mean_pTCorr = _mean_pTCorr;
    sig_pTCorr = _sig_pTCorr;
    flag_pTCorr = true;
    pTrand = new TRandom3();
};

double ioAjSparse::Aj() {
    return hopper[7];
};
bool ioAjSparse::fill(
        double EAbbc, double EAtpc, double TrigEt, double ZDCx, double Vz, double leadPt, double matchPt){
    hopper[0] = EAbbc;
    hopper[1] = EAtpc;
    hopper[2] = TrigEt;
    hopper[3] = ZDCx;
    hopper[4] = Vz;
    hopper[5] = leadPt;
    hopper[6] = matchPt;
    hopper[7] = (leadPt - matchPt)/(leadPt+matchPt);
    data->Fill(hopper,weight);
    return (leadPt >= 20 && matchPt >= 10);
};
bool ioAjSparse::fill(double* _in, vector<PseudoJet>& jets) {
    if (jets.size() < 2) return false;
    double lead_phi = jets[0].phi();
    double matchPt = -1.;
    for (unsigned int i{1}; i<jets.size(); ++i) {
        auto& jet = jets[i];
        if (io_absDphi(lead_phi,jet.phi())>dphi_min) {
            for (int k=0; k<5; ++k) hopper[k] = _in[k];
            hopper[5] = jets[0].perp();
            hopper[6] = jet.perp();
            double pT_0 = hopper[5];
            double pT_1 = hopper[6];
            if (flag_pTCorr) {
                pT_0 -= _in[1] * jet_area * (pTrand->Gaus(mean_pTCorr, sig_pTCorr));
                pT_1 -= _in[1] * jet_area * (pTrand->Gaus(mean_pTCorr, sig_pTCorr));
            }
            /* } else { */
            hopper[7] = (pT_0-pT_1)/(pT_0+pT_1);
            /* } */
            data->Fill(hopper,weight);

            return (hopper[5] >= 20 && hopper[6] >= 10);
        }
    }
    return false;
};

void ioAjSparse::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 7) throw std::runtime_error(
        Form("fatal: error in ioAjSparse, axis(%i) not valid, but by <7",
        i_axis)
    );
    data ->GetAxis(i_axis)->SetRange(i0, i1);
};

TH1D* ioAjSparse::hg_axis (int i_axis, double norm){
    TH1D* hg;
    TH1D* _hg = (TH1D*) data->Projection(i_axis,"E");
    if (bins != nullptr) {
        hg = (TH1D*) _hg->Rebin(*bins, ioUniqueName(), *bins);
        delete _hg;
        if (scaleByBinWidth) io_scaleByBinWidth(hg);
    } else {
        hg = _hg;
        hg->SetName(ioUniqueName());
    }
    if (norm != 0.) {
        hg->Scale(1./norm);
    }
    return hg;
};

/* ioAjSparse_dPhi::ioAjSparse_dPhi(THnSparseD* _data ) : data { _data } {}; */
/* ioAjSparse_dPhi::ioAjSparse_dPhi(const char* bin_file, const char* tag, double _rec_match) : */
/*     recoil_match {_rec_match} { */
/*     dphi_min = IO_pi-recoil_match; */

/*     TString s_tag = tag; */
/*     /1* int i_bins = (s_tag.Contains("_10")) ? 10 : 3; *1/ */
/*     ioBinVec bin_EAbbc { bin_file, "EAbbc_10bin" }; */
/*     ioBinVec bin_EAtpc { bin_file, "EAtpc_10bin" }; //"EAtpc_3bin" }; */
/*     ioBinVec bin_TrigEt    {{ 0.,0.,30, 30.}}; */

/*     ioBinVec info_ZDCx { bin_file, "zdcX_mu_sigma_206" }; */
/*     ioBinVec bin_ZDCx  {{ info_ZDCx[6], info_ZDCx[2], info_ZDCx[3], info_ZDCx[7] }}; */

/*     ioBinVec info_vz   { bin_file, "vz_mu_sigma_206" }; */
/*     ioBinVec bin_vz    {{ info_vz  [6], info_vz  [2], info_vz  [3], info_vz  [7] }}; */

/*     ioBinVec bin_leadPt  {{ 0., 0., 70, 70.  }}; */
/*     ioBinVec bin_matchPt {{ 0., 0., 70, 70.  }}; */
/*     ioBinVec bin_AJ      {{ 0., 0., 200., 1. }}; */
/*     ioBinVec bin_dPhi    {{ 0., 0., 100., IO_twopi }}; */

/*     // get the ZDCx bins: */
/*     nbins[0] = bin_EAbbc; */
/*     nbins[1] = bin_EAtpc; */
/*     nbins[2] = bin_TrigEt; */
/*     nbins[3] = bin_ZDCx; */
/*     nbins[4] = bin_vz; */
/*     nbins[5] = bin_leadPt; */ 
/*     nbins[6] = bin_matchPt; */
/*     nbins[7] = bin_AJ; */
/*     nbins[8] = bin_dPhi; */

/*     data = new THnSparseD(Form("data_%s",tag), */
/*             "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz;leadPt;matchPt;AJ;#Delta#phi", */
/*             9, nbins, NULL, NULL); */
/*     data->SetBinEdges(0,bin_EAbbc); */
/*     data->SetBinEdges(1,bin_EAtpc); */
/*     data->SetBinEdges(2,bin_TrigEt); */
/*     data->SetBinEdges(3,bin_ZDCx); */
/*     data->SetBinEdges(4,bin_vz); */
/*     data->SetBinEdges(5,bin_leadPt); */
/*     data->SetBinEdges(6,bin_matchPt); */
/*     data->SetBinEdges(7,bin_AJ); */
/*     data->SetBinEdges(8,bin_dPhi); */
/*     data->Sumw2(); */
/* }; */
/* void ioAjSparse_dPhi::write() { */ 
/*     cout << " entries for " << data->GetName() << ": " << data->GetEntries() << endl; */
/*     data->Write(); */
/* }; */

/* void ioAjSparse_dPhi::set_pTCorr ( double _mean_pTCorr, double _sig_pTCorr ) { */
/*     mean_pTCorr = _mean_pTCorr; */
/*     sig_pTCorr = _sig_pTCorr; */
/*     flag_pTCorr = true; */
/*     pTrand = new TRandom3(); */
/* }; */

/* double ioAjSparse_dPhi::Aj() { */
/*     return hopper[7]; */
/* }; */
/* bool ioAjSparse_dPhi::fill( */
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
/* bool ioAjSparse_dPhi::fill(double* _in, vector<PseudoJet>& jets) { */
/*     if (jets.size() < 2) return false; */
/*     double lead_phi = jets[0].phi(); */
/*     double matchPt = -1.; */
/*     for (unsigned int i{1}; i<jets.size(); ++i) { */
/*         auto& jet = jets[i]; */
/*         double dPhi {io_dphi(lead_phi,jet.phi())}; */
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
/*             if (dPhi<0) dPhi += IO_twopi; */
/*             hopper[8] = dPhi; */
/*             data->Fill(hopper,weight); */
/*             return (hopper[5] >= 10 && hopper[6] >= 10); */
/*         } */
/*     } */
/*     return false; */
/* }; */

/* void ioAjSparse_dPhi::range_axes (int i_axis, int i0, int i1) { */
/*     if (i_axis > 7) throw std::runtime_error( */
/*         Form("fatal: error in ioAjSparse_dPhi, axis(%i) not valid, but by <7", */
/*         i_axis) */
/*     ); */
/*     data ->GetAxis(i_axis)->SetRange(i0, i1); */
/* }; */

/* TH1D* ioAjSparse_dPhi::hg_axis (int i_axis, double norm){ */
/*     TH1D* hg; */
/*     TH1D* _hg = (TH1D*) data->Projection(i_axis,"E"); */
/*     if (bins != nullptr) { */
/*         hg = (TH1D*) _hg->Rebin(*bins, ioUniqueName(), *bins); */
/*         delete _hg; */
/*         if (scaleByBinWidth) io_scaleByBinWidth(hg); */
/*     } else { */
/*         hg = _hg; */
/*         hg->SetName(ioUniqueName()); */
/*     } */
/*     if (norm != 0.) { */
/*         hg->Scale(1./norm); */
/*     } */
/*     return hg; */
/* }; */

/* TH2D* ioAjSparse_dPhi::hg2_axis (int xAxis, int yAxis, double norm){ */
/*     TH2D* hg = (TH2D*) data->Projection(yAxis, xAxis); */
/*     if (norm) hg->Scale(norm); */
/*     return hg; */
/* }; */


