#include "ioTHnSparse.h"
#include "io_fnc.h"
#include "ioClass.h"

ioJetSpectraSparse::ioJetSpectraSparse(THnSparseD* _data_jet, THnSparseD* _data_trig) :
    data_trig { _data_trig }, data_jet { _data_jet } {};
ioJetSpectraSparse::ioJetSpectraSparse(const char* bin_file, const char* tag) {
    TString s_tag = tag;
    int i_bins = (s_tag.Contains("_10")) ? 10 : 3;
    ioBinVec bin_EAbbc { bin_file, Form("EAbbc_%ibin",i_bins) };
    ioBinVec bin_EAtpc { bin_file, Form("EAtpc_%ibin",i_bins) }; //"EAtpc_3bin" };
    ioBinVec bin_TrigEt    {{ 0.,0.,30,30.}};

    ioBinVec info_ZDCx { bin_file, "zdcX_mu_sigma_206" };
    ioBinVec bin_ZDCx  {{ info_ZDCx[6], info_ZDCx[2], info_ZDCx[3], info_ZDCx[7] }};

    ioBinVec info_vz   { bin_file, "vz_mu_sigma_206" };
    ioBinVec bin_vz    {{ info_vz  [6], info_vz  [2], info_vz  [3], info_vz  [7] }};

    ioBinVec bin_JetPt {{ 0., 0., 70, 70.}};
    ioBinVec bin_absDphi {{ 0., 0., 64., IO_pi }};

    // get the ZDCx bins:

    nbins[0] = bin_EAbbc;
    nbins[1] = bin_EAtpc;
    nbins[2] = bin_TrigEt;
    nbins[3] = bin_ZDCx;
    nbins[4] = bin_vz;
    nbins[5] = bin_JetPt;
    nbins[6] = bin_absDphi;

    data_trig = new THnSparseD(Form("data_trig%s",tag),
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz",
            5, nbins, NULL, NULL);
    data_trig->SetBinEdges(0,bin_EAbbc);
    data_trig->SetBinEdges(1,bin_EAtpc);
    data_trig->SetBinEdges(2,bin_TrigEt);
    data_trig->SetBinEdges(3,bin_ZDCx);
    data_trig->SetBinEdges(4,bin_vz);

    data_jet = new THnSparseD(Form("data_jet%s",tag),
            "jets;EAbbc;EAtpc;TrigEt;ZDCx;Vz;Jet #it{p}_{T}",
            7, nbins, NULL, NULL);
    data_jet->SetBinEdges(0,bin_EAbbc);
    data_jet->SetBinEdges(1,bin_EAtpc);
    data_jet->SetBinEdges(2,bin_TrigEt);
    data_jet->SetBinEdges(3,bin_ZDCx);
    data_jet->SetBinEdges(4,bin_vz);
    data_jet->SetBinEdges(5,bin_JetPt);
    data_jet->SetBinEdges(6,bin_absDphi);

    data_jet->Sumw2();
};
void ioJetSpectraSparse::write() { 
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
    data_trig->Fill(hopper);
};
void ioJetSpectraSparse::fill_jetpt_absDphi(double jetpt, double absDphi) {
    hopper[5] = jetpt;
    hopper[6] = absDphi;
    data_jet->Fill(hopper);
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
        io_scaleByBinWidth(hg);
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
ioAjSparse::ioAjSparse(const char* bin_file, const char* tag, double _rec_match) :
    recoil_match {_rec_match} {
    dphi_min = IO_pi-recoil_match;

    TString s_tag = tag;
    /* int i_bins = (s_tag.Contains("_10")) ? 10 : 3; */
    ioBinVec bin_EAbbc { bin_file, "EAbbc_10bin" };
    ioBinVec bin_EAtpc { bin_file, "EAtpc_10bin" }; //"EAtpc_3bin" };
    ioBinVec bin_TrigEt    {{ 0.,0.,30,30.}};

    ioBinVec info_ZDCx { bin_file, "zdcX_mu_sigma_206" };
    ioBinVec bin_ZDCx  {{ info_ZDCx[6], info_ZDCx[2], info_ZDCx[3], info_ZDCx[7] }};

    ioBinVec info_vz   { bin_file, "vz_mu_sigma_206" };
    ioBinVec bin_vz    {{ info_vz  [6], info_vz  [2], info_vz  [3], info_vz  [7] }};

    ioBinVec bin_leadPt  {{ 0., 0., 70, 70.  }};
    ioBinVec bin_matchPt {{ 0., 0., 70, 70.  }};
    ioBinVec bin_AJ      {{ 0., 0., 200., 1. }};

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
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz;leadPt;matchPt;AJ",
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
void ioAjSparse::write() { 
    data->Write();
};
void ioAjSparse::fill(
        double EAbbc, double EAtpc, double TrigEt, double ZDCx, double Vz, double leadPt, double matchPt){
    hopper[0] = EAbbc;
    hopper[1] = EAtpc;
    hopper[2] = TrigEt;
    hopper[3] = ZDCx;
    hopper[4] = Vz;
    hopper[5] = leadPt;
    hopper[6] = matchPt;
    hopper[7] = (leadPt - matchPt)/(leadPt+matchPt);
    data->Fill(hopper);
};
void ioAjSparse::fill(double* _in, vector<PseudoJet>& jets) {
    if (jets.size() < 2) return;
    double lead_phi = jets[0].phi();
    double matchPt = -1.;
    for (unsigned int i{1}; i<jets.size(); ++i) {
        auto& jet = jets[i];
        if (io_absDphi(lead_phi,jet.phi())>dphi_min) {
            for (int k=0; k<5; ++k) hopper[k] = _in[k];
            hopper[5] = jets[0].perp();
            hopper[6] = jet.perp();
            hopper[7] = (hopper[5]-hopper[6])/(hopper[5]+hopper[6]);
            data->Fill(hopper);
            break;
        }
    }
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
        io_scaleByBinWidth(hg);
    } else {
        hg = _hg;
        hg->SetName(ioUniqueName());
    }
    if (norm != 0.) {
        hg->Scale(1./norm);
    }
    return hg;
};




ioTrackSparse::ioTrackSparse(THnSparseD* _data_track, THnSparseD* _data_trig) :
    data_trig { _data_trig }, data_track { _data_track } {};
ioTrackSparse::ioTrackSparse(const char* bin_file, const char* tag) {
    TString s_tag = tag;
    int i_bins = (s_tag.Contains("_10")) ? 10 : 3;
    ioBinVec bin_EAbbc { bin_file, Form("EAbbc_%ibin",i_bins) };
    ioBinVec bin_EAtpc { bin_file, Form("EAtpc_%ibin",i_bins) }; //"EAtpc_3bin" };
    ioBinVec bin_TrigEt    {{ 0.,0.,30,30.}};

    ioBinVec info_ZDCx { bin_file, "zdcX_mu_sigma_206" };
    ioBinVec bin_ZDCx  {{ info_ZDCx[6], info_ZDCx[2], info_ZDCx[3], info_ZDCx[7] }};

    ioBinVec info_vz   { bin_file, "vz_mu_sigma_206" };
    ioBinVec bin_vz    {{ info_vz  [6], info_vz  [2], info_vz  [3], info_vz  [7] }};

    ioBinVec bin_trackPt { bin_file, "trackpt_resolution" };
    ioBinVec bin_absDphi {{ 0., 0., 64., IO_pi }};

    // get the ZDCx bins:

    nbins[0] = bin_EAbbc;
    nbins[1] = bin_EAtpc;
    nbins[2] = bin_TrigEt;
    nbins[3] = bin_ZDCx;
    nbins[4] = bin_vz;
    nbins[5] = bin_trackPt;
    nbins[6] = bin_absDphi;

    data_trig = new THnSparseD(Form("data_trig%s",tag),
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz",
            5, nbins, NULL, NULL);
    data_trig->SetBinEdges(0,bin_EAbbc);
    data_trig->SetBinEdges(1,bin_EAtpc);
    data_trig->SetBinEdges(2,bin_TrigEt);
    data_trig->SetBinEdges(3,bin_ZDCx);
    data_trig->SetBinEdges(4,bin_vz);

    data_track = new THnSparseD(Form("data_track%s",tag),
            "tracks;EAbbc;EAtpc;TrigEt;ZDCx;Vz;track #it{p}_{T}",
            7, nbins, NULL, NULL);
    data_track->SetBinEdges(0,bin_EAbbc);
    data_track->SetBinEdges(1,bin_EAtpc);
    data_track->SetBinEdges(2,bin_TrigEt);
    data_track->SetBinEdges(3,bin_ZDCx);
    data_track->SetBinEdges(4,bin_vz);
    data_track->SetBinEdges(5,bin_trackPt);
    data_track->SetBinEdges(6,bin_absDphi);

    data_track->Sumw2();
};
void ioTrackSparse::ioTrackSparse::write() { 
    data_trig->Write();
    data_track->Write();
};
void ioTrackSparse::fill_trig(
        double EAbbc, double EAtpc, double TrigEt, double ZDCx, double Vz){
    hopper[0] = EAbbc;
    hopper[1] = EAtpc;
    hopper[2] = TrigEt;
    hopper[3] = ZDCx;
    hopper[4] = Vz;
    data_trig->Fill(hopper);
};
void ioTrackSparse::fill_trackpt_absDphi(double trackpt, double absDphi) {
    hopper[5] = trackpt;
    hopper[6] = absDphi;
    data_track->Fill(hopper);
};
void ioTrackSparse::range_axes (int i_axis, int i0, int i1) {
    if (i_axis > 6) throw std::runtime_error(
        Form("fatal: error in ioTrackSparse, axis(%i) not valid, but by <7",
        i_axis)
    );

    n_triggers = -1.;
    data_track ->GetAxis(i_axis)->SetRange(i0, i1);
    if (i_axis < 5) data_trig->GetAxis(i_axis)->SetRange(i0, i1);
};

TH1D* ioTrackSparse::hg_axis (int i_axis, double norm, bool use_track_data){
    TH1D* hg;
    TH1D* _hg = (TH1D*) (use_track_data ? data_track : data_trig)->Projection(i_axis,"E");
    if (bins != nullptr) {
        hg = (TH1D*) _hg->Rebin(*bins, ioUniqueName(), *bins);
        delete _hg;
        io_scaleByBinWidth(hg);
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
