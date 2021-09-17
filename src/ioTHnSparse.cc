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

TH1D* ioJetSpectraSparse::hg_axis (int i_axis, double norm){
    TH1D* hg;
    TH1D* _hg = (TH1D*) data_jet->Projection(i_axis,"E");
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
