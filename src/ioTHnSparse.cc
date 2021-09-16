#include "ioTHnSparse.h"
#include "io_fnc.h"
#include "ioClass.h"

ioJetSpectraSparse::ioJetSpectraSparse(const char* bin_file, const char* tag) {
    ioBinVec bin_EAbbc { bin_file, "EAbbc_3bin" };
    ioBinVec bin_EAtpc { bin_file, "EAtpc_3bin" };
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

    data_trig = new THnSparseD(Form("data_trigs%s",tag),
            "triggers;EAbbc;EAtpc;TrigEt;ZDCx;Vz",
            5, nbins, NULL, NULL);
    data_trig->SetBinEdges(0,bin_EAbbc);
    data_trig->SetBinEdges(1,bin_EAtpc);
    data_trig->SetBinEdges(2,bin_TrigEt);
    data_trig->SetBinEdges(3,bin_ZDCx);
    data_trig->SetBinEdges(4,bin_vz);

    data_jet = new THnSparseD(Form("data_jet%s",tag),
            "jets;EAbbc;EAtpc;TrigEt;ZDCx;Vz;Jet #it{p}_{T}",
            6, nbins, NULL, NULL);
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
void ioJetSpectraSparse::range_EAbbc (int i0, int i1) {
    data_jet->GetAxis(0)->SetRange(i0, i1);
};
void ioJetSpectraSparse::range_EAtpc (int i0, int i1) {
    data_jet->GetAxis(1)->SetRange(i0, i1);
};
void ioJetSpectraSparse::range_TrigEt (int i0, int i1) {
    data_jet->GetAxis(2)->SetRange(i0, i1);
};
void ioJetSpectraSparse::range_ZDCx (int i0, int i1) {
    data_jet->GetAxis(3)->SetRange(i0, i1);
};
void ioJetSpectraSparse::range_Vz (int i0, int i1) {
    data_jet->GetAxis(4)->SetRange(i0, i1);
};
void ioJetSpectraSparse::range_JetPt (int i0, int i1) {
    data_jet->GetAxis(5)->SetRange(i0, i1);
};
void ioJetSpectraSparse::range_absDphi (int i0, int i1) {
    data_jet->GetAxis(6)->SetRange(i0, i1);
};
void ioJetSpectraSparse::range8_absDphi (int i) {
    int i0 = 8*i;
    int i1 = i0+7;
    data_jet->GetAxis(6)->SetRange(i0, i1);
};

TH1D* ioJetSpectraSparse::hg_EAbbc (){
    TH1D* hg = (TH1D*) data_jet->Projection(0,"E");
    hg->SetName(ioUniqueName());
    return hg;
};
TH1D* ioJetSpectraSparse::hg_EAtpc (){
    TH1D* hg = (TH1D*) data_jet->Projection(1,"E");
    hg->SetName(ioUniqueName());
    return hg;
};
TH1D* ioJetSpectraSparse::hg_TrigEt (){
    TH1D* hg = (TH1D*) data_jet->Projection(2,"E");
    hg->SetName(ioUniqueName());
    return hg;
};
TH1D* ioJetSpectraSparse::hg_ZDCx (){
    TH1D* hg = (TH1D*) data_jet->Projection(3,"E");
    hg->SetName(ioUniqueName());
    return hg;
};
TH1D* ioJetSpectraSparse::hg_Vz (){
    TH1D* hg = (TH1D*) data_jet->Projection(4,"E");
    hg->SetName(ioUniqueName());
    return hg;
};
TH1D* ioJetSpectraSparse::hg_JetPt (){
    TH1D* hg = (TH1D*) data_jet->Projection(5,"E");
    hg->SetName(ioUniqueName());
    return hg;
};
TH1D* ioJetSpectraSparse::hg_JetPt64 (int i0, int i1){
    range_absDphi(i0,i1); 
    TH1D* hg = (TH1D*) data_jet->Projection(5,"E");
    hg->SetName(ioUniqueName());
    return hg;
};
TH1D* ioJetSpectraSparse::hg_JetPt8 (int i) {
    range8_absDphi(i); 
    TH1D* hg = (TH1D*) data_jet->Projection(5,"E");
    hg->SetName(ioUniqueName());
    return hg;
};
