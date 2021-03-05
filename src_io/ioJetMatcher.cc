#include "ioJetMatcher.h"
#include <iostream>
#include "io_fnc.h"
#include "TH1D.h"
#include "TH2D.h"

double* ioEdges_pAuJet_prelim_13bins(){
    return  new double[14]{0.,1.,2.,3.,4.,5.,6.,8.,10.,15.,20.,25.,36.,50.};
};

// jet_match_floats
ioJetMatcher_float::ioJetMatcher_float(float _eta, float _phi, float _pT) :
    eta{_eta}, phi{_phi}, pT{_pT}, is_matched{false} {};
ioJetMatcher_float::ioJetMatcher_float() :
    eta{0.}, phi{0.}, pT{0.}, is_matched{false} {};

double ioJetMatcher_float::operator()(const ioJetMatcher_float& B, const float jetR2) {
    const float deta {eta-B.eta};
    const float dphi {io_dphi(phi, B.phi)};
    const float dist2 = deta*deta + dphi*dphi;
    if (dist2 == 0.) return 0.0001;
    else if (dist2 > jetR2) return 0.;
    else return dist2;
};

bool operator==(const ioJetMatcher_float& L, const ioJetMatcher_float R) 
    { return L.pT == R.pT; };

bool operator<(const ioJetMatcher_float& L, const ioJetMatcher_float R) 
    { return L.pT < R.pT; };

bool operator>(const ioJetMatcher_float& L, const ioJetMatcher_float R) 
    { return L.pT > R.pT; };

ioJetMatcher::ioJetMatcher(const char* _tag, int nbins, double lo_bin, double hi_bin, float _jet_R) {
    tag = _tag;
    /* cout << " making tag: " << _tag << endl; */
    /* double* edges = new double(nbins+1); */
    /* double wide  = (hi_bin-lo_bin)/nbins; */
    /* for (int i{0};i<=nbins;++i) edges[i] = lo_bin + i*wide; */
    hg_MC  = TH1D( Form("%s_MC",_tag),
            Form("%s MC jets (includes misses;#it{p}_{T};dN_{jets}/d#it{p}_{T}",_tag),
            nbins, lo_bin, hi_bin ); 
    hg_reco  = TH1D ( Form("%s_reco",_tag),
            Form("%s reco jets (includes fakes;#it{p}_{T};dN_{jets}/d#it{p}_{T}",_tag),
            nbins, lo_bin, hi_bin ); 
    hg_response  = TH2D ( Form("%s_responseM",_tag),
            Form("%s response matrix for matched jets;"
                 "#it{p}_{T} reco(-nstructed);#it{p}_{T} MC",_tag),
            nbins, lo_bin, hi_bin, nbins, lo_bin, hi_bin ); 
	jet_R2 = _jet_R*_jet_R;
    data_MC.clear();
    data_reco.clear();

    /* cout << " hg_MC: " << hg_MC.GetXaxis()->GetNbins() << endl; */
    /* for (int i{1}; i<= hg_MC.GetXaxis()->GetNbins(); ++i) */ 
    /*     cout << " ("<<i<<") " << hg_MC.GetXaxis()->GetBinCenter(i) ; */
    /* cout << endl; */

    /* cout << " hg_response: " << hg_response.GetXaxis()->GetNbins() << endl; */
    /* for (int i{1}; i<= hg_response.GetXaxis()->GetNbins(); ++i) */ 
    /*     cout << " ("<<i<<") " << hg_response.GetXaxis()->GetBinCenter(i) ; */
    /* cout << endl; */
    /* cout << " hg_response: " << hg_response.GetYaxis()->GetNbins() << endl; */
    /* for (int i{1}; i<= hg_response.GetYaxis()->GetNbins(); ++i) */ 
    /*     cout << " ("<<i<<") " << hg_response.GetYaxis()->GetBinCenter(i) ; */
    /* cout << endl; */
};

ioJetMatcher::ioJetMatcher(const char* _tag, int nbins, double* edges, float _jet_R) :
    hg_MC { Form("%s_MC",_tag),
            Form("%s MC jets (includes misses;#it{p}_{T};dN_{jets}/d#it{p}_{T}",_tag),
            nbins,edges },
    hg_reco { Form("%s_reco",_tag),
            Form("%s reco jets (includes fakes;#it{p}_{T};dN_{jets}/d#it{p}_{T}",_tag),
            nbins,edges },
    hg_response { Form("%s_responseM",_tag),
            Form("%s response matrix for matched jets;"
                 "#it{p}_{T} reco(-nstructed);#it{p}_{T} MC",_tag),
            nbins,edges,nbins,edges },
	jet_R2{_jet_R*_jet_R},
    data_MC{}, data_reco{}, tag{_tag}
{};
    /* cout << " making tag: " << tag << endl; */
    /* cout << " nbins: " << nbins << endl; */
    /* for (int i{0}; i<=nbins; ++i) cout << " ("<<i<<") " << edges[i]; */
    /* cout << endl; */

    /* cout << " hg_MC: " << hg_MC.GetXaxis()->GetNbins() << endl; */
    /* for (int i{1}; i<= hg_MC.GetXaxis()->GetNbins(); ++i) */ 
        /* cout << " ("<<i<<") " << hg_MC.GetXaxis()->GetBinCenter(i) ; */
    /* cout << endl; */

    /* cout << " hg_response: " << hg_response.GetXaxis()->GetNbins() << endl; */
    /* for (int i{1}; i<= hg_response.GetXaxis()->GetNbins(); ++i) */ 
        /* cout << " ("<<i<<") " << hg_response.GetXaxis()->GetBinCenter(i) ; */
    /* cout << endl; */
    /* cout << " hg_response: " << hg_response.GetYaxis()->GetNbins() << endl; */
    /* for (int i{1}; i<= hg_response.GetYaxis()->GetNbins(); ++i) */ 
        /* cout << " ("<<i<<") " << hg_response.GetYaxis()->GetBinCenter(i) ; */
    /* cout << endl; */
/* }; */

void ioJetMatcher::addjet_MC(float eta, float phi, float pT) {
    data_MC.push_back({eta,phi,pT});
};

void ioJetMatcher::addjet_reco(float eta, float phi, float pT) {
    data_reco.push_back({eta,phi,pT});
};

void ioJetMatcher::do_matching_highfirst(double W) {
    //sort into high to low pT order
    /* std::sort(data_MC.begin(), data_MC.end(), std::greater<ioJetMatcher_float>()); */
    /* std::sort(data_reco.begin(), data_reco.end(), std::greater<ioJetMatcher_float>()); */

    for (auto& MC : data_MC) {
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (reco.is_matched) continue;
            if (MC(reco,0.16)) {
                reco.is_matched = true;
                hg_MC.Fill(MC.pT,W);
                hg_reco.Fill(reco.pT,W);
                hg_response.Fill(reco.pT, MC.pT,W);
                found_match = true;
                break;
            }
        }
        if (!found_match) hg_MC.Fill(MC.pT,W); // enter a miss
    }
    for (auto& reco : data_reco) {
        if (!reco.is_matched) hg_reco.Fill(reco.pT,W); // enter a fake
    }
    reset();
    // algorithm:
    // starting with the highest to lowest MC jet -- try and match to the highest to lowest reco jets
    // and keep first match

};
void ioJetMatcher::reset() {
    data_MC.clear();
    data_reco.clear();
};
void ioJetMatcher::write(bool with_miss_fakes, bool scale_by_bin_width, bool add_unified2D) {
    /* cout << " z0 " << endl; */
    if (scale_by_bin_width) {
        io_scaleByBinWidth(&hg_MC);
        io_scaleByBinWidth(&hg_reco);
        io_scaleByBinWidth(&hg_response);
    }
    hg_MC.Write();
    hg_reco.Write();
    hg_response.Write();

    if (with_miss_fakes) {
        TH1D* miss = (TH1D*) hg_response.ProjectionY(Form("%s_miss",tag.c_str()));
        miss->Scale(-1.);
        miss->Add(&hg_MC);
        miss->Write();

        TH1D* fakes = (TH1D*) hg_response.ProjectionX(Form("%s_fakes",tag.c_str()));
        fakes->Scale(-1.);
        fakes->Add(&hg_reco);
        fakes->Write();

    /* cout << " z10 " << endl; */
        if (add_unified2D) { // make a TH2D histogram is misses and face in a new, 
                             // dummy bin to the left and bottem of the TH2D.
            TAxis *x_axis = hg_response.GetXaxis();
            TAxis *y_axis = hg_response.GetYaxis();

            int x_nbins = x_axis->GetNbins();
            int y_nbins = y_axis->GetNbins();

            double *x_edges = new double [ x_nbins+2 ];
            double *y_edges = new double [ y_nbins+2 ];
    /* TH1D* __1unified = new TH1D( "__1name__","title;x;y", 13, 0., 50.); */

    /* cout << " z11 " << endl; */
            x_edges[0] = x_axis->GetBinLowEdge(1) - x_axis->GetBinWidth(1);
            for (int i{1}; i<=x_nbins; ++i) x_edges[i] = x_axis->GetBinLowEdge(i);
            x_edges[x_nbins+1] = x_axis->GetBinUpEdge(x_nbins);

            y_edges[0] = y_axis->GetBinLowEdge(1) - y_axis->GetBinWidth(1);
            for (int i{1}; i<=y_nbins; ++i) y_edges[i] = y_axis->GetBinLowEdge(i);
            y_edges[y_nbins+1] = y_axis->GetBinUpEdge(y_nbins);

            TH2D unified { Form("%s_responseFakeMiss",tag.c_str()), 
                Form("%s;Reconstructed (first row is fake jets);Pythia (first column is misses)",
                        tag.c_str()), x_nbins+1, x_edges, y_nbins+1, y_edges };
        /* TH1D* __unified = new TH1D( "name__","title;x;y", 13, 0., 50.); */

            for (int ix{1}; ix <=x_nbins; ++ix) {
                unified.SetBinContent(ix,1, fakes->GetBinContent(ix));
                unified.SetBinError  (ix,1, fakes->GetBinError(ix));
            }
            for (int iy{1}; iy <=y_nbins; ++iy) {
                unified.SetBinContent(1,iy, miss->GetBinContent(iy));
                unified.SetBinError  (1,iy, miss->GetBinError(iy));
            }
            for (int ix{1}; ix <=x_nbins; ++ix)
            for (int iy{1}; iy <=y_nbins; ++iy) {
                unified.SetBinContent(ix+1,iy+1, hg_response.GetBinContent(ix,iy));
                unified.SetBinError  (ix+1,iy+1, hg_response.GetBinError  (ix,iy));
            }
            unified.Write();
        }
    }
};

