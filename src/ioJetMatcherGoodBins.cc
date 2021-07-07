#include "ioJetMatcherGoodBins.h"
/* #include "ioJetMatcher.h" */
#include <iostream>
#include "io_fnc.h"
#include "TH1D.h"
#include "TH2D.h"

ioJetMatcherGoodBins::ioJetMatcherGoodBins (
            const char* _name, 
            const char* bin_file, 
            const char* tag_M, 
            const char* tag_T, 
            const char* goodbin_file,
            ioOptMap opt,
            ioOptMap dict
    ) :
    name {_name}, 
    _rand{},
    options { dict + opt},
    nXsec  { options["nXsec"] },
    jet_R2 { (double)(options["R"].val()) * (double)(options["R"].val()) },
    hg_pthb_cnt { Form("pthg_cnt_%s",_name),"Counter;#hat{#it{p}}_{T}-bin;N_{events}",
        nXsec, -0.5, (double) nXsec-0.5 },
    hg_pthb_cnt_A { Form("pthg_cnt_%s_A",_name),"Counter;#hat{#it{p}}_{T}-bin;N_{events}",
        nXsec, -0.5, (double) nXsec-0.5 },
    hg_pthb_cnt_B { Form("pthg_cnt_%s_B",_name),"Counter;#hat{#it{p}}_{T}-bin;N_{events}",
        nXsec, -0.5, (double) nXsec-0.5 },
    debug{ (bool)options["debug"] },
    ratio_AtoB { options["ratio_AtoB"].val() },
    hg1_R2_match {
        Form("hg1_R2_match_%s",_name), 
        "Matched jets; #sqrt((#Delta#phi)^2+(#delta#eta)^2)",
        100, 0., 0.2 
    },
    limit_to_goodbins {strcmp(goodbin_file,"")!=0}, 
    in_T_bounds{ bin_file, tag_T },
    in_M_bounds{ bin_file, tag_M }

{
    cout << " debug: " << debug << endl;
    if (debug) {
        cout << " nXsec : " << nXsec << endl;
        cout << " jet_R2 : " << jet_R2 << endl;
        cout << " jet_R2 : " << jet_R2 << endl;
    }
    /* cout << " STARTED 0" << endl; */
    ioBinVec bins_M { bin_file, tag_M };
    /* cout << " STARTED 1" << endl; */
    ioBinVec bins_T { bin_file, tag_T };
    /* cout << " STARTED 2" << endl; */
    for (int i{0}; i<nXsec; ++i) {
        v_miss.push_back({Form("miss_%i__%s",i,_name),
            Form("Misses for #hat{#it{p}}_{T}-bin %i;#it{p}_{T}-miss;N",i),
            bins_T, bins_T});
        A_miss.push_back({Form("A_miss_%i__%s",i,_name),
            Form("Misses for #hat{#it{p}}_{T}-bin %i;#it{p}_{T}-miss;N",i),
            bins_T, bins_T});
        B_miss.push_back({Form("B_miss_%i__%s",i,_name),
            Form("Misses for #hat{#it{p}}_{T}-bin %i;#it{p}_{T}-miss;N",i),
            bins_T, bins_T});
        v_T.push_back({Form("T_%i__%s",i,_name),
            Form("Truth jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T}-true;N",i),
            bins_T, bins_T});
        v_M.push_back({Form("M_%i__%s",i,_name),
            Form("Measured jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T};N",i),
            bins_M, bins_M});
        v_fake.push_back({Form("fake_%i__%s",i,_name),
            Form("Fake jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T};N",i),
            bins_M, bins_M});
        A_fake.push_back({Form("A_fake_%i__%s",i,_name),
            Form("Fake jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T};N",i),
            bins_M, bins_M});
        B_fake.push_back({Form("B_fake_%i__%s",i,_name),
            Form("Fake jets for #hat{#it{p}}_{T}-bin %i;#it{p}_{T};N",i),
            bins_M, bins_M});
        v_match.push_back({Form("match_%i__%s",i,_name),
            Form("Matched jets for #hat{#it{p}}_{T}-bin %i;"
                "#it{p}_{T}-measured;#it{p}_{T}-truth;",i),
            bins_M, bins_M, bins_T, bins_T});
        A_match.push_back({Form("A_match_%i__%s",i,_name),
            Form("Matched jets for #hat{#it{p}}_{T}-bin %i;"
                "#it{p}_{T}-measured;#it{p}_{T}-truth;",i),
            bins_M, bins_M, bins_T, bins_T});
        B_match.push_back({Form("B_match_%i__%s",i,_name),
            Form("Matched jets for #hat{#it{p}}_{T}-bin %i;"
                "#it{p}_{T}-measured;#it{p}_{T}-truth;",i),
            bins_M, bins_M, bins_T, bins_T});
    }
    if (limit_to_goodbins) {
        ioBinVec gbins_T { goodbin_file, options["tag_goodbins_T"].c_str() };
        ioBinVec gbins_M { goodbin_file, options["tag_goodbins_M"].c_str() };

        hg_goodmatch = new TH2D("hg_goodbins__local",";Measured;Truth",
                gbins_M, gbins_M, gbins_T, gbins_T );
        for (int i{0}; i<nXsec; ++i) {
            goodmatchbins.push_back({goodbin_file, Form("match_pthb_%i",i)});
        }

        hg_goodmiss = new TH1D("hg_goodmiss__local",";Measured;Truth",
                gbins_T, gbins_T );
        for (int i{0}; i<nXsec; ++i) {
            goodmissbins.push_back({goodbin_file, Form("miss_pthb_%i",i)});
        }

        hg_goodfake = new TH1D("hg_goodfake__local",";Measured;Truth",
                gbins_M, gbins_M );
        for (int i{0}; i<nXsec; ++i) {
            goodfakebins.push_back({goodbin_file, Form("miss_pthb_%i",i)});
        }
    }
    /* cout << " STARTED 4" << endl; */
}
void ioJetMatcherGoodBins::addjet_MC(float eta, float phi, float pT) {
    if (in_T_bounds(pT)) data_MC.push_back({eta,phi,pT});
};
void ioJetMatcherGoodBins::addjet_reco(float eta, float phi, float pT) {
    if (in_M_bounds(pT)) data_reco.push_back({eta,phi,pT});
};

bool ioJetMatcherGoodBins::do_matching(int iXsecBin) {
    /* cout << " do_matching" << endl; */
    vector<double> fakes;
    vector<double> misses;
    vector<pair<double,double>> matches;
    vector<double> v_delta_R2;

    for (auto& MC : data_MC) {
        bool found_match {false};
        for (auto& reco : data_reco) {
            if (reco.is_matched) continue;
            double delta_R2 { MC(reco,jet_R2) };
            if (delta_R2 != 0) {
                // cut the event if matched pair is out of bounds
                if (limit_to_goodbins) {
                    int ibin { hg_goodmatch->FindBin(reco.pT, MC.pT) };
                    if (!goodmatchbins[iXsecBin](ibin)) return clear_MCreco();
                }
                found_match = true;
                reco.is_matched = true;
                matches.push_back({reco.pT, MC.pT});
                v_delta_R2.push_back(delta_R2);
                break;
            }
        }
        if (!found_match) { 
            if (limit_to_goodbins) {
                int ibin { hg_goodmiss->FindBin(MC.pT) };
                if (!goodmissbins[iXsecBin](ibin)) return clear_MCreco();
            }
            misses.push_back(MC.pT); 
        }
    }
    for (auto& reco : data_reco) {
        if (!reco.is_matched) {
            if (limit_to_goodbins) {
                int ibin { hg_goodfake->FindBin(reco.pT) };
                if (!goodfakebins[iXsecBin](ibin)) return clear_MCreco();
            }
            fakes.push_back(reco.pT);
        }
    }
    /* cout << " GOOD " << endl; */
    clear_MCreco();

    // now that is it not an outlier event, fill the data
    hg_pthb_cnt.Fill((double)iXsecBin);
    bool fillA = _rand.Uniform() < ratio_AtoB;
    if (fillA) hg_pthb_cnt_A.Fill((double)iXsecBin);
    else       hg_pthb_cnt_B.Fill((double)iXsecBin);
    
    int n{0};
    for (auto& match : matches) {
        /* cout << " auto " << match.first << " " << match.second << endl; */
        v_match[iXsecBin].Fill(match.first, match.second);
        if (fillA) A_match[iXsecBin].Fill(match.first, match.second);
        else       B_match[iXsecBin].Fill(match.first, match.second);
        v_T[iXsecBin].Fill(match.second);
        v_M[iXsecBin].Fill(match.first);
        hg1_R2_match.Fill(v_delta_R2[n]);
        ++n;
    }
    for (auto& miss : misses) {
        v_miss[iXsecBin].Fill(miss);
        v_T[iXsecBin].Fill(miss);
        if (fillA) A_miss[iXsecBin].Fill(miss);
        else       B_miss[iXsecBin].Fill(miss);
    }
    for (auto& fake : fakes) {
        v_fake[iXsecBin].Fill(fake);
        v_M[iXsecBin].Fill(fake);
        if (fillA) A_fake[iXsecBin].Fill(fake);
        else       B_fake[iXsecBin].Fill(fake);
    };
    return false;
};
bool ioJetMatcherGoodBins::clear_MCreco() {
    /* cout << " Clearing MCreco " << endl; */
    data_MC.clear();
    data_reco.clear();
    return true;
};
void ioJetMatcherGoodBins::write() {
    hg_pthb_cnt.Write();
    hg_pthb_cnt_A.Write();
    hg_pthb_cnt_B.Write();
    hg1_R2_match.Write();

    for (auto& h : v_T)      h.Write();
    for (auto& h : v_M)      h.Write();
    for (auto& h : v_fake)   h.Write();
    for (auto& h : v_miss)   h.Write();
    for (auto& h : v_match)  h.Write();

    for (auto& h : A_fake)   h.Write();
    for (auto& h : A_miss)   h.Write();
    for (auto& h : A_match)  h.Write();

    for (auto& h : B_fake)   h.Write();
    for (auto& h : B_miss)   h.Write();
    for (auto& h : B_match)  h.Write();

};
