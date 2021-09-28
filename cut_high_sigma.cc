TH2D* cut_high_sigmaX(TH2D* hg, double n_sigma=5.) {
    // cut all bin content above n_sigma from mean

    if (n_sigma==0) return hg;
    TAxis* y_axis = hg->GetYaxis();
    TAxis* x_axis = hg->GetXaxis();
    double pre_sum = hg->Integral();
    for (int y=1;y<y_axis->GetNbins();++y) {
        TH1D* proj = hg->ProjectionX(ioUniqueName(),y,y);
        double mu = proj->GetMean();
        double sigma = proj->GetStdDev();
        int i_bin = x_axis->FindBin(mu+n_sigma*sigma+8);
        for (int i=i_bin+1;i<=x_axis->GetNbins();++i) {
            if (hg->GetBinContent(i,y)!=0) {
                /* cout << " THIS BIN: " << x_axis->GetBinCenter(i) << ":" << y_axis->GetBinCenter(y) << "  with cutoff: " << mu+n_sigma*sigma << "  content: " << hg->GetBinContent(i,y) << endl; */
                hg->SetBinContent(i,y,0.);
                hg->SetBinError(i,y,0.);
            }
        }

    }
    double post_sum = hg->Integral();
    cout << " percent cut: " << 100.*(pre_sum-post_sum)/pre_sum << endl;
    return hg;
};
