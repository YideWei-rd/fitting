void fix_Fitting() {
    using namespace RooFit;
    TFile *f = TFile::Open("histo.root");
    TH1D *h = (TH1D*)f->Get("InvMass");
    RooRealVar x("x", "Z mass", 60, 120);
    RooDataHist data("data", "dataset with Z mass", x, h);
    int totEvents = data.sumEntries();

    RooRealVar mean("mean", "mean", 91, 85, 95);
    RooRealVar sigmaL("sigmaL", "left width", 2.5, 0.5, 5.0);
    RooRealVar sigmaR("sigmaR", "right width", 3.5, 0.5, 5.0);
    RooRealVar alphaL("alphaL", "alphaL", 1.5, 0.1, 10.0);
    RooRealVar nL("nL", "nL", 2.0, 0.1, 20.0);
    RooRealVar alphaR("alphaR", "alphaR", 2.0, 0.1, 10.0);
    RooRealVar nR("nR", "nR", 3.0, 0.1, 20.0);

    RooCrystalBall cb("cb", "Double-sided Crystal Ball", x, mean, sigmaL, sigmaR, alphaL, nL, alphaR, nR);

//    RooRealVar slope("slope", "slope", -0.0378565);
//    slope.setConstant(true);
//    RooExponential expo("expo", "expo", x, slope);

    RooRealVar c0("c0", "c0", -0.939952);// -0.825192 Chebychev line
    c0.setConstant(true);                // -0.939952, 0.317173 Chebychev second order
    RooRealVar c1("c1", "c1", 0.317173); // -0.0378565  exponential slope
    c1.setConstant(true);
    RooChebychev cheb_bkg("cheb_bkg", "Chebychev background", x, RooArgList(c0, c1));

    RooRealVar nsig("nsig", "#signal events", 0.9 * totEvents, 0, totEvents);
    RooRealVar nbkg("nbkg", "#background events", 0.1 * totEvents, 0, totEvents);

    RooAddPdf model("model", "Signal + Background", RooArgList(cb, cheb_bkg), RooArgList(nsig, nbkg));
    RooAbsReal* nll = model.createNLL(data, Extended(kTRUE));

    RooFitResult* fitRes = model.fitTo(data, Save(), Extended(kTRUE));
    fitRes->Print("v");

    RooPlot* frame = x.frame(Title("Crystal ball + Chebychev(second order)"));
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(cheb_bkg), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(cb), LineStyle(kDotted), LineColor(kGreen));

    TCanvas* c = new TCanvas("c", "c", 800, 600);
    frame->Draw();

    int nFloatingParams = fitRes->floatParsFinal().getSize();
    double chi2 = frame->chiSquare(nFloatingParams);
    int ndf = frame->residHist()->GetN() - nFloatingParams;
    TString tex1 = TString::Format("Entries: %d", totEvents);
    TString tex2 = TString::Format("Sig yield: %.0f #pm %.0f", nsig.getVal(), nsig.getError());
    TString tex3 = TString::Format("Bkg yield: %.0f #pm %.0f", nbkg.getVal(), nbkg.getError());
    TString tex4 = TString::Format("#chi^{2} / NDF = %.2f / %d", chi2*ndf, ndf);

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.DrawLatex(0.6,0.8,tex1);
    latex.DrawLatex(0.6,0.75,tex2);
    latex.DrawLatex(0.6,0.7,tex3);
    latex.DrawLatex(0.6,0.65,tex4);

    c->SaveAs("Zmass.png");

    x.setRange("signal_region", 80, 100);
    RooAbsReal* sig_frac = cb.createIntegral(x, NormSet(x), Range("signal_region"));
    RooAbsReal* bkg_frac = cheb_bkg.createIntegral(x, NormSet(x), Range("signal_region"));
    double sig_frac_val = sig_frac->getVal();
    double bkg_frac_val = bkg_frac->getVal();
    double nsig_subrange = nsig.getVal() * sig_frac_val;
    double nbkg_subrange = nbkg.getVal() * bkg_frac_val;

    double frac_sig_err = sig_frac->getPropagatedError(*fitRes);
    double frac_bkg_err = bkg_frac->getPropagatedError(*fitRes);

    // Use the standard error propagation formula for a product: N_sub = N_tot * f
    // (err_N_sub / N_sub)^2 = (err_N_tot / N_tot)^2 + (err_f / f)^2
    double nsig_subrange_err = nsig_subrange * sqrt( pow(nsig.getError() / nsig.getVal(), 2) + pow(frac_sig_err / sig_frac_val, 2) );
    double nbkg_subrange_err = nbkg_subrange * sqrt( pow(nbkg.getError() / nbkg.getVal(), 2) + pow(frac_bkg_err / bkg_frac_val, 2) );

    cout << "----------------------------------------------------" << endl;
    cout << "Calculation for mass range [80, 100] GeV:" << endl;
    cout << "----------------------------------------------------" << endl;
    cout << "Fraction of signal in sub-range: " << sig_frac_val << " +/- " << frac_sig_err << endl;
    cout << "Fraction of background in sub-range: " << bkg_frac_val << " +/- " << frac_bkg_err << endl;
    cout << "----------------------------------------------------" << endl;
    cout << "Signal yield in sub-range: " << nsig_subrange << " +/- " << nsig_subrange_err << " events" << endl;
    cout << "Background yield in sub-range: " << nbkg_subrange << " +/- " << nbkg_subrange_err << " events" << endl;
    cout << "----------------------------------------------------" << endl;

    RooPlot* frame_nll = nsig.frame(Title("NLL scan of nsig"), Range(nsig.getVal() - 5*nsig.getError(), nsig.getVal() + 5*nsig.getError()));
    nll->plotOn(frame_nll, ShiftToZero());
    frame_nll->Draw();
    c->SaveAs("NLL.png");
}
