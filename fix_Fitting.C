void fix_Fitting() {
    using namespace RooFit;
    TFile *f = TFile::Open("histo.root");
    TH1D *h = (TH1D*)f->Get("InvMass");
    RooRealVar x("x", "Z mass", 60, 120);
    RooDataHist data("data", "dataset with Z mass", x, h);
    
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

    RooRealVar nsig("nsig", "#signal events", 1000, 0, 1e6);
    RooRealVar nbkg("nbkg", "#background events", 500, 0, 1e6);
    
    RooAddPdf model("model", "Signal + Background", RooArgList(cb, cheb_bkg), RooArgList(nsig, nbkg));
    
    RooFitResult* fitRes = model.fitTo(data, Save(), Extended(kTRUE));
    fitRes->Print("v");
    
    RooPlot* frame = x.frame(Title("Crystal ball + Chebychev(second order)"));
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components(cheb_bkg), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components(cb), LineStyle(kDotted), LineColor(kGreen));
    
    TCanvas* c = new TCanvas("c_roofit", "RooFit Zmass Fit", 800, 600);
    frame->Draw();
    //c->SaveAs("Zfit_Exponential.png");
    std::cout << "Entries: " << data.sumEntries() << std::endl;    
    std::cout << "Signal yield: " << nsig.getVal() << " ± " << nsig.getError() << std::endl;
    std::cout << "Background yield: " << nbkg.getVal() << " ± " << nbkg.getError() << std::endl;
}
