void create_pt_weight() {
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TFile *file_lo = TFile::Open("test_lo/histos_lo.root");
    TFile *file_nlo = TFile::Open("test_nlo/histos_nlo.root");

    std::string var = "Z_pt";

    TH1D* hist_lo_orig = (TH1D*)file_lo->Get(var.c_str());
    TH1D* hist_nlo_orig = (TH1D*)file_nlo->Get(var.c_str());

    // --- Create the Ratio Histogram (The "pT weight") ---
    TH1D* h_ratio = (TH1D*)hist_lo_orig->Clone("h_ratio");
    if(h_ratio->GetSumw2N()==0) h_ratio->Sumw2();

    TH1D* hist_nlo_scaled_for_ratio = (TH1D*)hist_nlo_orig->Clone("hist_nlo_scaled_for_ratio");
    if(hist_nlo_scaled_for_ratio->GetSumw2N()==0) hist_nlo_scaled_for_ratio->Sumw2();

    if (hist_nlo_scaled_for_ratio->Integral() > 0) {
        hist_nlo_scaled_for_ratio->Scale(hist_lo_orig->Integral() / hist_nlo_scaled_for_ratio->Integral());
    }
    h_ratio->Divide(hist_lo_orig, hist_nlo_scaled_for_ratio, 1.0, 1.0, "");

    TFile *pt_weight = new TFile("z_pt_weight.root", "RECREATE");
    h_ratio->SetName("pt_weight");
    h_ratio->SetTitle("pT Weight for Z boson");
    h_ratio->Write();
    pt_weight->Close();
    delete pt_weight;
    std::cout << "pT weight saved to z_pt_weight.root" << std::endl;

    file_lo->Close();
    file_nlo->Close();
    delete file_lo;
    delete file_nlo;
}
