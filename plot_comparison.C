void plot_comparison() {
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TFile *file_lo = TFile::Open("test_lo/histos_lo.root");
    TFile *file_nlo = TFile::Open("test_output/Z_plots.root");

    std::vector<std::string> variables = {"InvMass", "Z_pt", "Z_eta", "Z_phi"};
    std::vector<std::string> x_titles = {
        "Invariant Mass M_{#mu#mu} [GeV]",
        "Z p_{T} [GeV]",
        "Z #eta",
        "Z #phi"
    };

    for (size_t i = 0; i < variables.size(); ++i) {
        std::string var = variables[i];
        TH1D* hist_lo = (TH1D*)file_lo->Get(var.c_str());
        TH1D* hist_nlo = (TH1D*)file_nlo->Get(var.c_str());

        TH1D* hist_lo_norm = (TH1D*)hist_lo->Clone(("hist_lo_norm_" + var).c_str());
        TH1D* hist_nlo_norm = (TH1D*)hist_nlo->Clone(("hist_nlo_norm_" + var).c_str());
        if (hist_lo_norm->Integral() > 0) hist_lo_norm->Scale(1.0 / hist_lo_norm->Integral());
        if (hist_nlo_norm->Integral() > 0) hist_nlo_norm->Scale(1.0 / hist_nlo_norm->Integral());

        TCanvas* c = new TCanvas(("c_" + var).c_str(), var.c_str(), 1000, 800);
        TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        pad1->SetBottomMargin(0.02);
        pad1->Draw();
        pad1->cd();

        // --- Style setup ---
        hist_lo_norm->SetLineColor(kBlue + 1);
        hist_lo_norm->SetLineWidth(2);
        hist_lo_norm->SetFillColorAlpha(kBlue, 0.35);
        hist_nlo_norm->SetLineColor(kBlack);
        hist_nlo_norm->SetMarkerColor(kBlack);
        hist_nlo_norm->SetMarkerStyle(20);
        hist_nlo_norm->SetMarkerSize(0.8);

        hist_lo_norm->GetYaxis()->SetTitle("Normalized Events");
        hist_lo_norm->GetYaxis()->SetTitleSize(0.05);
        hist_lo_norm->GetYaxis()->SetTitleOffset(0.8);
        hist_lo_norm->GetYaxis()->SetLabelSize(0.04);
        hist_lo_norm->GetXaxis()->SetLabelSize(0);

        if (var.find("pt") != std::string::npos) {
            pad1->SetLogy(1);
            hist_lo_norm->SetMaximum(hist_lo_norm->GetMaximum() * 10.0);
            hist_lo_norm->SetMinimum(1e-5);
        } else {
            pad1->SetLogy(0);
            hist_lo_norm->SetMaximum(1.35 * std::max(hist_lo_norm->GetMaximum(), hist_nlo_norm->GetMaximum()));
            hist_lo_norm->SetMinimum(0);
        }

        hist_lo_norm->Draw("HIST");
        hist_nlo_norm->Draw("E P SAME");

        TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
        legend->AddEntry(hist_lo_norm, "Z #rightarrow #mu#mu (LO)", "lf");
        legend->AddEntry(hist_nlo_norm, "Z #rightarrow ll (NLO)", "ep");
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->Draw();

        c->cd();
        TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.3);
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();

        TH1D* h_ratio = (TH1D*)hist_nlo_norm->Clone(("h_ratio_" + var).c_str());
        h_ratio->Divide(hist_lo_norm);

        h_ratio->SetMarkerStyle(20);
        h_ratio->SetLineColor(kBlack);
        h_ratio->GetYaxis()->SetTitle("NLO / LO");
        h_ratio->GetXaxis()->SetTitle(x_titles[i].c_str());
        h_ratio->GetYaxis()->SetRangeUser(0.0, 2.0);
        h_ratio->SetTitle("");

        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.12);
        h_ratio->GetYaxis()->SetTitleOffset(0.3);
        h_ratio->GetYaxis()->SetLabelSize(0.1);
        h_ratio->GetYaxis()->CenterTitle();
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetTitleOffset(1.0);
        h_ratio->GetXaxis()->SetLabelSize(0.1);

        h_ratio->Draw("ep");
        TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0,
                                h_ratio->GetXaxis()->GetXmax(), 1.0);
        line->SetLineStyle(2);
        line->SetLineColor(kGray + 2);
        line->Draw();

        c->SaveAs(TString::Format("verification_final_%s.png", var.c_str()));
        delete c;
    }

    file_lo->Close();
    file_nlo->Close();
}

