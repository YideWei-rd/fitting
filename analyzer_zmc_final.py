import ROOT
import os
import sys
import time

PT_WEIGHT_FILE = "z_pt_weight.root"
PT_WEIGHT_HISTNAME = "pt_weight"
TREENAME = "rootTupleTreeVeryLoose/tree"

start_time = time.time()
ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(16)

file_list_name = sys.argv[1] if len(sys.argv) > 1 else "file_list.txt"
print(f"Using file list: {file_list_name}")

with open(file_list_name, "r") as f:
    file_list = [line.strip() for line in f.read().splitlines() if line.strip()]

ch = ROOT.TChain(TREENAME)
for fname in file_list:
    ch.Add(fname)

df = ROOT.RDataFrame(ch)
total = df.Count().GetValue()
total_h = ROOT.TH1D("total_h", "total_h", 1, 0, 1)
total_h.SetBinContent(1, total)
print(f"Total events: {total}")
df = df.Filter("ROOT::VecOps::Sum(GenPromptParticlePt > 0) == 1", "Keep only gen-level Z->mumu events")
total_Z = df.Count().GetValue()
print(f"Events after gen-level Z->mumu selection: {total_Z}")
# --- Include the existing physics analyzer header ---
ROOT.gInterpreter.Declare('#include "analyzer_zmc.h"')

# Note: ROOT/Python handles string conversion, but using .encode('utf-8') is safer for char* arguments.
ROOT.init_pt_hist(PT_WEIGHT_FILE.encode('utf-8'), PT_WEIGHT_HISTNAME.encode('utf-8'))


df = df.Define("genZpt", "getGenZpt(GenPromptParticlePt)")
df = df.Define("ptw", "getPtWeight(genZpt)")
df = df.Define("weight_total", "ptw*((GenWeight[0] > 0) ? 1.0 : -1.0)")

sum_ptw = df.Sum("weight_total").GetValue()
print(f"Sum of event weights (ptw only) = {sum_ptw:.6g}")

h_ptw = df.Histo1D(("h_ptw","ptw distribution; ptw; entries",100,0,5), "ptw")
h_ptw_hist = h_ptw.GetValue()
print(f"ptw histogram: entries={h_ptw_hist.GetEntries():.0f}, integral={h_ptw_hist.Integral():.6g}")

h_genZpt = df.Histo1D(("h_genZpt_rew","genZpt reweighted; genZpt; sum(ptw)",50,0,200),"genZpt","weight_total")
h_genZpt_hist = h_genZpt.GetValue()
print(f"genZpt reweighted: integral={h_genZpt_rew_hist.Integral():.6g}")

cutflow = ROOT.TH1D("cutflow", "Weighted cutflow (ptw only)", 10, 0, 10)
cutflow.GetXaxis().SetBinLabel(1, "total")
cutflow.SetBinContent(1, sum_ptw)

step_names = []
step_cum_eff = []
step_step_eff = []

df1 = df.Define("rejectedVertices",
                "create_listOfRejectedVertices(IsGoodRecoVertex, SumPt, VectorSumPt, VLTrackNByVert, \
                 Cluster1Size, Cluster2Size, Cluster1Chi2, Cluster1DoF, Cluster2Chi2, Cluster2DoF, TotalChi2, TotalDoF, \
                 VLMuonPt, VLMuonEta, VLMuonPhi, MuonVertInd)")
df1 = df1.Define("rejectedVertices_size", "rejectedVertices.size()")
df1 = df1.Filter("rejectedVertices_size > 0")

w_cut1 = df1.Sum("weight_total").GetValue()
step_names.append("rejectedVertices")
step_cum_eff.append(w_cut1 / sum_ptw if sum_ptw != 0 else 0)
step_step_eff.append(w_cut1 / sum_ptw if sum_ptw != 0 else 0)
cutflow.GetXaxis().SetBinLabel(2, "rejectedVertices")
cutflow.SetBinContent(2, w_cut1)

df2 = df1.Define("isGoodMuon", "select_muonsFromGoodVertices(rejectedVertices, MuonVertInd)")
df2 = df2.Define("isGoodMuon_size", "std::count(isGoodMuon.begin(), isGoodMuon.end(), true)")
df2 = df2.Filter("isGoodMuon_size > 1")

w_cut2 = df2.Sum("weight_total").GetValue()
step_names.append("goodMuon")
step_cum_eff.append(w_cut2 / sum_ptw if sum_ptw != 0 else 0)
step_step_eff.append(w_cut2 / w_cut1 if w_cut1 != 0 else 0)
cutflow.GetXaxis().SetBinLabel(3, "goodMuon")
cutflow.SetBinContent(3, w_cut2)

df3 = df2.Define("goodMuonPtIsoEta", "select_muons_ptIsoEta(isGoodMuon,VLMuonPt,10.,MuonIso03,0.1,VLMuonEta,2.5)")
df3 = df3.Define("goodMuonPtIsoEta_size", "std::count(goodMuonPtIsoEta.begin(), goodMuonPtIsoEta.end(), true)")
df3 = df3.Filter("goodMuonPtIsoEta_size > 1")

w_cut3 = df3.Sum("weight_total").GetValue()
step_names.append("ptIsoEta")
step_cum_eff.append(w_cut3 / sum_ptw if sum_ptw != 0 else 0)
step_step_eff.append(w_cut3 / w_cut2 if w_cut2 != 0 else 0)
cutflow.GetXaxis().SetBinLabel(4, "ptIsoEta")
cutflow.SetBinContent(4, w_cut3)

df4 = df3.Define("muonpairs", "create_pairsOfGoodMuons(goodMuonPtIsoEta)")
df4 = df4.Define("pairsMuonsSamePV", "create_pairsMuonsSamePV(muonpairs, MuonVertInd)")
df4 = df4.Define("pairsMuonsSamePV_size", "pairsMuonsSamePV.size()")
df4 = df4.Filter("pairsMuonsSamePV_size > 0")

w_cut4 = df4.Sum("weight_total").GetValue()
step_names.append("pairsSamePV")
step_cum_eff.append(w_cut4 / sum_ptw if sum_ptw != 0 else 0)
step_step_eff.append(w_cut4 / w_cut3 if w_cut3 != 0 else 0)
cutflow.GetXaxis().SetBinLabel(5, "pairsSamePV")
cutflow.SetBinContent(5, w_cut4)

df5 = df4.Define("pairsMuonsPtAsymmetry", "create_pairsMuonsPtAsymmetry(pairsMuonsSamePV,VLMuonPt,20.)")
df5 = df5.Define("pairsMuonsPtAsymmetry_size", "pairsMuonsPtAsymmetry.size()")
df5 = df5.Filter("pairsMuonsPtAsymmetry_size > 0")

w_cut5 = df5.Sum("weight_total").GetValue()
step_names.append("ptAsymmetry")
step_cum_eff.append(w_cut5 / sum_ptw if sum_ptw != 0 else 0)
step_step_eff.append(w_cut5 / w_cut4 if w_cut4 != 0 else 0)
cutflow.GetXaxis().SetBinLabel(6, "ptAsymmetry")
cutflow.SetBinContent(6, w_cut5)

df6 = df5.Define("SelectedPairsOpCharge", "idx_0f_trackpairs_opcharge(pairsMuonsPtAsymmetry, VLMuonCharge)")
df6 = df6.Define("NumOfPairsTracksOpCharge", "SelectedPairsOpCharge.size()")
df6 = df6.Filter("NumOfPairsTracksOpCharge > 0")

w_cut6 = df6.Sum("weight_total").GetValue()
step_names.append("oppositeCharge")
step_cum_eff.append(w_cut6 / sum_ptw if sum_ptw != 0 else 0)
step_step_eff.append(w_cut6 / w_cut5 if w_cut5 != 0 else 0)
cutflow.GetXaxis().SetBinLabel(7, "oppositeCharge")
cutflow.SetBinContent(7, w_cut6)

df7 = df6.Define("InvMassVector", "track_inv_mass(SelectedPairsOpCharge,VLMuonPt,VLMuonEta,VLMuonPhi)")
df7 = df7.Define("InvMass", "InvMassVector[0].second")
df7 = df7.Filter("InvMass >= 80 && InvMass <= 100")

w_cut7 = df7.Sum("weight_total").GetValue()
step_names.append("InvMass")
step_cum_eff.append(w_cut7 / sum_ptw if sum_ptw != 0 else 0)
step_step_eff.append(w_cut7 / w_cut6 if w_cut6 != 0 else 0)
cutflow.GetXaxis().SetBinLabel(8, "InvMass")
cutflow.SetBinContent(8, w_cut7)

df7 = df7.Define("muon_indices", "InvMassVector[0].first") \
         .Define("z_p4", "get_z_p4(muon_indices, VLMuonPt, VLMuonEta, VLMuonPhi)") \
         .Define("Z_pt", "(float)z_p4.Pt()") \
         .Define("Z_eta", "(float)z_p4.Eta()") \
         .Define("Z_phi", "(float)z_p4.Phi()") \

# --- histograms ---
h_ptw   = df7.Histo1D(("ptw",   "pT weight distribution;pT weight;Events", 100, 0, 5), "ptw")
h_genZ  = df7.Histo1D(("genZpt_rew", "Generator-level Z p_{T};p_{T} [GeV];Weighted Events", 50, 0, 200), "genZpt", "weight_total")
h_m     = df7.Histo1D(("InvMass", "Invariant Mass;M_{#mu#mu} [GeV];Events", 60, 60, 120), "InvMass", "weight_total")
h_pt    = df7.Histo1D(("Z_pt", "Z boson p_{T};p_{T} [GeV];Events", 100, 0, 200), "Z_pt", "weight_total")
h_eta   = df7.Histo1D(("Z_eta", "Z boson #eta;#eta;Events", 60, -5, 5), "Z_eta", "weight_total")
h_phi   = df7.Histo1D(("Z_phi", "Z boson #phi;#phi;Events", 64, -3.2, 3.2), "Z_phi", "weight_total")

output_path = "test_output"
os.makedirs(output_path, exist_ok=True)
rootfile_path = os.path.join(output_path, "Z_plots.root")
f_out = ROOT.TFile(rootfile_path, "RECREATE")

for h in [total_h, cutflow, h_ptw, h_genZpt, h_m, h_pt, h_eta, h_phi]:
    h.Write()


f_out.Close()
print(f"Wrote output to {rootfile_path}")

print("\n***** Weighted Cutflow (ptw only) *****")
for i, name in enumerate(step_names):
    step_eff = step_step_eff[i]
    cum_eff = step_cum_eff[i]
    print(f"{name:20s} | step-eff = {step_eff:.6f} | cum-eff = {cum_eff:.6f}")

end_time = time.time()
print(f"\nExecution time: {end_time - start_time:.1f} s")
