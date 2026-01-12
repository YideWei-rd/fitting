import ROOT
import os
import sys
import time

start_time = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(16)

treename = "rootTupleTreeVeryLoose/tree"
file_list_name = sys.argv[1] if len(sys.argv) > 1 else "file_list.txt"
print(f"Using file list: {file_list_name}")

with open(file_list_name, "r") as f:
    file_list = f.read().splitlines()

ch = ROOT.TChain(treename)
for tree in file_list:
    ch.Add(tree)

df = ROOT.RDataFrame(ch)

ROOT.gInterpreter.Declare('#include "analyzer_zmc.h"')

sums_ptrs = {} 

# Filter 0: Gen Level
df = df.Filter("ROOT::VecOps::Sum(GenPromptParticlePt > 0) == 1", "Keep only gen-level Z->mumu events")
df = df.Define("weight", "GenWeight")

sums_ptrs["total"] = df.Sum("weight")

# --- CUT 1: reject vertices ---
df1 = df.Define("rejectedVertices",
                "create_listOfRejectedVertices(IsGoodRecoVertex, SumPt, VectorSumPt, VLTrackNByVert, \
                 Cluster1Size, Cluster2Size, Cluster1Chi2, Cluster1DoF, Cluster2Chi2, Cluster2DoF, TotalChi2, TotalDoF, \
                 VLMuonPt, VLMuonEta, VLMuonPhi, MuonVertInd)") \
        .Define("numAcceptedVertices","std::count(rejectedVertices.begin(), rejectedVertices.end(), false)") \
        .Filter("numAcceptedVertices > 0")

sums_ptrs["rejectedVertices"] = df1.Sum("weight")

df1 = df1.Define("AcceptedVertexSumPt", "get_accepted_vertices_sumpt(rejectedVertices, SumPt)")

h_accepted_sumpt_ptr = df1.Histo1D(
    ("AcceptedVertexSumPt", "Sum p_{T} of Accepted Vertices;Sum p_{T} [GeV];Weighted Entries", 100, 0, 500), 
    "AcceptedVertexSumPt", 
    "weight"
)

# --- CUT 2: good muons ---
df2 = df1.Define("isGoodMuon","select_muonsFromGoodVertices(rejectedVertices, MuonVertInd)") \
         .Define("isGoodMuon_size","std::count(isGoodMuon.begin(), isGoodMuon.end(), true)") \
         .Filter("isGoodMuon_size > 1")

sums_ptrs["goodMuon"] = df2.Sum("weight")

# --- CUT 3: pt, iso, eta ---
df3 = df2.Define("goodMuonPtIsoEta","select_muons_ptIsoEta(isGoodMuon,VLMuonPt,10.,MuonIso03,0.1,VLMuonEta,2.5)") \
         .Define("goodMuonPtIsoEta_size","std::count(goodMuonPtIsoEta.begin(), goodMuonPtIsoEta.end(), true)") \
         .Filter("goodMuonPtIsoEta_size > 1")

sums_ptrs["ptIsoEta"] = df3.Sum("weight")

# --- CUT 4: muon pairs ---
df4 = df3.Define("muonpairs","create_pairsOfGoodMuons(goodMuonPtIsoEta)") \
         .Define("pairsMuonsSamePV","create_pairsMuonsSamePV(muonpairs, MuonVertInd)") \
         .Define("pairsMuonsSamePV_size","pairsMuonsSamePV.size()") \
         .Filter("pairsMuonsSamePV_size > 0")

sums_ptrs["pairsSamePV"] = df4.Sum("weight")

# --- CUT 5: pt asymmetry ---
df5 = df4.Define("pairsMuonsPtAsymmetry","create_pairsMuonsPtAsymmetry(pairsMuonsSamePV,VLMuonPt,20.)") \
         .Define("pairsMuonsPtAsymmetry_size","pairsMuonsPtAsymmetry.size()") \
         .Filter("pairsMuonsPtAsymmetry_size > 0")

sums_ptrs["ptAsymmetry"] = df5.Sum("weight")

# --- CUT 6: opposite charge ---
df6 = df5.Define("SelectedPairsOpCharge","idx_0f_trackpairs_opcharge(pairsMuonsPtAsymmetry, VLMuonCharge)") \
         .Define("NumOfPairsTracksOpCharge","SelectedPairsOpCharge.size()") \
         .Filter("NumOfPairsTracksOpCharge > 0")

sums_ptrs["oppositeCharge"] = df6.Sum("weight")

# --- CUT 7: InvMass window ---
df7 = df6.Define("InvMassVector","track_inv_mass(SelectedPairsOpCharge,VLMuonPt,VLMuonEta,VLMuonPhi)") \
         .Define("InvMass","InvMassVector[0].second") \
         .Filter("InvMass >= 80 && InvMass <= 100")

sums_ptrs["InvMass"] = df7.Sum("weight")

invmass_h_ptr = df7.Histo1D(("InvMass","InvMass",60,60,120),"InvMass","weight")

print("Starting event loop...")

invmass_h = invmass_h_ptr.GetValue() 

print("Event loop finished.")

step_names = ["total", "rejectedVertices", "goodMuon", "ptIsoEta", 
              "pairsSamePV", "ptAsymmetry", "oppositeCharge", "InvMass"]

cutflow = ROOT.TH1D("cutflow", "Weighted cutflow", 10, 0, 10)

prev_w = 0
sumw_total = sums_ptrs["total"].GetValue() 
print(f"Total Sum of GenWeights: {sumw_total}")

print("\n***** Weighted Cutflow *****")
for i, name in enumerate(step_names):
    current_w = sums_ptrs[name].GetValue()
    
    cutflow.GetXaxis().SetBinLabel(i+1, name)
    cutflow.SetBinContent(i+1, current_w)
    
    if i == 0:
        step_eff = 1.0
        cum_eff = 1.0
        prev_w = current_w
    else:
        cum_eff = current_w / sumw_total if sumw_total > 0 else 0
        step_eff = current_w / prev_w if prev_w > 0 else 0
        prev_w = current_w
        
    print(f"{name:20s} | step-eff = {step_eff:.4f} | cum-eff = {cum_eff:.4f}")

# --- Output ---
output_path = "output_dir"
os.makedirs(output_path, exist_ok=True)
rootfile_path = os.path.join(output_path, "histos_nlo_GenWeight.root")
file_out = ROOT.TFile(rootfile_path, "RECREATE")

cutflow.Write()
invmass_h.Write()
h_accepted_sumpt_ptr.Write()
file_out.Close()

end_time = time.time()
print(f"\nExecution time: {end_time - start_time:.1f} s")
