import ROOT
import os
import time
import sys
start_time = time.time()
# Enable batch mode so graphics are not shown (but saved as pdf)
ROOT.gROOT.SetBatch(True)
# Enable multi-threading
#ROOT.ROOT.EnableImplicitMT(4)

treename = "rootTupleTreeVeryLoose/tree"

file_list_name = sys.argv[1] if len(sys.argv) > 1 else "file_list.txt"
print(f"Using file list: {file_list_name}")

with open(file_list_name, "r") as f:
    my_list = f.read().splitlines()

ch = ROOT.TChain(treename)
for tree in my_list:
    ch.Add(tree)

df = ROOT.RDataFrame(ch)
total = df.Count().GetValue()
total_h = ROOT.TH1D("total_h", "total_h", 1, 0, 1)
total_h.SetBinContent(1, total)

# define custom fuctions
ROOT.gInterpreter.Declare('#include "analyzer.h"')

# Request cut-flow report
report = df.Report()

df = df.Define("sumpt_size","SumPt.size()")
total_vertices = df.Sum("sumpt_size")

# create list of trigger objects
df = df.Define("TriggerPt", "combineTriggers(HLTMu10p5IP3p5Part012345Pt, HLTMu12IP6Part01234Pt, HLTMu7IP4Part01234Pt, HLTMu8IP3Part012345Pt, \
               HLTMu8IP5Part01234Pt, HLTMu8IP6Part01234Pt, HLTMu8p5IP3p5Part012345Pt, HLTMu9IP0Part0Pt, HLTMu9IP3Part0Pt, \
               HLTMu9IP4Part01234Pt, HLTMu9IP5Part01234Pt, HLTMu9IP6Part012345Pt)")
df = df.Define("TriggerEta", "combineTriggers(HLTMu10p5IP3p5Part012345Eta, HLTMu12IP6Part01234Eta, HLTMu7IP4Part01234Eta, HLTMu8IP3Part012345Eta, \
               HLTMu8IP5Part01234Eta, HLTMu8IP6Part01234Eta, HLTMu8p5IP3p5Part012345Eta, HLTMu9IP0Part0Eta, HLTMu9IP3Part0Eta, \
               HLTMu9IP4Part01234Eta, HLTMu9IP5Part01234Eta, HLTMu9IP6Part012345Eta)")
df = df.Define("TriggerPhi", "combineTriggers(HLTMu10p5IP3p5Part012345Phi, HLTMu12IP6Part01234Phi, HLTMu7IP4Part01234Phi, HLTMu8IP3Part012345Phi, \
               HLTMu8IP5Part01234Phi, HLTMu8IP6Part01234Phi, HLTMu8p5IP3p5Part012345Phi, HLTMu9IP0Part0Phi, HLTMu9IP3Part0Phi, \
               HLTMu9IP4Part01234Phi, HLTMu9IP5Part01234Phi, HLTMu9IP6Part012345Phi)")

df = df.Define("TriggerNum", "TriggerPt.size()")
df = df.Filter("TriggerNum > 0", "Filter events with at least one trigger")

df = df.Define("doTrigsMatch", "doTrigsMatch(TriggerPt, TriggerEta, TriggerPhi, VLMuonPt, VLMuonEta, VLMuonPhi, MuonVertInd, SumPt)")
df = df.Filter("doTrigsMatch", "One trigger must match to vertex")


# reject vertices
df = df.Define("rejectedVertices", "create_listOfRejectedVertices(IsGoodRecoVertex, SumPt, VectorSumPt, VLTrackNByVert, \
                                                                  Cluster1Size, Cluster2Size, Cluster1Chi2, Cluster1DoF, Cluster2Chi2, Cluster2DoF, TotalChi2, TotalDoF, \
                                                                  TriggerPt, TriggerEta, TriggerPhi, VLMuonPt, VLMuonEta, VLMuonPhi, MuonVertInd)")
df = df.Define("rejectedVertices_size","rejectedVertices.size()")
df = df.Filter("rejectedVertices_size > 0","rejectedVertices_size > 0")

df = df.Define("numAcceptedVertices","std::count(rejectedVertices.begin(), rejectedVertices.end(), false)")
accepted_vertices = df.Sum("numAcceptedVertices")


# Select Muons from Good Vertices
df = df.Define("isGoodMuon","select_muonsFromGoodVertices(rejectedVertices, MuonVertInd)")

df = df.Define("isGoodMuon_size","std::count(isGoodMuon.begin(), isGoodMuon.end(), true)")

df = df.Filter("isGoodMuon_size > 1", "isGoodMuon_size > 1")

# Selects muons with pt>10, iso<0.1, |eta|<2.5                  
df = df.Define("goodMuonPtIsoEta","select_muons_ptIsoEta(isGoodMuon,VLMuonPt, 10., MuonIso03, 0.1, VLMuonEta, 2.5)")
df = df.Define("goodMuonPtIsoEta_size","std::count(goodMuonPtIsoEta.begin(), goodMuonPtIsoEta.end(), true)")
df = df.Filter("goodMuonPtIsoEta_size > 1","goodMuonPtIsoEta_size > 1")

# note: add id cut VLMuonPOGIDBit passPOGIDCutBasedLoose

# Selects pairs of muons from save PV
df= df.Define("muonpairs","create_pairsOfGoodMuons(goodMuonPtIsoEta)")
df = df.Define("pairsMuonsSamePV","create_pairsMuonsSamePV(muonpairs, MuonVertInd)")
df = df.Define("pairsMuonsSamePV_size","pairsMuonsSamePV.size()")
df = df.Filter("pairsMuonsSamePV_size > 0","pairsMuonsSamePV > 0")

# Selects pairs of muons with one with pt>20, the other with pt>10      
df = df.Define("pairsMuonsPtAsymmetry","create_pairsMuonsPtAsymmetry(pairsMuonsSamePV, VLMuonPt, 20.)")
df = df.Define("pairsMuonsPtAsymmetry_size","pairsMuonsPtAsymmetry.size()")
df = df.Filter("pairsMuonsPtAsymmetry_size > 0","pairsMuonsPtAsymmetry_size > 0")

# Creates a vector called SelectedPairsOpCharge that holds the pairs of indices of tracks that have opposite charge
df = df.Define("SelectedPairsOpCharge","idx_0f_trackpairs_opcharge(pairsMuonsPtAsymmetry, VLMuonCharge)")
df = df.Define("NumOfPairsTracksOpCharge","SelectedPairsOpCharge.size()")
df = df.Define("SelectedPairsSameCharge","idx_0f_trackpairs_samecharge(pairsMuonsPtAsymmetry, VLMuonCharge)")
df = df.Define("NumOfPairsTracksSameCharge","SelectedPairsSameCharge.size()")
df1 = df.Filter("NumOfPairsTracksSameCharge > 0","num of pairs same charge > 0")
df = df.Filter("NumOfPairsTracksOpCharge > 0","num of pairs op charge > 0")


df1 = df1.Define("DimuonMassVec_Bkg", "track_inv_mass(SelectedPairsSpCharge, VLMuonPt, VLMuonEta, VLMuonPhi)")
df1 = df1.Define("DimuonMass_Bkg", "DimuonMassVec_Bkg[0].second")
mass_ss_hist = df1.Histo1D(("Mass_SameSignPairs", "Mass_bkg", 120, 0, 120), "DimuonMass_Bkg")

# Calculates the invariant mass of each pair of tracks that meet the requirements, sorts the values by descending order, stores in a vector called InvMassVector
df = df.Define("InvMassVector","track_inv_mass(SelectedPairsOpCharge, VLMuonPt, VLMuonEta, VLMuonPhi)")
df = df.Define("InvMass","InvMassVector[0].second")
invmass = df.Histo1D(("InvMass","InvMass", 120, 0, 120), "InvMass")
df = df.Filter("InvMass >= 80 && InvMass <= 100","80 <= InvMass <= 100")

h=[invmass, mass_ss_hist] # this is the list to hold our histos, for easier plotting later.

output_path = "test"

# Create the full output path if it doesn't exist
if not os.path.exists(output_path):
    os.makedirs(output_path)

rootfile_path = os.path.join(output_path, "histos.root")
rootfile = ROOT.TFile(rootfile_path, "RECREATE")
total_h.Write()
for histo in h:
    histo.Draw("HIST")
    histo.Write()
# Print report on cutflow 
print('\n*****Cutflow Report*****')
report.Print()

print(f"total_vertices: {total_vertices.GetValue()}")
print(f"accepted_vertices: {accepted_vertices.GetValue()}")

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")
exit()
