# Execute this when you login to hexfarm/lxplus (every login, only once)
# source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc12-opt/setup.sh
# ----
# Some tutorials:
# https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html
# https://root.cern.ch/doc/master/group__tutorial__dataframe.html 
# https://root.cern.ch/doc/master/df102__NanoAODDimuonAnalysis_8py.html (for analysis, histo making)
# https://root.cern.ch/doc/master/df103__NanoAODHiggsAnalysis_8py.html (for analysis, histo making)
# https://root.cern.ch/doc/master/df007__snapshot_8py.html (for making new trees from existing ones)
# ----
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
# df = ROOT.RDataFrame(treename, "/home/jth155/instanton1/CMSSW_10_6_29/src/Multilepton/RootTupleMaker/analysisTree.root") # for a single file, one can also do this.

# define custom fuctions
ROOT.gInterpreter.Declare('#include "analyzer_zmc.h"')

# Request cut-flow report
report = df.Report()


# reject vertices
df = df.Define("rejectedVertices", "create_listOfRejectedVertices(IsGoodRecoVertex, SumPt, VectorSumPt, VLTrackNByVert, \
                                                                  Cluster1Size, Cluster2Size, Cluster1Chi2, Cluster1DoF, Cluster2Chi2, Cluster2DoF, TotalChi2, TotalDoF, \
                                                                  VLMuonPt, VLMuonEta, VLMuonPhi, MuonVertInd)")
df = df.Define("rejectedVertices_size","rejectedVertices.size()")
df = df.Filter("rejectedVertices_size > 0","rejectedVertices_size > 0")



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
df = df.Filter("NumOfPairsTracksOpCharge > 0","num of pairs op charge > 0")




# Calculates the invariant mass of each pair of tracks that meet the requirements, sorts the values by descending order, stores in a vector called InvMassVector
df = df.Define("InvMassVector","track_inv_mass(SelectedPairsOpCharge, VLMuonPt, VLMuonEta, VLMuonPhi)")
df = df.Define("InvMass","InvMassVector[0].second")
df = df.Filter("InvMass >= 80 && InvMass <= 100","80 <=InvMass <= 100")

invmass = df.Histo1D(("InvMass","InvMass", 60, 60, 120), "InvMass")


h=[invmass] # this is the list to hold our histos, for easier plotting later.

#output_dir = "/cms/insta/jth155/results/"
#output_folder = "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/0000"
#output_folder = "test"
output_path = "test"

# Create the full output path if it doesn't exist
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Plot
#ROOT.gStyle.SetOptStat('emruo')
#ROOT.gStyle.SetTextFont(42)
#c = ROOT.TCanvas("c", "", 800, 700)
#c.SetLogx(0)
#c.SetLogy(0)
rootfile_path = os.path.join(output_path, "histos.root")
rootfile = ROOT.TFile(rootfile_path, "RECREATE")
total_h.Write()
for histo in h:
    histo.Draw("HIST")
    # c.SaveAs('{}.pdf'.format(histo.GetTitle()))
   # save_path = os.path.join(output_path, '{}.pdf'.format(histo.GetTitle()))
    #c.SaveAs(save_path)
    histo.Write()
# Print report on cutflow 
print('\n*****Cutflow Report*****')
report.Print()

# print(f"total_vertices: {total_vertices.GetValue()}")
# print(f"accepted_vertices: {accepted_vertices.GetValue()}")

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")
exit()
