#include <ROOT/RVec.hxx>
#include "TMath.h"
#include <TLorentzVector.h>
#include <vector>
#include <cmath>
#include <iostream>
#include "Math/Vector4D.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h" 
using namespace ROOT::VecOps;
using namespace ROOT::Math;

typedef std::vector<float>* vpF;
typedef std::vector<int>* vpI;

template <typename T>

float sum ( ROOT::VecOps::RVec<float> &x ){
  float sum=0;
  for (std::size_t i = 0; i != x.size(); ++i) {
    sum += x.at(i);
  }
  return sum;
}

int sum (const ROOT::VecOps::RVec<int> &x ){
  int sum=0;
  for (std::size_t i = 0; i != x.size(); ++i) {
    sum += x.at(i);
  }
  return sum;
}

float invMass( float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2 ){
  TLorentzVector a, b;
  a.SetPtEtaPhiM( pt1, eta1, phi1, mass1 );
  b.SetPtEtaPhiM( pt2, eta2, phi2, mass2 );
  float mass = (a+b).M();
  //
  return mass;
}


// std::pair<int, int> top_two_idx(const ROOT::VecOps::RVec<int> &vertIdx) {
//     std::vector<std::pair<float, int>> indexedNums;
//     for (int i = 0; i < vertIdx.size(); ++i) {
//         if (vertIdx[i] == 0) {
//             indexedNums.emplace_back(vertIdx[i], i);
//         }
//     }

//     return {indexedNums[0].second, indexedNums[1].second};
// }

std::pair<int, int> top_two_idx(const ROOT::VecOps::RVec<float> &nums) {
    std::vector<std::pair<float, int>> indexedNums;
    for (int i = 0; i < nums.size(); ++i) {
        indexedNums.emplace_back(nums[i], i);
    }

    std::sort(indexedNums.begin(), indexedNums.end(), std::greater<std::pair<float, int>>());

    return {indexedNums[0].second, indexedNums[1].second};
}

// vertex level filters

bool select_goodvertices(int i,
                         const ROOT::VecOps::RVec<int> &isgoodvertex) {

  if (isgoodvertex[i] == 0) return true;
  return false;
}

bool sumpt_filter(int i,
                  const ROOT::VecOps::RVec<float> &sumpt, 
                  const ROOT::VecOps::RVec<float> &vectorsumpt) {
  if (sumpt[i] > 20. && vectorsumpt[i]/sumpt[i] >= 0.4) return true;
  return false;
}

bool numTrackPerVertex_filter(int i,
                              const ROOT::VecOps::RVec<int> &trackNByVert) {
  if ((int) trackNByVert[i] < 3) return true;
  return false;
}

bool reVertexingFilter(int i, 
                       const ROOT::VecOps::RVec<int> &c1size, const ROOT::VecOps::RVec<int> &c2size, 
                       const ROOT::VecOps::RVec<float> &c1chi2, const ROOT::VecOps::RVec<int> &c1dof, 
                       const ROOT::VecOps::RVec<float> &c2chi2, const ROOT::VecOps::RVec<int> &c2dof, 
                       const ROOT::VecOps::RVec<float> &tchi2, const ROOT::VecOps::RVec<int> &tdof) {
  if(c1size[i]>=3 && c2size[i]>=3){
      if(tchi2[i]>0 && c1chi2[i]>0 && c2chi2[i]>0 && tdof[i]>0 && c1dof[i]>0 && c2dof[i]>0){
          float oneVert = tchi2[i] / tdof[i];
          float twoVert = (c1chi2[i] + c2chi2[i]) / (c1dof[i] + c2dof[i]);
          if(oneVert - twoVert >= 10*sqrt(2)) return true;
      }
  }
  return false;
}

void muonTriggerObjectReject(std::vector<bool>& out,
                            const ROOT::VecOps::RVec<float>& trigPt, 
                            const ROOT::VecOps::RVec<float>& trigEta, 
                            const ROOT::VecOps::RVec<float>& trigPhi, 
                            const ROOT::VecOps::RVec<float>& muonPt, 
                            const ROOT::VecOps::RVec<float>& muonEta, 
                            const ROOT::VecOps::RVec<float>& muonPhi, 
                            const ROOT::VecOps::RVec<int>& muonVert, 
                            const ROOT::VecOps::RVec<float>& sumPt){
  TLorentzVector muon = TLorentzVector();
  TLorentzVector trig = TLorentzVector();

  for(unsigned int i = 0; i < muonPt.size(); i++){
    if(muonVert[i] < 0){
        continue;
    }
    muon.SetPtEtaPhiM(1, muonEta[i], muonPhi[i], 1);

    for(unsigned int j = 0; j < trigPt.size(); j++){
      trig.SetPtEtaPhiM(1, trigEta[j], trigPhi[j], 1);

      // if(muon.DrEtaPhi(trig) < 0.01 && fabs(1 - muonPt[i]/trigPt[j]) <= 0.1){
      if(muon.DrEtaPhi(trig) < 0.01){
          out[muonVert[i]] = true;
          break;
      }
    }
  }

}

bool doTrigsMatch(const ROOT::VecOps::RVec<float>& trigPt, 
                  const ROOT::VecOps::RVec<float>& trigEta, 
                  const ROOT::VecOps::RVec<float>& trigPhi, 
                  const ROOT::VecOps::RVec<float>& muonPt, 
                  const ROOT::VecOps::RVec<float>& muonEta, 
                  const ROOT::VecOps::RVec<float>& muonPhi, 
                  const ROOT::VecOps::RVec<int>& muonVert,
                  const ROOT::VecOps::RVec<float> &sumpt) {
  std::vector<bool> rejectedVertices(sumpt.size());
  muonTriggerObjectReject(rejectedVertices, trigPt, trigEta, trigPhi, muonPt, muonEta, muonPhi, muonVert, sumpt);
  for (int i = 0; i < rejectedVertices.size(); i++) {
    if (rejectedVertices[i]) return true;
  }
  return false;
}

std::vector<bool> create_listOfRejectedVertices(const ROOT::VecOps::RVec<int> &isgoodvertex,
                                                const ROOT::VecOps::RVec<float> &sumpt, 
                                                const ROOT::VecOps::RVec<float> &vectorsumpt,
                                                const ROOT::VecOps::RVec<int> &trackNByVert,
                                                const ROOT::VecOps::RVec<int> &c1size, const ROOT::VecOps::RVec<int> &c2size, 
                                                const ROOT::VecOps::RVec<float> &c1chi2, const ROOT::VecOps::RVec<int> &c1dof, 
                                                const ROOT::VecOps::RVec<float> &c2chi2, const ROOT::VecOps::RVec<int> &c2dof, 
                                                const ROOT::VecOps::RVec<float> &tchi2, const ROOT::VecOps::RVec<int> &tdof,
                                                const ROOT::VecOps::RVec<float>& muonPt, 
                                                const ROOT::VecOps::RVec<float>& muonEta, 
                                                const ROOT::VecOps::RVec<float>& muonPhi, 
                                                const ROOT::VecOps::RVec<int>& muonVert) {
  std::vector<bool> rejectedVertices(sumpt.size());
   for (int i = 0; i < rejectedVertices.size(); i++) {
     if (select_goodvertices(i, isgoodvertex)) rejectedVertices[i] = true;
     else if (sumpt_filter(i, sumpt, vectorsumpt)) rejectedVertices[i] = true;
     else if (numTrackPerVertex_filter(i, trackNByVert)) rejectedVertices[i] = true;
     else if (reVertexingFilter(i, c1size, c2size, c1chi2, c1dof, c2chi2, c2dof, tchi2, tdof)) rejectedVertices[i] = true;
   }
  

  return rejectedVertices;

}

std::vector<int> testfunction(const ROOT::VecOps::RVec<int> &isgoodvertex) {
    std::vector<int> rejectedVertices(isgoodvertex.size());
    return rejectedVertices;
}

template <typename... Args>
ROOT::VecOps::RVec<float> combineTriggers(const Args&... triggers) {
    ROOT::VecOps::RVec<float> combined;
    
    // Compute total size for preallocation
    size_t total_size = (triggers.size() + ... + 0);
    combined.reserve(total_size);
    
    // Insert each trigger vector into the combined vector
    (combined.insert(combined.end(), triggers.begin(), triggers.end()), ...);

    return combined;
}


std::vector<bool> select_muonsFromGoodVertices(const std::vector<bool> &rejectedVertices, 
                                              const ROOT::VecOps::RVec<int> &muonvertind) {
  std::vector<bool> muons(muonvertind.size());
  for (int i = 0; i < muonvertind.size(); i++) {
    int pv = muonvertind[i];
    if (rejectedVertices[pv] == false) muons[i] = true;
  }
  return muons;
}

// muon level filters
std::vector<bool> select_muons_ptIsoEta(std::vector<bool> &isgoodmuon,
                                  const ROOT::VecOps::RVec<float> &pt, float ptMin, 
                                  const ROOT::VecOps::RVec<float> &iso, float isoMax, 
                                  const ROOT::VecOps::RVec<float> &eta, float etaMax) {
  for (int i = 0; i < isgoodmuon.size(); ++i) {
    if (!(isgoodmuon[i])) continue;
    bool ptCut = pt[i] > ptMin;
    bool isoCut = iso[i] < isoMax;
    // bool isoCut = iso[i] < (isoMax + 1.0)*pt[i];
    bool etaCut = fabs(eta[i]) < etaMax;    
    if (!ptCut || !isoCut || !etaCut) isgoodmuon[i] = false;
  }

  return isgoodmuon;
} 

std::vector<std::pair<int,int>> create_pairsOfGoodMuons(const std::vector<bool> &isgoodmuon) {
  std::vector<std::pair<int,int>> muonpairs;
  for (int i = 0; i < isgoodmuon.size()-1; i++) {
    if (!(isgoodmuon[i])) continue;
    for (int j = i+1; j < isgoodmuon.size(); j++) {
      if (isgoodmuon[j]) muonpairs.push_back({i,j});
    }
  }
  return muonpairs;
}

std::vector<std::pair<int,int>> create_pairsMuonsSamePV(const std::vector<std::pair<int,int>> &muonpairs, const ROOT::VecOps::RVec<int> &muonvertind) {
  std::vector<std::pair<int,int>> newMuonPairs;
  for (int i = 0; i < muonpairs.size(); ++i) {
    int idx1 = muonpairs[i].first;
    int idx2 = muonpairs[i].second;
    if (muonvertind[idx1] == muonvertind[idx2]) newMuonPairs.push_back(muonpairs[i]);
  }
  return newMuonPairs;
}

std::vector<std::pair<int,int>> create_pairsMuonsPtAsymmetry(const std::vector<std::pair<int,int>> &muonpairs, 
                                                              const ROOT::VecOps::RVec<float> &pt, float ptMin) {
  std::vector<std::pair<int,int>> newMuonPairs;
  for (int i = 0; i < muonpairs.size(); ++i) {
    int idx1 = muonpairs[i].first;
    if (pt[idx1] > ptMin) {
      newMuonPairs.push_back(muonpairs[i]);
    }
  }
  return newMuonPairs;
}

std::vector<std::pair<int,int>> idx_0f_trackpairs_opcharge(const std::vector<std::pair<int,int>> &trackPairs, const ROOT::VecOps::RVec<int> &charge) {
  std::vector<pair<int,int>> vec;
  for (int i = 0; i < trackPairs.size(); ++i) {
    int idx1 = trackPairs[i].first;
    int idx2 = trackPairs[i].second;
    if (charge[idx1] != charge[idx2]) {
      vec.push_back(trackPairs[i]);
    }
  }
  return vec;
}

std::vector<std::pair<int,int>> idx_0f_trackpairs_samecharge(const std::vector<std::pair<int,int>> &trackPairs, const ROOT::VecOps::RVec<int> &charge) {
  std::vector<pair<int,int>> vec;
  for (int i = 0; i < trackPairs.size(); ++i) {
    int idx1 = trackPairs[i].first;
    int idx2 = trackPairs[i].second;
    if (charge[idx1] == charge[idx2]) {
      vec.push_back(trackPairs[i]);
    }
  }
  return vec;
}

std::vector<std::pair<std::pair<int, int>, float>> track_inv_mass(const std::vector<std::pair<int,int>> &trackPairs, 
                                  const ROOT::VecOps::RVec<float> &pt, 
                                  const ROOT::VecOps::RVec<float> &eta, 
                                  const ROOT::VecOps::RVec<float> &phi) {
  std::vector<std::pair<std::pair<int, int>, float>> vec;
  for (int i = 0; i < trackPairs.size(); ++i) {
    int idx1 = trackPairs[i].first;
    int idx2 = trackPairs[i].second;
    float invmass = invMass(pt[idx1], eta[idx1], phi[idx1], 0.1, pt[idx2], eta[idx2], phi[idx2], 0.1);
    vec.push_back({trackPairs[i], invmass});
  }

  // sort by distance closest to 90 GeV
  std::sort(vec.begin(), vec.end(), [](const std::pair<std::pair<int, int>, float> &a, const std::pair<std::pair<int, int>, float> &b) {
    return fabs(a.second - 90.0) < fabs(b.second - 90.0);
  });

  return vec;
}

float calculate_pt_Z(float pt1, float phi1, float pt2, float phi2) {
  // Calculate the x and y components of the Z boson's transverse momentum
  float pxZ = pt1 * cos(phi1) + pt2 * cos(phi2);
  float pyZ = pt1 * sin(phi1) + pt2 * sin(phi2);

  // Calculate the magnitude of the transverse momentum
  float ptZ = sqrt(pxZ * pxZ + pyZ * pyZ);

  return ptZ;
}

float calculate_pt_Z_TLorentz(float pt1, float eta1, float phi1, float mass1, 
                              float pt2, float eta2, float phi2, float mass2) {
    // Create TLorentzVector for the first particle
    TLorentzVector vec1;
    vec1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);

    // Create TLorentzVector for the second particle
    TLorentzVector vec2;
    vec2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);

    // Sum the two TLorentzVectors to get the Z boson TLorentzVector
    TLorentzVector vecZ = vec1 + vec2;

    // Return the transverse momentum of the Z boson
    return vecZ.Pt();
}

float sum_pt_of_vertex_with_Z(const float &ZVertIdx, const ROOT::VecOps::RVec<int> &pvIdx, const ROOT::VecOps::RVec<float> &pt) {
  float sumpt = 0;
  for (int i = 0; i < pt.size(); ++i) {
    if (pvIdx[i] == ZVertIdx) {
      sumpt += pt[i];
    }
  }

  return sumpt;
}


std::vector<int> muonMatching(const ROOT::VecOps::RVec<int> &teta, 
                              const ROOT::VecOps::RVec<int> &tphi, 
                              const ROOT::VecOps::RVec<int> &meta, 
                              const ROOT::VecOps::RVec<int> &mphi,
                              const ROOT::VecOps::RVec<int> &mvert){
    std::vector<int> matchedMuonVerts;
    std::vector<int> matchedMuonInd;
    if(teta.size() == 0){
        matchedMuonVerts.emplace_back(-1);
        // return matchedMuonVerts;
        matchedMuonInd.emplace_back(-1);
        return matchedMuonInd;
    }
    //loop through trig objs
    for(unsigned int i = 0; i < teta.size(); i++){
        TLorentzVector trig = TLorentzVector(); 
        trig.SetPtEtaPhiM(1, teta[i], tphi[i], 1);
        float minDr = FLT_MAX;
        matchedMuonVerts.emplace_back(-1);
        matchedMuonInd.emplace_back(-1);
        //loop through muons
        for(unsigned int j = 0; j < mvert.size(); j++){
            TLorentzVector muon = TLorentzVector(); 
            muon.SetPtEtaPhiM(1, meta.at(j), mphi.at(j), 1);
            if(muon.DrEtaPhi(trig) < minDr){
                minDr = muon.DrEtaPhi(trig);
                //record vertex of matched muon 
                matchedMuonVerts.at(i) = mvert.at(j);
                // record index of matched muon
                matchedMuonInd.at(i) = j;
            }
        }
    }
    // return matchedMuonVerts;
    return matchedMuonInd;
}

int muonMatching2(const ROOT::VecOps::RVec<float>& trigPt, 
                                const ROOT::VecOps::RVec<float>& trigEta, 
                                const ROOT::VecOps::RVec<float>& trigPhi, 
                                const ROOT::VecOps::RVec<float>& muonPt, 
                                const ROOT::VecOps::RVec<float>& muonEta, 
                                const ROOT::VecOps::RVec<float>& muonPhi, 
                                const ROOT::VecOps::RVec<int>& muonVert){
  TLorentzVector muon = TLorentzVector();
  TLorentzVector trig = TLorentzVector();
  int matchedIdx;

  for(unsigned int i = 0; i < muonPt.size(); i++){
    if(muonVert[i] < 0){
        continue;
    }
    muon.SetPtEtaPhiM(1, muonEta[i], muonPhi[i], 1);

    for(unsigned int j = 0; j < trigPt.size(); j++){
      trig.SetPtEtaPhiM(1, trigEta[j], trigPhi[j], 1);

      // if(muon.DrEtaPhi(trig) < 0.01 && fabs(1 - muonPt[i]/trigPt[j]) <= 0.1){
      if(muon.DrEtaPhi(trig) < 0.01){
          matchedIdx = i;
          break;
      }
    }
  }

  return matchedIdx;
}


int numOfTrigMuons(std::vector<int> &trigMuInd, int mu1Idx, int mu2Idx) {
  int numOfTrigMu = 0;
  for (int i = 0; i < trigMuInd.size(); ++i) {
    if (trigMuInd[i] == mu1Idx || trigMuInd[i] == mu2Idx) {
      numOfTrigMu++;
    }
  }
  return numOfTrigMu;
}



//  Calculates the 4-vector of the Z boson from two muon indices.
PtEtaPhiMVector get_z_p4(
    const std::pair<unsigned int, unsigned int>& indices,
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& phi
) {
    // Muon mass in GeV/c^2
    const float muon_mass = 0.1056583745; 
    PtEtaPhiMVector p1(pt[indices.first], eta[indices.first], phi[indices.first], muon_mass);
    PtEtaPhiMVector p2(pt[indices.second], eta[indices.second], phi[indices.second], muon_mass);
    return p1 + p2;
}

// Global static variables for histogram access
static TFile* __pt_file__ = nullptr;
static TH1D* __pt_hist__ = nullptr;

// Initializes the global Z pT weighting histogram.
void init_pt_hist(const char* filename, const char* histname) {
    if(__pt_file__) return;
    __pt_file__ = TFile::Open(filename);
    if(!__pt_file__) {
        std::cerr << "init_pt_hist: ERROR - cannot open file: " << filename << std::endl;
        return;
    }
    __pt_hist__ = dynamic_cast<TH1D*>(__pt_file__->Get(histname));
    if(!__pt_hist__) {
        std::cerr << "init_pt_hist: ERROR - cannot find histogram '" << histname
                  << "' in file: " << filename << std::endl;
    } else {
        std::cout << "init_pt_hist: loaded histogram '" << histname
                  << "' (" << __pt_hist__->GetNbinsX() << " bins, x-range "
                  << __pt_hist__->GetXaxis()->GetXmin() << " - "
                  << __pt_hist__->GetXaxis()->GetXmax() << ")." << std::endl;
    }
}


 //Retrieves the weight for a given Z boson pT from the initialized histogram.
double getPtWeight(double zpt) {
    if(!__pt_hist__) return 1.0;
    if(zpt < 0.0) return 1.0;
    
    int bin = __pt_hist__->FindBin(zpt);
    int nbins = __pt_hist__->GetNbinsX();
    if(bin < 1) bin = 1;
    if(bin > nbins) bin = nbins;
    
    double w = __pt_hist__->GetBinContent(bin);
    if(std::isnan(w) || std::isinf(w)) return 1.0;
    
    return w;
}


 // Finds the maximum pT of generator-level prompt particles (ROOT RVec overload).
double getGenZpt(const ROOT::VecOps::RVec<float>& v) {
    if(v.empty()) return -1.0;
    double m = v[0];
    for(size_t i = 1; i < v.size(); ++i) {
        if(v[i] > m) m = v[i];
    }
    return m;
}

