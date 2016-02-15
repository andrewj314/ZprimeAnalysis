/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects and pairs 
   It prints a table following the cutFlow for rateEvt, relEff, and cumEff 

Need to specify
1. See Declare constants
*/
/////
//   To run: root -l Efficiency_cutflow_table_rateEvt_relEff_cumEff.cc+ 
/////
/////
//   Prepare Root and Roofit
/////
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
/////
//   Declare constants
/////
//Path - rootpla; Luminosity (for normalization of MC to data)
const string path     = "/uscms_data/d3/andrewj/CMSSW_7_4_7/src/Efficiencies/";
const string rootpla  = "OutTree.root";
const double Luminosity = 19600; //pb^-1
const bool noLumiNorm = true; //true means NO luminosity normalization done
const bool noPUcorr   = true; //true means NO PU corr done
const double cand1_acc[2]       = {20, 2.5};
const double cand2_acc[2]       = {20, 2.5};
const int cand1_numTotCut = 11; //It is different (-1) from what you could see in ZpmSelection.cc, but still correct
const int cand2_numTotCut = 12; //It is different (-1) from what you could see in ZpmSelection.cc, but still correct
const int pair_numTotCut  = 6;
const int numTotCut       = 4; //dr+cand1_numTotCut+cand2_numTotCut+pair_numTotCut 
//charge==-1 && cosDphi<-0.95 && met>=30 && pZetaMt-3.051*pZetaVisMt>-50 && cosDphiLMet<0.2 && nbjet==0
double cut_sig[pair_numTotCut] = {-1,-0.95,20,-50,0.2,0};
/////
//   Declare functions 
/////
TFile* Call_TFile();
void cutflow_eff(TTree* tree, int nentries);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
void ini_print(int nentries);
void measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight);
/////
//   Main function
/////
void test() {
 //Call tree 
 TFile* f = Call_TFile();
 cout << "File called\n";
 TTree* tree; f->GetObject("TNT/BOOM;1",tree);
 cout << "Tree initialized\n";
 int nentries = tree->GetEntries(); 
 cout << "nentries = " << nentries << "\n";
 cutflow_eff(tree, nentries);
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile() {
 string file_name = path+rootpla;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Efficiency for the object
/////
void cutflow_eff(TTree* tree, int nentries){ 
 //Define strings to print
 vector<string> varname(numTotCut); 
 vector<string> varcut(numTotCut);
 cout << "Initializing strings:\n";
 initialize_strings(varname,varcut);
 //Call variables 
 int num1Cut, num2Cut, njet, nbjet;
 double wgt_lumi, wgt_pu;

 vector<double> *Muon_isoSum , *Muon_isoCharParPt ;
 vector<double> *Muon_pt , *Muon_eta, *Muon_phi, *Muon_dz_pv, *Muon_energy, *Muon_iso;
 vector<double> *Muon_isoCharged, *Muon_isoNeutralHadron , *Muon_isoPhoton, *Muon_isoPU;
 vector<double> *Muon_charge, *Muon_chi2, *Muon_p, *Muon_matchedStat, *Muon_dxy_pv, *Muon_dxy_bs; 
 vector<double> *Muon_dz_bs, *Muon_dzError, *Muon_dxyError, *Muon_ndof, *Muon_vtx, *Muon_vty, *Muon_vtz; 
 vector<double> *Muon_track_pt, *Muon_track_ptError, *Muon_validHits, *Muon_validHitsInner, *Muon_TLayers; 
 vector<bool>   *Muon_loose, *Muon_medium, *Muon_tight, *Muon_soft, *Muon_isHighPt, *Muon_pf, *Muon_isGlobal;   
 vector<double> *Muon_dB, *Muon_besttrack_pt, *Muon_besttrack_ptError, *Muon_tunePBestTrack_pt;
 vector<double> *Muon_isTrackerMuon, *Muon_POGisGood;
 vector<double> *Muon_chi2LocalPosition, *Muon_trkKink, *Muon_segmentCompatibility, *Muon_validFraction, *Muon_pixelLayersWithMeasurement, *Muon_qualityhighPurity;
 vector<double> *Muon_tunePBestTrackType;
 vector<double> *Muon_track_PCAx_bs, *Muon_track_PCAy_bs, *Muon_track_PCAz_bs;
 vector<double> *Muon_track_PCAx_pv, *Muon_track_PCAy_pv, *Muon_track_PCAz_pv;
 vector<double> *Muon_trackFitErrorMatrix_00, *Muon_trackFitErrorMatrix_01, *Muon_trackFitErrorMatrix_02, *Muon_trackFitErrorMatrix_11;
 vector<double> *Muon_trackFitErrorMatrix_12, *Muon_trackFitErrorMatrix_22;

 vector<double> *patElectron_pt , *patElectron_eta, *patElectron_phi, *patElectron_energy, *patElectron_charge;
 vector<double> *patElectron_gsfTrack_ndof, *patElectron_gsfTrack_dxy_pv, *patElectron_dxyError, *patElectron_gsfTrack_normChi2; 
 vector<double> *patElectron_gsfTrack_dz_pv, *patElectron_gsfTrack_dz_bs; 
 vector<double> *patElectron_gsfTrack_vtx, *patElectron_gsfTrack_vty, *patElectron_gsfTrack_vtz;
 vector<double> *patElectron_gsfTrack_dxy_bs, *isoChargedHadrons_, *isoNeutralHadrons_, *isoPhotons_, *isoPU_;
 vector<double> *patElectron_gsfTrack_PCAx_bs, *patElectron_gsfTrack_PCAy_bs, *patElectron_gsfTrack_PCAz_bs;
 vector<double> *patElectron_gsfTrack_PCAx_pv, *patElectron_gsfTrack_PCAy_pv, *patElectron_gsfTrack_PCAz_pv;
 vector<int>    *passVetoId_, *passLooseId_, *passMediumId_, *passTightId_, *passHEEPId_, *passConversionVeto_, *expectedMissingInnerHits;
 vector<double> *patElectron_gsfTrackFitErrorMatrix_00, *patElectron_gsfTrackFitErrorMatrix_01, *patElectron_gsfTrackFitErrorMatrix_02;
 vector<double> *patElectron_gsfTrackFitErrorMatrix_11, *patElectron_gsfTrackFitErrorMatrix_12, *patElectron_gsfTrackFitErrorMatrix_22;

 vector<double> *Tau_eta, *Tau_phi, *Tau_pt, *Tau_energy, *Tau_charge, *Tau_chargedIsoPtSum, *Tau_neutralIsoPtSum, *Tau_puCorrPtSum ;
 vector<double> *Tau_leadChargedCandPt, *Tau_leadChargedCandCharge, *Tau_leadChargedCandEta, *Tau_leadChargedCandPhi, *Tau_nProngs;
 vector<double> *Tau_leadChargedCandChi2, *Tau_leadChargedCandValidHits, *Tau_leadChargedCandDxy_pv, *Tau_leadChargedCandDxy_bs, *Tau_leadChargedCandDz_bs;
 vector<double> *Tau_leadChargedCandDz_pv, *Tau_leadChargedCandDzError, *Tau_leadChargedCandDxyError, *Tau_leadChargedCandNdof, *Tau_leadChargedCandVtx;
 vector<double> *Tau_leadChargedCandVty, *Tau_leadChargedCandVtz, *Tau_leadChargedCandTrack_pt, *Tau_leadChargedCandTrack_ptError;
 vector<double> *Tau_leadChargedCandTrack_PCAx_bs, *Tau_leadChargedCandTrack_PCAy_bs, *Tau_leadChargedCandTrack_PCAz_bs;
 vector<double> *Tau_leadChargedCandTrack_PCAx_pv, *Tau_leadChargedCandTrack_PCAy_pv, *Tau_leadChargedCandTrack_PCAz_pv;
 vector<double> *Tau_leadChargedCandTrackFitErrorMatrix_00, *Tau_leadChargedCandTrackFitErrorMatrix_01, *Tau_leadChargedCandTrackFitErrorMatrix_02;
 vector<double> *Tau_leadChargedCandTrackFitErrorMatrix_11, *Tau_leadChargedCandTrackFitErrorMatrix_12, *Tau_leadChargedCandTrackFitErrorMatrix_22;
 vector<double> *Tau_defaultDxy, *Tau_defaultDxyError, *Tau_defaultDxySig, *Tau_defaultFlightLengthSig, *Tau_default_PCAx_pv, *Tau_default_PCAy_pv, *Tau_default_PCAz_pv;
 vector<double> *Tau_defaultFlightLengthX, *Tau_defaultFlightLengthY, *Tau_defaultFlightLengthZ;
 vector <int>   *Tau_decayModeFinding, *Tau_decayModeFindingOldDMs, *Tau_decayModeFindingNewDMs;
 vector <int>   *Tau_byLooseCombinedIsolationDeltaBetaCorr, *Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
 vector <int>   *Tau_byMediumCombinedIsolationDeltaBetaCorr, *Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
 vector <int>   *Tau_byTightCombinedIsolationDeltaBetaCorr, *Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
 vector <int>   *Tau_byLooseIsolationMVA3newDMwLT, *Tau_byLooseIsolationMVA3newDMwoLT, *Tau_byMediumIsolationMVA3newDMwLT;
 vector <int>   *Tau_byMediumIsolationMVA3newDMwoLT, *Tau_byLooseIsolationMva3oldDMwLT, *Tau_byLooseIsolationMVA3oldDMwoLT;
 vector <int>   *Tau_byMediumIsolationMva3oldDMwLT, *Tau_byMediumIsolationMVA3oldDMwoLT, *Tau_byTightIsolationMVA3newDMwLT;
 vector <int>   *Tau_byTightIsolationMVA3newDMwoLT, *Tau_byTightIsolationMva3oldDMwLT, *Tau_byTightIsolationMVA3oldDMwoLT;
 vector <int>   *Tau_againstMuonLoose2, *Tau_againstMuonLoose3, *Tau_againstMuonTight2, *Tau_againstMuonTight3;
 vector <int>   *Tau_againstElectronMVALooseMVA5, *Tau_againstElectronMVAMediumMVA5, *Tau_againstElectronMVATightMVA5;
 vector <int>   *Tau_byVLooseCombinedIsolationDeltaBetaCorr;

 Muon_isoSum = 0;
 Muon_isoCharParPt = 0;
 Muon_pt = 0;
 Muon_eta = 0;
 Muon_phi = 0;
 Muon_dz_pv = 0;
 Muon_energy = 0;
 Muon_iso = 0;
 Muon_isoCharged = 0; 
 Muon_isoNeutralHadron = 0; 
 Muon_isoPhoton = 0; 
 Muon_isoPU = 0;
 Muon_charge = 0; 
 Muon_chi2 = 0; 
 Muon_p = 0; 
 Muon_matchedStat = 0; 
 Muon_dxy_pv = 0; 
 Muon_dxy_bs = 0;
 Muon_dz_bs = 0; 
 Muon_dzError = 0; 
 Muon_dxyError = 0; 
 Muon_ndof = 0; 
 Muon_vtx = 0; 
 Muon_vty = 0; 
 Muon_vtz = 0;
 Muon_track_pt = 0; 
 Muon_track_ptError = 0; 
 Muon_validHits = 0; 
 Muon_validHitsInner = 0; 
 Muon_TLayers = 0;
 Muon_loose = 0; 
 Muon_medium = 0;
 Muon_tight = 0; 
 Muon_soft = 0; 
 Muon_isHighPt = 0; 
 Muon_pf = 0; 
 Muon_isGlobal = 0;
 Muon_dB = 0;
 Muon_besttrack_pt = 0; 
 Muon_besttrack_ptError = 0; 
 Muon_tunePBestTrack_pt = 0;
 Muon_isTrackerMuon = 0; 
 Muon_POGisGood = 0;
 Muon_chi2LocalPosition = 0; 
 Muon_trkKink = 0; 
 Muon_segmentCompatibility = 0; 
 Muon_validFraction = 0; 
 Muon_pixelLayersWithMeasurement = 0; 
 Muon_qualityhighPurity = 0;
 Muon_tunePBestTrackType = 0;
 Muon_track_PCAx_bs = 0; 
 Muon_track_PCAy_bs = 0; 
 Muon_track_PCAz_bs = 0;
 Muon_track_PCAx_pv = 0;
 Muon_track_PCAy_pv = 0; 
 Muon_track_PCAz_pv = 0;
 Muon_trackFitErrorMatrix_00 = 0; 
 Muon_trackFitErrorMatrix_01 = 0; 
 Muon_trackFitErrorMatrix_02 = 0; 
 Muon_trackFitErrorMatrix_11 = 0;
 Muon_trackFitErrorMatrix_12 = 0;
 Muon_trackFitErrorMatrix_22 = 0;

 patElectron_pt = 0;
 patElectron_eta = 0;
 patElectron_phi = 0;
 patElectron_energy = 0;
 patElectron_charge = 0;
 patElectron_gsfTrack_ndof = 0; 
 patElectron_gsfTrack_dxy_pv = 0;
 patElectron_dxyError = 0;
 patElectron_gsfTrack_normChi2 = 0;
 patElectron_gsfTrack_dz_pv = 0; 
 patElectron_gsfTrack_dz_bs = 0;
 patElectron_gsfTrack_vtx = 0; 
 patElectron_gsfTrack_vty = 0; 
 patElectron_gsfTrack_vtz = 0;
 patElectron_gsfTrack_dxy_bs = 0;
 isoChargedHadrons_ = 0; 
 isoNeutralHadrons_ = 0; 
 isoPhotons_ = 0; 
 isoPU_ = 0;
 patElectron_gsfTrack_PCAx_bs = 0; 
 patElectron_gsfTrack_PCAy_bs = 0; 
 patElectron_gsfTrack_PCAz_bs = 0;
 patElectron_gsfTrack_PCAx_pv = 0; 
 patElectron_gsfTrack_PCAy_pv = 0; 
 patElectron_gsfTrack_PCAz_pv = 0;
 passVetoId_ = 0; 
 passLooseId_ = 0; 
 passMediumId_ = 0; 
 passTightId_ = 0; 
 passHEEPId_ = 0; 
 passConversionVeto_ = 0; 
 expectedMissingInnerHits = 0;
 patElectron_gsfTrackFitErrorMatrix_00 = 0; 
 patElectron_gsfTrackFitErrorMatrix_01 = 0; 
 patElectron_gsfTrackFitErrorMatrix_02 = 0;
 patElectron_gsfTrackFitErrorMatrix_11 = 0; 
 patElectron_gsfTrackFitErrorMatrix_12 = 0; 
 patElectron_gsfTrackFitErrorMatrix_22 = 0;

 Tau_eta = 0; 
 Tau_phi = 0; 
 Tau_pt = 0; 
 Tau_energy = 0; 
 Tau_charge = 0;
 Tau_chargedIsoPtSum = 0; 
 Tau_neutralIsoPtSum = 0; 
 Tau_puCorrPtSum = 0;
 Tau_leadChargedCandPt = 0;
 Tau_leadChargedCandCharge = 0; 
 Tau_leadChargedCandEta = 0; 
 Tau_leadChargedCandPhi = 0; 
 Tau_nProngs = 0;
 Tau_leadChargedCandChi2 = 0; 
 Tau_leadChargedCandValidHits = 0; 
 Tau_leadChargedCandDxy_pv = 0; 
 Tau_leadChargedCandDxy_bs = 0; 
 Tau_leadChargedCandDz_bs = 0;
 Tau_leadChargedCandDz_pv = 0; 
 Tau_leadChargedCandDzError = 0; 
 Tau_leadChargedCandDxyError = 0; 
 Tau_leadChargedCandNdof = 0;
 Tau_leadChargedCandVtx = 0;
 Tau_leadChargedCandVty = 0;
 Tau_leadChargedCandVtz = 0; 
 Tau_leadChargedCandTrack_pt = 0;
 Tau_leadChargedCandTrack_ptError = 0;
 Tau_leadChargedCandTrack_PCAx_bs = 0;
 Tau_leadChargedCandTrack_PCAy_bs = 0; 
 Tau_leadChargedCandTrack_PCAz_bs = 0;
 Tau_leadChargedCandTrack_PCAx_pv = 0; 
 Tau_leadChargedCandTrack_PCAy_pv = 0; 
 Tau_leadChargedCandTrack_PCAz_pv = 0;
 Tau_leadChargedCandTrackFitErrorMatrix_00 = 0; 
 Tau_leadChargedCandTrackFitErrorMatrix_01 = 0; 
 Tau_leadChargedCandTrackFitErrorMatrix_02 = 0;
 Tau_leadChargedCandTrackFitErrorMatrix_11 = 0; 
 Tau_leadChargedCandTrackFitErrorMatrix_12 = 0;
 Tau_leadChargedCandTrackFitErrorMatrix_22 = 0;
 Tau_defaultDxy = 0; 
 Tau_defaultDxyError = 0; 
 Tau_defaultDxySig = 0; 
 Tau_defaultFlightLengthSig = 0;
 Tau_default_PCAx_pv = 0; 
 Tau_default_PCAy_pv = 0; 
 Tau_default_PCAz_pv = 0;
 Tau_defaultFlightLengthX = 0; 
 Tau_defaultFlightLengthY = 0; 
 Tau_defaultFlightLengthZ = 0;
 Tau_decayModeFinding = 0; 
 Tau_decayModeFindingOldDMs = 0; 
 Tau_decayModeFindingNewDMs = 0;
 Tau_byLooseCombinedIsolationDeltaBetaCorr = 0; 
 Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
 Tau_byMediumCombinedIsolationDeltaBetaCorr = 0; 
 Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
 Tau_byTightCombinedIsolationDeltaBetaCorr = 0; 
 Tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
 Tau_byLooseIsolationMVA3newDMwLT = 0; 
 Tau_byLooseIsolationMVA3newDMwoLT = 0; 
 Tau_byMediumIsolationMVA3newDMwLT = 0;
 Tau_byMediumIsolationMVA3newDMwoLT = 0; 
 Tau_byLooseIsolationMva3oldDMwLT = 0; 
 Tau_byLooseIsolationMVA3oldDMwoLT = 0;
 Tau_byMediumIsolationMva3oldDMwLT = 0; 
 Tau_byMediumIsolationMVA3oldDMwoLT = 0; 
 Tau_byTightIsolationMVA3newDMwLT = 0;
 Tau_byTightIsolationMVA3newDMwoLT = 0;
 Tau_byTightIsolationMva3oldDMwLT = 0; 
 Tau_byTightIsolationMVA3oldDMwoLT = 0;
 Tau_againstMuonLoose2 = 0; 
 Tau_againstMuonLoose3 = 0; 
 Tau_againstMuonTight2 = 0; 
 Tau_againstMuonTight3 = 0;
 Tau_againstElectronMVALooseMVA5 = 0; 
 Tau_againstElectronMVAMediumMVA5 = 0; 
 Tau_againstElectronMVATightMVA5 = 0;
 Tau_byVLooseCombinedIsolationDeltaBetaCorr = 0;



 TBranch *b_Muon_pt = 0;
 TBranch *b_Muon_eta = 0;
 TBranch *b_Muon_isoSum = 0;
 TBranch *b_Muon_isoCharParPt = 0;
 TBranch *b_Muon_phi = 0;
 TBranch *b_Muon_dz_pv = 0;
 TBranch *b_Muon_energy = 0;
 TBranch *b_Muon_iso = 0;
 TBranch *b_Muon_isoCharged = 0; 
 TBranch *b_Muon_isoNeutralHadron = 0; 
 TBranch *b_Muon_isoPhoton = 0; 
 TBranch *b_Muon_isoPU = 0;
 TBranch *b_Muon_charge = 0; 
 TBranch *b_Muon_chi2 = 0; 
 TBranch *b_Muon_p = 0; 
 TBranch *b_Muon_matchedStat = 0; 
 TBranch *b_Muon_dxy_pv = 0; 
 TBranch *b_Muon_dxy_bs = 0;
 TBranch *b_Muon_dz_bs = 0; 
 TBranch *b_Muon_dzError = 0; 
 TBranch *b_Muon_dxyError = 0; 
 TBranch *b_Muon_ndof = 0; 
 TBranch *b_Muon_vtx = 0; 
 TBranch *b_Muon_vty = 0; 
 TBranch *b_Muon_vtz = 0;
 TBranch *b_Muon_track_pt = 0; 
 TBranch *b_Muon_track_ptError = 0; 
 TBranch *b_Muon_validHits = 0; 
 TBranch *b_Muon_validHitsInner = 0; 
 TBranch *b_Muon_TLayers = 0;
 TBranch *b_Muon_loose = 0; 
 TBranch *b_Muon_medium = 0;
 TBranch *b_Muon_tight = 0; 
 TBranch *b_Muon_soft = 0; 
 TBranch *b_Muon_isHighPt = 0; 
 TBranch *b_Muon_pf = 0; 
 TBranch *b_Muon_isGlobal = 0;
 TBranch *b_Muon_dB = 0;
 TBranch *b_Muon_besttrack_pt = 0; 
 TBranch *b_Muon_besttrack_ptError = 0; 
 TBranch *b_Muon_tunePBestTrack_pt = 0;
 TBranch *b_Muon_isTrackerMuon = 0; 
 TBranch *b_Muon_POGisGood = 0;
 TBranch *b_Muon_chi2LocalPosition = 0; 
 TBranch *b_Muon_trkKink = 0; 
 TBranch *b_Muon_segmentCompatibility = 0; 
 TBranch *b_Muon_validFraction = 0; 
 TBranch *b_Muon_pixelLayersWithMeasurement = 0; 
 TBranch *b_Muon_qualityhighPurity = 0;
 TBranch *b_Muon_tunePBestTrackType = 0;
 TBranch *b_Muon_track_PCAx_bs = 0; 
 TBranch *b_Muon_track_PCAy_bs = 0; 
 TBranch *b_Muon_track_PCAz_bs = 0;
 TBranch *b_Muon_track_PCAx_pv = 0;
 TBranch *b_Muon_track_PCAy_pv = 0; 
 TBranch *b_Muon_track_PCAz_pv = 0;
 TBranch *b_Muon_trackFitErrorMatrix_00 = 0; 
 TBranch *b_Muon_trackFitErrorMatrix_01 = 0; 
 TBranch *b_Muon_trackFitErrorMatrix_02 = 0; 
 TBranch *b_Muon_trackFitErrorMatrix_11 = 0;
 TBranch *b_Muon_trackFitErrorMatrix_12 = 0;
 TBranch *b_Muon_trackFitErrorMatrix_22 = 0;

 TBranch *b_patElectron_pt = 0;
 TBranch *b_patElectron_eta = 0;
 TBranch *b_patElectron_phi = 0;
 TBranch *b_patElectron_energy = 0;
 TBranch *b_patElectron_charge = 0;
 TBranch *b_patElectron_gsfTrack_ndof = 0;
 TBranch *b_patElectron_gsfTrack_dxy_pv = 0;
 TBranch *b_patElectron_dxyError = 0;
 TBranch *b_patElectron_gsfTrack_normChi2 = 0;
 TBranch *b_patElectron_gsfTrack_dz_pv = 0;
 TBranch *b_patElectron_gsfTrack_dz_bs = 0;
 TBranch *b_patElectron_gsfTrack_vtx = 0;
 TBranch *b_patElectron_gsfTrack_vty = 0;
 TBranch *b_patElectron_gsfTrack_vtz = 0;
 TBranch *b_patElectron_gsfTrack_dxy_bs = 0;
 TBranch *b_isoChargedHadrons_ = 0;
 TBranch *b_isoNeutralHadrons_ = 0;
 TBranch *b_isoPhotons_ = 0;
 TBranch *b_isoPU_ = 0;
 TBranch *b_patElectron_gsfTrack_PCAx_bs = 0;
 TBranch *b_patElectron_gsfTrack_PCAy_bs = 0;
 TBranch *b_patElectron_gsfTrack_PCAz_bs = 0;
 TBranch *b_patElectron_gsfTrack_PCAx_pv = 0;
 TBranch *b_patElectron_gsfTrack_PCAy_pv = 0;
 TBranch *b_patElectron_gsfTrack_PCAz_pv = 0;
 TBranch *b_passVetoId_ = 0;
 TBranch *b_passLooseId_ = 0;
 TBranch *b_passMediumId_ = 0;
 TBranch *b_passTightId_ = 0;
 TBranch *b_passHEEPId_ = 0;
 TBranch *b_passConversionVeto_ = 0;
 TBranch *b_expectedMissingInnerHits = 0;
 TBranch *b_patElectron_gsfTrackFitErrorMatrix_00 = 0;
 TBranch *b_patElectron_gsfTrackFitErrorMatrix_01 = 0;
 TBranch *b_patElectron_gsfTrackFitErrorMatrix_02 = 0;
 TBranch *b_patElectron_gsfTrackFitErrorMatrix_11 = 0;
 TBranch *b_patElectron_gsfTrackFitErrorMatrix_12 = 0;
 TBranch *b_patElectron_gsfTrackFitErrorMatrix_22 = 0;

 TBranch *b_Tau_eta = 0; 
 TBranch *b_Tau_phi = 0; 
 TBranch *b_Tau_pt = 0; 
 TBranch *b_Tau_energy = 0; 
 TBranch *b_Tau_charge = 0;
 TBranch *b_Tau_chargedIsoPtSum = 0; 
 TBranch *b_Tau_neutralIsoPtSum = 0; 
 TBranch *b_Tau_puCorrPtSum = 0;
 TBranch *b_Tau_leadChargedCandPt = 0;
 TBranch *b_Tau_leadChargedCandCharge = 0; 
 TBranch *b_Tau_leadChargedCandEta = 0; 
 TBranch *b_Tau_leadChargedCandPhi = 0; 
 TBranch *b_Tau_nProngs = 0;
 TBranch *b_Tau_leadChargedCandChi2 = 0; 
 TBranch *b_Tau_leadChargedCandValidHits = 0; 
 TBranch *b_Tau_leadChargedCandDxy_pv = 0; 
 TBranch *b_Tau_leadChargedCandDxy_bs = 0; 
 TBranch *b_Tau_leadChargedCandDz_bs = 0;
 TBranch *b_Tau_leadChargedCandDz_pv = 0; 
 TBranch *b_Tau_leadChargedCandDzError = 0; 
 TBranch *b_Tau_leadChargedCandDxyError = 0; 
 TBranch *b_Tau_leadChargedCandNdof = 0;
 TBranch *b_Tau_leadChargedCandVtx = 0;
 TBranch *b_Tau_leadChargedCandVty = 0;
 TBranch *b_Tau_leadChargedCandVtz = 0; 
 TBranch *b_Tau_leadChargedCandTrack_pt = 0;
 TBranch *b_Tau_leadChargedCandTrack_ptError = 0;
 TBranch *b_Tau_leadChargedCandTrack_PCAx_bs = 0;
 TBranch *b_Tau_leadChargedCandTrack_PCAy_bs = 0; 
 TBranch *b_Tau_leadChargedCandTrack_PCAz_bs = 0;
 TBranch *b_Tau_leadChargedCandTrack_PCAx_pv = 0; 
 TBranch *b_Tau_leadChargedCandTrack_PCAy_pv = 0; 
 TBranch *b_Tau_leadChargedCandTrack_PCAz_pv = 0;
 TBranch *b_Tau_leadChargedCandTrackFitErrorMatrix_00 = 0; 
 TBranch *b_Tau_leadChargedCandTrackFitErrorMatrix_01 = 0; 
 TBranch *b_Tau_leadChargedCandTrackFitErrorMatrix_02 = 0;
 TBranch *b_Tau_leadChargedCandTrackFitErrorMatrix_11 = 0; 
 TBranch *b_Tau_leadChargedCandTrackFitErrorMatrix_12 = 0;
 TBranch *b_Tau_leadChargedCandTrackFitErrorMatrix_22 = 0;
 TBranch *b_Tau_defaultDxy = 0; 
 TBranch *b_Tau_defaultDxyError = 0; 
 TBranch *b_Tau_defaultDxySig = 0; 
 TBranch *b_Tau_defaultFlightLengthSig = 0;
 TBranch *b_Tau_default_PCAx_pv = 0; 
 TBranch *b_Tau_default_PCAy_pv = 0; 
 TBranch *b_Tau_default_PCAz_pv = 0;
 TBranch *b_Tau_defaultFlightLengthX = 0; 
 TBranch *b_Tau_defaultFlightLengthY = 0; 
 TBranch *b_Tau_defaultFlightLengthZ = 0;
 TBranch *b_Tau_decayModeFinding = 0; 
 TBranch *b_Tau_decayModeFindingOldDMs = 0; 
 TBranch *b_Tau_decayModeFindingNewDMs = 0;
 TBranch *b_Tau_byLooseCombinedIsolationDeltaBetaCorr = 0; 
 TBranch *b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
 TBranch *b_Tau_byMediumCombinedIsolationDeltaBetaCorr = 0; 
 TBranch *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
 TBranch *b_Tau_byTightCombinedIsolationDeltaBetaCorr = 0; 
 TBranch *b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
 TBranch *b_Tau_byLooseIsolationMVA3newDMwLT = 0; 
 TBranch *b_Tau_byLooseIsolationMVA3newDMwoLT = 0; 
 TBranch *b_Tau_byMediumIsolationMVA3newDMwLT = 0;
 TBranch *b_Tau_byMediumIsolationMVA3newDMwoLT = 0; 
 TBranch *b_Tau_byLooseIsolationMva3oldDMwLT = 0; 
 TBranch *b_Tau_byLooseIsolationMVA3oldDMwoLT = 0;
 TBranch *b_Tau_byMediumIsolationMva3oldDMwLT = 0; 
 TBranch *b_Tau_byMediumIsolationMVA3oldDMwoLT = 0; 
 TBranch *b_Tau_byTightIsolationMVA3newDMwLT = 0;
 TBranch *b_Tau_byTightIsolationMVA3newDMwoLT = 0;
 TBranch *b_Tau_byTightIsolationMva3oldDMwLT = 0; 
 TBranch *b_Tau_byTightIsolationMVA3oldDMwoLT = 0;
 TBranch *b_Tau_againstMuonLoose2 = 0; 
 TBranch *b_Tau_againstMuonLoose3 = 0; 
 TBranch *b_Tau_againstMuonTight2 = 0; 
 TBranch *b_Tau_againstMuonTight3 = 0;
 TBranch *b_Tau_againstElectronMVALooseMVA5 = 0; 
 TBranch *b_Tau_againstElectronMVAMediumMVA5 = 0; 
 TBranch *b_Tau_againstElectronMVATightMVA5 = 0;
 TBranch *b_Tau_byVLooseCombinedIsolationDeltaBetaCorr = 0;



 tree->SetBranchAddress("Muon_pt",&Muon_pt,&b_Muon_pt);
 tree->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
/*
 tree->SetBranchAddress("Muon_isoSum",&Muon_isoSum,&b_Muon_isoSum);
 tree->SetBranchAddress("Muon_isoCharParPt",&Muon_isoCharParPt,&b_Muon_isoCharParPt);
 tree->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
 tree->SetBranchAddress("Muon_dz_pv",&Muon_dz_pv,&b_Muon_dz_pv);
 tree->SetBranchAddress("Muon_energy",&Muon_energy,&b_Muon_energy);
 tree->SetBranchAddress("Muon_iso",&Muon_iso,&b_Muon_iso);
 tree->SetBranchAddress("Muon_isoCharged",&Muon_isoCharged,&b_Muon_isoCharged); 
 tree->SetBranchAddress("Muon_isoNeutralHadron",&Muon_isoNeutralHadron,&b_Muon_isoNeutralHadron); 
 tree->SetBranchAddress("Muon_isoPhoton",&Muon_isoPhoton,&b_Muon_isoPhoton); 
 tree->SetBranchAddress("Muon_isoPU",&Muon_isoPU,&b_Muon_isoPU);
 tree->SetBranchAddress("Muon_charge",&Muon_charge,&b_Muon_charge); 
 tree->SetBranchAddress("Muon_chi2",&Muon_chi2,&b_Muon_chi2); 
 tree->SetBranchAddress("Muon_p",&Muon_p,&b_Muon_p); 
 tree->SetBranchAddress("Muon_matchedStat",&Muon_matchedStat,&b_Muon_matchedStat); 
 tree->SetBranchAddress("Muon_dxy_pv",&Muon_dxy_pv,&b_Muon_dxy_pv); 
 tree->SetBranchAddress("Muon_dxy_bs",&Muon_dxy_bs,&b_Muon_dxy_bs);
 tree->SetBranchAddress("Muon_dz_bs",&Muon_dz_bs,&b_Muon_dz_bs); 
 tree->SetBranchAddress("Muon_dzError",&Muon_dzError,&b_Muon_dzError); 
 tree->SetBranchAddress("Muon_dxyError",&Muon_dxyError,&b_Muon_dxyError); 
 tree->SetBranchAddress("Muon_ndof",&Muon_ndof,&b_Muon_ndof); 
 tree->SetBranchAddress("Muon_vtx",&Muon_vtx,&b_Muon_vtx); 
 tree->SetBranchAddress("Muon_vty",&Muon_vty,&b_Muon_vty); 
 tree->SetBranchAddress("Muon_vtz",&Muon_vtz,&b_Muon_vtz);
 tree->SetBranchAddress("Muon_track_pt",&Muon_track_pt,&b_Muon_track_pt); 
 tree->SetBranchAddress("Muon_track_ptError",&Muon_track_ptError,&b_Muon_track_ptError); 
 tree->SetBranchAddress("Muon_validHits",&Muon_validHits,&b_Muon_validHits); 
 tree->SetBranchAddress("Muon_validHitsInner",&Muon_validHitsInner,&b_Muon_validHitsInner); 
 tree->SetBranchAddress("Muon_TLayers",&Muon_TLayers,&b_Muon_TLayers);
 tree->SetBranchAddress("Muon_loose",&Muon_loose,&b_Muon_loose); 
 tree->SetBranchAddress("Muon_medium",&Muon_medium,&b_Muon_medium);
 tree->SetBranchAddress("Muon_tight",&Muon_tight,&b_Muon_tight); 
 tree->SetBranchAddress("Muon_soft",&Muon_soft,&b_Muon_soft); 
 tree->SetBranchAddress("Muon_isHighPt",&Muon_isHighPt,&b_Muon_isHighPt); 
 tree->SetBranchAddress("Muon_pf",&Muon_pf,&b_Muon_pf); 
 tree->SetBranchAddress("Muon_isGlobal",&Muon_isGlobal,&b_Muon_isGlobal);
 tree->SetBranchAddress("Muon_dB",&Muon_dB,&b_Muon_dB);
 tree->SetBranchAddress("Muon_besttrack_pt",&Muon_besttrack_pt,&b_Muon_besttrack_pt); 
 tree->SetBranchAddress("Muon_besttrack_ptError",&Muon_besttrack_ptError,&b_Muon_besttrack_ptError); 
 tree->SetBranchAddress("Muon_tunePBestTrack_pt",&Muon_tunePBestTrack_pt,&b_Muon_tunePBestTrack_pt);
 tree->SetBranchAddress("Muon_isTrackerMuon",&Muon_isTrackerMuon,&b_Muon_isTrackerMuon); 
 tree->SetBranchAddress("Muon_POGisGood",&Muon_POGisGood,&b_Muon_POGisGood);
 tree->SetBranchAddress("Muon_chi2LocalPosition",&Muon_chi2LocalPosition,&b_Muon_chi2LocalPosition); 
 tree->SetBranchAddress("Muon_trkKink",&Muon_trkKink,&b_Muon_trkKink); 
 tree->SetBranchAddress("Muon_segmentCompatibility",&Muon_segmentCompatibility,&b_Muon_segmentCompatibility); 
 tree->SetBranchAddress("Muon_validFraction",&Muon_validFraction,&b_Muon_validFraction); 
 tree->SetBranchAddress("Muon_pixelLayersWithMeasurement",&Muon_pixelLayersWithMeasurement,&b_Muon_pixelLayersWithMeasurement); 
 tree->SetBranchAddress("Muon_qualityhighPurity",&Muon_qualityhighPurity,&b_Muon_qualityhighPurity);
 tree->SetBranchAddress("Muon_tunePBestTrackType",&Muon_tunePBestTrackType,&b_Muon_tunePBestTrackType);
 tree->SetBranchAddress("Muon_track_PCAx_bs",&Muon_track_PCAx_bs,&b_Muon_track_PCAx_bs); 
 tree->SetBranchAddress("Muon_track_PCAy_bs",&Muon_track_PCAy_bs,&b_Muon_track_PCAy_bs); 
 tree->SetBranchAddress("Muon_track_PCAz_bs",&Muon_track_PCAz_bs,&b_Muon_track_PCAz_bs);
 tree->SetBranchAddress("Muon_track_PCAx_pv",&Muon_track_PCAx_pv,&b_Muon_track_PCAx_pv);
 tree->SetBranchAddress("Muon_track_PCAy_pv",&Muon_track_PCAy_pv,&b_Muon_track_PCAy_pv); 
 tree->SetBranchAddress("Muon_track_PCAz_pv",&Muon_track_PCAz_pv,&b_Muon_track_PCAz_pv);
 tree->SetBranchAddress("Muon_trackFitErrorMatrix_00",&Muon_trackFitErrorMatrix_00,&b_Muon_trackFitErrorMatrix_00); 
 tree->SetBranchAddress("Muon_trackFitErrorMatrix_01",&Muon_trackFitErrorMatrix_01,&b_Muon_trackFitErrorMatrix_01); 
 tree->SetBranchAddress("Muon_trackFitErrorMatrix_02",&Muon_trackFitErrorMatrix_02,&b_Muon_trackFitErrorMatrix_02); 
 tree->SetBranchAddress("Muon_trackFitErrorMatrix_11",&Muon_trackFitErrorMatrix_11,&b_Muon_trackFitErrorMatrix_11);
 tree->SetBranchAddress("Muon_trackFitErrorMatrix_12",&Muon_trackFitErrorMatrix_12,&b_Muon_trackFitErrorMatrix_12);
 tree->SetBranchAddress("Muon_trackFitErrorMatrix_22",&Muon_trackFitErrorMatrix_22,&b_Muon_trackFitErrorMatrix_22);

*/



 tree->SetBranchAddress("patElectron_pt",&patElectron_pt,&b_patElectron_pt);
 tree->SetBranchAddress("patElectron_eta",&patElectron_eta,&b_patElectron_eta);
/*
 tree->SetBranchAddress("patElectron_phi",&patElectron_phi,&b_patElectron_phi);
 tree->SetBranchAddress("patElectron_energy",&patElectron_energy,&b_patElectron_energy);
 tree->SetBranchAddress("patElectron_charge",&patElectron_charge,&b_patElectron_charge);
 tree->SetBranchAddress("patElectron_gsfTrack_ndof",&patElectron_gsfTrack_ndof,&b_patElectron_gsfTrack_ndof);
 tree->SetBranchAddress("patElectron_gsfTrack_dxy_pv",&patElectron_gsfTrack_dxy_pv,&b_patElectron_gsfTrack_dxy_pv);
 tree->SetBranchAddress("patElectron_dxyError",&patElectron_dxyError,&b_patElectron_dxyError);
 tree->SetBranchAddress("patElectron_gsfTrack_normChi2",&patElectron_gsfTrack_normChi2,&b_patElectron_gsfTrack_normChi2);
 tree->SetBranchAddress("patElectron_gsfTrack_dz_pv",&patElectron_gsfTrack_dz_pv,&b_patElectron_gsfTrack_dz_pv);
 tree->SetBranchAddress("patElectron_gsfTrack_dz_bs",&patElectron_gsfTrack_dz_bs,&b_patElectron_gsfTrack_dz_bs);
 tree->SetBranchAddress("patElectron_gsfTrack_vtx",&patElectron_gsfTrack_vtx,&b_patElectron_gsfTrack_vtx);
 tree->SetBranchAddress("patElectron_gsfTrack_vty",&patElectron_gsfTrack_vty,&b_patElectron_gsfTrack_vty);
 tree->SetBranchAddress("patElectron_gsfTrack_vtz",&patElectron_gsfTrack_vtz,&b_patElectron_gsfTrack_vtz);
 tree->SetBranchAddress("patElectron_gsfTrack_dxy_bs",&patElectron_gsfTrack_dxy_bs,&b_patElectron_gsfTrack_dxy_bs);
 tree->SetBranchAddress("isoChargedHadrons_",&isoChargedHadrons_,&b_isoChargedHadrons_);
 tree->SetBranchAddress("isoNeutralHadrons_",&isoNeutralHadrons_,&b_isoNeutralHadrons_);
 tree->SetBranchAddress("isoPhotons_",&isoPhotons_,&b_isoPhotons_);
 tree->SetBranchAddress("isoPU_",&isoPU_,&b_isoPU_);
 tree->SetBranchAddress("patElectron_gsfTrack_PCAx_bs",&patElectron_gsfTrack_PCAx_bs,&b_patElectron_gsfTrack_PCAx_bs);
 tree->SetBranchAddress("patElectron_gsfTrack_PCAy_bs",&patElectron_gsfTrack_PCAy_bs,&b_patElectron_gsfTrack_PCAy_bs);
 tree->SetBranchAddress("patElectron_gsfTrack_PCAz_bs",&patElectron_gsfTrack_PCAz_bs,&b_patElectron_gsfTrack_PCAz_bs);
 tree->SetBranchAddress("patElectron_gsfTrack_PCAx_pv",&patElectron_gsfTrack_PCAx_pv,&b_patElectron_gsfTrack_PCAx_pv);
 tree->SetBranchAddress("patElectron_gsfTrack_PCAy_pv",&patElectron_gsfTrack_PCAy_pv,&b_patElectron_gsfTrack_PCAy_pv);
 tree->SetBranchAddress("patElectron_gsfTrack_PCAz_pv",&patElectron_gsfTrack_PCAz_pv,&b_patElectron_gsfTrack_PCAz_pv);
 tree->SetBranchAddress("passVetoId_",&passVetoId_,&b_passVetoId_);
 tree->SetBranchAddress("passLooseId_",&passLooseId_,&b_passLooseId_);
 tree->SetBranchAddress("passMediumId_",&passMediumId_,&b_passMediumId_);
 tree->SetBranchAddress("passTightId_",&passTightId_,&b_passTightId_);
 tree->SetBranchAddress("passHEEPId_",&passHEEPId_,&b_passHEEPId_);
 tree->SetBranchAddress("passConversionVeto_",&passConversionVeto_,&b_passConversionVeto_);
 tree->SetBranchAddress("expectedMissingInnerHits",&expectedMissingInnerHits,&b_expectedMissingInnerHits);
 tree->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_00",&patElectron_gsfTrackFitErrorMatrix_00,&b_patElectron_gsfTrackFitErrorMatrix_00);
 tree->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_01",&patElectron_gsfTrackFitErrorMatrix_01,&b_patElectron_gsfTrackFitErrorMatrix_01);
 tree->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_02",&patElectron_gsfTrackFitErrorMatrix_02,&b_patElectron_gsfTrackFitErrorMatrix_02);
 tree->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_11",&patElectron_gsfTrackFitErrorMatrix_11,&b_patElectron_gsfTrackFitErrorMatrix_11);
 tree->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_12",&patElectron_gsfTrackFitErrorMatrix_12,&b_patElectron_gsfTrackFitErrorMatrix_12);
 tree->SetBranchAddress("patElectron_gsfTrackFitErrorMatrix_22",&patElectron_gsfTrackFitErrorMatrix_22,&b_patElectron_gsfTrackFitErrorMatrix_22); 
*/


 tree->SetBranchAddress("Tau_eta",&Tau_eta,&b_Tau_eta); 
 tree->SetBranchAddress("Tau_phi",&Tau_phi,&b_Tau_phi); 
 tree->SetBranchAddress("Tau_pt",&Tau_pt,&b_Tau_pt); 
/*
 tree->SetBranchAddress("Tau_energy",&Tau_energy,&b_Tau_energy); 
 tree->SetBranchAddress("Tau_charge",&Tau_charge,&b_Tau_charge);
 tree->SetBranchAddress("Tau_chargedIsoPtSum",&Tau_chargedIsoPtSum,&b_Tau_chargedIsoPtSum); 
 tree->SetBranchAddress("Tau_neutralIsoPtSum",&Tau_neutralIsoPtSum,&b_Tau_neutralIsoPtSum); 
 tree->SetBranchAddress("Tau_puCorrPtSum",&Tau_puCorrPtSum,&b_Tau_puCorrPtSum);
 tree->SetBranchAddress("Tau_leadChargedCandPt",&Tau_leadChargedCandPt,&b_Tau_leadChargedCandPt);
 tree->SetBranchAddress("Tau_leadChargedCandCharge",&Tau_leadChargedCandCharge,&b_Tau_leadChargedCandCharge); 
 tree->SetBranchAddress("Tau_leadChargedCandEta",&Tau_leadChargedCandEta,&b_Tau_leadChargedCandEta); 
 tree->SetBranchAddress("Tau_leadChargedCandPhi",&Tau_leadChargedCandPhi,&b_Tau_leadChargedCandPhi); 
 tree->SetBranchAddress("Tau_nProngs",&Tau_nProngs,&b_Tau_nProngs);
 tree->SetBranchAddress("Tau_leadChargedCandChi2",&Tau_leadChargedCandChi2,&b_Tau_leadChargedCandChi2); 
 tree->SetBranchAddress("Tau_leadChargedCandValidHits",&Tau_leadChargedCandValidHits,&b_Tau_leadChargedCandValidHits); 
 tree->SetBranchAddress("Tau_leadChargedCandDxy_pv",&Tau_leadChargedCandDxy_pv,&b_Tau_leadChargedCandDxy_pv); 
 tree->SetBranchAddress("Tau_leadChargedCandDxy_bs",&Tau_leadChargedCandDxy_bs,&b_Tau_leadChargedCandDxy_bs); 
 tree->SetBranchAddress("Tau_leadChargedCandDz_bs",&Tau_leadChargedCandDz_bs,&b_Tau_leadChargedCandDz_bs);
 tree->SetBranchAddress("Tau_leadChargedCandDz_pv",&Tau_leadChargedCandDz_pv,&b_Tau_leadChargedCandDz_pv); 
 tree->SetBranchAddress("Tau_leadChargedCandDzError",&Tau_leadChargedCandDzError,&b_Tau_leadChargedCandDzError); 
 tree->SetBranchAddress("Tau_leadChargedCandDxyError",&Tau_leadChargedCandDxyError,&b_Tau_leadChargedCandDxyError); 
 tree->SetBranchAddress("Tau_leadChargedCandNdof",&Tau_leadChargedCandNdof,&b_Tau_leadChargedCandNdof);
 tree->SetBranchAddress("Tau_leadChargedCandVtx",&Tau_leadChargedCandVtx,&b_Tau_leadChargedCandVtx);
 tree->SetBranchAddress("Tau_leadChargedCandVty",&Tau_leadChargedCandVty,&b_Tau_leadChargedCandVty);
 tree->SetBranchAddress("Tau_leadChargedCandVtz",&Tau_leadChargedCandVtz,&b_Tau_leadChargedCandVtz); 
 tree->SetBranchAddress("Tau_leadChargedCandTrack_pt",&Tau_leadChargedCandTrack_pt,&b_Tau_leadChargedCandTrack_pt);
 tree->SetBranchAddress("Tau_leadChargedCandTrack_ptError",&Tau_leadChargedCandTrack_ptError,&b_Tau_leadChargedCandTrack_ptError);
 tree->SetBranchAddress("Tau_leadChargedCandTrack_PCAx_bs",&Tau_leadChargedCandTrack_PCAx_bs,&b_Tau_leadChargedCandTrack_PCAx_bs);
 tree->SetBranchAddress("Tau_leadChargedCandTrack_PCAy_bs",&Tau_leadChargedCandTrack_PCAy_bs,&b_Tau_leadChargedCandTrack_PCAy_bs); 
 tree->SetBranchAddress("Tau_leadChargedCandTrack_PCAz_bs",&Tau_leadChargedCandTrack_PCAz_bs,&b_Tau_leadChargedCandTrack_PCAz_bs);
 tree->SetBranchAddress("Tau_leadChargedCandTrack_PCAx_pv",&Tau_leadChargedCandTrack_PCAx_pv,&b_Tau_leadChargedCandTrack_PCAx_pv); 
 tree->SetBranchAddress("Tau_leadChargedCandTrack_PCAy_pv",&Tau_leadChargedCandTrack_PCAy_pv,&b_Tau_leadChargedCandTrack_PCAy_pv); 
 tree->SetBranchAddress("Tau_leadChargedCandTrack_PCAz_pv",&Tau_leadChargedCandTrack_PCAz_pv,&b_Tau_leadChargedCandTrack_PCAz_pv);
 tree->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_00",&Tau_leadChargedCandTrackFitErrorMatrix_00,&b_Tau_leadChargedCandTrackFitErrorMatrix_00); 
 tree->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_01",&Tau_leadChargedCandTrackFitErrorMatrix_01,&b_Tau_leadChargedCandTrackFitErrorMatrix_01); 
 tree->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_02",&Tau_leadChargedCandTrackFitErrorMatrix_02,&b_Tau_leadChargedCandTrackFitErrorMatrix_02);
 tree->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_11",&Tau_leadChargedCandTrackFitErrorMatrix_11,&b_Tau_leadChargedCandTrackFitErrorMatrix_11); 
 tree->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_12",&Tau_leadChargedCandTrackFitErrorMatrix_12,&b_Tau_leadChargedCandTrackFitErrorMatrix_12);
 tree->SetBranchAddress("Tau_leadChargedCandTrackFitErrorMatrix_22",&Tau_leadChargedCandTrackFitErrorMatrix_22,&b_Tau_leadChargedCandTrackFitErrorMatrix_22);
 tree->SetBranchAddress("Tau_defaultDxy",&Tau_defaultDxy,&b_Tau_defaultDxy); 
 tree->SetBranchAddress("Tau_defaultDxyError",&Tau_defaultDxyError,&b_Tau_defaultDxyError); 
 tree->SetBranchAddress("Tau_defaultDxySig",&Tau_defaultDxySig,&b_Tau_defaultDxySig); 
 tree->SetBranchAddress("Tau_defaultFlightLengthSig",&Tau_defaultFlightLengthSig,&b_Tau_defaultFlightLengthSig);
 tree->SetBranchAddress("Tau_default_PCAx_pv",&Tau_default_PCAx_pv,&b_Tau_default_PCAx_pv); 
 tree->SetBranchAddress("Tau_default_PCAy_pv",&Tau_default_PCAy_pv,&b_Tau_default_PCAy_pv); 
 tree->SetBranchAddress("Tau_default_PCAz_pv",&Tau_default_PCAz_pv,&b_Tau_default_PCAz_pv);
 tree->SetBranchAddress("Tau_defaultFlightLengthX",&Tau_defaultFlightLengthX,&b_Tau_defaultFlightLengthX); 
 tree->SetBranchAddress("Tau_defaultFlightLengthY",&Tau_defaultFlightLengthY,&b_Tau_defaultFlightLengthY); 
 tree->SetBranchAddress("Tau_defaultFlightLengthZ",&Tau_defaultFlightLengthZ,&b_Tau_defaultFlightLengthZ);
 tree->SetBranchAddress("Tau_decayModeFinding",&Tau_decayModeFinding,&b_Tau_decayModeFinding); 
 tree->SetBranchAddress("Tau_decayModeFindingOldDMs",&Tau_decayModeFindingOldDMs,&b_Tau_decayModeFindingOldDMs); 
 tree->SetBranchAddress("Tau_decayModeFindingNewDMs",&Tau_decayModeFindingNewDMs,&b_Tau_decayModeFindingNewDMs);
 tree->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr",&Tau_byLooseCombinedIsolationDeltaBetaCorr,&b_Tau_byLooseCombinedIsolationDeltaBetaCorr); 
 tree->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits",&Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits,&b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
 tree->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr",&Tau_byMediumCombinedIsolationDeltaBetaCorr,&b_Tau_byMediumCombinedIsolationDeltaBetaCorr); 
 tree->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits",&Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits,&b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
 tree->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr",&Tau_byTightCombinedIsolationDeltaBetaCorr,&b_Tau_byTightCombinedIsolationDeltaBetaCorr); 
 tree->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr3Hits",&Tau_byTightCombinedIsolationDeltaBetaCorr3Hits,&b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
 tree->SetBranchAddress("Tau_byLooseIsolationMVA3newDMwLT",&Tau_byLooseIsolationMVA3newDMwLT,&b_Tau_byLooseIsolationMVA3newDMwLT); 
 tree->SetBranchAddress("Tau_byLooseIsolationMVA3newDMwoLT",&Tau_byLooseIsolationMVA3newDMwoLT,&b_Tau_byLooseIsolationMVA3newDMwoLT); 
 tree->SetBranchAddress("Tau_byMediumIsolationMVA3newDMwLT",&Tau_byMediumIsolationMVA3newDMwLT,&b_Tau_byMediumIsolationMVA3newDMwLT);
 tree->SetBranchAddress("Tau_byMediumIsolationMVA3newDMwoLT",&Tau_byMediumIsolationMVA3newDMwoLT,&b_Tau_byMediumIsolationMVA3newDMwoLT); 
 tree->SetBranchAddress("Tau_byLooseIsolationMva3oldDMwLT",&Tau_byLooseIsolationMva3oldDMwLT,&b_Tau_byLooseIsolationMva3oldDMwLT); 
 tree->SetBranchAddress("Tau_byLooseIsolationMVA3oldDMwoLT",&Tau_byLooseIsolationMVA3oldDMwoLT,&b_Tau_byLooseIsolationMVA3oldDMwoLT);
 tree->SetBranchAddress("Tau_byMediumIsolationMva3oldDMwLT",&Tau_byMediumIsolationMva3oldDMwLT,&b_Tau_byMediumIsolationMva3oldDMwLT); 
 tree->SetBranchAddress("Tau_byMediumIsolationMVA3oldDMwoLT",&Tau_byMediumIsolationMVA3oldDMwoLT,&b_Tau_byMediumIsolationMVA3oldDMwoLT); 
 tree->SetBranchAddress("Tau_byTightIsolationMVA3newDMwLT",&Tau_byTightIsolationMVA3newDMwLT,&b_Tau_byTightIsolationMVA3newDMwLT);
 tree->SetBranchAddress("Tau_byTightIsolationMVA3newDMwoLT",&Tau_byTightIsolationMVA3newDMwoLT,&b_Tau_byTightIsolationMVA3newDMwoLT);
 tree->SetBranchAddress("Tau_byTightIsolationMva3oldDMwLT",&Tau_byTightIsolationMva3oldDMwLT,&b_Tau_byTightIsolationMva3oldDMwLT); 
 tree->SetBranchAddress("Tau_byTightIsolationMVA3oldDMwoLT",&Tau_byLooseIsolationMVA3oldDMwoLT,&b_Tau_byLooseIsolationMVA3oldDMwoLT);
 tree->SetBranchAddress("Tau_againstMuonLoose2",&Tau_againstMuonLoose2,&b_Tau_againstMuonLoose2); 
 tree->SetBranchAddress("Tau_againstMuonLoose3",&Tau_againstMuonLoose3,&b_Tau_againstMuonLoose3); 
 tree->SetBranchAddress("Tau_againstMuonTight2",&Tau_againstMuonTight2,&b_Tau_againstMuonTight2); 
 tree->SetBranchAddress("Tau_againstMuonTight3",&Tau_againstMuonTight3,&b_Tau_againstMuonTight3);
 tree->SetBranchAddress("Tau_againstElectronMVALooseMVA5",&Tau_againstElectronMVALooseMVA5,&b_Tau_againstElectronMVALooseMVA5); 
 tree->SetBranchAddress("Tau_againstElectronMVAMediumMVA5",&Tau_againstElectronMVAMediumMVA5,&b_Tau_againstElectronMVALooseMVA5); 
 tree->SetBranchAddress("Tau_againstElectronMVATightMVA5",&Tau_againstElectronMVATightMVA5,&b_Tau_againstElectronMVATightMVA5);
 tree->SetBranchAddress("Tau_byVLooseCombinedIsolationDeltaBetaCorr",&Tau_byVLooseCombinedIsolationDeltaBetaCorr,&b_Tau_byVLooseCombinedIsolationDeltaBetaCorr);




*/


 cout << "Cut sequence:\n";
 //Take entries for each variable of cut
 double cut_eff[numTotCut];
 for(int iniz=0; iniz<numTotCut; iniz++) cut_eff[iniz] = 0;
 double weight = 0;
 //All entries

 for(int i=0; i<nentries; i++){
  Long64_t tentry = tree->LoadTree(i);
  b_Muon_pt->GetEntry(tentry);
  b_Muon_eta->GetEntry(tentry);
  

  b_patElectron_pt->GetEntry(tentry);
  b_patElectron_eta->GetEntry(tentry);
 
  b_Tau_pt->GetEntry(tentry);
  
  Tau_pt->size() > 1 ? cout << "Tau 1 pt: " << Tau_pt->at(0) << " Tau 2 pt: " << Tau_pt->at(1) << "\n" : cout << "No tau_h pair for event " << i << "\n";



 
  //cout << "Muon Pt->size = " << Muon_pt->size() << "\n";
  //cout << "Muon eta size = " << Muon_eta->size() << "\n";
  //cout << "electron pt size = " << patElectron_pt->size() << "\n";
  //cout << "electron eta size = " << patElectron_eta->size() << "\n";

  //cout << "Muon Pt: " << Muon_pt->at(0) << "\n";
  //cout << "Muon eta: " << Muon_eta->at(0) << "\n";
  if(patElectron_pt->size() > 0){
   cout << "Electron Pt: " << patElectron_pt->at(0) << "\n";
   cout << "Electron eta: " << patElectron_eta->at(0) << "\n";
  }
  if(Muon_pt->size() > 0 && patElectron_pt->size() > 0){
    if(noPUcorr) wgt_pu = 1.;
    if(Muon_pt->at(0) >= cand1_acc[0]){
      cut_eff[0] += wgt_pu;
      if(Muon_eta->at(0) < cand1_acc[1]){
        cut_eff[1] += wgt_pu;
        if(patElectron_pt->at(0) >= cand2_acc[0]){
          cut_eff[2] += wgt_pu;
          if(patElectron_eta->at(0) < cand2_acc[1]){
            cut_eff[3] += wgt_pu;
            cout << "Passing Muon pt: " << Muon_pt->at(0) << "\n";
            cout << "Passing Electron Et: " << patElectron_pt->at(0) << "\n";
          }
        }
      }      
    }
  }
  weight = wgt_lumi;
 
 }//End all entries

 weight = weight*Luminosity/100.;
 //Evaluate efficiencies
 measure_efficiencies(nentries, cut_eff, varname, varcut, weight);
}
/////
//   Initialize the strings 
/////
void initialize_strings(vector<string> &varname, vector<string> &varcut){
  varname[0] = "Muon Pt:"; varname[1] = "Muon eta:";
  varname[2] = "Electron Et:"; varname[3] = "Electron SC eta:";
 //varname[0] = "$\\mu$ Pt"; varname[1] = "$\\mu$ $|\\eta|$ ";
 //varname[2] = "e Et "; varname[3] = "e $|\\eta_{SC}|$ ";
}
/////
//   Print initial strings
/////
void ini_print(int nentries){
 //\def\tablepagesize{\fontsize{7.5pt}{6pt}\selectfont} 
 //cout<<"\\vskip-5pt"<<endl;
 //cout<<"{\\tablepagesize"<<endl;
 //cout<<"\\begin{table}[htdp]"<<endl;
 //cout<<"\\begin{center}"<<endl;
 //cout<<"\\begin{tabular}{|c|c|c|c|}"<<endl;
 //cout<<"\\hline"<<endl;
 //cout<<"\\rowcolor{black} {\\color{white!90}Selection} & {\\color{white!90}Evt} & {\\color{white!90}Rel Eff (\\%)} & {\\color{white!90}Cum Eff (\\%)} \\\\"<<endl;
 //cout<<"{\\cellcolor{white!90}"<<"NoCuts} & "<<"{\\cellcolor{white!90}"<<nentries<<"} & {\\cellcolor{white!90} 100$\\pm$0} & {\\cellcolor{white!90}100$\\pm$0} \\\\"<<endl;
}
/////
//   Measure efficiencies of cuts
/////
void measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight){
 cout << "Measuring Efficiencies:\n";
 ini_print(nentries);
 double eff_Rel, err_eff_Rel, eff_Cum, err_eff_Cum;
 cout<<setiosflags(ios::fixed)<<setprecision(2);
 for(int i=0; i<numTotCut; i++){
  if(i==0){
   eff_Rel = (double) cut_eff[i]/nentries;
   err_eff_Rel = sqrt(eff_Rel*(1-eff_Rel)/(double)nentries);
  }else{
   eff_Rel = (double) cut_eff[i]/cut_eff[i-1];
   err_eff_Rel = sqrt(eff_Rel*(1-eff_Rel)/(double)cut_eff[i-1]);
  }
  eff_Cum = (double)cut_eff[i]/nentries;
  err_eff_Cum = sqrt(eff_Cum*(1-eff_Cum)/(double)nentries);
  //cout<<"\\hline"<<endl;
  if(noLumiNorm){
    cout << varname[i] << varcut[i] << " " << cut_eff[i] << " eff_Rel: " << eff_Rel*100 << " +/- " << err_eff_Rel*100 << " eff_Cum: " << eff_Cum*100 << " +/- " << err_eff_Cum*100 <<"\n";  
//cout<<"{\\cellcolor{white!90}"<<varname[i].c_str()<<varcut[i].c_str()<<"} & {\\cellcolor{white!90}"<<cut_eff[i]<<"} & {\\cellcolor{white!90}"<<eff_Rel*100<<"$\\pm$"<<err_eff_Rel*100<<"} & {\\cellcolor{white!90}"<<eff_Cum*100<<"$\\pm$"<<err_eff_Cum*100<<"}\\\\"<<endl;
  }else{
    cout << varname[i] << varcut[i] << " " << cut_eff[i]*weight << " eff_Rel: " << eff_Rel*100 << " +/- " << err_eff_Rel*100 << " eff_Cum: " << eff_Cum*100 << " +/- " << err_eff_Cum*100 <<"\n";   
//cout<<"{\\cellcolor{white!90}"<<varname[i].c_str()<<varcut[i].c_str()<<"} & {\\cellcolor{white!90}"<<cut_eff[i]*weight<<"} & {\\cellcolor{white!90}"<<eff_Rel*100<<"$\\pm$"<<err_eff_Rel*100<<"} & {\\cellcolor{white!90}"<<eff_Cum*100<<"$\\pm$"<<err_eff_Cum*100<<"}\\\\"<<endl;
  }
 }
 //cout<<"\\hline"<<endl;
 //cout<<"\\end{tabular}"<<endl;
 //cout<<"\\end{center}"<<endl;
 //cout<<"\\end{table}"<<endl;
 //cout<<"}"<<endl;
}
