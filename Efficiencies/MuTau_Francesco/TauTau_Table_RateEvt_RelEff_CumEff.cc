/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects and events 
   It prints a table following the cutFlow for rateEvt or relEff or cumEff, for all the files specified

Need to specify
1. See Declare constants
*/
/////
//   To run: root -l TauTau_Table_RateEvt_RelEff_CumEff.cc
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
#include <TLorentzVector.h>
#include <TMath.h>
#include <Math/SVector.h>
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "DataFormats/Math/interface/deltaR.h"
using namespace ROOT::Math;
/////
//   Declare constants
/////
//Path and root file
const string path           = "/uscms_data/d3/andrewj/CMSSW_7_4_7/src/Efficiencies/MuTau_Francesco/";
const char *samples[]       = {"Zprime"};
//Corrections
const double Luminosity     = 19600; //pb^-1
const bool noLumiNorm       = true; //true means NO luminosity normalization done
const bool noPUcorr         = true; //true means NO PU corr done
//Choose among relEff, cumEff, or rateEvt 
const bool relEff           = true;
const bool cumEff           = false;
const bool rateEvt          = false; 
//Acceptance and object selection
const double cand1_acc[2]   = {45, 2.1};
const double cand2_acc[2]   = {45, 2.1};
const int numtot_cuts       = 18;
//Signal selection
//charge==-1 && cosDphi<-0.95 && met>=30 && pZetaMt-3.051*pZetaVisMt>-50 && cosDphiLMet<0.2 && nbjet==0
//const int sr_cuts           = 6;
//double sigreg_cuts[sr_cuts] = {-1,-0.95,20,-50,0.2,0};
const double SETPRECISION = 3;
//For histos
const int const_bin = 200; const double const_from = 0; const double const_to = 1000;



/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
void cutflow_eff(TTree* tree, int nentries, int pos, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
void measure_efficiencies(int nentries, double cut_eff[], double weight, int pos, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]);
void print_table(vector<string> &rootplas, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]);
pair<bool, TLorentzVector> matchTauToGen(TLorentzVector p4Vector, TTree* tree, int nentry);
pair<bool, pair<TLorentzVector,TLorentzVector>> genHadronicTauPairExists(TTree* tree, int nentry);
pair<bool,TLorentzVector> leptonTau_matching(TTree* tree, int entry, int candid, TLorentzVector reco_cand);



TH1F* theFirstTauPt = new TH1F("theFirstTauPt", "pT of the first gen vis tau", const_bin, const_from, const_to);
TH1F* theSecondTauPt = new TH1F("theSecondTauPt", "pT of the second gen vis tau", 100, 0, 400);
TH1F* theFirstMatchedTauPt = new TH1F("theFirstMatchedTauPt", "pT of the first matched gen vis tau", const_bin, const_from, const_to);
TH1F* theSecondMatchedTauPt = new TH1F("theSecondMatchedTauPt", "pT of the second matched gen vis tau", 100, 0, 400);
TH1F* theFirstGenTauEta = new TH1F("theFirstGenTauEta", "eta of the first gen vis tau", 100, -2.5, 2.5);
TH1F* theSecondGenTauEta = new TH1F("theSecondGenTauEta", "eta of the second gen vis tau", 100, -2.5, 2.5);
TH1F* theFirstRecoTauEta = new TH1F("theFirstRecoTauEta", "eta of the first reco vis tau", 100, -2.5, 2.5);
TH1F* theSecondRecoTauEta = new TH1F("theSecondRecoTauEta", "eta of the second reco vis tau", 100, -2.5, 2.5);
TH1F* theFirstGenTauPhi = new TH1F("theFirstGenTauPhi", "phi of the first gen vis tau", 100, -4, 4);
TH1F* theSecondGenTauPhi = new TH1F("theSecondGenTauPhi", "phi of the second gen vis tau", 100, -4, 4);
TH1F* theFirstRecoTauPhi = new TH1F("theFirstRecoTauPhi", "phi of the first reco vis tau", 100, -4, 4);
TH1F* theSecondRecoTauPhi = new TH1F("theSecondRecoTauPhi", "phi of the second reco vis tau", 100, -4, 4);



/////
//   Main function
/////
void TauTau_Table_RateEvt_RelEff_CumEff() {
 //Run over all samples 
 vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
 const uint rootplas_size = rootplas.size();
 //Matrix with values of rel eff
 double sample_cuts[rootplas_size][numtot_cuts];
 double sample_cuts_err[rootplas_size][numtot_cuts];
 for(uint i=0; i<rootplas_size; i++){
  //Call tree  
  TFile* f = Call_TFile(rootplas[i]);
  cout << "Initializing tree";
  TTree* tree; f->GetObject("TNT/BOOM;1",tree);
  int nentries = tree->GetEntries(); 
  cout << "nentries = " << nentries << "\n";
  //Take values 
  cutflow_eff(tree,nentries,i,sample_cuts,sample_cuts_err);
 }
 //Print table with rel eff values
 print_table(rootplas,sample_cuts,sample_cuts_err);
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla){
 string file_name = path+rootpla+".root";
 cout<<"file_name "<<file_name<<endl;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Efficiency for the object
/////
void cutflow_eff(TTree* tree, int nentries, int pos, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]){ 
 /////
 //   Get variables
 /////
 //Tau
 vector<double> *Tau_pt, *Tau_eta, *Tau_phi, *Tau_energy, *Tau_leadChargedCandDz_pv;
 vector<int> *Tau_decayModeFindingNewDMs, *Tau_againstElectronMVALooseMVA5, *Tau_againstMuonLoose3, *Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
 Tau_pt = 0;
 Tau_eta = 0;
 Tau_phi = 0;
 Tau_energy = 0;
 Tau_leadChargedCandDz_pv = 0;
 Tau_decayModeFindingNewDMs = 0;
 Tau_againstElectronMVALooseMVA5 = 0;
 Tau_againstMuonLoose3 = 0;
 Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
 TBranch *b_Tau_pt = 0;
 TBranch *b_Tau_eta = 0;
 TBranch *b_Tau_phi = 0;
 TBranch *b_Tau_energy = 0;
 TBranch *b_Tau_leadChargedCandDz_pv = 0;
 TBranch *b_Tau_decayModeFindingNewDMs = 0;
 TBranch *b_Tau_againstElectronMVALooseMVA5 = 0;
 TBranch *b_Tau_againstMuonLoose3 = 0;
 TBranch *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0; 
 tree->SetBranchAddress("Tau_pt",&Tau_pt,&b_Tau_pt);
 tree->SetBranchAddress("Tau_eta",&Tau_eta,&b_Tau_eta);
 tree->SetBranchAddress("Tau_phi",&Tau_phi,&b_Tau_phi);
 tree->SetBranchAddress("Tau_energy",&Tau_energy,&b_Tau_energy);
 //tree->SetBranchAddress("Tau_leadChargedCandDz_pv",&Tau_leadChargedCandDz_pv,&b_Tau_leadChargedCandDz_pv);
 tree->SetBranchAddress("Tau_decayModeFindingNewDMs",&Tau_decayModeFindingNewDMs,&b_Tau_decayModeFindingNewDMs);
 tree->SetBranchAddress("Tau_againstElectronMVALooseMVA5",&Tau_againstElectronMVALooseMVA5,&b_Tau_againstElectronMVALooseMVA5);
 tree->SetBranchAddress("Tau_againstMuonLoose3",&Tau_againstMuonLoose3,&b_Tau_againstMuonLoose3);
 tree->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits",&Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits,&b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
 /////
 //   Start the cut sequence
 /////
 cout << "Cut sequence:\n";
 //Take entries for each variable of cut
 double cut_eff[numtot_cuts];
 for(int iniz=0; iniz<numtot_cuts; iniz++) cut_eff[iniz] = 0;
 double weight = 0;
 //All entries
 for(int i=0; i<nentries; i++){
  Long64_t tentry = tree->LoadTree(i);
 
  //Tau
  b_Tau_pt->GetEntry(tentry);
  b_Tau_eta->GetEntry(tentry);
  b_Tau_phi->GetEntry(tentry);
  b_Tau_energy->GetEntry(tentry);
 // b_Tau_leadChargedCandDz_pv->GetEntry(tentry);
  b_Tau_decayModeFindingNewDMs->GetEntry(tentry);
  b_Tau_againstElectronMVALooseMVA5->GetEntry(tentry);
  b_Tau_againstMuonLoose3->GetEntry(tentry);
  b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->GetEntry(tentry);
  //Get values
  double wgt_pu   = 1.;
  double wgt_lumi = 1.;
  if(noPUcorr) wgt_pu = 1.;


  if(genHadronicTauPairExists(tree, i).first){
    cut_eff[0] += wgt_pu;
    if(Tau_pt->size() >= 2){
      cut_eff[1] += wgt_pu;
      bool twoMatchedTaus = false;
      for(uint j = 0; j < Tau_pt->size(); j++){
        for(uint k = j+1; k < Tau_pt->size(); k++){
          TLorentzVector Tau1Vector;
	  Tau1Vector.SetPtEtaPhiE(Tau_pt->at(0), Tau_eta->at(0), Tau_phi->at(0), Tau_energy->at(0));
          TLorentzVector Tau2Vector;
	  Tau2Vector.SetPtEtaPhiE(Tau_pt->at(1), Tau_eta->at(1), Tau_phi->at(1), Tau_energy->at(1));
          if(leptonTau_matching(tree, i, 15, Tau1Vector).first && leptonTau_matching(tree, i, 15, Tau2Vector).first){
            twoMatchedTaus = true;
            theFirstMatchedTauPt->Fill(matchTauToGen(Tau1Vector, tree, i).second.Pt());
            theSecondMatchedTauPt->Fill(matchTauToGen(Tau2Vector, tree, i).second.Pt());
          }
        }
      }
      if(twoMatchedTaus){
        cut_eff[2] += wgt_pu;
        if(Tau_pt->at(0) >= cand1_acc[0]){
          cut_eff[3] += wgt_pu;
          if(Tau_eta->at(0) < cand1_acc[1]){
            cut_eff[4] += wgt_pu;
            if(Tau_pt->at(1) >= cand2_acc[0]){
              cut_eff[5] += wgt_pu;
              if(Tau_eta->at(1) < cand2_acc[1]){
                cut_eff[6] += wgt_pu;
                if(Tau_decayModeFindingNewDMs->at(0) > 0.5){
                  cut_eff[7] += wgt_pu;
                  if(/*Tau_leadChargedCandDz_pv->at(0) < 0.2*/1==1){
                    cut_eff[8] += wgt_pu;
                    if(Tau_decayModeFindingNewDMs->at(1) > 0.5){
                      cut_eff[9] += wgt_pu;
                      if(/*Tau_leadChargedCandDz_pv->at(1) < 0.2*/1==1){
                        cut_eff[10] += wgt_pu;
                        double DeltaR = deltaR(Tau_eta->at(0), Tau_phi->at(0), Tau_eta->at(1), Tau_phi->at(1));
                        if(DeltaR > 0.5){
                          cut_eff[11] += wgt_pu;
                          if(Tau_againstElectronMVALooseMVA5->at(0) > 0.5){
                            cut_eff[12] += wgt_pu;
                            if(Tau_againstMuonLoose3->at(0) > 0.5){
                              cut_eff[13] += wgt_pu;
                              if(Tau_againstElectronMVALooseMVA5->at(1) > 0.5){
                                cut_eff[14] += wgt_pu;
                                if(Tau_againstMuonLoose3->at(1) > 0.5){
                                  cut_eff[15] += wgt_pu;
                                  if(Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(0) < 1.0){
                                    cut_eff[16] += wgt_pu;
                                    if(Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(1) < 1.0){
                                      cut_eff[17] += wgt_pu;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  weight = wgt_lumi;
 }//End all entries
 //To be changed
 weight = weight*Luminosity/100.;
 /////
 //   Define strings to print
 /////
 vector<string> varname(numtot_cuts); 
 vector<string> varcut(numtot_cuts);
 initialize_strings(varname,varcut);
 cout << "Strings initialized:\n";
 //Evaluate efficiencies
 measure_efficiencies(nentries, cut_eff, weight, pos, sample_cuts, sample_cuts_err);
}


pair<bool, pair<TLorentzVector, TLorentzVector>> genHadronicTauPairExists(TTree* tree, int nentry){

  vector<double> *Gen_pt, *Gen_eta, *Gen_phi, *Gen_energy;
  vector<int> *Gen_status, *Gen_pdg_id, *Gen_BmotherIndex;
  Gen_pt = 0;
  Gen_eta = 0;
  Gen_phi = 0;
  Gen_energy = 0;
  Gen_status = 0;
  Gen_pdg_id = 0;
  Gen_BmotherIndex = 0;
  TBranch *b_Gen_pt = 0;
  TBranch *b_Gen_eta = 0;
  TBranch *b_Gen_phi = 0;
  TBranch *b_Gen_energy = 0;
  TBranch *b_Gen_status = 0;
  TBranch *b_Gen_pdg_id = 0;
  TBranch *b_Gen_BmotherIndex = 0;

  tree->SetBranchAddress("Gen_pt",&Gen_pt,&b_Gen_pt);
  tree->SetBranchAddress("Gen_eta",&Gen_eta,&b_Gen_eta);
  tree->SetBranchAddress("Gen_phi",&Gen_phi,&b_Gen_phi);
  tree->SetBranchAddress("Gen_energy",&Gen_energy,&b_Gen_energy);
  tree->SetBranchAddress("Gen_status",&Gen_status,&b_Gen_status);
  tree->SetBranchAddress("Gen_pdg_id",&Gen_pdg_id,&b_Gen_pdg_id);
  tree->SetBranchAddress("Gen_BmotherIndex",&Gen_BmotherIndex,&b_Gen_BmotherIndex);


  Long64_t tentry = tree->LoadTree(nentry);
  b_Gen_pt->GetEntry(tentry);
  b_Gen_eta->GetEntry(tentry);
  b_Gen_phi->GetEntry(tentry);
  b_Gen_energy->GetEntry(tentry);
  b_Gen_status->GetEntry(tentry);
  b_Gen_pdg_id->GetEntry(tentry);
  b_Gen_BmotherIndex->GetEntry(tentry);

  bool genHadTauPair = false;
  bool IsItAHadronicDecay; // boolean used to specify whether a gen tau lepton decays hadronically
  pair<TLorentzVector, TLorentzVector> genTau4Momenta;
  TLorentzVector theFirstGenTau(0,0,0,0); // initialize the 4-momentum vector of the vis gen tau matched to the reco tau
  genTau4Momenta.first = theFirstGenTau;
  TLorentzVector theSecondGenTau(0,0,0,0);
  genTau4Momenta.second = theSecondGenTau;
  vector<int> tempTauIndexVector; // vector which contains the index (in the gen collection) of the gen level taus (before decaying)
  tempTauIndexVector.clear(); // clear any previous declaration from memory
  vector<bool> IsItAHadronicDecayVector; // vector of booleans which contain the information about whether a gen tau lepton decays hadronically
  IsItAHadronicDecayVector.clear(); // clear any previous declaration from memory
  TLorentzVector theNeutrino(0,0,0,0);
  vector<TLorentzVector> theNeutrinoVector; //vector of gen neutrinos used to subtract from gen taus in order to get visible momentum
  theNeutrinoVector.clear();

  //---Loop over gen particles to find the tau neutrinos and then store the index of each tau neutrino's mother (a tau).
  for(uint jj = 0; jj < Gen_pt->size(); jj++) {
    if(Gen_BmotherIndex->at(jj) >= 0 && (abs(Gen_pdg_id->at(jj)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(jj))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(jj)) == 2)
    /*&& Gen_pt->at(Gen_BmotherIndex->at(jj)) > 45 && Gen_eta->at(Gen_BmotherIndex->at(jj)) < 2.1*/ ) {
      tempTauIndexVector.push_back(Gen_BmotherIndex->at(jj));
      theNeutrino.SetPtEtaPhiE(Gen_pt->at(jj), Gen_eta->at(jj), Gen_phi->at(jj), Gen_energy->at(jj));
      theNeutrinoVector.push_back(theNeutrino);
    }
  }

  //require exactly 2 gen taus
  if(tempTauIndexVector.size() == 2) {
    //cout << "entry " << nentry << "has two taus\n";
    //---Loop over the gen taus and determine whether it decays hadronically.
    for(uint jjj = 0; jjj < tempTauIndexVector.size(); jjj++) {
      IsItAHadronicDecay = true;
      for(uint jj = 0; jj < Gen_pt->size(); jj++) {
        if(Gen_BmotherIndex->at(jj) >= 0 &&  ((abs(Gen_pdg_id->at(jj)) == 12) || (abs(Gen_pdg_id->at(jj)) == 14)) && (Gen_BmotherIndex->at(jj) == tempTauIndexVector.at(jjj)) ) {
          IsItAHadronicDecay = false; // it is not a hadronic tau decay since it decayed to a electron/muon neutrino
        }
      }
      IsItAHadronicDecayVector.push_back(IsItAHadronicDecay);
    }

    for(uint i = 0; i < tempTauIndexVector.size(); i++){
      for(uint j = i+1; j < tempTauIndexVector.size(); j++){
        if(IsItAHadronicDecayVector.at(i) && IsItAHadronicDecayVector.at(j) /*&& Gen_BmotherIndex->at(tempTauIndexVector.at(i)) == Gen_BmotherIndex->at(tempTauIndexVector.at(j))*/){
          genHadTauPair = true;
          theFirstGenTau.SetPtEtaPhiE(Gen_pt->at(tempTauIndexVector.at(i)), Gen_eta->at(tempTauIndexVector.at(i)), Gen_phi->at(tempTauIndexVector.at(i)), Gen_energy->at(tempTauIndexVector.at(i)));
          theFirstGenTau = theFirstGenTau - theNeutrinoVector.at(i);
          theSecondGenTau.SetPtEtaPhiE(Gen_pt->at(tempTauIndexVector.at(j)), Gen_eta->at(tempTauIndexVector.at(j)), Gen_phi->at(tempTauIndexVector.at(j)), Gen_energy->at(tempTauIndexVector.at(j)));
          theSecondGenTau = theSecondGenTau - theNeutrinoVector.at(j);
          genTau4Momenta.first = theFirstGenTau;
          genTau4Momenta.second = theSecondGenTau;
        }
      }
    }
    pair<bool, pair<TLorentzVector, TLorentzVector>> twoTauInformation(genHadTauPair, genTau4Momenta);
    return twoTauInformation;
  }
  else{
    pair<bool, pair<TLorentzVector, TLorentzVector>> twoTauInformation(genHadTauPair, genTau4Momenta);
    return twoTauInformation;
  }
}



pair<bool, TLorentzVector> matchTauToGen(TLorentzVector p4Vector, TTree* tree, int nentry) {

  vector<double> *Gen_pt, *Gen_eta, *Gen_phi, *Gen_energy;
  vector<int> *Gen_status, *Gen_pdg_id, *Gen_BmotherIndex;
  Gen_pt = 0;
  Gen_eta = 0;
  Gen_phi = 0;
  Gen_energy = 0;
  Gen_status = 0;
  Gen_pdg_id = 0;
  Gen_BmotherIndex = 0;
  TBranch *b_Gen_pt = 0;
  TBranch *b_Gen_eta = 0;
  TBranch *b_Gen_phi = 0;
  TBranch *b_Gen_energy = 0;
  TBranch *b_Gen_status = 0;
  TBranch *b_Gen_pdg_id = 0;
  TBranch *b_Gen_BmotherIndex = 0;

  tree->SetBranchAddress("Gen_pt",&Gen_pt,&b_Gen_pt);
  tree->SetBranchAddress("Gen_eta",&Gen_eta,&b_Gen_eta);
  tree->SetBranchAddress("Gen_phi",&Gen_phi,&b_Gen_phi);
  tree->SetBranchAddress("Gen_energy",&Gen_energy,&b_Gen_energy);
  tree->SetBranchAddress("Gen_status",&Gen_status,&b_Gen_status);
  tree->SetBranchAddress("Gen_pdg_id",&Gen_pdg_id,&b_Gen_pdg_id);
  tree->SetBranchAddress("Gen_BmotherIndex",&Gen_BmotherIndex,&b_Gen_BmotherIndex);

  Long64_t tentry = tree->LoadTree(nentry);
  b_Gen_pt->GetEntry(tentry);
  b_Gen_eta->GetEntry(tentry);
  b_Gen_phi->GetEntry(tentry);
  b_Gen_energy->GetEntry(tentry);
  b_Gen_status->GetEntry(tentry);
  b_Gen_pdg_id->GetEntry(tentry);
  b_Gen_BmotherIndex->GetEntry(tentry);

  //Do Gen-Tau_h matching

  bool isGenMatched = false; // by default, the reco-gen tau matching boolean is set to false
  bool IsItAHadronicDecay; // boolean used to specify whether a gen tau lepton decays hadronically
  TLorentzVector theGenObject(0,0,0,0); // initialize the 4-momentum vector of the vis gen tau matched to the reco tau
  TLorentzVector theNeutrinoObject(0,0,0,0); // initialize the 4-momentum vector of the tau neutrino from the tau decay
  TLorentzVector theRecoTau(0,0,0,0); // initialize the 4-momentum vector of the reco tau to match
  vector<bool> IsItAHadronicDecayVector; // vector of booleans which contain the information about whether a gen tau lepton decays hadronically
  IsItAHadronicDecayVector.clear(); // clear any previous declaration from memory
  vector<int> tempTauIndexVector; // vector which contains the index (in the gen collection) of the gen level taus (before decaying)
  tempTauIndexVector.clear(); // clear any previous declaration from memory
  vector<TLorentzVector> tempNeutrinoMomentumVector; // vector of lorentz 4-momentum vectors for each tau neutrino from the tau decay
  tempNeutrinoMomentumVector.clear(); // clear any previous declaration from memory

  //---Loop over gen particles to find the tau neutrinos and then store the index of each tau neutrino's mother (a tau).
  //---Also store the tau neutrino's 4-momentum vector in order to calculate the visible tau 4-momentum at a later point.


  for(int jj = 0; jj < Gen_pt->size(); jj++) {
    if(Gen_BmotherIndex->at(jj) >= 0 && (abs(Gen_pdg_id->at(jj)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(jj))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(jj)) == 2) ) {
      tempTauIndexVector.push_back(Gen_BmotherIndex->at(jj));
      theNeutrinoObject.SetPtEtaPhiE(Gen_pt->at(jj), Gen_eta->at(jj), Gen_phi->at(jj), Gen_energy->at(jj));
      tempNeutrinoMomentumVector.push_back(theNeutrinoObject);
    }
  }

  //---Perform matching only if there is at least one gen tau in the event
  if(tempTauIndexVector.size() < 2 ) cout << "Size of tau index vector = " << tempTauIndexVector.size() << "\n";
  if(tempTauIndexVector.size() > 0) {
    //---Loop over the gen taus and determine whether it decays hadronically.

    for(int jjj = 0; jjj < tempTauIndexVector.size(); jjj++) {
      IsItAHadronicDecay = true;
      for(int jj = 0; jj < Gen_pt->size(); jj++) {
        if( ((abs(Gen_pdg_id->at(jj)) == 12) || (abs(Gen_pdg_id->at(jj)) == 14)) && (Gen_BmotherIndex->at(jj) == tempTauIndexVector.at(jjj)) ) {
          IsItAHadronicDecay = false; // it is not a hadronic tau decay since it decayed to a electron/muon neutrino
        }
      }
      IsItAHadronicDecayVector.push_back(IsItAHadronicDecay);
    }
    //---Loop over the gen taus and calculate the 4-momentum of the visible products (i.e. subtract the 4-momentum of the tau neutrino)

    for(int jjj = 0; jjj < tempTauIndexVector.size(); jjj++) {
      for(int jj = 0; jj < Gen_pt->size(); jj++) {
        if(jj == tempTauIndexVector.at(jjj)) {
          theGenObject.SetPtEtaPhiE(Gen_pt->at(jj), Gen_eta->at(jj), Gen_phi->at(jj), Gen_energy->at(jj)); // 4-momentum of the gen tau
          theGenObject = theGenObject - tempNeutrinoMomentumVector.at(jjj); // subtract the 4-momentum of the tau neutrino (visible tau)
          if( (IsItAHadronicDecayVector.at(jjj)) && p4Vector.DeltaR(theGenObject) <= 0.3){
            isGenMatched = true;
          }
        }
      }
    }
    pair<bool, TLorentzVector> GenMatchedInformation(isGenMatched,theGenObject);
    return GenMatchedInformation;
  } else {
    pair<bool, TLorentzVector> GenMatchedInformation(isGenMatched,theGenObject);
    return GenMatchedInformation;
  }

}



pair<bool,TLorentzVector> leptonTau_matching(TTree* tree, int entry, int candid, TLorentzVector reco_cand){
 bool matched = false; TLorentzVector gen_cand(0.,0.,0.,0.);
 vector<double> *Gen_pdg_id, *Gen_motherpdg_id, *Gen_pt, *Gen_eta, *Gen_phi, *Gen_energy;
 Gen_pdg_id = 0;
 Gen_motherpdg_id = 0;
 Gen_pt = 0;
 Gen_eta = 0;
 Gen_phi = 0;
 Gen_energy = 0;
 TBranch *b_Gen_pdg_id = 0;
 TBranch *b_Gen_motherpdg_id = 0;
 TBranch *b_Gen_pt = 0;
 TBranch *b_Gen_eta = 0;
 TBranch *b_Gen_phi = 0;
 TBranch *b_Gen_energy = 0;
 tree->SetBranchAddress("Gen_pdg_id",&Gen_pdg_id,&b_Gen_pdg_id);
 tree->SetBranchAddress("Gen_motherpdg_id",&Gen_motherpdg_id,&b_Gen_motherpdg_id);
 tree->SetBranchAddress("Gen_pt",&Gen_pt,&b_Gen_pt);
 tree->SetBranchAddress("Gen_eta",&Gen_eta,&b_Gen_eta);
 tree->SetBranchAddress("Gen_phi",&Gen_phi,&b_Gen_phi);
 tree->SetBranchAddress("Gen_energy",&Gen_energy,&b_Gen_energy);
 vector<int> *Gen_BmotherIndex;
 vector<double> *Gen_status;
 Gen_BmotherIndex = 0;
 Gen_status = 0;
 TBranch* b_Gen_BmotherIndex = 0;
 TBranch* b_Gen_status = 0;
 tree->SetBranchAddress("Gen_BmotherIndex",&Gen_BmotherIndex,&b_Gen_BmotherIndex);
 tree->SetBranchAddress("Gen_status",&Gen_status,&b_Gen_status);
 Long64_t tentry = tree->LoadTree(entry);
 b_Gen_pdg_id->GetEntry(tentry);
 b_Gen_motherpdg_id->GetEntry(tentry);
 b_Gen_BmotherIndex->GetEntry(tentry);
 b_Gen_status->GetEntry(tentry);
 b_Gen_pt->GetEntry(tentry);
 b_Gen_eta->GetEntry(tentry);
 b_Gen_phi->GetEntry(tentry);
 b_Gen_energy->GetEntry(tentry);
 //Loop over gen candidates to look for matching
 for(uint gp=0; gp<Gen_pdg_id->size(); gp++){
   //Take a tau
  if(!(abs(Gen_pdg_id->at(gp))==15 && Gen_status->at(gp)==2)) continue;
  gen_cand.SetPtEtaPhiE(Gen_pt->at(gp),Gen_eta->at(gp),Gen_phi->at(gp),Gen_energy->at(gp));
  TLorentzVector gen_neu(0.,0.,0.,0.);
  //Look tau decays
  int etau    = 0;
  int mutau   = 0;
  int tauhtau = 0;
  for(int gpd=0; gpd<int(Gen_pdg_id->size()); gpd++){
   if(!(Gen_BmotherIndex->at(gpd)==gp)) continue;
   if(abs(Gen_pdg_id->at(gpd))==11 || abs(Gen_pdg_id->at(gpd))==12 || abs(Gen_pdg_id->at(gpd))==16) etau    += 1;
   if(abs(Gen_pdg_id->at(gpd))==13 || abs(Gen_pdg_id->at(gpd))==14 || abs(Gen_pdg_id->at(gpd))==16) mutau   += 1;
   if(abs(Gen_pdg_id->at(gpd))==12 || abs(Gen_pdg_id->at(gpd))==14 || abs(Gen_pdg_id->at(gpd))==16){
    tauhtau += 1;
    gen_neu.SetPtEtaPhiE(Gen_pt->at(gpd),Gen_eta->at(gpd),Gen_phi->at(gpd),Gen_energy->at(gpd));
    gen_cand = gen_cand - gen_neu;
   }
  }
  if(VectorUtil::DeltaR(reco_cand,gen_cand)<0.3){
   if((candid==11 && etau==3) || (candid==13 && mutau==3) || (candid==15 && tauhtau==1)){
    matched = true;
    break;
   }
  }
 }
 pair<bool,TLorentzVector> thepair = make_pair(matched,gen_cand);
 return thepair;
}




/////
//   Initialize the strings 
/////
void initialize_strings(vector<string> &varname, vector<string> &varcut){
 varname[0] = "Gen_tau_h_pair_exists";
 varname[1] = "At least 2 reco tau exist\t"; varname[2] = "diTau gen matching\t";
 varname[3] = "Tau1 Pt > 45 \t\t"; varname[4] = "Tau1 |eta| < 2.1 \t\t";
 varname[5] = "Tau2 Pt > 45 \t\t\t"; varname[6] = "Tau2 |eta| < 2.1 \t\t";
 varname[7] = "Tau1 DMF new DMs\t"; varname[8] = "Tau1 dz < 0.2 \t\t\t";
 varname[9] = "Tau2 DMF new DMs\t"; varname[10] = "Tau2 dz < 0.2 \t\t\t";
 varname[11] = "TauTau DeltaR > 0.5 \t\t"; varname[12] = "Tau1 anti-e MVA5\t";
 varname[13] = "Tau1 anti-mu loose\t"; varname[14] = "Tau2 anti-e MVA5\t";
 varname[15] = "Tau2 anti-mu loose\t"; varname[16] = "Tau1 Iso DB 3 hits\t";
 varname[17] = "Tau2 Iso DB 3 hits\t";
}
/////
//   Measure efficiencies of cuts
/////
void measure_efficiencies(int nentries, double cut_eff[], double weight, int pos, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]){
 double eff_Rel, err_eff_Rel, eff_Cum, err_eff_Cum;
 for(int i=0; i<numtot_cuts; i++){
  //Rel Eff
  if(i==0){
   eff_Rel = (double) cut_eff[i]/nentries;
   err_eff_Rel = sqrt(eff_Rel*(1-eff_Rel)/(double)nentries);
  }else{
   eff_Rel = (double) cut_eff[i]/cut_eff[i-1];
   err_eff_Rel = sqrt(eff_Rel*(1-eff_Rel)/(double)cut_eff[i-1]);
  }
  //Cum Eff
  eff_Cum = (double)cut_eff[i]/nentries;
  err_eff_Cum = sqrt(eff_Cum*(1-eff_Cum)/(double)nentries);
  //Take values
  if(relEff){
   sample_cuts[pos][i]     = eff_Rel; 
   sample_cuts_err[pos][i] = err_eff_Rel; 
  }
  if(cumEff){
   sample_cuts[pos][i]     = eff_Cum;
   sample_cuts_err[pos][i] = err_eff_Cum;
  }
  if(rateEvt){
   if(noLumiNorm) weight = 1.;
   sample_cuts[pos][i]     = cut_eff[i]*weight;
   double err_rateEvt      = sqrt(cut_eff[i])*weight;
   sample_cuts_err[pos][i] = err_rateEvt;
  }
 }
}
/////
//   Print table of rel eff for all samples
/////
void print_table(vector<string> &rootplas, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]){
 //Table ini part
 //\def\tablepagesize{\fontsize{9pt}{8pt}\selectfont}
 cout<<"{\\tablepagesize"<<endl;
 cout<<"\\begin{table}[htdp]"<<endl;
 cout<<"\\begin{center}"<<endl;
 cout<<"\\begin{tabular}{|c|";
 for(uint i=0; i<rootplas.size(); i++){
  cout<<"c|";
  if(i==rootplas.size()-1) cout<<"}"<<endl;
 }
 cout<<"\\hline"<<endl;
 cout<<" Selection ";
 for(uint i=0; i<rootplas.size(); i++){
  cout<<" & "<<rootplas[i]<<"";
  if(i==rootplas.size()-1) cout<<"\\\\"<<endl;
 }
 //Table body
 cout<<setiosflags(ios::fixed)<<setprecision(SETPRECISION); 
 //Define strings to print
 vector<string> varname(numtot_cuts);
 vector<string> varcut(numtot_cuts);
 initialize_strings(varname,varcut);
 for(int j=0; j<numtot_cuts; j++) for(uint i=0; i<rootplas.size(); i++){
  if(i==0){
   cout<<"\\hline"<<endl;
   if(!rateEvt){
cout<<""<<varname[j].c_str()<<varcut[j].c_str()<<" & "<<sample_cuts[i][j]*100<<"$\\pm$"<<sample_cuts_err[i][j]*100<<"";
   }else{
cout<<""<<varname[j].c_str()<<varcut[j].c_str()<<" & "<<sample_cuts[i][j]<<"$\\pm$"<<sample_cuts_err[i][j]<<"";
   }
  }else{
   if(!rateEvt){
    cout<<" & "<<sample_cuts[i][j]*100<<"$\\pm$"<<sample_cuts_err[i][j]*100<<""; 
   }else{
    cout<<" & "<<sample_cuts[i][j]<<"$\\pm$"<<sample_cuts_err[i][j]<<"";
   }
  }
  if(i==rootplas.size()-1) cout<<"\\\\"<<endl;
 } 
 //Table end part
 cout<<"\\hline"<<endl;
 cout<<"\\end{tabular}"<<endl;
 cout<<"\\end{center}"<<endl;
 cout<<"\\end{table}"<<endl;
 cout<<"}"<<endl;
}
