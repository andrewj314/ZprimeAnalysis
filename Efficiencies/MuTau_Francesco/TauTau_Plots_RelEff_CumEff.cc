/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects and events 
   It draws a plot following the cutFlow 

Need to specify
1. See Declare constants
*/
/////
//   To run: root -l TauTau_Plots_RelEff_CumEff.cc+
/////
/////
//   Prepare Root and Roofit
/////
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
/////
//   Declare constants
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
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
TH1F* cutflow_eff(TTree* tree, int nentries, int pos, string rootpla);
bool genHadronicTauPairExists(TTree* tree, int nentry);
bool genMuTauPairExists(TTree* tree, int nentry);

pair<bool, TLorentzVector> matchTauToGen(TLorentzVector p4Vector, TTree* tree, int nentry);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
TH1F* measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight, int pos, string rootplas);
void setTDRStyle();
//Path and root file
const string path           = "/uscms_data/d3/andrewj/CMSSW_7_4_7/src/Efficiencies/MuTau_Francesco/";//"/home/francescoromeovb/Shared_win_ubu_/Work/RunII_HighMassDiTau/Efficiency/";
const char *samples[]       = {"Zprime"};// {"SusyHTT"};
//Corrections
const double Luminosity     = 19600; //pb^-1
const bool noLumiNorm       = true; //true means NO luminosity normalization done
const bool noPUcorr         = true; //true means NO PU corr done
//Choose among relEff, cumEff, or rateEvt 
const bool relEff           = true;
const bool cumEff           = false;
//Acceptance and object selection
const double cand1_acc[2]   = {45, 2.1};
const double cand2_acc[2]   = {45, 2.1};
const int numtot_cuts       = 18;
//Signal selection
//charge==-1 && cosDphi<-0.95 && met>=30 && pZetaMt-3.051*pZetaVisMt>-50 && cosDphiLMet<0.2 && nbjet==0
//const int sr_cuts           = 6;
//double sigreg_cuts[sr_cuts] = {-1,-0.95,20,-50,0.2,0};
const double SETPRECISION = 3;
/////
//   Main function
/////
void TauTau_Plots_RelEff_CumEff() {
 setTDRStyle();
 //For Rel Eff
 TCanvas* c1 = new TCanvas("Eff_rel","Eff_rel",100,100,1300,900); 
 c1->cd();
 TLegend *leg = new TLegend(0.2, 0.3, 0.4, 0.75);
 leg->SetHeader("Samples");
 leg->SetBorderSize(0);
 string nameFile = "RelEff_";
 //Run over all samples 
 vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
 for(uint i=0; i<rootplas.size(); i++){
  //Call tree  
  cout << "call TFile\n";
  TFile* f = Call_TFile(rootplas[i]);
  TTree* tree; f->GetObject("TNT/BOOM;1",tree);
  int nentries = tree->GetEntries(); 
  //Take plot
  cout << "Instantiate TH1F\n";
  TH1F* hRelEff = new TH1F("test",rootplas[i].c_str(),numtot_cuts,0,numtot_cuts);
  cout << "Cut Flow\n";
  hRelEff = cutflow_eff(tree,nentries,i,rootplas[i]);
  cout << "Cut Flow passed\n";

  if(i==0){
   hRelEff->Draw("PE1");
  }else{
   hRelEff->Draw("PE1same");
  }
  leg->AddEntry(hRelEff,rootplas[i].c_str(),"LP");
  if(i+1!=rootplas.size()){nameFile = nameFile+rootplas[i]+"_";}else{nameFile = nameFile+rootplas[i];} 
 }
 //leg->Draw(); 
 //nameFile = nameFile+".pdf";
 if(relEff) nameFile = "RelEffAll.pdf";
 if(cumEff) nameFile = "CumEffAll.pdf";
 c1->SaveAs(nameFile.c_str());

 
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla){
 string file_name = path+rootpla+".root";
 cout<<"File_name "<<file_name<<endl;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Efficiency for the object
/////
TH1F* cutflow_eff(TTree* tree, int nentries, int pos, string rootpla){
 cout << "Define TH1\n"; 
 TH1F *hCurrRelEff = new TH1F(rootpla.c_str(),rootpla.c_str(),numtot_cuts,0,numtot_cuts);
 cout << "TH1 defined\n";
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
// tree->SetBranchAddress("Tau_leadChargedCandDz_pv",&Tau_leadChargedCandDz_pv,&b_Tau_leadChargedCandDz_pv);
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
  //b_Tau_leadChargedCandDz_pv->GetEntry(tentry);
  b_Tau_decayModeFindingNewDMs->GetEntry(tentry);
  b_Tau_againstElectronMVALooseMVA5->GetEntry(tentry);
  b_Tau_againstMuonLoose3->GetEntry(tentry);
  b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->GetEntry(tentry);
  //get values
  double wgt_pu   = 1.;
  double wgt_lumi = 1.;
  if(noPUcorr) wgt_pu = 1.;

  if(genHadronicTauPairExists(tree, i)){
    cut_eff[0] += wgt_pu;
    if(Tau_pt->size() >= 2){
      cut_eff[1] += wgt_pu;
      bool twoMatchedTaus = false;
      for(uint j = 0; j < Tau_pt->size(); j++){
        for(uint k = j+1; k < Tau_pt->size(); k++){
          TLorentzVector Tau1Vector;
          Tau1Vector.SetPtEtaPhiE(Tau_pt->at(j), Tau_eta->at(j), Tau_phi->at(j), Tau_energy->at(j));
          TLorentzVector Tau2Vector;
          Tau2Vector.SetPtEtaPhiE(Tau_pt->at(k), Tau_eta->at(k), Tau_phi->at(k), Tau_energy->at(k));
          if(matchTauToGen(Tau1Vector, tree, i).first && matchTauToGen(Tau2Vector, tree, i).first){
            twoMatchedTaus = true;
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
                  if(/*Tau_leadChargedCandDz_pv->at(0) < 0.2*/true){
                    cut_eff[8] += wgt_pu;
                    if(Tau_decayModeFindingNewDMs->at(1) > 0.5){
                      cut_eff[9] += wgt_pu;
                      if(/*Tau_leadChargedCandDz_pv->at(1) < 0.2*/true){
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
 hCurrRelEff = measure_efficiencies(nentries, cut_eff, varname, varcut, weight, pos, rootpla);
 return hCurrRelEff;
}



bool genMuTauPairExists(TTree* tree, int nentry){

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

  bool genMuTauPair = false;
  bool IsItAHadronicDecay; // boolean used to specify whether a gen tau lepton decays hadronically
  bool IsItAMuonicDecay; // boolean used to specify whether a gen tau decays into a muon
  vector<int> tempTauIndexVector; // vector which contains the index (in the gen collection) of the gen level taus (before decaying)
  tempTauIndexVector.clear(); // clear any previous declaration from memory
  vector<bool> IsItAHadronicDecayVector; // vector of booleans which contain the information about whether a gen tau lepton decays hadronically
  IsItAHadronicDecayVector.clear(); // clear any previous declaration from memory
  vector<bool> IsItAMuonicDecayVector; //vector of booleans which contains the information about whether a gen tau decays into a muon
  IsItAMuonicDecayVector.clear();
  
  //---Loop over gen particles to find the tau neutrinos and then store the index of each tau neutrino's mother (a tau).
  for(uint jj = 0; jj < Gen_pt->size(); jj++) {
    if(Gen_BmotherIndex->at(jj) >= 0 && (abs(Gen_pdg_id->at(jj)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(jj))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(jj)) == 2) ) {
      tempTauIndexVector.push_back(Gen_BmotherIndex->at(jj));
    }
  }


  //---Perform matching only if there is at least one gen tau in the event
  if(tempTauIndexVector.size() > 1) {
    //---Loop over the gen taus and determine whether it decays hadronically.
    for(uint jjj = 0; jjj < tempTauIndexVector.size(); jjj++) {
      IsItAHadronicDecay = true;
      IsItAMuonicDecay = false;
      for(uint jj = 0; jj < Gen_pt->size(); jj++) {
        if(Gen_BmotherIndex->at(jj) >= 0 &&  ((abs(Gen_pdg_id->at(jj)) == 12) || (abs(Gen_pdg_id->at(jj)) == 14)) && (Gen_BmotherIndex->at(jj) == tempTauIndexVector.at(jjj)) ) {
          IsItAHadronicDecay = false; // it is not a hadronic tau decay since it decayed to a electron/muon neutrino
        }
        if(Gen_BmotherIndex->at(jj) >= 0 && abs(Gen_pdg_id->at(jj)) == 13 && Gen_BmotherIndex->at(jj) == tempTauIndexVector.at(jjj)){
	  IsItAMuonicDecay = true; // it is a muonic tau decay since it decayed into a muon
	}
      }
      IsItAHadronicDecayVector.push_back(IsItAHadronicDecay);
      IsItAMuonicDecayVector.push_back(IsItAMuonicDecay);
    }
    
    for(uint i = 0; i < tempTauIndexVector.size(); i++){
      for(uint j = i+1; j < tempTauIndexVector.size(); j++){
        if(((IsItAHadronicDecayVector.at(i) && IsItAMuonicDecayVector.at(j)) || (IsItAMuonicDecayVector.at(i) && IsItAHadronicDecayVector.at(j))) && abs(Gen_pdg_id->at(Gen_BmotherIndex->at(tempTauIndexVector.at(i)))) == abs(Gen_pdg_id->at(Gen_BmotherIndex->at(tempTauIndexVector.at(j))))){
	  genMuTauPair = true;
	} 
      }
    }
    
  } 

  return genMuTauPair;

}



bool genHadronicTauPairExists(TTree* tree, int nentry){

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
  vector<int> tempTauIndexVector; // vector which contains the index (in the gen collection) of the gen level taus (before decaying)
  tempTauIndexVector.clear(); // clear any previous declaration from memory
  vector<bool> IsItAHadronicDecayVector; // vector of booleans which contain the information about whether a gen tau lepton decays hadronically
  IsItAHadronicDecayVector.clear(); // clear any previous declaration from memory

  
  //---Loop over gen particles to find the tau neutrinos and then store the index of each tau neutrino's mother (a tau).
  for(uint jj = 0; jj < Gen_pt->size(); jj++) {
    if(Gen_BmotherIndex->at(jj) >= 0 && (abs(Gen_pdg_id->at(jj)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(jj))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(jj)) == 2) ) {
      tempTauIndexVector.push_back(Gen_BmotherIndex->at(jj));
    }
  }


  //---Perform matching only if there is at least one gen tau in the event
  if(tempTauIndexVector.size() > 1) {
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
        if(IsItAHadronicDecayVector.at(i) && IsItAHadronicDecayVector.at(j) && abs(Gen_pdg_id->at(Gen_BmotherIndex->at(tempTauIndexVector.at(i)))) == abs(Gen_pdg_id->at(Gen_BmotherIndex->at(tempTauIndexVector.at(j))))){
	  genHadTauPair = true;
	} 
      }
    }
    
  } 

  return genHadTauPair;

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
  vector<bool> IsItAHadronicDecayVector; // vector of booleans which contain the information about whether a gen tau lepton decays hadronically
  IsItAHadronicDecayVector.clear(); // clear any previous declaration from memory
  vector<uint> tempTauIndexVector; // vector which contains the index (in the gen collection) of the gen level taus (before decaying)
  tempTauIndexVector.clear(); // clear any previous declaration from memory
  vector<TLorentzVector> tempNeutrinoMomentumVector; // vector of lorentz 4-momentum vectors for each tau neutrino from the tau decay
  tempNeutrinoMomentumVector.clear(); // clear any previous declaration from memory

  //---Loop over gen particles to find the tau neutrinos and then store the index of each tau neutrino's mother (a tau).
  //---Also store the tau neutrino's 4-momentum vector in order to calculate the visible tau 4-momentum at a later point.
  for(uint jj = 0; jj < Gen_pt->size(); jj++) {
    if( (abs(Gen_pdg_id->at(jj)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(jj))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(jj)) == 2) ) {
      tempTauIndexVector.push_back(Gen_BmotherIndex->at(jj));
      theNeutrinoObject.SetPtEtaPhiE(Gen_pt->at(jj), Gen_eta->at(jj), Gen_phi->at(jj), Gen_energy->at(jj));
      tempNeutrinoMomentumVector.push_back(theNeutrinoObject);
    }
  }

  //---Perform matching only if there is at least one gen tau in the event
  if(tempTauIndexVector.size() > 0) {
    //---Loop over the gen taus and determine whether it decays hadronically.
    for(uint jjj = 0; jjj < tempTauIndexVector.size(); jjj++) {
      IsItAHadronicDecay = true;
      for(uint jj = 0; jj < Gen_pt->size(); jj++) {
        if( ((abs(Gen_pdg_id->at(jj)) == 12) || (abs(Gen_pdg_id->at(jj)) == 14)) && (Gen_BmotherIndex->at(jj) == tempTauIndexVector.at(jjj)) ) {
          IsItAHadronicDecay = false; // it is not a hadronic tau decay since it decayed to a electron/muon neutrino
        }
      }
      IsItAHadronicDecayVector.push_back(IsItAHadronicDecay);
    }
    //---Loop over the gen taus and calculate the 4-momentum of the visible products (i.e. subtract the 4-momentum of the tau neutrino)
    for(uint jjj = 0; jjj < tempTauIndexVector.size(); jjj++) {
      for(uint jj = 0; jj < Gen_pt->size(); jj++) {
        if(jj == tempTauIndexVector.at(jjj)) {
          theGenObject.SetPtEtaPhiE(Gen_pt->at(jj), Gen_eta->at(jj), Gen_phi->at(jj), Gen_energy->at(jj)); // 4-momentum of the gen tau
          theGenObject = theGenObject - tempNeutrinoMomentumVector.at(jjj); // subtract the 4-momentum of the tau neutrino (visible tau)
          cout << "DeltaR = " << p4Vector.DeltaR(theGenObject) << "\n";
          if( (IsItAHadronicDecayVector.at(jjj)) && (p4Vector.DeltaR(theGenObject) <= 0.5) ) {isGenMatched = true;}
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


/////
//   Initialize the strings 
/////
void initialize_strings(vector<string> &varname, vector<string> &varcut){
 varname[0] = "Gen_tau_h_pair_exists";
 varname[1] = "Reco_tau_h_pair_exists\t"; varname[2] = "diTau gen matching\t";
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
TH1F* measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight, int pos, string rootpla){
 TH1F *hCurrRelEff = new TH1F(rootpla.c_str(),rootpla.c_str(),numtot_cuts,0,numtot_cuts);
 hCurrRelEff->SetTitle(0);
 if(relEff) hCurrRelEff->GetYaxis()->SetTitle("Relative Efficiency");
 if(cumEff) hCurrRelEff->GetYaxis()->SetTitle("Cumulative Efficiency");
 hCurrRelEff->LabelsOption("v");
 hCurrRelEff->SetMaximum(1.1);
 hCurrRelEff->SetMinimum(0);
 if(rootpla=="data"){
  hCurrRelEff->SetMarkerStyle(20);
  hCurrRelEff->SetLineColor(1);
  hCurrRelEff->SetMarkerColor(1);
 }else if(rootpla=="SSM1000"){
  hCurrRelEff->SetMarkerStyle(21);
  hCurrRelEff->SetLineColor(1);
  hCurrRelEff->SetMarkerColor(1);
 }else{
  hCurrRelEff->SetMarkerStyle(pos+20);
  hCurrRelEff->SetLineColor(pos+2);
  hCurrRelEff->SetMarkerColor(pos+2);
 }
 double eff_Rel, err_eff_Rel, eff_Cum, err_eff_Cum;
 cout<<setiosflags(ios::fixed)<<setprecision(2);
 cout<<"Dataset is: "<<rootpla<<endl;
 for(int i=0; i<numtot_cuts; i++){
  if(i==0){
   eff_Rel = (double) cut_eff[i]/nentries;
   err_eff_Rel = sqrt(eff_Rel*(1-eff_Rel)/(double)nentries);
  }else{
   eff_Rel = (double) cut_eff[i]/cut_eff[i-1];
   err_eff_Rel = sqrt(eff_Rel*(1-eff_Rel)/(double)cut_eff[i-1]);
  }
  eff_Cum = (double)cut_eff[i]/nentries;
  err_eff_Cum = sqrt(eff_Cum*(1-eff_Cum)/(double)nentries);
  cout<<"\\hline"<<endl;
  if(noLumiNorm){
   cout<<"{\\cellcolor{white!90}"<<varname[i].c_str()<<varcut[i].c_str()<<"} & {\\cellcolor{white!90}"<<cut_eff[i]<<"} & {\\cellcolor{white!90}"<<eff_Rel*100<<"$\\pm$"<<err_eff_Rel*100<<"} & {\\cellcolor{white!90}"<<eff_Cum*100<<"$\\pm$"<<err_eff_Cum*100<<"}\\\\"<<endl;
  }else{
   cout<<"{\\cellcolor{white!90}"<<varname[i].c_str()<<varcut[i].c_str()<<"} & {\\cellcolor{white!90}"<<cut_eff[i]*weight<<"} & {\\cellcolor{white!90}"<<eff_Rel*100<<"$\\pm$"<<err_eff_Rel*100<<"} & {\\cellcolor{white!90}"<<eff_Cum*100<<"$\\pm$"<<err_eff_Cum*100<<"}\\\\"<<endl;
  }
  if(relEff){
   hCurrRelEff->SetBinContent(i+1,eff_Rel);
   hCurrRelEff->SetBinError(i+1,err_eff_Rel);
  }
  if(cumEff){
   hCurrRelEff->SetBinContent(i+1,eff_Cum);
   hCurrRelEff->SetBinError(i+1,err_eff_Cum);
  }
  string binLabel = varname[i]+varcut[i];
  hCurrRelEff->GetXaxis()->SetBinLabel(i+1,binLabel.c_str());
  hCurrRelEff->LabelsOption("v");
 }
 return hCurrRelEff;
}
/////
//   Set setTDRStyle_modified (from link https://twiki.cern.ch/twiki/pub/CMS/TRK10001/setTDRStyle_modified.C)
/////
void setTDRStyle(){
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);
  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(0);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

//  tdrStyle->SetEndErrorSize(0);
  tdrStyle->SetErrorX(0.5);
//  tdrStyle->SetErrorMarker(20);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1111);
  tdrStyle->SetFitFormat("5.4g");
  //tdrStyle->SetFuncColor(1);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.025);
  tdrStyle->SetPadBottomMargin(0.275);
  tdrStyle->SetPadLeftMargin(0.09);
  tdrStyle->SetPadRightMargin(0.025);

  // For the Global title:

  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.7);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);
  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}
