/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects  
   It makes plots for rel eff following cut flow 
   For event selection you can choose between rel eff plots or N-1 eff (only integral or vs event quantities)
   Special eff plots can also be done  

Need to specify
1. See Declare constants

To do
- //This part depends on how you define the cut flow
  //v1 = Mu Pt, v2 = Mu eta
  ... May improve it 
Notes
Depending on what you are studying (objects, sigreg, speceff)
you must pay attention of what you plot (pt,eta,phi,massVis,nvert).
E.g. you can plot massVis if you do not have 2 candidates
*/
/////
//   To run: root -l TauTau_Plots_EffpT_Effeta_Effpu.cc+ 
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
#include <TMath.h>
#include <Math/SVector.h>
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/normalizedPhi.h"
#include <TLorentzVector.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <TROOT.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TEnv.h>
#include <TChain.h>


using namespace ROOT::Math;
/////
//   Declare constants
/////
//Path and root file
const string path       = "/uscms_data/d3/andrewj/CMSSW_7_4_7/src/ZprimeSamples/ZprimeToTauTau_M_2000/";
const string rootpla    = "Zprime.root";
//Corrections
const double Luminosity = 19600; //pb^-1
const bool noLumiNorm   = true; //true means NO luminosity normalization done
const bool noPUcorr     = true; //true means NO PU corr done
//Plots
const bool save_plots = true; 
//Choose gen pt,eta or reco pt,eta
const bool gen_pteta  = false;
//Choose Pt aut Eta aut NVert and comments constants for others
//Pt
const bool pt_plots   = true;
const int const_bin = 10; const double const_from = 0; const double const_to = 1000;
//const string titleXaxis = "Gen_pT (GeV/c)"; const int lmcol = 4; //Blue
const string titleXaxis = "Reco_pT (GeV/c)"; const int lmcol = 2; //Red
const string part = "recopT_"; //genpT or recopT
//Eta
const bool eta_plots = false;
//const int const_bin  = 31; const double const_from = -3.1; const double const_to = 3.1;
//const string titleXaxis = "Gen_#eta"; const int lmcol = 4; //Red
//const string part = "geneta_";
//Phi
const bool phi_plots = false;
//const int const_bin  = 33; const double const_from = -3.3; const double const_to = 3.3;
//const string titleXaxis = "Gen_#phi"; const int lmcol = 4; //Blue
//const string part = "genphi_";
//NVert
const bool nvert_plots = false;
//const int const_bin = 50; const double const_from = 0; const double const_to = 50;
//const string titleXaxis = "Vertices"; const int lmcol = 4; //Black
//const string part = "nvert_";
//MassVis
const bool massvis_plots = false;
//const int const_bin = 20; const double const_from = 0; const double const_to = 1500;
//const string titleXaxis = "massVis(e,#mu,E_{T}^{miss}) (GeV/c^{2})"; const int lmcol = 4; //Black
//const string part = "massvis_";
//Acceptance and object selection
const double cand1_acc[2]   = {45, 2.1};
const double cand2_acc[2]   = {45, 2.1};
const int numtot_cuts       = 22;
//Signal selection
//charge==-1 && cosDphi<-0.95 && met>=30 && pZetaMt-3.051*pZetaVisMt>-50 && cosDphiLMet<0.2 && nbjet==0
//const int sr_cuts           = 6;
//double sigreg_cuts[sr_cuts] = {-1,-0.95,20,-50,0.2,0};
/////
//   Declare functions 
/////
TFile* Call_TFile();
void cutflow_eff_obj(TTree* tree, int nentries);
//void cutflow_eff_sr(TTree* tree, int nentries);
//void cutflow_eff_spec(TTree* tree, int nentries);
pair<bool, TLorentzVector> matchTauToGen(TLorentzVector p4Vector, TTree* tree, int nentry);
pair<bool, pair<TLorentzVector,TLorentzVector>> genHadronicTauPairExists(TTree* tree, int nentry);
double visibleMass(TLorentzVector obj1, TLorentzVector obj2);
double visiblePlusMETMass(TLorentzVector obj1, TLorentzVector obj2, TLorentzVector METVector);
double visiblePlusDeltaPtMass(TLorentzVector obj1, TLorentzVector obj2);
double svFitMass(TLorentzVector obj1, TLorentzVector obj2, TLorentzVector METVector);
double svFitDeltaPtMass(TLorentzVector obj1, TLorentzVector obj2);
double pZeta(TLorentzVector obj1, TLorentzVector obj2, TLorentzVector METVector);
double pZetaVis(TLorentzVector obj1, TLorentzVector obj2);
int nBTags(TTree *tree, int nentry);
bool passBJetCuts(int nobj);

void draw_eff(TTree* tree, int nentries, vector<vector<double> > cut_eff);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
void hist_ratio(TH1* hpn, TH1* hpd, string varname, string varcut, int num);
void setTDRStyle();

TH1* theFirstTauPt = new TH1F("theFirstTauPt", "pT of the first gen vis tau", const_bin, const_from, const_to);
TH1* theFirstMatchedTauPt = new TH1F("theFirstMatchedTauPt", "pT of the first matched gen vis tau", const_bin, const_from, const_to);

TH1F* theSecondTauPt = new TH1F("theSecondTauPt", "pT of the second gen vis tau", 100, 0, 400);
TH1F* theSecondMatchedTauPt = new TH1F("theSecondMatchedTauPt", "pT of the second matched gen vis tau", 100, 0, 400);

TH1I* tauMultiplicity = new TH1I("tauMultiplicity", "Number of reco taus per event", 10, 0, 10);


/////
//   Main function
/////
void TauTau_Plots_EffpT_Effeta_Effpu(){
 setTDRStyle();
 //Call tree 
 TFile* f = Call_TFile();
 TTree* tree; 
 f->GetObject("BOOM",tree);
 int nentries = tree->GetEntries(); 
 cout << "nentries = " << nentries;
 cutflow_eff_obj(tree, nentries);
 //if(sigreg)  cutflow_eff_sr(tree, nentries);
 //if(speceff) cutflow_eff_spec(tree, nentries);
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
void cutflow_eff_obj(TTree* tree, int nentries){ 
 /////
 //   Get variables
 /////



 TH1* hpn = new TH1F("hpn","hpn",const_bin,const_from,const_to);
 hpn->Sumw2();
 TH1* hpd = new TH1F("hpd","hpd",const_bin,const_from,const_to);
 hpd->Sumw2();

 TH1* DMF_n  = new TH1F("DMF_n","Tau DMF",const_bin,const_from,const_to);
 DMF_n->Sumw2();
 TH1* DMF_d = new TH1F("DMF_d","hpd",const_bin,const_from,const_to);
 DMF_d->Sumw2();

 TH1* DMFNewDMs_n  = new TH1F("DMFNewDMs_n","Tau DMF",const_bin,const_from,const_to);
 DMFNewDMs_n->Sumw2();
 TH1* DMFNewDMs_d = new TH1F("DMFNewDMs_d","hpd",const_bin,const_from,const_to);
 DMFNewDMs_d->Sumw2();

 TH1* againstMuonLoose3_n = new TH1F("againstMuonLoose3_n","againstMuonLoose3_n",const_bin,const_from,const_to);
 againstMuonLoose3_n->Sumw2();
 TH1* againstMuonLoose3_d = new TH1F("againstMuonLoose3_d","againstMuonLoose3_d",const_bin,const_from,const_to);
 againstMuonLoose3_d->Sumw2();

 TH1* againstMuonTight3_n = new TH1F("againstMuonTight3_n","againstMuonTight3_n",const_bin,const_from,const_to);
 againstMuonTight3_n->Sumw2();
 TH1* againstMuonTight3_d = new TH1F("againstMuonTight3_d","againstMuonTight3_d",const_bin,const_from,const_to);
 againstMuonTight3_d->Sumw2();

 TH1* againstElectronLoose_n = new TH1F("againstElectronLoose_n","againstElectronLoose_n",const_bin,const_from,const_to);
 againstElectronLoose_n->Sumw2();
 TH1* againstElectronLoose_d = new TH1F("againstElectronLoose_d","againstElectronLoose_d",const_bin,const_from,const_to);
 againstElectronLoose_d->Sumw2();

 TH1* againstElectronMedium_n = new TH1F("againstElectronMedium_n","againstElectronMedium_n",const_bin,const_from,const_to);
 againstElectronMedium_n->Sumw2();
 TH1* againstElectronMedium_d = new TH1F("againstElectronMedium_d","againstElectronMedium_d",const_bin,const_from,const_to);
 againstElectronMedium_d->Sumw2();

 TH1* looseCombinedIsoDB3Hits_n = new TH1F("looseCombinedIsoDB3Hits_n","looseCombinedIsoDB3Hits_n",const_bin,const_from,const_to);
 looseCombinedIsoDB3Hits_n->Sumw2();
 TH1* looseCombinedIsoDB3Hits_d = new TH1F("looseCombinedIsoDB3Hits_d","looseCombinedIsoDB3Hits_d",const_bin,const_from,const_to);
 looseCombinedIsoDB3Hits_d->Sumw2();

 TH1* mediumCombinedIsoDB3Hits_n = new TH1F("mediumCombinedIsoDB3Hits_n","mediumCombinedIsoDB3Hits_n",const_bin,const_from,const_to);
 mediumCombinedIsoDB3Hits_n->Sumw2();
 TH1* mediumCombinedIsoDB3Hits_d = new TH1F("mediumCombinedIsoDB3Hits_d","mediumCombinedIsoDB3Hits_d",const_bin,const_from,const_to);
 mediumCombinedIsoDB3Hits_d->Sumw2();

 TH1* tightCombinedIsoDB3Hits_n = new TH1F("tightCombinedIsoDB3Hits_n","tightCombinedIsoDB3Hits_n",const_bin,const_from,const_to);
 tightCombinedIsoDB3Hits_n->Sumw2();
 TH1* tightCombinedIsoDB3Hits_d = new TH1F("tightCombinedIsoDB3Hits_d","tightCombinedIsoDB3Hits_d",const_bin,const_from,const_to);
 tightCombinedIsoDB3Hits_d->Sumw2();

 

 TH1D* visMassHisto = new TH1D("visMassHisto", "visible Mass", 50, 0, 4000);
 TH1D* visPlusMetMassHisto = new TH1D("visPlusMetMassHisto", "visible +MET Mass", 50, 0, 4000);
 TH1D* visPlusDeltaPtMassHisto = new TH1D("visPlusDeltaPtMassHisto", "visible + #Delta p_{T} Mass", 50, 0, 4000);


 //Tau
 vector<double> *Tau_pt, *Tau_eta, *Tau_phi, *Tau_energy, *Tau_leadChargedCandDz_pv, *Tau_leadChargedCandDxy_bs, *Tau_leadChargedCandDxyError;
 vector<int> *Tau_decayModeFinding, *Tau_decayModeFindingNewDMs, *Tau_againstElectronMVALooseMVA5, *Tau_againstElectronMVAMediumMVA5;
 vector<int> *Tau_againstMuonLoose3, *Tau_againstMuonTight3, *Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits; 
 vector<int> *Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, *Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
 vector<int> *Tau_nProngs;
 vector<int> *Tau_charge;
 Tau_pt = 0;
 Tau_eta = 0;
 Tau_phi = 0;
 Tau_energy = 0;
 Tau_leadChargedCandDz_pv = 0;
 Tau_leadChargedCandDxy_bs = 0;
 Tau_leadChargedCandDxyError = 0;
 Tau_decayModeFinding = 0;
 Tau_decayModeFindingNewDMs = 0;
 Tau_againstElectronMVALooseMVA5 = 0;
 Tau_againstElectronMVAMediumMVA5 = 0;
 Tau_againstMuonLoose3 = 0;
 Tau_againstMuonTight3 = 0;
 Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
 Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
 Tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
 Tau_nProngs = 0;
 Tau_decayModeFinding = 0;
 Tau_charge = 0;
 TBranch *b_Tau_pt = 0;
 TBranch *b_Tau_eta = 0;
 TBranch *b_Tau_phi = 0;
 TBranch *b_Tau_energy = 0;
 TBranch *b_Tau_leadChargedCandDz_pv = 0;
 TBranch *b_Tau_leadChargedCandDxy_bs = 0;
 TBranch *b_Tau_leadChargedCandDxyError = 0;
 TBranch *b_Tau_decayModeFindingNewDMs = 0;
 TBranch *b_Tau_againstElectronMVALooseMVA5 = 0;
 TBranch *b_Tau_againstElectronMVAMediumMVA5 = 0;
 TBranch *b_Tau_againstMuonLoose3 = 0;
 TBranch *b_Tau_againstMuonTight3 = 0;
 TBranch *b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0; 
 TBranch *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0; 
 TBranch *b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0; 
 TBranch *b_Tau_nProngs = 0;
 TBranch *b_Tau_decayModeFinding = 0;
 TBranch *b_Tau_charge = 0;
 tree->SetBranchAddress("Tau_pt",&Tau_pt,&b_Tau_pt);
 tree->SetBranchAddress("Tau_eta",&Tau_eta,&b_Tau_eta);
 tree->SetBranchAddress("Tau_phi",&Tau_phi,&b_Tau_phi);
 tree->SetBranchAddress("Tau_energy",&Tau_energy,&b_Tau_energy);
 tree->SetBranchAddress("Tau_leadChargedCandDz_pv",&Tau_leadChargedCandDz_pv,&b_Tau_leadChargedCandDz_pv);
 tree->SetBranchAddress("Tau_leadChargedCandDxy_bs",&Tau_leadChargedCandDxy_bs,&b_Tau_leadChargedCandDxy_bs);
 tree->SetBranchAddress("Tau_leadChargedCandDxyError",&Tau_leadChargedCandDxyError,&b_Tau_leadChargedCandDxyError);
 tree->SetBranchAddress("Tau_decayModeFindingNewDMs",&Tau_decayModeFindingNewDMs,&b_Tau_decayModeFindingNewDMs);
 tree->SetBranchAddress("Tau_againstElectronMVALooseMVA5",&Tau_againstElectronMVALooseMVA5,&b_Tau_againstElectronMVALooseMVA5);
 tree->SetBranchAddress("Tau_againstElectronMVAMediumMVA5",&Tau_againstElectronMVAMediumMVA5,&b_Tau_againstElectronMVAMediumMVA5);
 tree->SetBranchAddress("Tau_againstMuonLoose3",&Tau_againstMuonLoose3,&b_Tau_againstMuonLoose3);
 tree->SetBranchAddress("Tau_againstMuonTight3",&Tau_againstMuonTight3,&b_Tau_againstMuonTight3);
 tree->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits",&Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits,&b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
 tree->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits",&Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits,&b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
 tree->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr3Hits",&Tau_byTightCombinedIsolationDeltaBetaCorr3Hits,&b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
 tree->SetBranchAddress("Tau_nProngs",&Tau_nProngs,&b_Tau_nProngs);
 tree->SetBranchAddress("Tau_decayModeFinding",&Tau_decayModeFinding,&b_Tau_decayModeFinding);
 tree->SetBranchAddress("Tau_charge",&Tau_charge,&b_Tau_charge);
 //MET
 Double_t Met_type1PF_px;
 Double_t Met_type1PF_py;
 Double_t Met_type1PF_pz;
 Double_t Met_type1PF_sumEt;

 
 TBranch *b_Met_type1PF_px;  //!
 TBranch *b_Met_type1PF_py;  //!
 TBranch *b_Met_type1PF_pz;  //!
 TBranch *b_Met_type1PF_sumEt;  //!
 tree->SetBranchAddress("Met_type1PF_px",&Met_type1PF_px,&b_Met_type1PF_px);
 tree->SetBranchAddress("Met_type1PF_py",&Met_type1PF_py,&b_Met_type1PF_py);
 tree->SetBranchAddress("Met_type1PF_pz",&Met_type1PF_pz,&b_Met_type1PF_pz);
 tree->SetBranchAddress("Met_type1PF_sumEt",&Met_type1PF_sumEt,&b_Met_type1PF_sumEt);

 vector<double> *Gen_pt;
 Gen_pt = 0;
 TBranch *b_Gen_pt = 0;
 tree->SetBranchAddress("Gen_pt",&Gen_pt,&b_Gen_pt);

 /////
 //   Take the efficiency conditions
 /////
 //Take entries for each variable of cut
 vector<vector<double> > cut_eff;
 for(int iniX=0; iniX<numtot_cuts; iniX++){
  vector<double>  cutentries;
  for(int iniY=0; iniY<nentries; iniY++){
   cutentries.push_back(0);
  }
  cut_eff.push_back(cutentries);
 }
 double weight = 0;
 //All entries
 for(int i=0; i<nentries; i++){
  Long64_t tentry = tree->LoadTree(i);
 //Tau
  b_Tau_pt->GetEntry(tentry);
  b_Tau_eta->GetEntry(tentry);
  b_Tau_phi->GetEntry(tentry);
  b_Tau_energy->GetEntry(tentry);
  b_Tau_leadChargedCandDz_pv->GetEntry(tentry);
  b_Tau_leadChargedCandDxy_bs->GetEntry(tentry);
  b_Tau_leadChargedCandDxyError->GetEntry(tentry);
  b_Tau_decayModeFindingNewDMs->GetEntry(tentry);
  b_Tau_againstElectronMVALooseMVA5->GetEntry(tentry);
  b_Tau_againstElectronMVAMediumMVA5->GetEntry(tentry);
  b_Tau_againstMuonLoose3->GetEntry(tentry);
  b_Tau_againstMuonTight3->GetEntry(tentry);
  b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->GetEntry(tentry);
  b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->GetEntry(tentry);
  b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits->GetEntry(tentry);
  b_Tau_nProngs->GetEntry(tentry);
  b_Tau_decayModeFinding->GetEntry(tentry);
  b_Tau_charge->GetEntry(tentry);
 //MET
  b_Met_type1PF_px->GetEntry(tentry);
  b_Met_type1PF_py->GetEntry(tentry);
  b_Met_type1PF_pz->GetEntry(tentry);
  b_Met_type1PF_sumEt->GetEntry(tentry);

  //Fill Multiplicity Histo
  tauMultiplicity->Fill(Tau_pt->size());


  //Get efficiencies
  int matchedTaus = 0;
  if(genHadronicTauPairExists(tree, i).first && Tau_pt->size() >= 2){
    cout << "Tau pt = " << Tau_pt->at(0) << "\n";
    for(int j = 0; j < 2; j++){
      TLorentzVector recoCand;
      recoCand.SetPtEtaPhiE(Tau_pt->at(j), Tau_eta->at(j), Tau_phi->at(j), Tau_energy->at(j));
      if(matchTauToGen(recoCand, tree, i).first) matchedTaus++;
    }
  }

  if(Tau_pt->size() >= 2){//matchedTaus==2)
   cut_eff[0][i] = 1;
   if(Tau_pt->at(0) > cand2_acc[0]){
    cut_eff[1][i] = 1; 
    if(fabs(Tau_eta->at(0)) < cand2_acc[1]){
     cut_eff[2][i] = 1;
     DMF_d->Fill(Tau_pt->at(0));
     DMFNewDMs_d->Fill(Tau_pt->at(0));
     if(Tau_decayModeFinding->at(0) > 0.5) DMF_n->Fill(Tau_pt->at(0));
     if(Tau_decayModeFindingNewDMs->at(0) > 0.5){
      DMFNewDMs_n->Fill(Tau_pt->at(0));
      cut_eff[3][i] = 1;
      if(Tau_leadChargedCandDz_pv->at(0) < 0.2){
       cut_eff[4][i] = 1;
       if(Tau_againstElectronMVALooseMVA5->at(0) > 0.5){
        cut_eff[5][i] = 1;
        if(Tau_againstMuonTight3->at(0) > 0.5){
         cut_eff[6][i] = 1;
         if(Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(0) > 0.5){
          cut_eff[7][i] = 1;
          if(Tau_pt->at(1) > cand2_acc[0]){
           cut_eff[8][i] = 1;
           if(fabs(Tau_eta->at(1)) < cand2_acc[1]){
            cut_eff[9][i] = 1; 
            if(Tau_decayModeFindingNewDMs->at(1) > 0.5){
             cut_eff[10][i] = 1;
	     if(Tau_leadChargedCandDz_pv->at(1) < 0.2){
	      cut_eff[11][i] = 1;
              if(Tau_againstElectronMVALooseMVA5->at(1) > 0.5){
               cut_eff[12][i] = 1;
               if(Tau_againstMuonTight3->at(1) > 0.5){
                cut_eff[13][i] = 1;
                if(Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(1) > 0.5){
                 cut_eff[14][i] = 1;
		 double DeltaR = deltaR(Tau_eta->at(0), Tau_phi->at(0), Tau_eta->at(1), Tau_phi->at(1));
                 if(DeltaR > 0.5){
		  cut_eff[15][i] = 1;
                  if(Tau_charge->at(0)*Tau_charge->at(1) < 0){
                   cut_eff[16][i] = 1;
                   float cosDPhi = cos(TMath::Abs(normalizedPhi(Tau_phi->at(0) - Tau_phi->at(1))));
		   if(cosDPhi < -0.95){
		    cut_eff[17][i] = 1;
                    if(Met_type1PF_sumEt > 20){
		     cut_eff[18][i] = 1;
                     TLorentzVector Tau1Vector;
		     TLorentzVector Tau2Vector;
		     TLorentzVector METVector(Met_type1PF_px, Met_type1PF_py, Met_type1PF_pz, Met_type1PF_sumEt);
		     Tau1Vector.SetPtEtaPhiE(Tau_pt->at(0), Tau_eta->at(0), Tau_phi->at(0), Tau_energy->at(0));
		     Tau2Vector.SetPtEtaPhiE(Tau_pt->at(1), Tau_eta->at(1), Tau_phi->at(1), Tau_energy->at(1));
		     if(pZeta(Tau1Vector, Tau2Vector, METVector) - 3.1*pZetaVis(Tau1Vector, Tau2Vector) > -50){
		      cut_eff[19][i] = 1;
                      if(nBTags >= 0){
		       cut_eff[20][i] = 1;
		       double Tau1Ip = fabs(Tau_leadChargedCandDxy_bs->at(0));
		       double Tau2Ip = fabs(Tau_leadChargedCandDxy_bs->at(1));
		       double Tau1IpErr = fabs(Tau_leadChargedCandDxyError->at(0));
		       double Tau2IpErr = fabs(Tau_leadChargedCandDxyError->at(1));
		       double combinedIpOverErr = TMath::Sqrt((Tau1Ip*Tau1Ip + 2*Tau1Ip*Tau2Ip + Tau2Ip*Tau2Ip)/(Tau1IpErr*Tau1IpErr + Tau2IpErr*Tau2IpErr));
		       if(combinedIpOverErr >= 2){
	 	        cut_eff[21][i] = 1;
                        cout << "Tau pt = " << Tau_pt->at(0) << "\n";
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
     }
    }          
   }  
  }
 }//End all entries
 
 //tauMultiplicity->Draw();

 TCanvas *c1 = new TCanvas();
 c1->cd();
 TGraphAsymmErrors* theDMFRatio = new TGraphAsymmErrors(DMF_n,DMF_d);
 theDMFRatio->BayesDivide(DMF_n, DMF_d);
 theDMFRatio->SetTitle("Tau DMF");
 theDMFRatio->SetLineColor(lmcol);
 theDMFRatio->SetMarkerColor(lmcol);
 string titleYaxis;
  titleYaxis = "#epsilon(DMF)";
 theDMFRatio->GetXaxis()->SetTitle("Reco Tau p_{T}");
 theDMFRatio->GetYaxis()->SetTitle(titleYaxis.c_str());
 theDMFRatio->GetYaxis()->SetRangeUser(0,1.25);
 theDMFRatio->Draw();


/*
 TCanvas *c5 = new TCanvas();
 c5->cd();
 TGraphAsymmErrors* theRatio = new TGraphAsymmErrors(hpn,hpd);
 theRatio->BayesDivide(hpn,hpd);
 theRatio->SetTitle("Matching Efficiency vs Gen pT");
 theRatio->SetLineColor(lmcol);
 theRatio->SetMarkerColor(lmcol);
 string titleYaxis;
  titleYaxis = "#epsilon(Matching)";
 theRatio->GetXaxis()->SetTitle("Gen Tau p_{T}");
 theRatio->GetYaxis()->SetTitle(titleYaxis.c_str());
 theRatio->GetYaxis()->SetRangeUser(0,1.25);
  theRatio->Draw();
 c5->SaveAs("MatchingEff.pdf");

 visMassHisto->SetLineColor(kRed);
 visPlusMetMassHisto->SetLineColor(kBlue);
 visPlusDeltaPtMassHisto->SetLineColor(kGreen);

 TCanvas *c6 = new TCanvas();
 c6->cd();
 visMassHisto->Draw();
 visPlusMetMassHisto->Draw("same");
 visPlusDeltaPtMassHisto->Draw("same");
*/


 //draw_eff(tree,nentries,cut_eff);
}
/////
//   Draw efficiency plots 
/////


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
    if(Gen_BmotherIndex->at(jj) >= 0 && (abs(Gen_pdg_id->at(jj)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(jj))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(jj)) == 2)  && Gen_pt->at(Gen_BmotherIndex->at(jj)) > 45 && Gen_eta->at(Gen_BmotherIndex->at(jj)) < 2.1 ) {
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
        //if(IsItAHadronicDecayVector.at(i) && IsItAHadronicDecayVector.at(j) && abs(Gen_pdg_id->at(Gen_BmotherIndex->at(tempTauIndexVector.at(i)))) == abs(Gen_pdg_id->at(Gen_BmotherIndex->at(tempTauIndexVector.at(j))))){
        if(IsItAHadronicDecayVector.at(i) && IsItAHadronicDecayVector.at(j) && Gen_BmotherIndex->at(tempTauIndexVector.at(i)) == Gen_BmotherIndex->at(tempTauIndexVector.at(j))){
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
    if( (abs(Gen_pdg_id->at(jj)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(jj))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(jj)) == 2) ) {
      tempTauIndexVector.push_back(Gen_BmotherIndex->at(jj));
      theNeutrinoObject.SetPtEtaPhiE(Gen_pt->at(jj), Gen_eta->at(jj), Gen_phi->at(jj), Gen_energy->at(jj));
      tempNeutrinoMomentumVector.push_back(theNeutrinoObject);
    }
  }

  //---Perform matching only if there is at least one gen tau in the event

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
          float phi_1 = p4Vector.Phi();
          float phi_2 = theGenObject.Phi();
	  float eta_1 = p4Vector.Eta();
	  float eta_2 = theGenObject.Eta();
	  auto dp=std::abs(phi_1 - phi_2); if (dp>(float)(M_PI)) dp-=(float)(2*M_PI);
	  float dR = std::sqrt((eta_1-eta_2)*(eta_1-eta_2) + dp*dp);
          if( (IsItAHadronicDecayVector.at(jjj)) && /*(p4Vector.DeltaR(theGenObject)*/dR <= 0.3)  isGenMatched = true;
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


double pZeta(TLorentzVector obj1, TLorentzVector obj2, TLorentzVector METVector){
 double zetaX = cos(obj1.Phi()) + cos(obj2.Phi());
 double zetaY = sin(obj1.Phi()) + sin(obj2.Phi());
 double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
 if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
 double visPx = obj1.Px() + obj2.Px();
 double visPy = obj1.Py() + obj2.Py();
 double px = visPx + METVector.Px();
 double py = visPy + METVector.Py();
 double pZeta = px*zetaX + py*zetaY;
 return pZeta;
}

double pZetaVis(TLorentzVector obj1, TLorentzVector obj2){
 double zetaX = cos(obj1.Phi()) + cos(obj2.Phi());
 double zetaY = sin(obj1.Phi()) + sin(obj2.Phi());
 double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
 if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
 double visPx = obj1.Px() + obj2.Px();
 double visPy = obj1.Py() + obj2.Py();
 double pZetaVis = visPx*zetaX + visPy*zetaY;
 return pZetaVis;
}


int nBtags(TTree* tree, int nentry){
 vector<double> *Jet_pt;
 Jet_pt = 0;
 TBranch *b_Jet_pt = 0;
 tree->SetBranchAddress("Jet_pt",&Jet_pt,&b_Jet_pt);
 Long64_t tentry = tree->LoadTree(nentry);
 b_Jet_pt->GetEntry(tentry);
 int nBJets = 0;
 for(int i = 0; i < Jet_pt->size(); i++){
  if(!passBJetCuts(i)) continue;
  nBJets++;
 }
 return nBJets;
}

bool passBJetCuts(int nobj){
 return true;
}


double visibleMass(TLorentzVector obj1, TLorentzVector obj2){
  TLorentzVector theVector = obj1 + obj2;
  return theVector.M();
}

double visiblePlusMETMass(TLorentzVector obj1, TLorentzVector obj2, TLorentzVector METVector){
  double px = obj1.Px() + obj2.Px() + METVector.Px();
  double py = obj1.Py() + obj2.Py() + METVector.Py();
  double pz = obj1.Pz() + obj2.Pz();
  double e = obj1.Energy() + obj2.Energy() + TMath::Sqrt(METVector.Px()*METVector.Px() + METVector.Py()*METVector.Py());
  TLorentzVector theVector(px, py, pz, e);
  return theVector.M();
}

double visiblePlusDeltaPtMass(TLorentzVector obj1, TLorentzVector obj2){
  double deltaPtx = -(obj1.Px() + obj2.Px());
  double deltaPty = -(obj1.Py() + obj2.Py());
  double px = obj1.Px() + obj2.Px() + deltaPtx;
  double py = obj1.Py() + obj2.Py() + deltaPty;
  double pz = obj1.Pz() + obj2.Pz();
  double e = obj1.Energy() + obj2.Energy() + TMath::Sqrt(deltaPtx*deltaPtx + deltaPty*deltaPty);
  TLorentzVector theVector(px, py, pz, e);
  return theVector.M();
}


void draw_eff(TTree* tree, int nentries, vector<vector<double> > cut_eff){
 //Define strings to print
 vector<string> varname(numtot_cuts);
 vector<string> varcut(numtot_cuts);
 initialize_strings(varname,varcut);
 //Get pt,eta,phi,nvert
//Tau
 vector<double> *Tau_pt, *Tau_eta, *Tau_phi;
 Tau_pt  = 0;
 Tau_eta = 0;
 Tau_phi = 0;
 TBranch *b_Tau_pt = 0;
 TBranch *b_Tau_eta = 0;
 TBranch *b_Tau_phi = 0;
 tree->SetBranchAddress("Tau_pt",&Tau_pt,&b_Tau_pt);
 tree->SetBranchAddress("Tau_eta",&Tau_eta,&b_Tau_eta);
 tree->SetBranchAddress("Tau_phi",&Tau_phi,&b_Tau_phi);
 //NVertices
 Int_t bestVertices;
 bestVertices = 0;
 TBranch *b_bestVertices;
 //tree->SetBranchAddress("bestVertices", &bestVertices, &b_bestVertices);
 //Draw eff
 for(int v=0; v<numtot_cuts-1; v++){
  bool isTau  = true;
  int taupos = -1;
  //This part depends on how you define the cut flow
  if(v<7 || v>13) taupos = 0;//Lead tau
  else    taupos = 1;//Sub-lead tau

  //Prepare histograms
  TH1* hpn = new TH1F("hpn","hpn",const_bin,const_from,const_to);
  hpn->Sumw2();
  TH1* hpd = new TH1F("hpd","hpd",const_bin,const_from,const_to);
  hpd->Sumw2();
  //Fill them
  double weight = 1;
  double wgt_pu = 1;
  for(int i=0; i<nentries; i++){
   Long64_t tentry = tree->LoadTree(i);
  //Tau
   b_Tau_pt->GetEntry(tentry);
   b_Tau_eta->GetEntry(tentry);
   b_Tau_phi->GetEntry(tentry);
   //nVerts
   //b_bestVertices->GetEntry(tentry);
   if(noLumiNorm) weight = 1.;
   if(noPUcorr) wgt_pu = 1.;
   if(isTau){
    if(cut_eff[v][i]==1){
     if(pt_plots)    hpd->Fill(Tau_pt->at(taupos),weight*wgt_pu);
     if(eta_plots)   hpd->Fill(Tau_eta->at(taupos),weight*wgt_pu);
     if(phi_plots)   hpd->Fill(Tau_phi->at(taupos),weight*wgt_pu);
     if(nvert_plots) hpd->Fill(bestVertices,weight*wgt_pu);
    }
    if(cut_eff[v+1][i]==1){
     if(pt_plots)    hpn->Fill(Tau_pt->at(taupos),weight*wgt_pu);
     if(eta_plots)   hpn->Fill(Tau_eta->at(taupos),weight*wgt_pu);
     if(phi_plots)   hpn->Fill(Tau_phi->at(taupos),weight*wgt_pu);
     if(nvert_plots) hpn->Fill(bestVertices,weight*wgt_pu);
    }
   }
  }
  //Draw plots with asym err
  hist_ratio(hpn,hpd,varname[v+1],varcut[v+1],v);
  delete hpn;
  delete hpd;
 }
}
/////
//   Initialize the strings 
/////
void initialize_strings(vector<string> &varname, vector<string> &varcut){
 varname[0]  = "2 reco tau_{h} (gen matched)";
 varname[1]  = "Lead tau Pt > 20";    varname[2]  = "Lead tau |eta| < 2.3";
 varname[3]  = "Lead tau DMF";        varname[4]  = "Lead tau dz < 0.2";
 varname[5]  = "Lead tau AgEle";      varname[6]  = "Lead tau AgMu";       
 varname[7]  = "Lead tau Iso";        varname[8]  = "SubLead tau Pt > 20";
 varname[9]  = "SubLead tau |eta| < 2.3"; 
 varname[10]  = "SubLead tau DMF";    varname[11] = "SubLead tau dz < 0.2";
 varname[12] = "SubLead tau AgEle";   varname[13] = "SubLead tau AgMu";   
 varname[14] = "SubLead tau Iso";     varname[15] = "#DeltaR(tau1,tau2)";
 varname[16] = "Tau1Tau2 Charge < 0"; varname[17] = "Cos(#Delta#phi Tau1, Tau2) < -0.95";
 varname[18] = "ME_{T} > 20";	      varname[19] = "PZeta - 3.1PZeta_vis < -50";
 varname[20] = "nBtags = 0";	      varname[21] = "Combined Ip/Err > 2";
}
/////
//   Evaluate ratio of histograms with asymm error
/////
void hist_ratio(TH1* hpn, TH1* hpd, string varname, string varcut, int num){
 //TCanvas* c1 = new TCanvas(varname.c_str(),varname.c_str(),100,100,480,400);
 TCanvas* c1 = new TCanvas(varname.c_str(),varname.c_str(),100,100,700,500);
 TGraphAsymmErrors* hpratio = new TGraphAsymmErrors(hpn,hpd);
 hpratio->BayesDivide(hpn,hpd);
 hpratio->SetLineColor(lmcol);
 hpratio->SetMarkerColor(lmcol);
 hpratio->SetTitle(0);
 hpratio->GetXaxis()->SetTitle(titleXaxis.c_str());
 char bin_size_c[1000]; float bin_size_f = (fabs(const_from-const_to)/const_bin); sprintf(bin_size_c,"%.2f",bin_size_f);
 string titleYaxis;
 //if(sigreg && !releffsig){
 // titleYaxis = "#epsilon(N-1 "+varname+varcut+")";//"+(string) bin_size_c;
 //}else{
  titleYaxis = "#epsilon("+varname+varcut+")";//"+(string) bin_size_c;
 //}
 hpratio->GetYaxis()->SetTitle(titleYaxis.c_str());
 hpratio->GetYaxis()->SetRangeUser(0,1.25);
 hpratio->Draw("PAZ9");
 c1->Update();
 char nums[1000]; sprintf(nums,"%d",num);
 string namefile = part+"_"+varname+".pdf";
 if(save_plots) c1->SaveAs(namefile.c_str());
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
  tdrStyle->SetErrorX(0.);
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
  tdrStyle->SetOptStat(1111); // To display the mean and RMS:   SetOptStat("mr");
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
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);

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
  tdrStyle->SetTitleSize(0.06, "X");
  tdrStyle->SetTitleSize(0.075, "Y");
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
