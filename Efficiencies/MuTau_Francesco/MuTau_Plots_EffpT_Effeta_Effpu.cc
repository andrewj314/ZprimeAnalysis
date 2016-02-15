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
//   To run: root -l MuTau_Plots_EffpT_Effeta_Effpu.cc+ 
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
using namespace ROOT::Math;
/////
//   Declare constants
/////
//Path and root file
const string path       = "/uscms_data/d3/andrewj/CMSSW_7_4_7/src/Efficiencies/MuTau_Francesco/";
const string rootpla    = "OutTree.root";
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
const int const_bin = 200; const double const_from = 0; const double const_to = 1000;
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
const double cand1_acc[2]   = {18, 2.1};
const double cand2_acc[2]   = {20, 2.3};
const int numtot_cuts       = 15;
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
void draw_eff(TTree* tree, int nentries, vector<vector<double> > cut_eff);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
void hist_ratio(TH1* hpn, TH1* hpd, string varname, string varcut, int num);
void setTDRStyle();
/////
//   Main function
/////
void MuTau_Plots_EffpT_Effeta_Effpu(){
 setTDRStyle();
 //Call tree 
 TFile* f = Call_TFile();
 TTree* tree; f->GetObject("TNT/BOOM;3",tree);
 int nentries = tree->GetEntries(); 
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
 //Muon
 vector<double> *Muon_pt, *Muon_eta, *Muon_phi, *Muon_energy, *Muon_dxy, *Muon_dz, *Muon_isoCharged, *Muon_isoNeutralHadron, *Muon_isoPhoton, *Muon_isoPU, *Muon_isMediumMuon, *Muon_isTrackerMuon, *Muon_charge;
 vector<int> *Muon_isglobal, *Muon_pf;  
 Muon_pt = 0;
 Muon_eta = 0;
 Muon_phi = 0;
 Muon_energy = 0;
 Muon_dxy = 0;
 Muon_dz = 0;
 Muon_isoCharged = 0;
 Muon_isoNeutralHadron = 0;
 Muon_isoPhoton = 0;
 Muon_isoPU = 0;
 Muon_charge = 0;
 Muon_isMediumMuon = 0;
 Muon_isTrackerMuon = 0;
 Muon_isglobal = 0;
 Muon_pf = 0;
 TBranch *b_Muon_pt = 0;
 TBranch *b_Muon_eta = 0;
 TBranch *b_Muon_phi = 0;
 TBranch *b_Muon_energy = 0;
 TBranch *b_Muon_dxy = 0;
 TBranch *b_Muon_dz = 0;
 TBranch *b_Muon_isoCharged = 0;
 TBranch *b_Muon_isoNeutralHadron = 0;
 TBranch *b_Muon_isoPhoton = 0;
 TBranch *b_Muon_isoPU = 0;
 TBranch *b_Muon_charge = 0;
 TBranch *b_Muon_isMediumMuon = 0;
 TBranch *b_Muon_isglobal = 0;
 TBranch *b_Muon_isTrackerMuon = 0;
 TBranch *b_Muon_pf = 0;
 tree->SetBranchAddress("Muon_pt",&Muon_pt,&b_Muon_pt);
 tree->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
 tree->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
 tree->SetBranchAddress("Muon_energy",&Muon_energy,&b_Muon_energy);
 tree->SetBranchAddress("Muon_dxy",&Muon_dxy,&b_Muon_dxy);
 tree->SetBranchAddress("Muon_dz",&Muon_dz,&b_Muon_dz);
 tree->SetBranchAddress("Muon_isoCharged",&Muon_isoCharged,&b_Muon_isoCharged);
 tree->SetBranchAddress("Muon_isoNeutralHadron",&Muon_isoNeutralHadron,&b_Muon_isoNeutralHadron);
 tree->SetBranchAddress("Muon_isoPhoton",&Muon_isoPhoton,&b_Muon_isoPhoton);
 tree->SetBranchAddress("Muon_isoPU",&Muon_isoPU,&b_Muon_isoPU);
 tree->SetBranchAddress("Muon_charge",&Muon_charge,&b_Muon_charge);
 tree->SetBranchAddress("Muon_isMediumMuon",&Muon_isMediumMuon,&b_Muon_isMediumMuon);
 tree->SetBranchAddress("Muon_isTrackerMuon",&Muon_isTrackerMuon,&b_Muon_isTrackerMuon);
 tree->SetBranchAddress("Muon_isglobal",&Muon_isglobal,&b_Muon_isglobal);
 tree->SetBranchAddress("Muon_pf",&Muon_pf,&b_Muon_pf);
 //Tau
 vector<double> *Tau_pt, *Tau_eta, *Tau_phi, *Tau_energy, *Tau_leadChargedCand_dz;
 vector<int> *Tau_decayModeFindingNewDMs, *Tau_againstElectronMVALooseMVA5, *Tau_againstMuonTight3, *Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
 Tau_pt = 0;
 Tau_eta = 0;
 Tau_phi = 0;
 Tau_energy = 0;
 Tau_leadChargedCand_dz = 0;
 Tau_decayModeFindingNewDMs = 0;
 Tau_againstElectronMVALooseMVA5 = 0;
 Tau_againstMuonTight3 = 0;
 Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
 TBranch *b_Tau_pt = 0;
 TBranch *b_Tau_eta = 0;
 TBranch *b_Tau_phi = 0;
 TBranch *b_Tau_energy = 0;
 TBranch *b_Tau_leadChargedCand_dz = 0;
 TBranch *b_Tau_decayModeFindingNewDMs = 0;
 TBranch *b_Tau_againstElectronMVALooseMVA5 = 0;
 TBranch *b_Tau_againstMuonTight3 = 0;
 TBranch *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0; 
 tree->SetBranchAddress("Tau_pt",&Tau_pt,&b_Tau_pt);
 tree->SetBranchAddress("Tau_eta",&Tau_eta,&b_Tau_eta);
 tree->SetBranchAddress("Tau_phi",&Tau_phi,&b_Tau_phi);
 tree->SetBranchAddress("Tau_energy",&Tau_energy,&b_Tau_energy);
 tree->SetBranchAddress("Tau_leadChargedCand_dz",&Tau_leadChargedCand_dz,&b_Tau_leadChargedCand_dz);
 tree->SetBranchAddress("Tau_decayModeFindingNewDMs",&Tau_decayModeFindingNewDMs,&b_Tau_decayModeFindingNewDMs);
 tree->SetBranchAddress("Tau_againstElectronMVALooseMVA5",&Tau_againstElectronMVALooseMVA5,&b_Tau_againstElectronMVALooseMVA5);
 tree->SetBranchAddress("Tau_againstMuonTight3",&Tau_againstMuonTight3,&b_Tau_againstMuonTight3);
 tree->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits",&Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits,&b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
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
  //Muon
  b_Muon_pt->GetEntry(tentry);
  b_Muon_eta->GetEntry(tentry);
  b_Muon_phi->GetEntry(tentry);
  b_Muon_energy->GetEntry(tentry);
  b_Muon_dxy->GetEntry(tentry);
  b_Muon_dz->GetEntry(tentry);
  b_Muon_isoCharged->GetEntry(tentry);
  b_Muon_isoNeutralHadron->GetEntry(tentry);
  b_Muon_isoPhoton->GetEntry(tentry);
  b_Muon_isoPU->GetEntry(tentry);
  b_Muon_charge->GetEntry(tentry);
  b_Muon_isMediumMuon->GetEntry(tentry);
  b_Muon_isglobal->GetEntry(tentry);
  b_Muon_isTrackerMuon->GetEntry(tentry);
  b_Muon_pf->GetEntry(tentry);
  //Tau
  b_Tau_pt->GetEntry(tentry);
  b_Tau_eta->GetEntry(tentry);
  b_Tau_phi->GetEntry(tentry);
  b_Tau_energy->GetEntry(tentry);
  b_Tau_leadChargedCand_dz->GetEntry(tentry);
  b_Tau_decayModeFindingNewDMs->GetEntry(tentry);
  b_Tau_againstElectronMVALooseMVA5->GetEntry(tentry);
  b_Tau_againstMuonTight3->GetEntry(tentry);
  b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->GetEntry(tentry);
  //Get efficiencies
  if(Muon_pt->size() > 0 && Tau_pt->size() > 0){
    cut_eff[0][i] = 1;
    if(Muon_pt->at(0) > cand1_acc[0]){
      cut_eff[1][i] = 1;
      if(fabs(Muon_eta->at(0)) < cand1_acc[1]){
        cut_eff[2][i] = 1;
        if(Tau_pt->at(0) > cand2_acc[0]){
          cut_eff[3][i] = 1;
          if(fabs(Tau_eta->at(0)) < cand2_acc[1]){
            cut_eff[4][i] = 1;
            if(Muon_dxy->at(0) < 0.045 && Muon_dz->at(0) < 0.2){
              cut_eff[5][i] = 1;
              if(Muon_isMediumMuon->at(0) > 0){
                cut_eff[6][i] = 1;
                if(Tau_decayModeFindingNewDMs->at(0) > 0.5){
                  cut_eff[7][i] = 1;
                  if(Tau_leadChargedCand_dz->at(0) < 0.2){
                    cut_eff[8][i] = 1;
                    PtEtaPhiEVector PtEtaPhiEMu0(Muon_pt->at(0),Muon_eta->at(0),Muon_phi->at(0),Muon_energy->at(0));
                    PtEtaPhiEVector PtEtaPhiETau0(Tau_pt->at(0),Tau_eta->at(0),Tau_phi->at(0),Tau_energy->at(0));
		    double Mu0Tau0_DeltaR = VectorUtil::DeltaR(PtEtaPhiEMu0,PtEtaPhiETau0);
		    if(Mu0Tau0_DeltaR > 0.5){
		      cut_eff[9][i] = 1;
		      double Muon_Iso;
		      if(Muon_isoCharged->size() > 0 && Muon_isoNeutralHadron->size() > 0 && Muon_isoPhoton->size() > 0 && Muon_isoPU->size() > 0){
                        Muon_Iso = (Muon_isoCharged->at(0) + max(Muon_isoNeutralHadron->at(0) + Muon_isoPhoton->at(0) - 0.5*Muon_isoPU->at(0), 0.0))/Muon_pt->at(0);
		      }else{
                        Muon_Iso = 10;
                      }
		      if(Muon_Iso < 0.1){
		        cut_eff[10][i] = 1;
		        if(Tau_againstElectronMVALooseMVA5->at(0) > 0.5){
			  cut_eff[11][i] = 1;
			  if(Tau_againstMuonTight3->at(0) > 0.5){
			    cut_eff[12][i] = 1;
			    if(Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(0) < 1.5){
			      cut_eff[13][i] = 1;
			      //Di-Muon Veto
			      if(Muon_pt->size() < 2){
  				cut_eff[14][i] = 1;
			      }else{
                                PtEtaPhiEVector PtEtaPhiEMu1(Muon_pt->at(1),Muon_eta->at(1),Muon_phi->at(1),Muon_energy->at(1));
		                double Mu0Mu1_DeltaR = VectorUtil::DeltaR(PtEtaPhiEMu0,PtEtaPhiEMu1);
				double Muon_Iso_1 = (Muon_isoCharged->at(1) + max(Muon_isoNeutralHadron->at(1) + Muon_isoPhoton->at(1) - 0.5*Muon_isoPU->at(1), 0.0))/Muon_pt->at(1);
                                bool secondLep  = false;
                                if(Muon_pt->at(1) > 15 && fabs(Muon_eta->at(1)) < 2.4 
                                   && Muon_isglobal->at(1) && Muon_isTrackerMuon->at(1) && Muon_pf->at(1) 
                                   && Muon_dxy->at(1) < 0.045 && Muon_dz->at(1) < 0.2 && Muon_Iso_1 < 0.3 
                                   && Mu0Mu1_DeltaR>0.15 && Muon_charge->at(0)*Muon_charge->at(1)==-1
                                  )
                                    secondLep = true;
			        if(!secondLep){
			          cut_eff[14][i] = 1;
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
 draw_eff(tree,nentries,cut_eff);
}
/////
//   Draw efficiency plots 
/////
void draw_eff(TTree* tree, int nentries, vector<vector<double> > cut_eff){
 //Define strings to print
 vector<string> varname(numtot_cuts);
 vector<string> varcut(numtot_cuts);
 initialize_strings(varname,varcut);
 //Get pt,eta,phi,nvert
 //Muon
 vector<double> *Muon_pt, *Muon_eta, *Muon_phi;
 Muon_pt = 0;
 Muon_eta = 0;
 Muon_phi = 0;
 TBranch *b_Muon_pt = 0;
 TBranch *b_Muon_eta = 0;
 TBranch *b_Muon_phi = 0;
 tree->SetBranchAddress("Muon_pt",&Muon_pt,&b_Muon_pt);
 tree->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
 tree->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
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
 tree->SetBranchAddress("bestVertices", &bestVertices, &b_bestVertices);
 //Draw eff
 for(int v=0; v<numtot_cuts-1; v++){
  bool isMu  = false;
  bool isTau = false;
  //This part depends on how you define the cut flow
  //v1 = Mu Pt, v2 = Mu eta
  if(0<=v && v<=1)   isMu  = true;
  //v3 = Tau Pt, v4 = Tau eta
  if(2<=v && v<=3)   isTau = true;
  //v5 = Mu dxy&dz, v6 = Mu ID
  if(4<=v && v<=5)   isMu  = true;
  //v7 = Tau DMF, v8 Tau dz
  if(6<=v && v<=7)   isTau = true;
  //v9 = MuTauDr, v10 MuIso
  if(8<=v && v<=9)   isMu  = true;
  //v11 = TauAntiE, v12 = TauAntiMu, v13 = TauIso
  if(10<=v && v<=12) isTau = true; 
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
   //Muon
   b_Muon_pt->GetEntry(tentry);
   b_Muon_eta->GetEntry(tentry);
   b_Muon_phi->GetEntry(tentry);
   //Tau
   b_Tau_pt->GetEntry(tentry);
   b_Tau_eta->GetEntry(tentry);
   b_Tau_phi->GetEntry(tentry);
   //nVerts
   b_bestVertices->GetEntry(tentry);
   if(noLumiNorm) weight = 1.;
   if(noPUcorr) wgt_pu = 1.;
   if(isMu){
    if(cut_eff[v][i]==1){
     if(pt_plots)    hpd->Fill(Muon_pt->at(0),weight*wgt_pu);
     if(eta_plots)   hpd->Fill(Muon_eta->at(0),weight*wgt_pu);
     if(phi_plots)   hpd->Fill(Muon_phi->at(0),weight*wgt_pu);
     if(nvert_plots) hpd->Fill(bestVertices,weight*wgt_pu);
    }
    if(cut_eff[v+1][i]==1){
     if(pt_plots)    hpn->Fill(Muon_pt->at(0),weight*wgt_pu);
     if(eta_plots)   hpn->Fill(Muon_eta->at(0),weight*wgt_pu);
     if(phi_plots)   hpn->Fill(Muon_phi->at(0),weight*wgt_pu);
     if(nvert_plots) hpn->Fill(bestVertices,weight*wgt_pu);
    }
   }
   if(isTau){
    if(cut_eff[v][i]==1){
     if(pt_plots)    hpd->Fill(Tau_pt->at(0),weight*wgt_pu);
     if(eta_plots)   hpd->Fill(Tau_eta->at(0),weight*wgt_pu);
     if(phi_plots)   hpd->Fill(Tau_phi->at(0),weight*wgt_pu);
     if(nvert_plots) hpd->Fill(bestVertices,weight*wgt_pu);
    }
    if(cut_eff[v+1][i]==1){
     if(pt_plots)    hpn->Fill(Tau_pt->at(0),weight*wgt_pu);
     if(eta_plots)   hpn->Fill(Tau_eta->at(0),weight*wgt_pu);
     if(phi_plots)   hpn->Fill(Tau_phi->at(0),weight*wgt_pu);
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
 varname[0] = ">= 1 Mu and 1 Tau\t";
 varname[1] = "Mu Pt > 18 \t\t"; varname[2] = "Mu |eta| < 2.1 \t\t";
 varname[3] = "Tau Pt > 20 \t\t\t"; varname[4] = "Tau |eta| < 2.3 \t\t";
 varname[5] = "Mu track dxy < 0.045 & dz < 0.2\t"; varname[6] = "Mu MediumID\t\t";
 varname[7] = "Tau DMF new DMs\t"; varname[8] = "Tau dz < 0.2 \t\t\t";
 varname[9] = "MuTau DeltaR > 0.5 \t\t"; varname[10] = "Mu Iso < 0.1\t\t";
 varname[11] = "Tau anti-e MVA5\t"; varname[12] = "Tau anti-mu tight\t";
 varname[13] = "Tau Iso DB 3hits < 1.5"; varname[14] = "Di-mu veto\t\t";
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
 string namefile = part+(string)nums+".pdf";
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
