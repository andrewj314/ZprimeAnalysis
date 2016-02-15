/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects  
   It makes plots for rel eff following the cut flow 
   For event selection you can choose between rel eff plots or N-1 eff (only integral or vs event quantities)
   Special eff plots can also be done  

Need to specify
1. See Declare constants

To do
- //This part depends on how you define the cut flow
  ... May need to change it according to the efficiency you want to study 
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
#include <TLorentzVector.h>
#include "TROOT.h"
#include "TSystem.h"
#include <utility>
#include <vector>
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
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
const bool pt_plots = true;
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
const double cand1_acc[2]   = {18, 2.1};
const double cand2_acc[2]   = {20, 2.3};
const int numtot_cuts       = 13;
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
int get_ditau_br(TTree* tree, int entry);
pair<bool,TLorentzVector> leptonTau_matching(TTree* tree, int entry, int candid, TLorentzVector reco_cand);
void setTDRStyle();
/////
//   Main function
/////
void TauTau_Plots_EffpT_Effeta_Effpu_Francesco(){
 setTDRStyle();
 //Call tree 
 TFile* f = Call_TFile();
 TTree* tree; f->GetObject("BOOM",tree);
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
 //Tau
 vector<double> *Tau_pt, *Tau_eta, *Tau_phi, *Tau_energy;//, *Tau_leadChargedCand_dz;
 vector<int> *Tau_decayModeFindingNewDMs, *Tau_againstElectronMVALooseMVA5, *Tau_againstMuonTight3, *Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
 Tau_pt = 0;
 Tau_eta = 0;
 Tau_phi = 0;
 Tau_energy = 0;
 //Tau_leadChargedCand_dz = 0;
 Tau_decayModeFindingNewDMs = 0;
 Tau_againstElectronMVALooseMVA5 = 0;
 Tau_againstMuonTight3 = 0;
 Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
 TBranch *b_Tau_pt = 0;
 TBranch *b_Tau_eta = 0;
 TBranch *b_Tau_phi = 0;
 TBranch *b_Tau_energy = 0;
 //TBranch *b_Tau_leadChargedCand_dz = 0;
 TBranch *b_Tau_decayModeFindingNewDMs = 0;
 TBranch *b_Tau_againstElectronMVALooseMVA5 = 0;
 TBranch *b_Tau_againstMuonTight3 = 0;
 TBranch *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0; 
 tree->SetBranchAddress("Tau_pt",&Tau_pt,&b_Tau_pt);
 tree->SetBranchAddress("Tau_eta",&Tau_eta,&b_Tau_eta);
 tree->SetBranchAddress("Tau_phi",&Tau_phi,&b_Tau_phi);
 tree->SetBranchAddress("Tau_energy",&Tau_energy,&b_Tau_energy);
 //tree->SetBranchAddress("Tau_leadChargedCand_dz",&Tau_leadChargedCand_dz,&b_Tau_leadChargedCand_dz);
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
  //Tau
  b_Tau_pt->GetEntry(tentry);
  b_Tau_eta->GetEntry(tentry);
  b_Tau_phi->GetEntry(tentry);
  b_Tau_energy->GetEntry(tentry);
  //b_Tau_leadChargedCand_dz->GetEntry(tentry);
  b_Tau_decayModeFindingNewDMs->GetEntry(tentry);
  b_Tau_againstElectronMVALooseMVA5->GetEntry(tentry);
  b_Tau_againstMuonTight3->GetEntry(tentry);
  b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->GetEntry(tentry);
  //Get the BR values
  int decaychannel = get_ditau_br(tree,i);
  //Find the num of matched candidates, among leading and subleading reco cands
  int matched_tau = 0;
  if(decaychannel==0 && Tau_pt->size()>=2){
   for(uint rc = 0; rc<2; rc++){
    TLorentzVector reco_cand(0.,0.,0.,0.);
    reco_cand.SetPtEtaPhiE(Tau_pt->at(rc),Tau_eta->at(rc),Tau_phi->at(rc),Tau_energy->at(rc));
    pair<bool,TLorentzVector> curr_pair = leptonTau_matching(tree,tentry,15,reco_cand);
    if(curr_pair.first) matched_tau++;
   }
  }
  //Get efficiencies
  if(matched_tau==2){
   cut_eff[0][i] = 1;
   if(Tau_pt->at(0) > cand2_acc[0]){
    cut_eff[1][i] = 1;
    if(fabs(Tau_eta->at(0)) < cand2_acc[1]){
     cut_eff[2][i] = 1;
     if(Tau_decayModeFindingNewDMs->at(0) > 0.5){
      cut_eff[3][i] = 1;
      if(Tau_againstElectronMVALooseMVA5->at(0) > 0.5){
       cut_eff[4][i] = 1;
       if(Tau_againstMuonTight3->at(0) > 0.5){
        cut_eff[5][i] = 1;
        if(Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(0) > 0.5){
         cut_eff[6][i] = 1;
         if(Tau_pt->at(1) > cand2_acc[0]){
          cut_eff[7][i] = 1;
          if(fabs(Tau_eta->at(1)) < cand2_acc[1]){
           cut_eff[8][i] = 1;
           if(Tau_decayModeFindingNewDMs->at(1) > 0.5){
            cut_eff[9][i] = 1;
            if(Tau_againstElectronMVALooseMVA5->at(1) > 0.5){
             cut_eff[10][i] = 1;
             if(Tau_againstMuonTight3->at(1) > 0.5){
              cut_eff[11][i] = 1;
              if(Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(1) > 0.5){
               cut_eff[12][i] = 1;
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
 //Int_t bestVertices;
 //bestVertices = 0;
 //TBranch *b_bestVertices;
 //tree->SetBranchAddress("bestVertices", &bestVertices, &b_bestVertices);
 //Draw eff
 for(int v=0; v<numtot_cuts-1; v++){//numtot_cuts-1; v++){
  //This part depends on how you define the cut flow
  bool isTau = true;
  int taupos = -1;
  if(v<6) taupos = 0;//Lead tau
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
     //if(nvert_plots) hpd->Fill(bestVertices,weight*wgt_pu);
    }
    if(cut_eff[v+1][i]==1){
     if(pt_plots)    hpn->Fill(Tau_pt->at(taupos),weight*wgt_pu);
     if(eta_plots)   hpn->Fill(Tau_eta->at(taupos),weight*wgt_pu);
     if(phi_plots)   hpn->Fill(Tau_phi->at(taupos),weight*wgt_pu);
     //if(nvert_plots) hpn->Fill(bestVertices,weight*wgt_pu);
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
 varname[3]  = "Lead tau DMF";        varname[4]  = "Lead tau AgEle";
 varname[5]  = "Lead tau AgMu";       varname[6]  = "Lead tau Iso";
 varname[7]  = "SubLead tau Pt > 20"; varname[8]  = "SubLead tau |eta| < 2.3";
 varname[9]  = "SubLead tau DMF";     varname[10] = "SubLead tau AgEle";
 varname[11] = "SubLead tau AgMu";    varname[12] = "SubLead tau Iso";
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
//   Tau info
/////
int get_ditau_br(TTree* tree, int entry){
 int decaychannel = -1;
 vector<double> *Gen_pdg_id, *Gen_motherpdg_id;
 Gen_pdg_id = 0;
 Gen_motherpdg_id = 0;
 TBranch *b_Gen_pdg_id = 0;
 TBranch *b_Gen_motherpdg_id = 0;
 tree->SetBranchAddress("Gen_pdg_id",&Gen_pdg_id,&b_Gen_pdg_id);
 tree->SetBranchAddress("Gen_motherpdg_id",&Gen_motherpdg_id,&b_Gen_motherpdg_id);
 vector<double> *Gen_BmotherIndex;
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
 //Search for tau decays of each tau
 vector<int> taudecays; taudecays.clear();
 for(uint gp=0; gp<Gen_pdg_id->size(); gp++){
  //Take a tau
  if(!(abs(Gen_pdg_id->at(gp))==15 && Gen_status->at(gp)==2)) continue;
  //Look at tau decays
  int etau    = 0;
  int mutau   = 0;
  int tauhtau = 0;
  for(uint gpd=0; gpd<Gen_pdg_id->size(); gpd++){
   if(!(Gen_BmotherIndex->at(gpd)==gp)) continue;
   if(abs(Gen_pdg_id->at(gpd))==11 || abs(Gen_pdg_id->at(gpd))==12 || abs(Gen_pdg_id->at(gpd))==16) etau    += 1;
   if(abs(Gen_pdg_id->at(gpd))==13 || abs(Gen_pdg_id->at(gpd))==14 || abs(Gen_pdg_id->at(gpd))==16) mutau   += 1;
   if(abs(Gen_pdg_id->at(gpd))==12 || abs(Gen_pdg_id->at(gpd))==14 || abs(Gen_pdg_id->at(gpd))==16) tauhtau += 1;
  }
  if(etau==3)    taudecays.push_back(1);
  if(mutau==3)   taudecays.push_back(2);
  if(tauhtau==1) taudecays.push_back(3);
 }
 if(taudecays.size()>=2){
  if(taudecays[0]==3  && taudecays[1]==3)                                          decaychannel = 0; //tauh tauh
  if((taudecays[0]==3 && taudecays[1]==2) || (taudecays[0]==2 && taudecays[1]==3)) decaychannel = 1; //tauh mu
  if((taudecays[0]==3 && taudecays[1]==1) || (taudecays[0]==1 && taudecays[1]==3)) decaychannel = 2; //tauh e
  if((taudecays[0]==2 && taudecays[1]==1) || (taudecays[0]==1 && taudecays[1]==2)) decaychannel = 3; //mu   e
  if(taudecays[0]==2  && taudecays[1]==2)                                          decaychannel = 4; //mu   mu 
  if(taudecays[0]==1  && taudecays[1]==1)                                          decaychannel = 5; //e    e  
 }
 return decaychannel;
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
