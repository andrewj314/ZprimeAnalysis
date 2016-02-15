/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects  
   It makes plots for rel eff following cut flow 
   For signal you can choose between rel eff plots or N-1 eff (only integral or vs nvert)
   Special eff plots can also be done  

Need to specify
1. See Declare constants

Notes
Depending on what you are studying (objects, sigreg, speceff)
you must pay attention of what you plot (pt,eta,phi,massVis,nvert).
E.g. you can plot massVis if you do not have 2 candidates
*/
/////
//   To run: root -l Efficiency_cutflow_cutrelEff_plots.cc+ 
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
const string path     = "/uscms_data/d3/andrewj/CMSSW_7_4_0_patch1/src/Efficiencies/"; //"/afs/cern.ch/user/f/fromeo/public/EXO/Efficiency/";
const string rootpla  = "Sel_emu_SSMToTauTau1250.root";
const double Luminosity = 19703.225;//19779.362;//19600; //pb^-1 
const bool noLumiNorm = true; //true means NO luminosity normalization done
const bool noPUcorr   = true; //true means NO PU corr done
const int cand1_acc       = 3;
const int cand2_acc       = 3;
const int cand1_numTotCut = 11; //It is different (-1) from what you could see in ZpmSelection.cc, but still correct
const int cand2_numTotCut = 12; //It is different (-1) from what you could see in ZpmSelection.cc, but still correct
const int numTotCut       = 24;//dR+cand1_numTotCut+cand2_numTotCut 
const int numsigCut       = 6;
const bool objects        = true;
const bool sigreg         = false; const bool releffsig = false; //false means N-1 eff for SR cuts
const bool speceff        = false;
const int numspeceffini   = 4;
const int  numspecefffin  = 7;
//Plots
const bool save_plots = true; 
//Choose gen pt,eta or reco pt,eta
const bool gen_pteta = true;
//Choose Pt aut Eta aut NVert and comments constants for others
//Pt
const bool pt_plots = false;
//const int const_bin = 100; const double const_from = 0; const double const_to = 1500;
//const string titleXaxis = "Gen_pT (GeV/c)"; const int lmcol = 4; //Blue
//const string titleXaxis = "Reco_pT (GeV/c)"; const int lmcol = 2; //Red
//const string part = "genpT_";
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
const bool nvert_plots = true;
const int const_bin = 50; const double const_from = 0; const double const_to = 50;
const string titleXaxis = "Vertices"; const int lmcol = 4; //Black
const string part = "nvert_";
//MassVis
const bool massvis_plots = false;
//const int const_bin = 20; const double const_from = 0; const double const_to = 1500;
//const string titleXaxis = "massVis(e,#mu,E_{T}^{miss}) (GeV/c^{2})"; const int lmcol = 4; //Black
//const string part = "massvis_";
/////
//   Declare functions 
/////
TFile* Call_TFile();
void cutflow_eff_obj(TTree* tree, int nentries);
void cutflow_eff_sr(TTree* tree, int nentries);
void cutflow_eff_spec(TTree* tree, int nentries);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
void initialize_strings2(vector<string> &varname, vector<string> &varcut);
void hist_ratio(TH1* hpn, TH1* hpd, string varname, string varcut, int num);
void setTDRStyle();
/////
//   Main function
/////
void Efficiency_cutflow_cutrelEff_plots() {
 setTDRStyle();
 //Call tree 
 TFile* f = Call_TFile();
 TTree* tree; f->GetObject("Selection/tree",tree);
 int nentries = tree->GetEntries(); 
 if(objects) cutflow_eff_obj(tree, nentries);
 if(sigreg)  cutflow_eff_sr(tree, nentries);
 if(speceff) cutflow_eff_spec(tree, nentries);
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
 //Define strings to print
 vector<string> varname(numTotCut+numsigCut); 
 vector<string> varcut(numTotCut+numsigCut);
 initialize_strings(varname,varcut);
 //Call variables 
 int num1Cut, num2Cut, nvert1, nvert2;
 tree->SetBranchAddress("cand1_numCut",&num1Cut);
 tree->SetBranchAddress("cand2_numCut",&num2Cut);
 tree->SetBranchAddress("cand1_nvert",&nvert1);
 tree->SetBranchAddress("cand2_nvert",&nvert2);
 double wgt_lumi, wgt_pu, pt1, eta1, phi1, pt2, eta2, phi2;
 tree->SetBranchAddress("wgt_lumi",&wgt_lumi);
 tree->SetBranchAddress("wgt_pu",&wgt_pu);
 if(gen_pteta){
  tree->SetBranchAddress("cand1_pt_gen",&pt1);
  tree->SetBranchAddress("cand1_eta_gen",&eta1);
  tree->SetBranchAddress("cand1_phi_gen",&phi1);
  tree->SetBranchAddress("cand2_pt_gen",&pt2);
  tree->SetBranchAddress("cand2_eta_gen",&eta2);
  tree->SetBranchAddress("cand2_phi_gen",&phi2);
 }else{
  tree->SetBranchAddress("cand1_pt",&pt1);
  tree->SetBranchAddress("cand1_eta",&eta1);
  tree->SetBranchAddress("cand1_phi",&phi1);
  tree->SetBranchAddress("cand2_pt",&pt2);
  tree->SetBranchAddress("cand2_eta",&eta2);
  tree->SetBranchAddress("cand2_phi",&phi2);
 }
 double charge, cosDphi, met, pZetaMt, pZetaVisMt, cosDphiLMet, massVis;
 tree->SetBranchAddress("charge",&charge);
 tree->SetBranchAddress("cosDphi",&cosDphi);
 tree->SetBranchAddress("met",&met);
 tree->SetBranchAddress("pZetaMt",&pZetaMt);
 tree->SetBranchAddress("pZetaVisMt",&pZetaVisMt);
 tree->SetBranchAddress("cosDphiLMet",&cosDphiLMet);
 tree->SetBranchAddress("massVis",&massVis);
 int pair_nvert, nbjet;
 tree->SetBranchAddress("pair_nvert",&pair_nvert);
 tree->SetBranchAddress("nbjet",&nbjet);
 //Make plots
 double weight = 0;
 //Acceptance
 //cand1 
 for(int v=1; v<cand1_acc; v++){ 
  //Prepare histograms
  TH1* hpn = new TH1F("hpn","hpn",const_bin,const_from,const_to);
  hpn->Sumw2();
  TH1* hpd = new TH1F("hpd","hpd",const_bin,const_from,const_to);
  hpd->Sumw2();
  //Fill them
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   weight = wgt_lumi;
   weight = weight*Luminosity/100.;
   if(noLumiNorm) weight = 1.;
   if(noPUcorr) wgt_pu = 1.;
   if(num1Cut>=v){
    if(pt_plots)    hpd->Fill(pt1,weight*wgt_pu);
    if(eta_plots)   hpd->Fill(eta1,weight*wgt_pu);
    if(phi_plots)   hpd->Fill(phi1,weight*wgt_pu);
    if(nvert_plots) hpd->Fill(nvert1,weight*wgt_pu);
   }
   if(num1Cut>=v+1){
    if(pt_plots)    hpn->Fill(pt1,weight*wgt_pu);
    if(eta_plots)   hpn->Fill(eta1,weight*wgt_pu);
    if(phi_plots)   hpn->Fill(phi1,weight*wgt_pu);
    if(nvert_plots) hpn->Fill(nvert1,weight*wgt_pu);
   }
  }
  //Draw plots with asym err
  hist_ratio(hpn,hpd,varname[v],varcut[v],v);
  delete hpn;
  delete hpd;
 }
 //cand2
 for(int v=1; v<cand2_acc; v++){
  //Prepare histograms
  TH1* hpn = new TH1F("hpn","hpn",const_bin,const_from,const_to);
  hpn->Sumw2();
  TH1* hpd = new TH1F("hpd","hpd",const_bin,const_from,const_to);
  hpd->Sumw2();
  //Fill them
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   weight = wgt_lumi;
   weight = weight*Luminosity/100.;
   if(noLumiNorm) weight = 1.;
   if(noPUcorr) wgt_pu = 1.;
   if(num2Cut>=v && num1Cut>=cand1_acc){
    if(pt_plots)    hpd->Fill(pt2,weight*wgt_pu);
    if(eta_plots)   hpd->Fill(eta2,weight*wgt_pu);
    if(phi_plots)   hpd->Fill(phi2,weight*wgt_pu);
    if(nvert_plots) hpd->Fill(nvert2,weight*wgt_pu);
   }
   if(num2Cut>=v+1 && num1Cut>=cand1_acc){
    if(pt_plots)    hpn->Fill(pt2,weight*wgt_pu);
    if(eta_plots)   hpn->Fill(eta2,weight*wgt_pu);
    if(phi_plots)   hpn->Fill(phi2,weight*wgt_pu);
    if(nvert_plots) hpn->Fill(nvert2,weight*wgt_pu);
   }
  }
  //Draw plots with asym err
  hist_ratio(hpn,hpd,varname[v+cand1_acc],varcut[v+cand1_acc],v+cand1_acc);
  delete hpn;
  delete hpd;
 }
 //Cand1
 for(int v=cand1_acc+1; v<cand1_numTotCut+1; v++){
  //Prepare histograms
  TH1* hpn = new TH1F("hpn","hpn",const_bin,const_from,const_to);
  hpn->Sumw2();
  TH1* hpd = new TH1F("hpd","hpd",const_bin,const_from,const_to);
  hpd->Sumw2();
  //Fill them
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   weight = wgt_lumi;
   weight = weight*Luminosity/100.;
   if(noLumiNorm) weight = 1.;
   if(noPUcorr) wgt_pu = 1.;
   if(num1Cut>=v){
    if(pt_plots)    hpd->Fill(pt1,weight*wgt_pu);
    if(eta_plots)   hpd->Fill(eta1,weight*wgt_pu);
    if(phi_plots)   hpd->Fill(phi1,weight*wgt_pu);
    if(nvert_plots) hpd->Fill(nvert1,weight*wgt_pu);
   }
   if(num1Cut>=v+1){
    if(pt_plots)    hpn->Fill(pt1,weight*wgt_pu);
    if(eta_plots)   hpn->Fill(eta1,weight*wgt_pu);
    if(phi_plots)   hpn->Fill(phi1,weight*wgt_pu);
    if(nvert_plots) hpn->Fill(nvert1,weight*wgt_pu);
   }
  }
  //Draw plots with asym err
  hist_ratio(hpn,hpd,varname[v+cand2_acc],varcut[v+cand2_acc],v+cand2_acc);
  delete hpn;
  delete hpd;
 } 
 //Cand2
 for(int v=cand2_acc+1; v<cand2_numTotCut+1; v++){
  //Prepare histograms
  TH1* hpn = new TH1F("hpn","hpn",const_bin,const_from,const_to);
  hpn->Sumw2();
  TH1* hpd = new TH1F("hpd","hpd",const_bin,const_from,const_to);
  hpd->Sumw2();
  //Fill them
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   weight = wgt_lumi;
   weight = weight*Luminosity/100.;
   if(noLumiNorm) weight = 1.;
   if(noPUcorr) wgt_pu = 1.;
   if(num2Cut>=v && num1Cut==cand1_numTotCut+1){
    if(pt_plots)    hpd->Fill(pt2,weight*wgt_pu);
    if(eta_plots)   hpd->Fill(eta2,weight*wgt_pu);
    if(phi_plots)   hpd->Fill(phi2,weight*wgt_pu);
    if(nvert_plots) hpd->Fill(nvert2,weight*wgt_pu);
   }
   if(num2Cut>=v+1 && num1Cut==cand1_numTotCut+1){
    if(pt_plots)    hpn->Fill(pt2,weight*wgt_pu);
    if(eta_plots)   hpn->Fill(eta2,weight*wgt_pu);
    if(phi_plots)   hpn->Fill(phi2,weight*wgt_pu);
    if(nvert_plots) hpn->Fill(nvert2,weight*wgt_pu);
   }
  }
  //Draw plots with asym err
  hist_ratio(hpn,hpd,varname[v+cand1_numTotCut],varcut[v+cand1_numTotCut],v+cand1_numTotCut);
  delete hpn;
  delete hpd;
 }
}
/////
//   Efficiency for the signal region 
/////
void cutflow_eff_sr(TTree* tree, int nentries){ 
 //Define strings to print
 vector<string> varname(numTotCut+numsigCut); 
 vector<string> varcut(numTotCut+numsigCut);
 initialize_strings(varname,varcut);
 //Call variables 
 int num1Cut, num2Cut, nvert1, nvert2;
 tree->SetBranchAddress("cand1_numCut",&num1Cut);
 tree->SetBranchAddress("cand2_numCut",&num2Cut);
 tree->SetBranchAddress("cand1_nvert",&nvert1);
 tree->SetBranchAddress("cand2_nvert",&nvert2);
 double wgt_lumi, wgt_pu, pt1, eta1, phi1, pt2, eta2, phi2;
 tree->SetBranchAddress("wgt_lumi",&wgt_lumi);
 tree->SetBranchAddress("wgt_pu",&wgt_pu);
 if(gen_pteta){
  tree->SetBranchAddress("cand1_pt_gen",&pt1);
  tree->SetBranchAddress("cand1_eta_gen",&eta1);
  tree->SetBranchAddress("cand1_phi_gen",&phi1);
  tree->SetBranchAddress("cand2_pt_gen",&pt2);
  tree->SetBranchAddress("cand2_eta_gen",&eta2);
  tree->SetBranchAddress("cand2_phi_gen",&phi2);
 }else{
  tree->SetBranchAddress("cand1_pt",&pt1);
  tree->SetBranchAddress("cand1_eta",&eta1);
  tree->SetBranchAddress("cand1_phi",&phi1);
  tree->SetBranchAddress("cand2_pt",&pt2);
  tree->SetBranchAddress("cand2_eta",&eta2);
  tree->SetBranchAddress("cand2_phi",&phi2);
 }
 double charge, cosDphi, met, pZetaMt, pZetaVisMt, cosDphiLMet, massVis;
 tree->SetBranchAddress("charge",&charge);
 tree->SetBranchAddress("cosDphi",&cosDphi);
 tree->SetBranchAddress("met",&met);
 tree->SetBranchAddress("pZetaMt",&pZetaMt);
 tree->SetBranchAddress("pZetaVisMt",&pZetaVisMt);
 tree->SetBranchAddress("cosDphiLMet",&cosDphiLMet);
 tree->SetBranchAddress("massVis",&massVis);
 int pair_nvert, nbjet;
 tree->SetBranchAddress("pair_nvert",&pair_nvert);
 tree->SetBranchAddress("nbjet",&nbjet);
 //Make plots
 double weight = 0;
 //Signal region
 for(int v=0; v<numsigCut; v++){
  //Prepare histograms
  TH1* hpn = new TH1F("hpn","hpn",const_bin,const_from,const_to);
  hpn->Sumw2();
  TH1* hpd = new TH1F("hpd","hpd",const_bin,const_from,const_to);
  hpd->Sumw2();
  //Fill them
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   weight = wgt_lumi;
   weight = weight*Luminosity/100.;
   if(noLumiNorm) weight = 1.;
   if(noPUcorr) wgt_pu = 1.;
   bool cond_d = false;
   //Relative eff
   if(releffsig){
    if(v==0 && pair_nvert!=-9999) cond_d = true;
    if(v==1 && pair_nvert!=-9999 && charge==-1) cond_d = true;
    if(v==2 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95) cond_d = true;
    if(v==3 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20) cond_d = true;
    if(v==4 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0) cond_d = true;  
    if(v==5 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50) cond_d = true;
   }else{
   //N-1 eff
    if(v==0 && pair_nvert!=-9999 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15)
    cond_d = true;
    if(v==1 && pair_nvert!=-9999 && charge==-1 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15)
    cond_d = true;
    if(v==2 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15)
    cond_d = true;
    if(v==3 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15)
    cond_d = true;
    if(v==4 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && cosDphiLMet<0.15)
    cond_d = true;
    if(v==5 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50)
    cond_d = true;
   }
   if(cond_d){
    //which one pt2 or pt1?
    //if(pt_plots)    hpd->Fill(pt2,weight*wgt_pu);
    //if(eta_plots)   hpd->Fill(eta2,weight*wgt_pu);
    //if(phi_plots)   hpd->Fill(phi2,weight*wgt_pu);
    if(nvert_plots)   hpd->Fill(pair_nvert,weight*wgt_pu);
    if(massvis_plots) hpd->Fill(massVis,weight*wgt_pu);
   }
   bool cond_n = false;
   //Relative eff 
   if(releffsig){
    if(v==0 && pair_nvert!=-9999 && charge==-1) cond_n = true;
    if(v==1 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95) cond_n = true;
    if(v==2 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20) cond_n = true;
    if(v==3 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0) cond_n = true;
    if(v==4 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50) cond_n = true;  
    if(v==5 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15) cond_n = true;
   }else{
   //N-1 eff
    if(pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15)
    cond_n = true;
   }
   if(cond_n){
    //which one pt2 or pt1?
    //if(pt_plots)    hpn->Fill(pt2,weight*wgt_pu);
    //if(eta_plots)   hpn->Fill(eta2,weight*wgt_pu);
    //if(phi_plots)   hpn->Fill(phi2,weight*wgt_pu);
    if(nvert_plots)   hpn->Fill(pair_nvert,weight*wgt_pu);
    if(massvis_plots) hpn->Fill(massVis,weight*wgt_pu);
   }
  }
  //Draw plots with asym err
  double eff_Nm1 = hpn->Integral()/hpd->Integral();
  double enden   = hpd->Integral();
  cout<<setiosflags(ios::fixed)<<setprecision(2);  
  cout<<setw(50)<<varname[v+numTotCut]<<setw(10)<<eff_Nm1*100<<setw(4)<<"$\\pm$"<<setw(4)<<sqrt(eff_Nm1*(1-eff_Nm1)/enden)*100<<endl;
  hist_ratio(hpn,hpd,varname[v+numTotCut],varcut[v+numTotCut],v+numTotCut);
  delete hpn;
  delete hpd;
 }
}
/////
//   Efficiency for the special selections
/////
void cutflow_eff_spec(TTree* tree, int nentries){ 
 //Define strings to print
 vector<string> varname(100); 
 vector<string> varcut(100);
 initialize_strings2(varname,varcut);
 //Call variables 
 int num1Cut, num2Cut, nvert1, nvert2;
 tree->SetBranchAddress("cand1_numCut",&num1Cut);
 tree->SetBranchAddress("cand2_numCut",&num2Cut);
 tree->SetBranchAddress("cand1_nvert",&nvert1);
 tree->SetBranchAddress("cand2_nvert",&nvert2);
 double wgt_lumi, wgt_pu, pt1, eta1, phi1, pt2, eta2, phi2;
 tree->SetBranchAddress("wgt_lumi",&wgt_lumi);
 tree->SetBranchAddress("wgt_pu",&wgt_pu);
 if(gen_pteta){
  tree->SetBranchAddress("cand1_pt_gen",&pt1);
  tree->SetBranchAddress("cand1_eta_gen",&eta1);
  tree->SetBranchAddress("cand1_phi_gen",&phi1);
  tree->SetBranchAddress("cand2_pt_gen",&pt2);
  tree->SetBranchAddress("cand2_eta_gen",&eta2);
  tree->SetBranchAddress("cand2_phi_gen",&phi2);
 }else{
  tree->SetBranchAddress("cand1_pt",&pt1);
  tree->SetBranchAddress("cand1_eta",&eta1);
  tree->SetBranchAddress("cand1_phi",&phi1);
  tree->SetBranchAddress("cand2_pt",&pt2);
  tree->SetBranchAddress("cand2_eta",&eta2);
  tree->SetBranchAddress("cand2_phi",&phi2);
 }
 double charge, cosDphi, met, pZetaMt, pZetaVisMt, cosDphiLMet, massVis, muiso;
 tree->SetBranchAddress("charge",&charge);
 tree->SetBranchAddress("cosDphi",&cosDphi);
 tree->SetBranchAddress("met",&met);
 tree->SetBranchAddress("pZetaMt",&pZetaMt);
 tree->SetBranchAddress("pZetaVisMt",&pZetaVisMt);
 tree->SetBranchAddress("cosDphiLMet",&cosDphiLMet);
 tree->SetBranchAddress("massVis",&massVis);
 tree->SetBranchAddress("muiso",&muiso);
 int pair_nvert, nbjet, isHsuEY, isshowerY, isecalisoY, isTrkPtIsoElY;
 tree->SetBranchAddress("pair_nvert",&pair_nvert);
 tree->SetBranchAddress("nbjet",&nbjet);
 tree->SetBranchAddress("isHsuEY",&isHsuEY);
 tree->SetBranchAddress("isshowerY",&isshowerY);
 tree->SetBranchAddress("isecalisoY",&isecalisoY);
 tree->SetBranchAddress("isTrkPtIsoElY",&isTrkPtIsoElY); 
 //Make plots
 double weight = 0;
 //Special selections
 for(int v=numspeceffini; v<numspecefffin; v++){
  //Prepare histograms
  TH1* hpn = new TH1F("hpn","hpn",const_bin,const_from,const_to);
  hpn->Sumw2();
  TH1* hpd = new TH1F("hpd","hpd",const_bin,const_from,const_to);
  hpd->Sumw2();
  //Fill them
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   weight = wgt_lumi;
   weight = weight*Luminosity/100.;
   if(noLumiNorm) weight = 1.;
   if(noPUcorr) wgt_pu = 1.;
   bool cond_d = false;
   if(v==0 && num1Cut==cand1_numTotCut+1) cond_d = true; //Overall ele  eff w.r.t. ele  acceptance cuts and muon     selection   
   if(v==1 && num2Cut==cand2_numTotCut+1) cond_d = true; //Overall muon eff w.r.t. muon acceptance cuts and electron selection
   if(v==2 && num2Cut==cand2_numTotCut+1) cond_d = true; //Overall muon, but iso eff w.r.t. muon acceptance cuts and electron selection
   if(v==3 && pair_nvert!=-9999) cond_d = true; //Mu iso eff w.r.t electron+muon (Need Sel_emu_SSMToTauTau1500_nomuiso.root) 
   if(v==4 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15) cond_d = true; //N-1 Mu iso  
   if(v==5 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15
      && isshowerY==1 && isHsuEY==1 && isTrkPtIsoElY==1) cond_d = true; //N-1 e Cal Iso
   if(v==6 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15
      && isshowerY==1 && isHsuEY==1 && isecalisoY==1) cond_d = true; //N-1 e TrkPt Iso
   if(v==7 && pair_nvert!=-9999) cond_d = true; //Overall topological cuts
   if(cond_d){
    if(pt_plots    && v==0) hpd->Fill(pt2,weight*wgt_pu);
    if(eta_plots   && v==0) hpd->Fill(eta2,weight*wgt_pu);
    if(phi_plots   && v==0) hpd->Fill(phi2,weight*wgt_pu);
    if(nvert_plots && v==0) hpd->Fill(nvert2,weight*wgt_pu);
    if(pt_plots    && 0<v && v<numspecefffin) hpd->Fill(pt1,weight*wgt_pu);
    if(eta_plots   && 0<v && v<numspecefffin) hpd->Fill(eta1,weight*wgt_pu);
    if(phi_plots   && 0<v && v<numspecefffin) hpd->Fill(phi1,weight*wgt_pu);
    if(nvert_plots && 0<v && v<numspecefffin) hpd->Fill(nvert1,weight*wgt_pu);
    if(massvis_plots) hpd->Fill(massVis,weight*wgt_pu);
   }
   bool cond_n = false;
   if(v==0 && pair_nvert!=-9999) cond_n = true;
   if(v==1 && pair_nvert!=-9999) cond_n = true;
   if(v==2 && num2Cut==cand2_numTotCut+1 && num1Cut>=cand1_numTotCut) cond_n = true;
   if(v==3 && pair_nvert!=-9999 && muiso<0.12) cond_n = true;
   if(v==4 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15 && muiso<0.12) cond_n = true; //N-1 Mu iso
   if(v==5 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15
      && isshowerY==1 && isHsuEY==1 && isecalisoY==1 && isTrkPtIsoElY==1) cond_n = true; //N-1 e Cal Iso
   if(v==6 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15
      && isshowerY==1 && isHsuEY==1 && isecalisoY==1 && isTrkPtIsoElY==1) cond_n = true; //N-1 e TrkPt Iso
   if(v==7 && pair_nvert!=-9999 && charge==-1 && cosDphi<-0.95 && met>=20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15) cond_n = true; //Overall topological cuts
   if(cond_n){
    if(pt_plots    && v==0) hpn->Fill(pt2,weight*wgt_pu);
    if(eta_plots   && v==0) hpn->Fill(eta2,weight*wgt_pu);
    if(phi_plots   && v==0) hpn->Fill(phi2,weight*wgt_pu);
    if(nvert_plots && v==0) hpn->Fill(nvert2,weight*wgt_pu);
    if(pt_plots    && 0<v && v<numspecefffin) hpn->Fill(pt1,weight*wgt_pu);
    if(eta_plots   && 0<v && v<numspecefffin) hpn->Fill(eta1,weight*wgt_pu);
    if(phi_plots   && 0<v && v<numspecefffin) hpn->Fill(phi1,weight*wgt_pu);
    if(nvert_plots && 0<v && v<numspecefffin) hpn->Fill(nvert1,weight*wgt_pu);
    if(massvis_plots) hpn->Fill(massVis,weight*wgt_pu);
   }
  }
  //Draw plots with asym err
  cout<<hpn->Integral()<<" "<<hpd->Integral()<<endl;
  double eff_Nm1 = hpn->Integral()/hpd->Integral();
  double enden   = hpd->Integral();
  cout<<setiosflags(ios::fixed)<<setprecision(2);
  cout<<setw(50)<<varname[v]<<setw(10)<<eff_Nm1*100<<setw(4)<<"$\\pm$"<<setw(4)<<sqrt(eff_Nm1*(1-eff_Nm1)/enden)*100<<endl;
  hist_ratio(hpn,hpd,varname[v],varcut[v],v);
  delete hpn;
  delete hpd;
 }
}
/////
//   Initialize the strings 
/////
void initialize_strings(vector<string> &varname, vector<string> &varcut){
 //Acceptance
 varname[0] = "#mu Pt"; varname[1] = "#mu |#eta| "; varname[2] = "#mu IsGlobal";
 varname[3] = "e Et"; varname[4] = "e |#eta_{SC}| "; varname[5] = "e IsEcalDriven";
 varname[6] = "#Delta R(e,#mu) ";
 varcut[0] = ">20 GeV"; varcut[1] = "<2.1"; varcut[2] = "";
 varcut[3] = ">20 GeV"; varcut[4] = "region"; varcut[5] = "";
 varcut[6] = ">0.3";
 //Muon
 varname[7] = "#mu |d_{xy}| ";  varname[8] = "#mu |d_{z}| "; varname[9] = "#mu #ChHits";
 varname[10] = "#mu # PixelHits";
 varname[11] = "#mu #MatchStat"; varname[12] = "#mu #TrkLayMeas";
 varname[13] = "#mu dpt/pt ";  varname[14] = "#mu Isol";
 varcut[7] = "<0.2"; varcut[8] = "<0.5";  varcut[9] = ">0"; varcut[10] = ">0";
 varcut[11] = ">1"; varcut[12] = ">8"; varcut[13] = "<0.3"; varcut[14] = "<0.12";
 //Electron
 varname[15] = "e #Delta#eta_{in} ";  varname[16] = "e #Delta#phi_{in} ";  varname[17] = "e H/E";
 varname[18] = "e (#sigma_{i#etai#eta}";
 varname[19] = "e ShowerShape ";  varname[20] = "e (EM+Had)Iso ";  varname[21] = "e TrkIsoTrkPt ";
 varname[22] = "e #LostHits "; varname[23] = "e |d_{xy}| ";
 varcut[15] = "<0.005 (<0.007)"; varcut[16] = "<0.06"; varcut[17] = "<0.05";
 varcut[18] = "<0.03)";
 varcut[19] = ""; varcut[20] = " "; varcut[21] = "<5";
 varcut[22] = "#leq1"; varcut[23] = "<0.02 (<0.05)";
 //Signal region
 varname[24] = "Ch_{#mu}*Ch_{e}"; varname[25] = "cos#Delta#phi(e,#mu)"; varname[26] = "E_{T}^{miss}";
 varname[27] = "#b-jet"; varname[28] = "CDF-#zeta"; varname[29] = "cos#Delta#phi(l_{lead},E_{T}^{miss})"; 
 varcut[24]  = "-1"; varcut[25] = "<-0.95"; varcut[26] = ">20";
 varcut[27] = "=0"; varcut[28]  = ""; varcut[29] = "<0.15"; 
}
/////
//   Initialize the strings 
/////
void initialize_strings2(vector<string> &varname, vector<string> &varcut){
 //Acceptance
 varname[0] = "HEEP e";  varname[1] = "High pT #mu + isolation"; varname[2] = "High pT #mu"; 
 varname[3] = "#mu Iso before SR"; varname[4] = "N-1 #mu Iso"; varname[5] = "N-1 e $Cal Iso$"; varname[6] = "N-1 e TrkIsoTrkPt";
 varname[7] = "Topol sel w.r.t obj sel";
 varcut[0] = ""; varcut[1] = ""; varcut[2] = "";
 varcut[3] = ""; varcut[4] = ""; varcut[5] = ""; varcut[6] = "";
 varcut[7] = "";
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
 if(sigreg && !releffsig){
  titleYaxis = "#epsilon(N-1 "+varname+varcut+")";//"+(string) bin_size_c;
 }else{
  titleYaxis = "#epsilon("+varname+varcut+")";//"+(string) bin_size_c;
 }
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
