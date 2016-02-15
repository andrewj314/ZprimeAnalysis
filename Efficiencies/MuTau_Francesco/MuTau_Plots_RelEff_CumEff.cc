/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects and events 
   It draws a plot following the cutFlow 

Need to specify
1. See Declare constants
*/
/////
//   To run: root -l MuTau_Plots_RelEff_CumEff.cc+
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
using namespace ROOT::Math;
/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
TH1F* cutflow_eff(TTree* tree, int nentries, int pos, string rootpla);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
TH1F* measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight, int pos, string rootplas);
void setTDRStyle();
//Path and root file
const string path           = "/uscms_data/d3/andrewj/CMSSW_7_4_7/src/Efficiencies/MuTau_Francesco/";//"/home/francescoromeovb/Shared_win_ubu_/Work/RunII_HighMassDiTau/Efficiency/";
const char *samples[]       = {"OutTree"};// {"SusyHTT"};
//Corrections
const double Luminosity     = 19600; //pb^-1
const bool noLumiNorm       = true; //true means NO luminosity normalization done
const bool noPUcorr         = true; //true means NO PU corr done
//Choose among relEff, cumEff, or rateEvt 
const bool relEff           = true;
const bool cumEff           = false;
//Acceptance and object selection
const double cand1_acc[2]   = {18, 2.1};
const double cand2_acc[2]   = {20, 2.3};
const int numtot_cuts       = 15;
//Signal selection
//charge==-1 && cosDphi<-0.95 && met>=30 && pZetaMt-3.051*pZetaVisMt>-50 && cosDphiLMet<0.2 && nbjet==0
//const int sr_cuts           = 6;
//double sigreg_cuts[sr_cuts] = {-1,-0.95,20,-50,0.2,0};
const double SETPRECISION = 3;
/////
//   Main function
/////
void MuTau_Plots_RelEff_CumEff() {
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
  TFile* f = Call_TFile(rootplas[i]);
  TTree* tree; f->GetObject("TNT/BOOM;3",tree);
  int nentries = tree->GetEntries(); 
  //Take plot
  TH1F* hRelEff = new TH1F(rootplas[i].c_str(),rootplas[i].c_str(),numtot_cuts,0,numtot_cuts);
  hRelEff = cutflow_eff(tree,nentries,i,rootplas[i]);
  if(i==0){
   hRelEff->Draw("PE1");
  }else{
   hRelEff->Draw("PE1same");
  }
  leg->AddEntry(hRelEff,rootplas[i].c_str(),"LP");
  if(i+1!=rootplas.size()){nameFile = nameFile+rootplas[i]+"_";}else{nameFile = nameFile+rootplas[i];} 
 }
 leg->Draw(); 
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
 TH1F *hCurrRelEff = new TH1F(rootpla.c_str(),rootpla.c_str(),numtot_cuts,0,numtot_cuts);
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
  //Get values
  double wgt_pu   = 1.;
  double wgt_lumi = 1.;
  if(noPUcorr) wgt_pu = 1.;
  if(Muon_pt->size() > 0 && Tau_pt->size() > 0){
    cut_eff[0] += wgt_pu;
    if(Muon_pt->at(0) > cand1_acc[0]){
      cut_eff[1] += wgt_pu;
      if(fabs(Muon_eta->at(0)) < cand1_acc[1]){
        cut_eff[2] += wgt_pu;
        if(Tau_pt->at(0) > cand2_acc[0]){
          cut_eff[3] += wgt_pu;
          if(fabs(Tau_eta->at(0)) < cand2_acc[1]){
            cut_eff[4] += wgt_pu;
            if(Muon_dxy->at(0) < 0.045 && Muon_dz->at(0) < 0.2){
              cut_eff[5] += wgt_pu;
              if(Muon_isMediumMuon->at(0) > 0){
                cut_eff[6] += wgt_pu;
                if(Tau_decayModeFindingNewDMs->at(0) > 0.5){
                  cut_eff[7] += wgt_pu;
                  if(Tau_leadChargedCand_dz->at(0) < 0.2){
                    cut_eff[8] += wgt_pu;
                    PtEtaPhiEVector PtEtaPhiEMu0(Muon_pt->at(0),Muon_eta->at(0),Muon_phi->at(0),Muon_energy->at(0));
                    PtEtaPhiEVector PtEtaPhiETau0(Tau_pt->at(0),Tau_eta->at(0),Tau_phi->at(0),Tau_energy->at(0));
		    double Mu0Tau0_DeltaR = VectorUtil::DeltaR(PtEtaPhiEMu0,PtEtaPhiETau0);
		    if(Mu0Tau0_DeltaR > 0.5){
		      cut_eff[9] += wgt_pu;
		      double Muon_Iso;
		      if(Muon_isoCharged->size() > 0 && Muon_isoNeutralHadron->size() > 0 && Muon_isoPhoton->size() > 0 && Muon_isoPU->size() > 0){
                        Muon_Iso = (Muon_isoCharged->at(0) + max(Muon_isoNeutralHadron->at(0) + Muon_isoPhoton->at(0) - 0.5*Muon_isoPU->at(0), 0.0))/Muon_pt->at(0);
		      }else{
                        Muon_Iso = 10;
                      }
		      if(Muon_Iso < 0.1){
		        cut_eff[10] += wgt_pu;
		        if(Tau_againstElectronMVALooseMVA5->at(0) > 0.5){
			  cut_eff[11] += wgt_pu;
			  if(Tau_againstMuonTight3->at(0) > 0.5){
			    cut_eff[12] += wgt_pu;
			    if(Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(0) < 1.5){
			      cut_eff[13] += wgt_pu;
			      //Di-Muon Veto
			      if(Muon_pt->size() < 2){
  				cut_eff[14] += wgt_pu;
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
			          cut_eff[14] += wgt_pu;
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
