/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects and events 
   It prints a table following the cutFlow for rateEvt or relEff or cumEff, for all the files specified

Need to specify
1. See Declare constants
*/
/////
//   To run: root -l MuTau_Table_RateEvt_RelEff_CumEff.cc
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
using namespace ROOT::Math;
/////
//   Declare constants
/////
//Path and root file
const string path           = "/uscms_data/d3/andrewj/CMSSW_7_4_7/src/Efficiencies/MuTau_Francesco/";
const char *samples[]       = {"OutTree"};
//Corrections
const double Luminosity     = 19600; //pb^-1
const bool noLumiNorm       = true; //true means NO luminosity normalization done
const bool noPUcorr         = true; //true means NO PU corr done
//Choose among relEff, cumEff, or rateEvt 
const bool relEff           = false;
const bool cumEff           = false;
const bool rateEvt          = true; 
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
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
void cutflow_eff(TTree* tree, int nentries, int pos, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
void measure_efficiencies(int nentries, double cut_eff[], double weight, int pos, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]);
void print_table(vector<string> &rootplas, double sample_cuts[][numtot_cuts], double sample_cuts_err[][numtot_cuts]);
/////
//   Main function
/////
void MuTau_Table_RateEvt_RelEff_CumEff() {
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
  TTree* tree; f->GetObject("TNT/BOOM;3",tree);
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
            if(Muon_dxy->at(0) < 100 /*0.045*/ && Muon_dz->at(0) < 100/*0.2*/){
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
 measure_efficiencies(nentries, cut_eff, weight, pos, sample_cuts, sample_cuts_err);
}
/////
//   Initialize the strings 
/////
void initialize_strings(vector<string> &varname, vector<string> &varcut){
 varname[0] = "$\\geq$ 1 Mu and 1 Tau\t";
 varname[1] = "Mu Pt $>$ 18 \t\t"; varname[2] = "Mu |eta| $<$ 2.1 \t\t";
 varname[3] = "Tau Pt $>$ 20 \t\t\t"; varname[4] = "Tau |eta| $<$ 2.3 \t\t";
 varname[5] = "Mu track dxy $<$ 0.045 \\& dz $<$ 0.2\t"; varname[6] = "Mu MediumID\t\t";
 varname[7] = "Tau DMF new DMs\t"; varname[8] = "Tau dz $<$ 0.2 \t\t\t";
 varname[9] = "MuTau DeltaR $>$ 0.5 \t\t"; varname[10] = "Mu Iso $<$ 0.1\t\t";
 varname[11] = "Tau anti-e MVA5\t"; varname[12] = "Tau anti-mu tight\t";
 varname[13] = "Tau Medium Iso DB 3hits $<$ 1.5"; varname[14] = "Di-mu veto\t\t";
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
