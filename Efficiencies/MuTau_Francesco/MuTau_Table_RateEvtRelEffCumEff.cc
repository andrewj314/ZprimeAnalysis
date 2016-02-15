/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects and events 
   It prints a table following the cutFlow for rateEvt, relEff, and cumEff 

Need to specify
1. See Declare constants

To Do
- Change weight = weight*Luminosity/100.; to properly normalise to luminosity
- Include neg gen weights for aMC@NLO samples
- When printing the event with what ever weights we are currently using the else of if(noLumiNorm){
  that is not very general
- we are not using varcut right now, need to decide if we wanto to keep it
*/
/////
//   To run: root -l MuTau_Table_RateEvtRelEffCumEff.cc+
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
const string path           = "/home/francescoromeovb/Shared_win_ubu_/Work/RunII_HighMassDiTau/Efficiency/";
const string rootpla        = "SusyHTT.root";
//Corrections
const double Luminosity     = 19600; //pb^-1
const bool noLumiNorm       = true; //true means NO luminosity normalization done
const bool noPUcorr         = true; //true means NO PU corr done
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
TFile* Call_TFile();
void cutflow_eff(TTree* tree, int nentries);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
void ini_print(int nentries);
void measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight);
/////
//   Main function
/////
void MuTau_Table_RateEvtRelEffCumEff(){
 TFile* f = Call_TFile();
 cout << "File called\n";
 TTree* tree; f->GetObject("BOOM",tree);
 cout << "Tree initialized\n";
 int nentries = tree->GetEntries(); 
 cout << "nentries = " << nentries << "\n";
 cutflow_eff(tree, nentries);
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(){
 string file_name = path+rootpla;
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Efficiency for the object
/////
void cutflow_eff(TTree* tree, int nentries){ 
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
 //Weights are usefull for considering gen weights (see aMC@NLO samples), pu, lumi...
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
 measure_efficiencies(nentries, cut_eff, varname, varcut, weight);
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
//   Print initial strings
/////
void ini_print(int nentries){
 //def\tablepagesize{\fontsize{7.5pt}{6pt}\selectfont} 
 cout<<"\\vskip-5pt"<<endl;
 cout<<"{\\tablepagesize"<<endl;
 cout<<"\\begin{table}[htdp]"<<endl;
 cout<<"\\begin{center}"<<endl;
 cout<<"\\begin{tabular}{|c|c|c|c|}"<<endl;
 cout<<"\\hline"<<endl;
 cout<<"\\rowcolor{black} {\\color{white!90}Selection} & {\\color{white!90}Evt} & {\\color{white!90}Rel Eff (\\%)} & {\\color{white!90}Cum Eff (\\%)} \\\\"<<endl;
 cout<<"{\\cellcolor{white!90}"<<"NoCuts} & "<<"{\\cellcolor{white!90}"<<nentries<<"} & {\\cellcolor{white!90} 100$\\pm$0} & {\\cellcolor{white!90}100$\\pm$0} \\\\"<<endl;
}
/////
//   Measure efficiencies of cuts
/////
void measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight){
 cout << "Measuring Efficiencies:\n";
 cout << "nEvents = " << nentries << "\n";
 ini_print(nentries);
 double eff_Rel, err_eff_Rel, eff_Cum, err_eff_Cum;
 cout<<setiosflags(ios::fixed)<<setprecision(SETPRECISION);
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
    //cout << varname[i] << varcut[i] << " surviving events: " << cut_eff[i] << " eff_Rel: " << eff_Rel*100 << " +/- " << err_eff_Rel*100 << " eff_Cum: " << eff_Cum*100 << " +/- " << err_eff_Cum*100 <<"\n";  
    cout<<"{\\cellcolor{white!90}"<<varname[i].c_str()<<varcut[i].c_str()<<"} & {\\cellcolor{white!90}"<<cut_eff[i]<<"} & {\\cellcolor{white!90}"<<eff_Rel*100<<"$\\pm$"<<err_eff_Rel*100<<"} & {\\cellcolor{white!90}"<<eff_Cum*100<<"$\\pm$"<<err_eff_Cum*100<<"}\\\\"<<endl;
  }else{
    //cout << varname[i] << varcut[i] << " surviving events: " << cut_eff[i]*weight << " eff_Rel: " << eff_Rel*100 << " +/- " << err_eff_Rel*100 << " eff_Cum: " << eff_Cum*100 << " +/- " << err_eff_Cum*100 <<"\n";   
    cout<<"{\\cellcolor{white!90}"<<varname[i].c_str()<<varcut[i].c_str()<<"} & {\\cellcolor{white!90}"<<cut_eff[i]*weight<<"} & {\\cellcolor{white!90}"<<eff_Rel*100<<"$\\pm$"<<err_eff_Rel*100<<"} & {\\cellcolor{white!90}"<<eff_Cum*100<<"$\\pm$"<<err_eff_Cum*100<<"}\\\\"<<endl;
  }
 }
 cout<<"\\hline"<<endl;
 cout<<"\\end{tabular}"<<endl;
 cout<<"\\end{center}"<<endl;
 cout<<"\\end{table}"<<endl;
 cout<<"}"<<endl;
}
  //Tau_pt->size() > 1 ? cout << "Tau 1 pt: " << Tau_pt->at(0) << " Tau 2 pt: " << Tau_pt->at(1) << "\n" : cout << "No tau_h pair for event " << i << "\n";
  //cout << "Muon Pt->size = " << Muon_pt->size() << "\n";
  //cout << "Muon eta size = " << Muon_eta->size() << "\n";
  //cout << "electron pt size = " << patElectron_pt->size() << "\n";
  //cout << "electron eta size = " << patElectron_eta->size() << "\n";
  //cout << "Muon Pt: " << Muon_pt->at(0) << "\n";
  //cout << "Muon eta: " << Muon_eta->at(0) << "\n";
