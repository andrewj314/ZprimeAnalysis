/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects and pairs 
   It prints a table following the cutFlow for rateEvt or relEff or cumEff, for all the files specified

Need to specify
1. See Declare constants
*/
/////
//   To run: root -l Efficiency_cutflow_relEff_cumEff_rateEvt_table.cc+
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
//Path - rootpla; Luminosity (for normalization of MC to data)
const string path     = "/uscms_data/d3/andrewj/CMSSW_7_4_0_patch1/src/Efficiencies/"; //"/afs/cern.ch/user/f/fromeo/public/EXO/Efficiency/";
const char *samples[] = {"SSMToTauTau1250"};
const string channel  = "emu"; const double Luminosity = 19703.225;//19779.362;//19600; //pb^-1
const bool noLumiNorm = false; //true means NO luminosity normalization done
const bool noPUcorr   = false; //true means NO PU corr done
const bool relEff     = false;
const bool cumEff     = false;
const bool rateEvt    = true; 
const int cand1_acc       = 3;
const int cand2_acc       = 3;
const int cand1_numTotCut = 11; //It is different (-1) from what you could see in ZpmSelection.cc, but still correct
const int cand2_numTotCut = 12; //It is different (-1) from what you could see in ZpmSelection.cc, but still correct
const int pair_numTotCut  = 6;
const int numTotCut       = 30; //dr+cand1_numTotCut+cand2_numTotCut+pair_numTotCut 
//charge==-1 && cosDphi<-0.95 && met>20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15
double cut_sig[pair_numTotCut] = {-1,-0.95,20,0,-50,0.15};
//Other
const double divide_num = 1; //Divide the printed number for divide_num (if e.g. printed number are too high, as it can be for rateEvt of data)
const double SETPRECISION = 3;
/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
void cutflow_eff(TTree* tree, int nentries, int pos, double sample_cuts[][numTotCut], double sample_cuts_err[][numTotCut]);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
void measure_efficiencies(int nentries, double cut_eff[], double weight, int pos, double sample_cuts[][numTotCut], double sample_cuts_err[][numTotCut]);
void print_table(vector<string> &rootplas, double sample_cuts[][numTotCut], double sample_cuts_err[][numTotCut]);
/////
//   Main function
/////
void Efficiency_cutflow_relEff_cumEff_rateEvt_table() {
 //Run over all samples 
 vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
 const uint rootplas_size = rootplas.size();
 //Matrix with values of rel eff
 double sample_cuts[rootplas_size][numTotCut];
 double sample_cuts_err[rootplas_size][numTotCut];
 for(uint i=0; i<rootplas_size; i++){
  //Call tree  
  TFile* f = Call_TFile(rootplas[i]);
  TTree* tree; f->GetObject("Selection/tree",tree);
  int nentries = tree->GetEntries(); 
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
 string file_name = path+"Sel_"+channel+"_"+rootpla+".root";
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Efficiency for the object
/////
void cutflow_eff(TTree* tree, int nentries, int pos, double sample_cuts[][numTotCut], double sample_cuts_err[][numTotCut]){ 
 //Call variables 
 int num1Cut, num2Cut, njet, nbjet;
 tree->SetBranchAddress("cand1_numCut",&num1Cut);
 tree->SetBranchAddress("cand2_numCut",&num2Cut);
 tree->SetBranchAddress("njet",&njet);
 tree->SetBranchAddress("nbjet",&nbjet);
 double wgt_lumi, wgt_pu, charge, dR, cosDphi, cosDphiLMet, pZetaVisMt, pZetaMt, pZeta, massT, massT2, massVis, met, jetHPT_pt, sumJetEt, EcalIso, HadrDepth1, EtEl, Kt6JetsRho, TrkPtIsoEl, Losthit, Gsfdxy, muiso;
 tree->SetBranchAddress("wgt_lumi",&wgt_lumi);
 tree->SetBranchAddress("wgt_pu",&wgt_pu);
 tree->SetBranchAddress("charge",&charge);
 tree->SetBranchAddress("dR",&dR);
 tree->SetBranchAddress("cosDphi",&cosDphi);
 tree->SetBranchAddress("cosDphiLMet",&cosDphiLMet);
 tree->SetBranchAddress("pZetaVisMt",&pZetaVisMt);
 tree->SetBranchAddress("pZetaMt",&pZetaMt);
 tree->SetBranchAddress("pZeta",&pZeta);
 tree->SetBranchAddress("massT",&massT);
 tree->SetBranchAddress("massT2",&massT2);
 tree->SetBranchAddress("massVis",&massVis);
 tree->SetBranchAddress("met",&met);
 tree->SetBranchAddress("jetHPT_pt",&jetHPT_pt);
 tree->SetBranchAddress("sumJetEt",&sumJetEt);
 tree->SetBranchAddress("EcalIso",&EcalIso);
 tree->SetBranchAddress("HadrDepth1",&HadrDepth1);
 tree->SetBranchAddress("EtEl",&EtEl);
 tree->SetBranchAddress("Kt6JetsRho",&Kt6JetsRho);
 tree->SetBranchAddress("TrkPtIsoEl",&TrkPtIsoEl);
 tree->SetBranchAddress("Losthit",&Losthit);
 tree->SetBranchAddress("Gsfdxy",&Gsfdxy);
 tree->SetBranchAddress("muiso",&muiso);
 //Take entries for each variable of cut
 double cut_eff[numTotCut];
 for(int iniz=0; iniz<numTotCut; iniz++) cut_eff[iniz] = 0;
 double weight = 0;
 //All entries
 for(int i=0; i<nentries; i++){
  tree->GetEntry(i);
  if(noPUcorr) wgt_pu = 1.;
  //Acceptance
  for(int j=0; j<cand1_acc; j++) if(j<num1Cut) cut_eff[j] += wgt_pu;//0,1,2
  for(int j=0; j<cand2_acc; j++) if(j<num2Cut && num1Cut>=cand1_acc) cut_eff[j+cand1_acc] += wgt_pu;//3,4,5 
  if(num1Cut>=cand1_acc+1 && num2Cut>=cand2_acc+1) cut_eff[cand1_acc+cand2_acc] += wgt_pu;//6
  //Cand1
  for(int j=cand1_acc+1; j<num1Cut; j++) if(j<num1Cut) cut_eff[j+cand2_acc] += wgt_pu;//7->14 (remember that num1Cut is cand1_numTotCut+1(dR))
  //Cand2
  for(int j=cand2_acc+1; j<num2Cut; j++) if(j<num2Cut && num1Cut==cand1_numTotCut+1) cut_eff[j+cand1_numTotCut] += wgt_pu;//15->23
  //Pair
  //charge==-1 && cosDphi<-0.95 && met>20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15
  if(num1Cut==cand1_numTotCut+1 && num2Cut==cand2_numTotCut+1){
   if(charge==cut_sig[0]){
    cut_eff[cand1_numTotCut+cand2_numTotCut+1] += wgt_pu;//24
    if(cosDphi<cut_sig[1]){
     cut_eff[cand1_numTotCut+cand2_numTotCut+2] += wgt_pu;//25
     if(met>=cut_sig[2]){
      cut_eff[cand1_numTotCut+cand2_numTotCut+3] += wgt_pu;//26
      if(nbjet==cut_sig[3]){
       cut_eff[cand1_numTotCut+cand2_numTotCut+4] += wgt_pu;//27
       if(pZetaMt-3.1*pZetaVisMt>cut_sig[4]){
        cut_eff[cand1_numTotCut+cand2_numTotCut+5] += wgt_pu;//28
        if(cosDphiLMet<cut_sig[5]){
         cut_eff[cand1_numTotCut+cand2_numTotCut+6] += wgt_pu;//29
        }
       }
      }
     }
    }
   }
  }
  weight = wgt_lumi;
 }//End all entries
 weight = weight*Luminosity/100.;
 //Evaluate efficiencies
 measure_efficiencies(nentries, cut_eff, weight, pos, sample_cuts, sample_cuts_err);
}
/////
//   Initialize the strings 
/////
void initialize_strings(vector<string> &varname, vector<string> &varcut){
 //Acceptance
 varname[0] = "$\\mu$ Pt"; varname[1] = "$\\mu$ $|\\eta|$ "; varname[2] = "$\\mu$ IsGlobal";
 varname[3] = "e Et"; varname[4] = "e $|\\eta_{SC}|$ "; varname[5] = "e IsEcalDriven";
 varname[6] = "$\\Delta R(e,\\mu)$ ";
 varcut[0] = "$>$20 GeV"; varcut[1] = "$<$2.1"; varcut[2] = "";
 varcut[3] = "$>$20 GeV"; varcut[4] = "$<$1.442 (1.56$<|\\eta_{SC}|<$2.5)"; varcut[5] = "";
 varcut[6] = "$>$0.3";
 //Muon
 varname[7] = "$\\mu$ $|d_{xy}|$ ";  varname[8] = "$\\mu$ $|d_{z}|$ "; varname[9] = "$\\mu$ Num Chamber Hits";
 varname[10] = "$\\mu$ Num Pixel Hits";
 varname[11] = "$\\mu$ Num Matched Stations"; varname[12] = "$\\mu$ Num Tracker Layers with Meas";
 varname[13] = "$\\mu$ $dpt/pt$ ";  varname[14] = "$\\mu$ Isolation";
 varcut[7] = "$<$0.2"; varcut[8] = "$<$0.5";  varcut[9] = "$>$0"; varcut[10] = "$>$0";
 varcut[11] = "$>$1"; varcut[12] = "$>$5"; varcut[13] = "$<$0.3"; varcut[14] = "$<$0.12";
 //Electron
 varname[15] = "e $\\Delta\\eta_{in}$ ";  varname[16] = "e $\\Delta\\phi_{in}$ ";  varname[17] = "e H/E";
 varname[18] = "e ($\\sigma_{i\\eta i\\eta}";
 varname[19] = "e $E^{2x5}/E^{5x5}>0.94 OR E^{1x5}/^{E5x5} > 0.83$ ";  varname[20] = "e $EM+HadDepth1Iso$ ";  varname[21] = "e TrkIsoTrkPt ";
 varname[22] = "e Num Lost Hits "; varname[23] = "e $|d_{xy}|$ ";
 varcut[15] = "$<$0.005 ($<$0.007)"; varcut[16] = "$<$0.06"; varcut[17] = "$<$0.05";
 varcut[18] = "<$0.03)";
 varcut[19] = ""; varcut[20] = "$<$2+0.03*Et+0.28$\\rho$ ($<$2.5+0.03*(Et-50)+0.28$\\rho$ (Et=50,Et<50; Et=Et,Et>50))"; varcut[21] = "$<$5";
 varcut[22] = "$\\leq$1"; varcut[23] = "$<$0.02 ($<$0.05)";
 //Pair
 //charge==-1 && cosDphi<-0.95 && met>=30 && pZetaMt-3.051*pZetaVisMt>-50 && cosDphiLMet<0.2 && nbjet==0
 varname[24] = "$Ch_{\\mu}*Ch_{e}$ "; varname[25] = "cos$\\Delta\\phi$ "; varname[26] = "${\\not E_T}$ ";
 varname[27] = "nbjet"; varname[28] = "$p\\zeta Mt-3.1*p\\zeta VisMt$ "; varname[29]= "cos$\\Delta\\phi(l,{\\not E_T})$ ";
 varcut[24]  = "=-1"; varcut[25] = "$<$-0.95"; varcut[26] = "$>$20 GeV/c";
 varcut[27] = "=0"; varcut[28]  = "$>$-50"; varcut[29] = "$<$0.15";
}
/////
//   Measure efficiencies of cuts
/////
void measure_efficiencies(int nentries, double cut_eff[], double weight, int pos, double sample_cuts[][numTotCut], double sample_cuts_err[][numTotCut]){
 double eff_Rel, err_eff_Rel, eff_Cum, err_eff_Cum;
 for(int i=0; i<numTotCut; i++){
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
   sample_cuts[pos][i]     = cut_eff[i]*weight/divide_num;
   double err_rateEvt      = sqrt(cut_eff[i])*weight;
   sample_cuts_err[pos][i] = err_rateEvt/divide_num;
  }
 }
}
/////
//   Print table of rel eff for all samples
/////
void print_table(vector<string> &rootplas, double sample_cuts[][numTotCut], double sample_cuts_err[][numTotCut]){
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
 vector<string> varname(numTotCut);
 vector<string> varcut(numTotCut);
 initialize_strings(varname,varcut);
 for(int j=0; j<numTotCut; j++) for(uint i=0; i<rootplas.size(); i++){
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
