/**
This Macro   
1. Studies efficiencies of cuts for the selection of objects and pairs 
   It draws a plot following the cutFlow 

Need to specify
1. See Declare constants
2. Check else if(rootpla=="SSM in case of signal
*/
/////
//   To run: root -l Efficiency_cutflow_relEff_cumEff_plots.cc+
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
const string path     = "/afs/cern.ch/user/f/fromeo/public/EXO/Efficiency/";
const char *samples[] = {"SSMToTauTau1250","TTJet"};
const string channel  = "emu"; const double Luminosity = 19703.225;//19779.362;//19600; //pb^-1
const bool noLumiNorm = true; //true means NO luminosity normalization done
const bool noPUcorr   = true; //true means NO PU corr done
const bool relEff     = true;
const bool cumEff     = false;
const int cand1_acc       = 3;
const int cand2_acc       = 3;
const int cand1_numTotCut = 11; //It is different (-1) from what you could see in ZpmSelection.cc, but still correct
const int cand2_numTotCut = 12; //It is different (-1) from what you could see in ZpmSelection.cc, but still correct
const int pair_numTotCut  = 6;
const int numTotCut       = 30; //dr+cand1_numTotCut+cand2_numTotCut+pair_numTotCut 
//charge==-1 && cosDphi<-0.95 && met>20 && nbjet==0 && pZetaMt-3.1*pZetaVisMt>-50 && cosDphiLMet<0.15
double cut_sig[pair_numTotCut] = {-1,-0.95,20,0,-50,0.15};
/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
TH1F* cutflow_eff(TTree* tree, int nentries, int pos, string rootpla);
void initialize_strings(vector<string> &varname, vector<string> &varcut);
TH1F* measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight, int pos, string rootplas);
void setTDRStyle();
/////
//   Main function
/////
void Efficiency_cutflow_relEff_cumEff_plots() {
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
  TTree* tree; f->GetObject("Selection/tree",tree);
  int nentries = tree->GetEntries(); 
  //Take plot
  TH1F* hRelEff = new TH1F(rootplas[i].c_str(),rootplas[i].c_str(),numTotCut,0,numTotCut);
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
 string file_name = path+"Sel_"+channel+"_"+rootpla+".root";
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
}
/////
//   Efficiency for the object
/////
TH1F* cutflow_eff(TTree* tree, int nentries, int pos, string rootpla){ 
 TH1F *hCurrRelEff = new TH1F(rootpla.c_str(),rootpla.c_str(),numTotCut,0,numTotCut);
 //Define strings to print
 vector<string> varname(numTotCut); 
 vector<string> varcut(numTotCut);
 initialize_strings(varname,varcut);
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
 hCurrRelEff = measure_efficiencies(nentries, cut_eff, varname, varcut, weight, pos, rootpla);
 return hCurrRelEff;
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
 varcut[11] = ">1"; varcut[12] = ">5"; varcut[13] = "<0.3"; varcut[14] = "<0.12";
 //Electron
 varname[15] = "e #Delta#eta_{in} ";  varname[16] = "e #Delta#phi_{in} ";  varname[17] = "e H/E";
 varname[18] = "e (#sigma_{i#etai#eta}"; 
 varname[19] = "e ShowerShape ";  varname[20] = "e CalIso ";  varname[21] = "e TrkIsoTrkPt ";
 varname[22] = "e #LostHits "; varname[23] = "e |d_{xy}| ";
 varcut[15] = "<0.005(7)"; varcut[16] = "<0.06"; varcut[17] = "<0.05";
 varcut[18] = "<0.03)";
 varcut[19] = ""; varcut[20] = " "; varcut[21] = "<5";
 varcut[22] = "#leq1"; varcut[23] = "<0.02(5)";
 //Pair
 //charge==-1 && cosDphi<-0.95 && met>=30 && pZetaMt-3.051*pZetaVisMt>-50 && cosDphiLMet<0.2 && nbjet==0
 varname[24] = "Ch_{#mu}*Ch_{e} "; varname[25] = "cos#Delta#phi "; varname[26] = "met ";
 varname[27] = "nbjet"; varname[28] = "p#zeta "; varname[29] = "cos#Delta#phi(l,met) ";
 varcut[24]  = "=-1";   varcut[25]  = "<-0.95";  varcut[26]  = ">20 GeV/c";
 varcut[27]  = "=0";    varcut[28]  = ">-50";    varcut[29]  = "<0.15";
}
/////
//   Measure efficiencies of cuts
/////
TH1F* measure_efficiencies(int nentries, double cut_eff[], vector<string> &varname, vector<string> &varcut, double weight, int pos, string rootpla){
 TH1F *hCurrRelEff = new TH1F(rootpla.c_str(),rootpla.c_str(),numTotCut,0,numTotCut);
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
 for(int i=0; i<numTotCut; i++){
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
