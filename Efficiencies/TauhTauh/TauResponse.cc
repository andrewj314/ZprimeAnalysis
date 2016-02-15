/**
This Macro   
1. Plots pair of variables in Profile in order the get the Tau response 

Need to specify
1. See Declare Constants
*/
/////
//   To run: root -l TauResponse.cc+
/////
/////
//   Prepare Root and Roofit
/////
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
using namespace std;
/////
//   Declare constants
/////
//Path - samples  
const string path      = "/uscms_data/d3/andrewj/CMSSW_7_4_15/src/ZprimeSamples/ZprimeSSM_M_2000/";
const char *samples[]  = {"OutTree_1"}; 
//For Plots
const bool   save_plots  = true;
const int    numVar           = 2;
const int    bin[numVar]      = {15, 15};
const double inRange[numVar]  = {0, 0};
const double endRange[numVar] = {1500, 0.1};
//const char  *titleaxes[] = {"Gen p_{T}^{vis}","(Reco p_{T}-Gen p_{T}^{vis})/Gen p_{T}^{vis}"};
const char  *titleaxes[] = {"Gen pT_{#tau^{Vis}_{had}}","(Reco pT_{#tau_{had}}-Gen pT_{#tau^{Vis}_{had}})/Gen pT_{#tau^{Vis}_{had}}"};
/////
//   Declare functions 
/////
TFile* Call_TFile(string rootpla);
void setTDRStyle();
/////
//   Main function
/////
void TauResponse(){
 setTDRStyle();
 vector<string> varTitleaxes(titleaxes, titleaxes + sizeof(titleaxes)/sizeof(titleaxes[0]));
 //Loop over samples
 vector<string> rootplas(samples, samples + sizeof(samples)/sizeof(samples[0]));
 for(uint i=0; i<rootplas.size(); i++){
  //Call tree and variables 
  TFile* f = Call_TFile(rootplas[i]);
  TTree* tree;
  f->GetObject("TNT/BOOM",tree);
  vector<double> *Gen_pt, *Tau_pt;
  vector<int> *Tau_decayModeFindingNewDMs, 
  TBranch *b_Gen_pt;   //!
  TBranch *b_Tau_pt;   //!


  tree->SetBranchAddress("Gen_pt",&Gen_pt,&b_Gen_pt);
  tree->SetBranchAddress("Tau_pt",&Tau_pt,&b_Tau_pt);
  //For plots
  TCanvas* c1 = new TCanvas(rootplas[i].c_str(),rootplas[i].c_str(),100,100,500,500);
  //Declare histo
  TProfile* h_var = new TProfile("hprof","hprof",bin[0],inRange[0],endRange[0]);
  h_var->SetTitle(0);
  h_var->SetMarkerStyle(20);
  h_var->SetMarkerColor(1);
  h_var->GetXaxis()->SetTitle(varTitleaxes[0].c_str());
  h_var->GetYaxis()->SetTitle(varTitleaxes[1].c_str());
  h_var->SetMaximum(endRange[1]);
  //Fill it
  double wratio = 0;
  for(int j=0; j<tree->GetEntries(); j++){
   tree->GetEntry(j);
   wratio = (recotaupt-genvistaupt)/genvistaupt;
   h_var->Fill(genvistaupt,wratio,1.);
  }   
  //Plot it 
  h_var->Draw("");
  //Line
  TLine* line1 = new TLine(0,0,1500,0);
  line1->SetLineColor(1);
  line1->SetLineWidth(3);
  line1->Draw("same");
  //Save plot
  string namefile = "TauResponse.pdf";
  if(save_plots) c1->SaveAs(namefile.c_str());
  delete tree;
 }
}
/////
//   Call TFile to be read
/////
TFile* Call_TFile(string rootpla) {
 string file_name = path+rootpla+".root";
 TFile* f = new TFile(file_name.c_str(),"update");
 return f;
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
  tdrStyle->SetOptFit(1);
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
  tdrStyle->SetOptStat(""); // To display the mean and RMS:   SetOptStat("mr");
  //tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatColor(kGray);
  tdrStyle->SetStatFont(42);

  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(0);
  tdrStyle->SetStatX(0.9); //Starting position on X axis
  tdrStyle->SetStatY(1.); //Starting position on Y axis
  tdrStyle->SetStatFontSize(0.025); //Vertical Size
  tdrStyle->SetStatW(0.15); //Horizontal size 
  // tdrStyle->SetStatStyle(Style_t style = 1001);

  // Margins:
  tdrStyle->SetPadTopMargin(0.025);
  tdrStyle->SetPadBottomMargin(0.125);
  tdrStyle->SetPadLeftMargin(0.165);
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
  tdrStyle->SetTitleYOffset(1.1);
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
  //tdrStyle->SetPalette(0,0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}
