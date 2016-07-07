
void LimitPlotter_overalid_brazil_TauTau() {
  
  //Load the files
  TFile * f1 = new TFile("TauTau.root");

  //Open the trees
  TTree *tree1 = (TTree*)f1->Get("limit");

  Double_t mh1;
  tree1->SetBranchAddress("mh", &mh1);

  Float_t quantileExpected1;
  tree1->SetBranchAddress("quantileExpected", &quantileExpected1);

  // Define the variables and branches
  Double_t limit1;

  tree1->SetBranchAddress("limit", &limit1);

  // Define the Canvas
  TCanvas *mass = new TCanvas("mass", "mass",81,80,500,604);
   mass->Range(-278.9672,-1.221595,3469.144,7.904693);
   mass->SetFillColor(0);
   mass->SetBorderMode(0);
   mass->SetBorderSize(2);
   mass->SetLogy();
   mass->SetTickx(1);
   mass->SetTicky(1);
   mass->SetLeftMargin(0.141129);
   mass->SetRightMargin(0.05846774);
   mass->SetFrameBorderMode(0);
   mass->SetFrameBorderMode(0); 
 
  // Define arrays needed for the limits
  const int val = 6;
  const int val2 = 6;
/*
  Double_t Zprime_signal_Scale_Factors[val] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  Double_t Zprime_signal_Masses[val]     = {500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0};
  Double_t VBG_signal_xSec[val]          = {9330.0, 468.0, 72.3, 17.3, 5.54, 1.2};
  Double_t Zprime_signal_Masses2[val2]   = {500.0, 1000.0, 1500.0, 2000.0, 2500., 3000.0};
  Double_t VBG_signal_xSec2[val2]        = {9330.0, 468.0, 72.3, 17.3, 5.54, 1.29};
  Double_t y_theory_up[val2]             = {9330.0*1.3, 468.0*1.3, 72.3*1.3, 17.3*1.3, 5.54*1.3, 1.29*1.3};
  Double_t y_theory_dn[val2]             = {9330.0*0.7, 468.0*0.7, 72.3*0.7, 17.3*0.7, 5.54*0.7, 1.29*0.7};
*/

  Double_t Zprime_signal_Scale_Factors[val] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  Double_t Zprime_signal_Masses[val]     = {500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0};
  Double_t VBG_signal_xSec[val]          = {9.330, 0.468, 0.0723, 0.0173, 0.00554, 0.00129};
  Double_t Zprime_signal_Masses2[val2]   = {500.0, 1000.0, 1500.0, 2000.0, 2500., 3000.0};
  Double_t VBG_signal_xSec2[val2]        = {9.330*1.3, 0.468*1.3, 0.0723*1.3, 0.0173*1.3, 0.00554*1.3, 0.00129*1.3};
  Double_t y_theory_up[val2]             = {9330.0*1.3, 468.0*1.3, 72.3*1.3, 17.3*1.3, 5.54*1.3, 1.29*1.3};
  Double_t y_theory_dn[val2]             = {9330.0*0.7, 468.0*0.7, 72.3*0.7, 17.3*0.7, 5.54*0.7, 1.29*0.7};

  Double_t Limit_Obs[val]     = {0};
  Double_t Limit_Exp_m2s[val] = {0};
  Double_t Limit_Exp_m1s[val] = {0};
  Double_t Limit_Exp[val]     = {0};
  Double_t Limit_Exp_p1s[val] = {0};
  Double_t Limit_Exp_p2s[val] = {0};

  int counter   = 0;
  int counter_1 = 0;
  int counter_2 = 0;
  int counter_3 = 0;
  int counter_4 = 0;

  // Gey entries
  Long64_t nentries1 = tree1->GetEntriesFast();

  for (Int_t jentry=0; jentry<nentries1;jentry++) {
    tree1->GetEntry(jentry);
    if (limit1 > 0.0001)
        {
         cout <<"quantileExpected1 "<<quantileExpected1<<endl; 
         if (mh1 == 500){
           if(quantileExpected1 < 0.17 && (quantileExpected1 > 0)) Limit_Exp_m1s[0] = limit1*VBG_signal_xSec[0];
           if(quantileExpected1 == 0.5) Limit_Exp[0] = limit1*VBG_signal_xSec[0];
           if(quantileExpected1 > 0.8 && (quantileExpected1 < 0.9)) Limit_Exp_p1s[0] = limit1*VBG_signal_xSec[0];
           if(quantileExpected1 > 0.024 && (quantileExpected1 < 0.027)) Limit_Exp_m2s[0] = limit1*VBG_signal_xSec[0];
           if(quantileExpected1 > 0.974 && (quantileExpected1 < 0.977)) Limit_Exp_p2s[0] = limit1*VBG_signal_xSec[0];
           if(quantileExpected1 == -1) Limit_Obs[0] = limit1*VBG_signal_xSec[0];
         }

         if (mh1 == 1000){
           if(quantileExpected1 < 0.17 && (quantileExpected1 > 0)) Limit_Exp_m1s[1] = limit1*VBG_signal_xSec[1];
           if(quantileExpected1 == 0.5) Limit_Exp[1] = limit1*VBG_signal_xSec[1];
           if(quantileExpected1 > 0.8 && (quantileExpected1 < 0.9)) Limit_Exp_p1s[1] = limit1*VBG_signal_xSec[1];
           if(quantileExpected1 > 0.024 && (quantileExpected1 < 0.026)) Limit_Exp_m2s[1] = limit1*VBG_signal_xSec[1];
           if(quantileExpected1 > 0.974 && (quantileExpected1 < 0.977)) Limit_Exp_p2s[1] = limit1*VBG_signal_xSec[1];
           if(quantileExpected1 == -1) Limit_Obs[1] = limit1*VBG_signal_xSec[1];

         }

         if (mh1 == 1500){
           if(quantileExpected1 < 0.17 && (quantileExpected1 > 0)) Limit_Exp_m1s[2] = limit1*VBG_signal_xSec[2];
           if(quantileExpected1 == 0.5) Limit_Exp[2] = limit1*VBG_signal_xSec[2];
           if(quantileExpected1 > 0.8 && (quantileExpected1 < 0.9)) Limit_Exp_p1s[2] = limit1*VBG_signal_xSec[2];
           if(quantileExpected1 > 0.024 && (quantileExpected1 < 0.026)) Limit_Exp_m2s[2] = limit1*VBG_signal_xSec[2];
           if(quantileExpected1 > 0.974 && (quantileExpected1 < 0.977)) Limit_Exp_p2s[2] = limit1*VBG_signal_xSec[2];
           if(quantileExpected1 == -1) Limit_Obs[2] = limit1*VBG_signal_xSec[2];
         }

         if (mh1 == 2000){
           if(quantileExpected1 < 0.17 && (quantileExpected1 > 0)) Limit_Exp_m1s[3] = limit1*VBG_signal_xSec[3];
           if(quantileExpected1 == 0.5) Limit_Exp[3] = limit1*VBG_signal_xSec[3];
           if(quantileExpected1 > 0.8 && (quantileExpected1 < 0.9)) Limit_Exp_p1s[3] = limit1*VBG_signal_xSec[3];
           if(quantileExpected1 > 0.024 && (quantileExpected1 < 0.026)) Limit_Exp_m2s[3] = limit1*VBG_signal_xSec[3];
           if(quantileExpected1 > 0.974 && (quantileExpected1 < 0.977)) Limit_Exp_p2s[3] = limit1*VBG_signal_xSec[3];
           if(quantileExpected1 == -1) Limit_Obs[3] = limit1*VBG_signal_xSec[3];
         }

         if (mh1 == 2500){
           if(quantileExpected1 < 0.17 && (quantileExpected1 > 0)) Limit_Exp_m1s[4] = limit1*VBG_signal_xSec[4];
           if(quantileExpected1 == 0.5) Limit_Exp[4] = limit1*VBG_signal_xSec[4];
           if(quantileExpected1 > 0.8 && (quantileExpected1 < 0.9)) Limit_Exp_p1s[4] = limit1*VBG_signal_xSec[4];
           if(quantileExpected1 > 0.024 && (quantileExpected1 < 0.026)) Limit_Exp_m2s[4] = limit1*VBG_signal_xSec[4];
           if(quantileExpected1 > 0.974 && (quantileExpected1 < 0.977)) Limit_Exp_p2s[4] = limit1*VBG_signal_xSec[4];
           if(quantileExpected1 == -1) Limit_Obs[4] = limit1*VBG_signal_xSec[4];
         }

         if (mh1 == 3000){
           if(quantileExpected1 < 0.17 && (quantileExpected1 > 0)) Limit_Exp_m1s[5] = limit1*VBG_signal_xSec[5];
           if(quantileExpected1 == 0.5) Limit_Exp[5] = limit1*VBG_signal_xSec[5];
           if(quantileExpected1 > 0.8 && (quantileExpected1 < 0.9)) Limit_Exp_p1s[5] = limit1*VBG_signal_xSec[5];
           if(quantileExpected1 > 0.024 && (quantileExpected1 < 0.026)) Limit_Exp_m2s[5] = limit1*VBG_signal_xSec[5];
           if(quantileExpected1 > 0.974 && (quantileExpected1 < 0.977)) Limit_Exp_p2s[5] = limit1*VBG_signal_xSec[5];
           if(quantileExpected1 == -1) Limit_Obs[5] = limit1*VBG_signal_xSec[5];
         }

        }
  
  }   

  TLegend *legendr = new TLegend(0.5460484,0.6334783,0.9230645,0.8852174,NULL,"brNDC");

   legendr->SetShadowColor(0);
   legendr->SetBorderSize(0);
   legendr->SetFillColor(0);

   //1 sigma
   TGraph *limit_1s_up = new TGraph(2*val+1);
   for (int i=0;i<val;i++)
     {
       limit_1s_up->SetPoint(i,Zprime_signal_Masses[i],Limit_Exp_p1s[i]);
       limit_1s_up->SetPoint(i+val,Zprime_signal_Masses[val-1-i],Limit_Exp_m1s[val-1-i]);
     }

   limit_1s_up->SetPoint(2*val,Zprime_signal_Masses[0],Limit_Exp_p1s[0]);
   limit_1s_up->SetLineStyle(2);
   limit_1s_up->SetFillColor(kYellow);

   // 2 sigma
   TGraph *limit_2s_up = new TGraph(2*val+1);
   for (int i=0;i<val;i++)
     {
       limit_2s_up->SetPoint(i,Zprime_signal_Masses[i],Limit_Exp_p2s[i]);
       limit_2s_up->SetPoint(i+val,Zprime_signal_Masses[val-1-i],Limit_Exp_m2s[val-1-i]);
     }
   limit_2s_up->SetPoint(2*val,Zprime_signal_Masses[0],Limit_Exp_p2s[0]);
   limit_2s_up->SetLineStyle(2);
   limit_2s_up->SetFillColor(kGreen);
   limit_2s_up->SetMinimum(0.0);
   limit_2s_up->SetMaximum(8.);
   limit_2s_up->SetTitle(0);
   limit_2s_up->GetXaxis()->SetTitle("m(Z') [GeV]");
   limit_2s_up->GetXaxis()->SetLabelFont(42);
   limit_2s_up->GetXaxis()->SetLabelSize(0.05);
   limit_2s_up->GetXaxis()->SetTitleSize(0.05);
   limit_2s_up->GetXaxis()->SetTitleOffset(0.9);
   limit_2s_up->GetXaxis()->SetTitleFont(42);

   limit_2s_up->GetYaxis()->SetTitle("#sigma(pp#rightarrow Z') x BR(Z'#rightarrow#tau#tau) [pb]");
   limit_2s_up->GetYaxis()->SetLabelFont(42);
   limit_2s_up->GetYaxis()->SetLabelSize(0.05);
   limit_2s_up->GetYaxis()->SetTitleSize(0.05);
   limit_2s_up->GetYaxis()->SetTitleOffset(1.29);
   limit_2s_up->GetYaxis()->SetTitleFont(42);   


   limit_2s_up->Draw("ALF2");
   limit_1s_up->Draw("LF2same");

   TGraph *limit_1s_dn_leg = new TGraph(val, Zprime_signal_Masses, Limit_Exp_m1s);
   limit_1s_dn_leg->SetLineWidth(2);
   limit_1s_dn_leg->SetLineStyle(4);

   //Expected
   TGraph *limit_exp = new TGraph(val,Zprime_signal_Masses,Limit_Exp);
   limit_exp->SetLineWidth(2);
   limit_exp->SetLineStyle(10);
   limit_exp->SetLineColor(1);
   limit_exp->Draw("C");

   TGraph *limit_exp_leg = new TGraph(val,Zprime_signal_Masses,Limit_Exp);
   limit_exp_leg->SetLineWidth(2);
   limit_exp_leg->SetLineStyle(10);

   for (int i = 0 ; i < val; i++)
     {
      cout<<"Mass "<<Zprime_signal_Masses[i]<<"  Limit "<<Limit_Exp[i]<<endl;
     }

   //Observed
   TGraph *limit_obs = new TGraph(val,Zprime_signal_Masses,Limit_Obs);
   limit_obs->SetLineWidth(3);
   limit_obs->SetLineColor(1);
   limit_obs->Draw("C");

   TGraph *limit_obs_leg = new TGraph(val,Zprime_signal_Masses,Limit_Obs);
   limit_obs_leg->SetLineWidth(3);

   // Theory
   TGraph *limit_theory = new TGraph(val2, Zprime_signal_Masses2, VBG_signal_xSec2);
   limit_theory->SetLineColor(kBlue);
   limit_theory->SetLineWidth(2);
   limit_theory->SetLineStyle(1);
   limit_theory->Draw("C");

   TGraph *limit_theory_up = new TGraph(val2, Zprime_signal_Masses2, y_theory_up);
   limit_theory_up->SetLineColor(kBlue);
   limit_theory_up->SetLineWidth(2);
   limit_theory_up->SetLineStyle(4);
   //limit_theory_up->Draw("sameC");

   TGraph *limit_theory_dn = new TGraph(val2, Zprime_signal_Masses2, y_theory_dn);
   limit_theory_dn->SetLineColor(kBlue);
   limit_theory_dn->SetLineWidth(2);
   limit_theory_dn->SetLineStyle(4);
   //limit_theory_dn->Draw("sameC");

   //Legend
   legendr->AddEntry(limit_obs_leg,"Observed","L");
   legendr->AddEntry(limit_1s_up,"Expected #pm 1#sigma","fl");
   legendr->AddEntry(limit_exp_leg, "Expected", "L"); 
   legendr->AddEntry(limit_2s_up,"Expected #pm 2#sigma","fl");
   legendr->AddEntry(limit_theory,"#sigma(LO*1.3)","L");
   legendr->SetFillStyle(0);
   legendr->Draw();

   TLatex *   tex = new TLatex(300,8.5,"CMS Preliminary");
   tex->SetTextSize(0.04347826);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(2200.,8.5,"2.2 fb^{-1} (13 TeV)");
   tex->SetTextSize(0.04347826);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(700,5.5,"#tau_{h}#tau_{h}");
   tex->SetTextSize(0.05847826);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(89.4153,230,"m( #tilde{#chi}^{#pm}_{1}) - m( #tilde{#chi}^{0}_{1}) = 50 GeV");
   tex->SetTextColor(50);
   tex->SetTextSize(0.03847826);
   tex->SetLineWidth(2);
   //tex->Draw();
      tex = new TLatex(89.4153,160,"m( #tilde{#chi}^{#pm}_{1}) - m( #tilde{#tau}) = 5 GeV");
   tex->SetTextColor(50);
   tex->SetTextSize(0.03847826);
   tex->SetLineWidth(2);
   //tex->Draw();
      tex = new TLatex(89.4153,75,"m( #tilde{#chi}^{0}_{1}) = 0 GeV");
   tex->SetTextColor(32);
   tex->SetTextSize(0.03847826);
   tex->SetLineWidth(2);
   //tex->Draw();
      tex = new TLatex(89.4153,112,"m( #tilde{#chi}^{#pm}_{1}) - m( #tilde{#tau}) = 5 GeV");
   tex->SetTextColor(32);
   tex->SetTextSize(0.03847826);
   tex->SetLineWidth(2);
   //tex->Draw();
      tex = new TLatex(89.4153,340,"#tilde{#chi}^{#pm}_{1} #rightarrow #tilde{#tau} #nu_{#tau} , #tilde{#chi}^{0}_{2} #rightarrow #tilde{#tau} #tau ");
   tex->SetTextSize(0.03847826);
   tex->SetLineWidth(2);
   //tex->Draw();
   //Save plot
   mass->SaveAs("Limit_VBF.pdf");


}

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
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  //tdrStyle->SetFuncColor(1);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("e"); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kGray);
  tdrStyle->SetStatFont(42);

  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(0);
  tdrStyle->SetStatX(1.); //Starting position on X axis
  tdrStyle->SetStatY(1.); //Starting position on Y axis
  tdrStyle->SetStatFontSize(0.025); //Vertical Size
  tdrStyle->SetStatW(0.15); //Horizontal size
  // tdrStyle->SetStatStyle(Style_t style = 1001);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.125);
  tdrStyle->SetPadLeftMargin(0.105);
  tdrStyle->SetPadRightMargin(0.1);

  // For the Global title:

  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.9);
  tdrStyle->SetTitleOffset(0.7, "Y"); // Another way to set the Offset

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

  tdrStyle->cd();


}
