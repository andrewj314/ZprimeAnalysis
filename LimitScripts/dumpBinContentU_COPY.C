#include <iostream> 
#include <fstream>

dumpBinContent(){
  TFile* f1 = new TFile("DIRECTORY/FILE.root");
  cout <<"DIRECTORY/FILE.root"<<endl;
  //TH1F* h_mjj = DiJetMass;
  TH1D* h_mjj = m_effective;
  //h_mjj->Rebin(2); 
  cout << "Nbins "<<h_mjj->GetXaxis()->GetNbins()<<endl; 
  ofstream outputFile;
  outputFile.open("output_DIRECTORY/yieldsU_FILE.txt") ;
  double bin_sum = 0.0;
  
  for (int i = 0; i < (h_mjj->GetXaxis()->GetNbins()+1); i++)
    {
      double error = 999.9;
      if (h_mjj->GetBinContent(i) > 0.0){
        double error = h_mjj->GetBinError(i)/h_mjj->GetBinContent(i);
        error = 1. + error;
        if (error > 10.){
         error = 999.9;
        }
      }
      outputFile <<error<< endl;      
/*
      if(i <= 9){bin_sum += h_mjj->GetBinContent(i);} 
      if (i == 9){outputFile <<bin_sum << endl; bin_sum = 0.0;}
      if((i > 9) && (i <= 20)) {outputFile <<h_mjj->GetBinContent(i) << endl;}
      if (i > 20){bin_sum += h_mjj->GetBinContent(i);}
      if (i == h_mjj->GetXaxis()->GetNbins()){outputFile <<bin_sum << endl;}
      float bin_value = h_mjj->GetBinContent(i);
*/
    }

  outputFile.close(); 
}