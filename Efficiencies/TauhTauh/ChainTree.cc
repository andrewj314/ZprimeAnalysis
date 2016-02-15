/**
This Macro
1. Uses a chain to merge the content of the trees that have been output by the analyzer and save it on a new rootpla 

Need to specify
0. See Declare constants
1. Do "voms-proxy-init --voms cms" if you read remote files
   Use lcg-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/fromeo/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
   or corresponding syntax you use a different tier, to redirect the list of files you want to process in a text file
2. Adjust the list of files (paht+names) to be readable by root
   0,$s/\/cms/root:\/\/xrootd-cms.infn.it\//g in the file specified in name_file
*/
/////
//   To run: root -l ChainTree.cc+  
/////
/////
//   Prepare Root and Roofit
/////
#include "TTree.h"
#include "TTreePlayer.h"
#include "TFile.h"
#include "TChain.h"
#include "TFileCollection.h"
#include <iostream>
using namespace std;
////
//   Declare constants
/////
const string tree_name    = "TNT/BOOM"; //The name of the tree you have defined in your rootplas
const string sample       = "TTJets";
const string name_file    = sample+".txt"; //List with the path of the rootfiles
const string name_rootple = sample+".root"; //The name of the new rootpla
/////
//   Main function
/////
void ChainTree(){
 //Put here the tree
 TTree *tree = new TTree("BOOM","BOOM");
 tree->SetMaxTreeSize(99000000000);
 //Create new file
 TFile *newfile = new TFile(name_rootple.c_str(),"recreate");
 newfile->cd();
 //Create chain, merge trees 
 TChain * chain      = new TChain(tree_name.c_str(),"");
 TFileCollection* fc = new TFileCollection("list", "list",name_file.c_str());
 chain->AddFileInfoList((TCollection*)fc->GetList());
 //Save it
 tree    = chain->CopyTree("");
 newfile = tree->GetCurrentFile();
 tree    = NULL;
 newfile->Write();
 newfile->Close();
 delete newfile; 
}
