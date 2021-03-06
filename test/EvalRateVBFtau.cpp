// evaluates rate  - compile with c++ -lm -o EvalRateVBF test/EvalRateVBF.cpp `root-config --glibs --cflags`
// and launch the executable 
#include <sstream>
#include <fstream>
#include <map>
#include <iostream>
#include <fstream>
#include <utility>
#include <regex>
#include <tuple>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "../interface/object.h"

using namespace std;

void appendFromFileList (TChain* chain, TString filename)
{
 
  std::ifstream infile(filename.Data());
  std::string line;
  while (std::getline(infile, line))
    {
      line = line.substr(0, line.find("#", 0)); // remove comments introduced by #
      while (line.find(" ") != std::string::npos) line = line.erase(line.find(" "), 1); // remove white spaces
      while (line.find("\n") != std::string::npos) line = line.erase(line.find("\n"), 1); // remove new line characters
      while (line.find("\r") != std::string::npos) line = line.erase(line.find("\r"), 1); // remove carriage return characters
      if (!line.empty()) // skip empty lines
	chain->Add(line.c_str());
    }
  return;
}


int main(int argc, char** argv){

 
  float nbStudiedRun =96.;
  float thisLumiRun = 1750E30;
  float scaleToLumi = 2.E34;
  float scale = 0.001*(nbStudiedRun*11245.6)*scaleToLumi/thisLumiRun;  
 

  TString LumiTarget = "scaleToLumi20E33";
  TString year = "2017";
  
  cout << "Scale factor: " << scale << endl;
 
  TString directory = "/data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesHighPU_2017_fill6194/";
  TString fileList = "fileLists/ZeroBias2017_HighPU.list";
  bool emulated  =false;
  bool reweight = false;
  
  
  ifstream PUFile("utils/PU_per_LS_fill6194_2017.txt"); //2017

  
  TTree * tInput;
  TChain * cInput;
  TChain * lInput; 



  
  if (emulated){
    cInput = new TChain ("l1UpgradeEmuTree/L1UpgradeTree");
    appendFromFileList(cInput, fileList);
    lInput = new TChain ("l1EventTree/L1EventTree");
    appendFromFileList(lInput, fileList);
  }else{
    cInput= new TChain ("l1UpgradeTree/L1UpgradeTree");
    appendFromFileList(cInput, fileList);
    lInput = new TChain ("l1EventTree/L1EventTree");
    appendFromFileList(lInput, fileList);

  }

  TString fOutNameVBF;


  if(emulated){
    fOutNameVBF = directory+"Emu_rateL1_VBF_"+LumiTarget+"_"+year+".root";

  }else{
    fOutNameVBF = directory+"rateL1_VBF_"+LumiTarget+"_"+year+"test.root";

  }
  
  
  lInput->SetMakeClass(1);
  cInput->SetMakeClass(1);  
  
  Int_t lumi ;
  UInt_t run ;

  UShort_t stage2_jetN ;
  std::vector<float> stage2_jetEt;
  std::vector<float> stage2_jetEta;
  std::vector<float> stage2_jetPhi;
  std::vector<short> stage2_jetBx;
  
  

  
  // set branch and variables
  TBranch *b_lumi ;
  TBranch *b_run ;

  
  TBranch *b_stage2_jetN ;
  TBranch *b_stage2_jetEt;
  TBranch *b_stage2_jetEta;
  TBranch *b_stage2_jetPhi;
  TBranch *b_stage2_jetBx;
  
  

  

    lInput ->SetBranchAddress("lumi", &lumi , &b_lumi );
    lInput ->SetBranchAddress("run", &run , &b_run );

    cInput ->SetBranchAddress("nJets", &stage2_jetN , &b_stage2_jetN);
    cInput ->SetBranchAddress("jetEta", &stage2_jetEta, &b_stage2_jetEta);
    cInput ->SetBranchAddress("jetPhi", &stage2_jetPhi, &b_stage2_jetPhi);
    cInput ->SetBranchAddress("jetEt", &stage2_jetEt, &b_stage2_jetEt);
    cInput ->SetBranchAddress("jetBx", &stage2_jetBx, &b_stage2_jetBx);
    
   
  TFile* fOutVBF = new TFile (Form("%s",fOutNameVBF.Data()), "recreate");

   fOutVBF->cd(); 
  //VBF

  TH2D* DiJet2D_Pass = new TH2D ("DiJet2D", "events firing, 2-dim (lead, mjj)", 800, 0, 800, 100,30, 130);
  TH2D* Ratio_DiJet2D = new TH2D ("Ratio_DiJet2D", "VBF 2-dim ratio", 800, 0, 800, 100, 30, 130);
  TH2D* Rate_DiJet2D = new TH2D ("Rate_DiJet2D", "VBF 2-dim rate", 800, 0, 800, 100, 30, 130);

  TH2D* DiJet2D_Sub_Pass = new TH2D ("DiJet2D_Sub", "events firing, 2-dim (lead, sublead)", 30, 30, 60, 100,30, 130);
  TH2D* Ratio_Sub_DiJet2D = new TH2D ("Ratio_Sub_DiJet2D", "VBF 2-dim ratio", 30, 30, 60, 100, 30, 130);
  TH2D* Rate_Sub_DiJet2D = new TH2D ("Rate_Sub_DiJet2D", "VBF 2-dim rate", 30, 30, 60, 100, 30, 130); 

  TH2D* DiJet2D_Pass_rej = new TH2D ("DiJet2D_rej", "events firing, 2-dim (lead, mjj), no TT28", 800, 0, 800, 100,30, 130);
  TH2D* Ratio_DiJet2D_rej = new TH2D ("Ratio_DiJet2D_rej", "VBF 2-dim ratio", 800, 0, 800, 100, 30, 130);
  TH2D* Rate_DiJet2D_rej = new TH2D ("Rate_DiJet2D_rej", "VBF 2-dim rate", 800, 0, 800, 100, 30, 130);

  TH2D* DiJet2D_Sub_Pass_rej = new TH2D ("DiJet2D_Sub_rej", "events firing, 2-dim (lead, sublead), no TT28", 30, 30, 60, 100,30, 130);
  TH2D* Ratio_Sub_DiJet2D_rej = new TH2D ("Ratio_Sub_DiJet2D_rej", "VBF 2-dim ratio", 30, 30, 60, 100, 30, 130);
  TH2D* Rate_Sub_DiJet2D_rej = new TH2D ("Rate_Sub_DiJet2D_rej", "VBF 2-dim rate", 30, 30, 60, 100, 30, 130); 
  
  TH1D* VBF_lead = new TH1D ("VBF_lead", "events firing, 1-dim (lead jet)", 100,30, 130);
  TH1D* Ratio_VBF_lead = new TH1D ("Ratio_VBF_lead", "VBF 1-dim ratio", 100,30, 130);
  TH1D* Rate_VBF_lead = new TH1D ("Rate_VBF_lead", "VBF 1-dim rate", 100,30, 130);
  
  TH1D* Mjj30 = new TH1D ("Mjj30", "", 100, 0, 800);
  TH1D* jetsRes = new TH1D ("jetsRes", "", 60, 0, 30);
  TH1D* Mjj35 = new TH1D ("Mjj35", "", 100, 0, 800);
  TH1D* Mjj40 = new TH1D ("Mjj40", "", 100, 0, 800);
  TH1D* Mjj45 = new TH1D ("Mjj45", "", 100, 0, 800);
  TH1D* Mjj50 = new TH1D ("Mjj50", "", 100, 0, 800);
  TH1D* Mjj55 = new TH1D ("Mjj55", "", 100, 0, 800);

   
 
  std::map<Int_t,Float_t> PU_per_LS;
  std::string str; 
  while (std::getline(PUFile, str))
    {
      TString temp(str);

      //temp.ReplaceAll("5412,283171,",""); //2016

             regex reg("6194,[0-9]+,");
          temp = regex_replace(str, reg, "");
      int pos_coma = temp.First(",");
      TString LS_str(temp,pos_coma);
      // cout<<LS_str<<endl;
      TString Replacing = LS_str ;
      Replacing += ",";
      temp.ReplaceAll(Replacing.Data(),"");
      TString PU_str = temp;
      // cout<<PU_str<<endl;
      std::istringstream ss_LS(LS_str.Data());
      Int_t LS ;
      ss_LS >> LS;
      std::istringstream ss_PU(PU_str.Data());
      Float_t PU ;
      ss_PU >> PU;     
      PU_per_LS.insert(std::pair<Int_t,Float_t>(LS , PU ));

    }
  // analyze data    
  long int nEvents = 0;

  int nEventsPass = 0;



  std::vector<object> jet30;   
  std::vector<object> jet30rej;   



  std::vector< std::tuple<double,int,int> > mjj_pass; //VBF
  std::vector< std::tuple<double,int,int,double, double> > mjj_pass_sortPt; //VBF
  std::vector< std::tuple<double,int,int> > mjj_rej; //VBF
  std::vector< std::tuple<double,int,int,double, double> > mjj_rej_sortPt; //VBF


  
  for (Long64_t iEv =0 ;iEv<30000; ++iEv){
    lumi = 0;			   
    run = 0;
    
    stage2_jetN =0;		   
    stage2_jetEt.clear(); 
    stage2_jetEta.clear();
    stage2_jetPhi.clear();
    int got = 0;
    
    lInput->GetEntry(iEv);
    got = cInput->GetEntry(iEv);
    
    if (got == 0) break;

   
  
     
    if (iEv%1000 == 0) cout << iEv << " / " << nEvents << endl;
       
    jet30.clear();
    jet30rej.clear();
    mjj_pass.clear();
    mjj_pass_sortPt.clear();
    mjj_rej.clear();
    mjj_rej_sortPt.clear();


    

    if(run>302674) continue;
    if(lumi<135 || lumi>264) continue; 
    
    if(PU_per_LS.find(lumi)==PU_per_LS.end()) continue;
    Float_t weight = 0;
    if (reweight){
      weight = PU_per_LS[56]/PU_per_LS[lumi];
    }else{
      weight =  1.;
    }

    nEventsPass ++;
    


    
   
    for (long int iL1 = 0; iL1 < stage2_jetEt.size(); iL1++){ //loop on jets
      // selections
      double jetPt  = stage2_jetEt.at(iL1);
      double jetEta  = fabs(stage2_jetEta.at(iL1));
      int jetBx = fabs(stage2_jetBx.at(iL1));
     
      if(jetPt>30. && jetBx == 0) {
	jet30.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	if(jetEta>2.7 && jetEta<3.0){
	  if(jetPt>60.) jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	}else{
	  jet30rej.push_back(object(stage2_jetEt.at(iL1),stage2_jetEta.at(iL1),stage2_jetPhi.at(iL1),-999)) ;
	}
      }
      
    }
    


    std::sort (jet30.begin(),jet30.end());
    std::sort (jet30rej.begin(),jet30rej.end());
    

    //VBF

       
    if (jet30.size() >= 2){
      for (int iJet = 0; iJet <jet30.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet30.size(); kJet++){      
	  if (kJet!=iJet) {
	    	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet30[iJet].Et(),
			      jet30[iJet].Eta(),
			      jet30[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet30[kJet].Et(),
			      jet30[kJet].Eta(),
			      jet30[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_pass.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    if(jetPair.M()>=620)    mjj_pass_sortPt.push_back(make_tuple(jetPair.M(),iJet,kJet,jet30[iJet].Et(),jet30[kJet].Et()));
	    double jetDiff = jet30[iJet].Et()-jet30[kJet].Et();
	    double jetSum  = jet30[iJet].Et()+jet30[kJet].Et();
	    jetsRes->Fill(jetSum/jetDiff);
	    Mjj30->Fill(jetPair.M());
	    if(jet30[kJet].Et()>35) Mjj35->Fill(jetPair.M());
	    if(jet30[kJet].Et()>40) Mjj40->Fill(jetPair.M());
	    if(jet30[kJet].Et()>45) Mjj45->Fill(jetPair.M());
	    if(jet30[kJet].Et()>50) Mjj50->Fill(jetPair.M());
	    if(jet30[kJet].Et()>55) Mjj55->Fill(jetPair.M());
	    
	  }
	  
	}
	
      }
      std::sort(mjj_pass.begin(),mjj_pass.end());
      std::sort(mjj_pass_sortPt.begin(),mjj_pass_sortPt.end(),SortMjjByJetThreshold);
      double Mjj_max = std::get<0>(*(mjj_pass.rbegin()));      
      DiJet2D_Pass -> Fill (Mjj_max,jet30[0].Et(),weight);        
      if(mjj_pass_sortPt.size()>0) {
	double subJet_min = std::get<4>(*(mjj_pass_sortPt.rbegin()));
	if (Mjj_max>620) DiJet2D_Sub_Pass -> Fill (subJet_min,jet30[0].Et(),weight);  
	if (subJet_min>35 && Mjj_max>620) VBF_lead->Fill(jet30[0].Et(),weight); 
      }else{
	DiJet2D_Sub_Pass->Fill(-1,-1);
      }
    } else {
      DiJet2D_Pass -> Fill (-1, -1);
      DiJet2D_Sub_Pass->Fill(-1,-1);
    }


    //VBF rejecting TT28
    
    if (jet30rej.size() >= 2){
      for (int iJet = 0; iJet <jet30rej.size(); iJet++){      
	for (int kJet = iJet+1; kJet <jet30rej.size(); kJet++){      
	  if (kJet!=iJet) {
	    
	    TLorentzVector ijet;
	    ijet.SetPtEtaPhiM(
			      jet30rej[iJet].Et(),
			      jet30rej[iJet].Eta(),
			      jet30rej[iJet].Phi(),
			      0.);
	    TLorentzVector kjet;
	    kjet.SetPtEtaPhiM(
			      jet30rej[kJet].Et(),
			      jet30rej[kJet].Eta(),
			      jet30rej[kJet].Phi(),
			      0.);
	    TLorentzVector jetPair = ijet+kjet;
	    mjj_rej.push_back(make_tuple(jetPair.M(),iJet,kJet));
	    if(jetPair.M()>=620) mjj_rej_sortPt.push_back(make_tuple(jetPair.M(),iJet,kJet,jet30rej[iJet].Et(),jet30rej[kJet].Et()));
	  }
	}
      }
      std::sort(mjj_rej.begin(),mjj_rej.end());
      std::sort(mjj_rej_sortPt.begin(),mjj_rej_sortPt.end(),SortMjjByJetThreshold);
      double Mjj_max =  std::get<0>(*(mjj_rej.rbegin()));
      DiJet2D_Pass_rej -> Fill ( Mjj_max,jet30rej[0].Et(),weight);  
      if(mjj_rej_sortPt.size()>0) {
	double subJet_min =  std::get<4>(*(mjj_rej_sortPt.rbegin()));
	DiJet2D_Sub_Pass_rej -> Fill (subJet_min,jet30rej[0].Et(),weight);  
      }else{
	DiJet2D_Sub_Pass_rej->Fill(-1,-1);
      }
    } else {
      DiJet2D_Pass_rej -> Fill (-1, -1);
      DiJet2D_Sub_Pass_rej->Fill(-1,-1);
    }  
      
  }

  // compute rate plots

  cout << endl;
  cout << "Computing rates..." << endl; 

  //VBF
  for (int i = 1; i <=DiJet2D_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiJet2D_Pass->GetNbinsY();j++ ){
      double binDiJet2D = 1.*(DiJet2D_Pass->Integral(i, DiJet2D_Pass->GetNbinsX()+1,j,DiJet2D_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      binDiJet2D *=scale;
      Rate_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      //rej
      binDiJet2D = 1.*(DiJet2D_Pass_rej->Integral(i, DiJet2D_Pass_rej->GetNbinsX()+1,j,DiJet2D_Pass_rej->GetNbinsY()+1))/nEventsPass;
      Ratio_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
      binDiJet2D *=scale;
      Rate_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
    }
  }
  
  for (int j = 1; j<= VBF_lead->GetNbinsX();j++ ){
      double binDiJet = 1.*(VBF_lead->Integral(j, VBF_lead->GetNbinsX()+1))/nEventsPass;
      Ratio_VBF_lead -> SetBinContent (j, binDiJet);        
      binDiJet *=scale;
      Rate_VBF_lead ->  SetBinContent (j, binDiJet);
  }
  
  for (int i = 1; i <=DiJet2D_Sub_Pass->GetNbinsX(); i++){
    for (int j = 1; j<= DiJet2D_Sub_Pass->GetNbinsY();j++ ){
      
      double binDiJet2D = 1.*(DiJet2D_Sub_Pass->Integral(i, DiJet2D_Sub_Pass->GetNbinsX()+1,j,DiJet2D_Sub_Pass->GetNbinsY()+1))/nEventsPass;
      Ratio_Sub_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      binDiJet2D *=scale;
      Rate_Sub_DiJet2D -> SetBinContent (i, j, binDiJet2D);        
      //rej
      binDiJet2D = 1.*(DiJet2D_Sub_Pass_rej->Integral(i, DiJet2D_Sub_Pass_rej->GetNbinsX()+1,j,DiJet2D_Sub_Pass_rej->GetNbinsY()+1))/nEventsPass;
      Ratio_Sub_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
      binDiJet2D *=scale;
      Rate_Sub_DiJet2D_rej -> SetBinContent (i, j, binDiJet2D);        
    }
  }
 
  int xbin = Rate_DiJet2D->GetXaxis()->FindBin(620.0);
  int ybin = Rate_DiJet2D->GetYaxis()->FindBin(90.0);
  cout<<" Rate VBFseed 90 30 620  "<<Rate_DiJet2D->GetBinContent(xbin,ybin)<<" kHz"<<endl;
  xbin = Rate_Sub_DiJet2D->GetXaxis()->FindBin(30.0);
  ybin = Rate_Sub_DiJet2D->GetYaxis()->FindBin(90.0);
  cout<<" Rate VBFseed 90 30 620  "<<Rate_Sub_DiJet2D->GetBinContent(xbin,ybin)<<" kHz"<<endl; 
  fOutVBF->cd();
  fOutVBF -> Write();
  cout <<"Output saved in "<<fOutNameVBF<<endl;
}

