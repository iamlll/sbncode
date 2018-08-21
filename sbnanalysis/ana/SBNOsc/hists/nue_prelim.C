#include <iostream>
#include <math.h>
#include <fstream>
#include "TCanvas.h"

using namespace std;

TLegend* CreateNamedLegend(string legTitle, THStack* stack, vector<string> titles, bool percent){ 
  char* cstrLeg = new char[legTitle.length()+1];
  strcpy(cstrLeg, legTitle.c_str());
  int totalEntries = 0;
  TLegend* leg = new TLegend(0.45,0.7,0.9,0.9,cstrLeg);
  vector<int> successes;
  for(auto obj: *stack->GetHists()){
    TH1D* hist = (TH1D*) obj;
    int entries=0;
    for(int bin=0; bin < hist->GetNbinsX(); bin++){
      entries += hist->GetBinContent(bin);
    }
    totalEntries += entries;
    successes.push_back(entries);
  }
  char buffer[20];
  int num=0; //num of hists
  if(percent==true){
    for(auto obj: *stack->GetHists()){
      TH1D* hist = (TH1D*) obj;
      double proportion = ((double)successes[num])/totalEntries*100;
      char* cstrEntry = new char [titles[num].length()+1];
      strcpy (cstrEntry, titles[num].c_str());
      sprintf(buffer, "%s%s%f%s",cstrEntry, " (", proportion, "%)");
      leg->AddEntry(hist, buffer, "f");
      num++;
    }
  }
  else{
    for(auto obj: *stack->GetHists()){
      TH1D* hist = (TH1D*) obj;
      char* cstrEntry = new char [titles[num].length()+1];
      strcpy (cstrEntry, titles[num].c_str());
      sprintf(buffer, "%s%s%d%s",cstrEntry, " (", successes[num], ")");
      leg->AddEntry(hist, buffer, "f");
      num++;
    }
  }
  return leg;
}  

void nue_prelim(){
  TFile* nue = new TFile("output.root");

  THStack* shAssns = (THStack*)nue->Get("showerE_stack");
  TH2D* nue_vs_reco = (TH2D*)nue->Get("nuE_vs_reco");
  THStack* nuereco_type_stack = (THStack*)nue->Get("nuereco_stack"); 
  THStack* prelimCuts = (THStack*)nue->Get("prelim_stack_nue");
  TH1D* dist_from_vertex = (TH1D*)nue->Get("dist_from_vertex");
  THStack* vertexDist_truth = (THStack*)nue->Get("truthVD_stack");
  THStack* vertexDist_reco = (THStack*)nue->Get("recoVD_stack");

  //Draw histograms
  TCanvas* shcanv = new TCanvas("shcanv","shcanv",800,800);
  shAssns->Draw();
  CreateNamedLegend("shower assns",shAssns,{"showers assoc. w/ #nu", "showers not assoc. w/ #nu"},false)->Draw();

  TCanvas* pcanv = new TCanvas("pcanv","pcanv",1000,1000);
  prelimCuts->Draw("nostack");  
  CreateNamedLegend("",prelimCuts,{"Total #nu interactions","Interactions in fid vol","true #nu_e CC interactions","#nu_e CC matched w/ shower","matched #nu_e CC w/ prim. e- shower"},false)->Draw();

  TCanvas* distcanv = new TCanvas("distcanv","distcanv",1500,750);
  distcanv->Divide(2,1);
  distcanv->cd(1);
  dist_from_vertex->Draw();
  distcanv->cd(2);
  vertexDist_truth->Draw();
  CreateNamedLegend("",vertexDist_truth,{"true #nu_e CC tr+sh","photon showers","cosmic rays"},false)->Draw();

  TCanvas *energycanv = new TCanvas("energycanv","energycanv",800,800);
  nue_vs_reco->Draw();

  TCanvas *nutypecanv = new TCanvas("nutypecanv","nutypecanv",800,800);
  nutypecanv->Divide(3,2);
  int count = 0;
  for(auto hist : nuereco_type_stack.GetHists()){
    TLegend* leg = new TLegend(0.75,0.45,0.95,0.6);
    nutypecanv->cd(count);
    hist->Draw("colz");
    if(count==0)leg->AddEntry(hist, "CCQE", "f");
    leg->Draw(); 
    count++;
  }
}
