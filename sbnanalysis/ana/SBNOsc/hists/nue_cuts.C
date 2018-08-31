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

void nue_cuts(){
  //TFile* nue = new TFile("nueoutput.root");
  TFile* nue = new TFile("output_SBNOsc_NueSelection.root"); 
  
  THStack* cut1 = (THStack*)nue->Get("cut1stack");
  THStack* cut1reco = (THStack*)nue->Get("cut1stack_reco");
  TCanvas* c1 = new TCanvas("c1","c1",2000,1000);
  c1->Divide(2,1);
  c1->cd(1);
  cut1->Draw("nostack");
  CreateNamedLegend("#nu_e CC candidates (matched)", cut1,{"nu_e CC interactions","\" w/ 1+ matched E_sh > 200 MeV","\" which are prim. e- sh"},false)->Draw();
  c1->cd(2);
  cut1reco->Draw("nostack");  
  CreateNamedLegend("Cut 1, reconstructed E_#nu", cut1reco, {"ev w/ 1+ matched E_sh > 200 MeV","\" which are prim. e- sh","empty"},false)->Draw();
  
  THStack* fig11 = (THStack*)nue->Get("fig11stack");
  TCanvas* canv11 = new TCanvas("canv11", "canv11", 1000, 1000);
  fig11->Draw();
  CreateNamedLegend("", fig11, {"#mu->#nue","K+->#nue","K0->#nue","NC single #gamma", "#nu_mu CC"},false)->Draw();
  
  
  THStack* cut2 = (THStack*)nue->Get("cut2stack");
  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  cut2->Draw();
  CreateNamedLegend("NC #gamma prod", cut2, {"event rejected","event passed"},true)->Draw();

  TCanvas* canv2a = new TCanvas("canv2a","2nd photon cut canv",1600,800);
  canv2a->Divide(2,1);
  canv2a->cd(1);
  THStack* pi2gamma = (THStack*)nue->Get("2ndPhotonCut_stack");
  pi2gamma->Draw("nostack");
  CreateNamedLegend("2nd photon cut", pi2gamma, {"all #gamma showers w/in 5 cm","#gamma showers w/ #pi0 parent","#gamma sh with E_#gamma1 > 200 MeV","#pi0->2#gamma","2nd #gamma with E_#gamma2 > 100 MeV","# 2ndary #gamma showers w/in fid. vol"},false)->Draw();
  canv2a->cd(2);
  THStack* cut2a = (THStack*)nue->Get("cut2a_stack");
  cut2a->Draw("nostack");
  CreateNamedLegend("2nd photon cut", cut2a, {"all #nu events","rejected events from 2nd photon cut"}, false)->Draw();

  TCanvas* canv2b = new TCanvas("canv2b","conversion gap cut canv", 1600,800);
  canv2b->Divide(2,1);
  canv2b->cd(1);
  THStack* visiblevertex = (THStack*)nue->Get("visiblevertexstack"); 
  visiblevertex->Draw("nostack");
  CreateNamedLegend("conversion gap cut",visiblevertex,{"# w/in 5 cm", "# from visible vertex", "# converting >3 cm away from vertex"},false)->Draw();
  canv2b->cd(2);
  THStack* cut2b = (THStack*)nue->Get("cut2b_stack");
  cut2b->Draw("nostack");
  CreateNamedLegend("Conversion gap cut", cut2b, {"All #nu events","visible #nu events","rejected from conversion gap cut"},false)->Draw();
  
  
  TCanvas* dEdxcanv_gamma = new TCanvas("dEdxcanv_gamma","gamma dEdx cut canv",800,800);
  THStack* dEdx_gamma = (THStack*)nue->Get("gammastack");
  dEdx_gamma->Draw("nostack");
  CreateNamedLegend("dEdx cut", dEdx_gamma, {"all #gamma showers w/in 5 cm","#gamma showers from unrejected #nu event","# #gamma showers passing dE/dx cut"}, false)->Draw();
  
  
  THStack* cut3 = (THStack*)nue->Get("cut3stack");
  TCanvas* c3 = new TCanvas("c3","c3",800,800);
  cut3->Draw();
  CreateNamedLegend("CC #nu_#mu", cut3, {"Failed #mu","Passed #mu"},true)->Draw();
  
  
  TH1D* dEdx = (TH1D*)nue->Get("dEdx");
  THStack* dEdx_stack = (THStack*)nue->Get("showerStack");
  TCanvas* dEdxcanv = new TCanvas("dEdxcanv","dEdxcanv",2000,1000);
  dEdxcanv->Divide(2,1);
  dEdxcanv->cd(1);
  dEdx->Draw();
  dEdxcanv->cd(2);
  dEdx_stack->Draw();
  CreateNamedLegend("dEdx",dEdx_stack,{"#gamma","e"},false)->Draw();
  
}  
