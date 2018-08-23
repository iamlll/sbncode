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

void DrawRectangle(double xbotleft, double ybotleft, double xtoprt, double ytoprt, Color_t color){
  TLine* left = new TLine(xbotleft, ybotleft, xbotleft, ytoprt); 
  TLine* right = new TLine(xtoprt, ybotleft, xtoprt, ytoprt);
  TLine* bottom = new TLine(xbotleft, ybotleft, xtoprt, ybotleft);
  TLine* top = new TLine(xbotleft, ytoprt, xtoprt, ytoprt);

  left->SetLineColor(color);
  right->SetLineColor(color);
  top->SetLineColor(color);
  bottom->SetLineColor(color);
  left->SetLineWidth(3);
  right->SetLineWidth(3);
  top->SetLineWidth(3);
  bottom->SetLineWidth(3);

  left->Draw();
  right->Draw();
  top->Draw();
  bottom->Draw();
}

void DrawNamedRect(double xbotleft, double ybotleft, double xtoprt, double ytoprt, Color_t color, std::string name, TLegend* leg){
  TLine* left = new TLine(xbotleft, ybotleft, xbotleft, ytoprt);
  TLine* right = new TLine(xtoprt, ybotleft, xtoprt, ytoprt);
  TLine* bottom = new TLine(xbotleft, ybotleft, xtoprt, ybotleft);
  TLine* top = new TLine(xbotleft, ytoprt, xtoprt, ytoprt);

  left->SetLineColor(color);
  right->SetLineColor(color);
  top->SetLineColor(color);
  bottom->SetLineColor(color);
  left->SetLineWidth(3);
  right->SetLineWidth(3);
  top->SetLineWidth(3);
  bottom->SetLineWidth(3);

  left->Draw();
  right->Draw();
  top->Draw();
  bottom->Draw();

  char* charname = new char [name.length()+1];
  strcpy (charname, name.c_str());
  leg->AddEntry(top,charname,"l");
}

void nue_XYZ(){
  //TFile* nue = new TFile("XYZoutput.root");
  TFile* nue = new TFile("output_SBNOsc_NueSelection.root");
   
  TH2D* xy = (TH2D*)nue->Get("nu_vtx_XY");
  TH2D* xz = (TH2D*)nue->Get("nu_vtx_XZ");
  TH2D* yz = (TH2D*)nue->Get("nu_vtx_YZ");
            
  //Draw histograms
  TLegend* legXY = new TLegend(0.75,0.45,0.95,0.6);
  TLegend* legYZ = new TLegend(0.75,0.45,0.95,0.6);
  TLegend* legXZ = new TLegend(0.75,0.45,0.95,0.6);
  TCanvas* nucanv = new TCanvas("nucanv","nucanv",2400,800);
  nucanv->Divide(3,1);
  nucanv->cd(1);
  xy->Draw();
  //TPC 0
  DrawRectangle(-199.15,-200,-2.65,200,kBlue);//act vol
  DrawRectangle(-199.15+25,-200+25,-2.65-25,200-25, kRed);//fid vol
  //TPC 1
  DrawNamedRect(2.65,-200,199.15,200,kBlue, "active vol", legXY);
  DrawNamedRect(2.65+25,-200+25,199.15-25,200-25, kRed, "fiducial vol", legXY);
  DrawNamedRect(-260.1,-271.15,260.1,271.15,kYellow, "bounding box", legXY);
  legXY->Draw();

  nucanv->cd(2);  
  yz->Draw();
  DrawNamedRect(-200,0,200,500,kBlue, "active vol", legYZ);
  DrawNamedRect(-200+25,25,200-25,500-25, kRed, "fiducial vol", legYZ);
  DrawNamedRect(-271.15,-143.1,271.15,559.6,kYellow, "bounding box", legYZ);
  legYZ->Draw();

  nucanv->cd(3);
  xz->Draw();
  //TPC 0
  DrawNamedRect(-199.15,0,-2.65,500,kBlue, "active vol", legXZ);
  DrawNamedRect(-199.15+25,25,-2.65-25,500-25, kRed, "fiducial vol", legXZ);
  //TPC 1
  DrawRectangle(2.65,0,199.15,500,kBlue);
  DrawRectangle(2.65+25,25,199.15-25,500-25, kRed);
  DrawNamedRect(-260.1,-143.1,260.1,559.6,kYellow, "bounding box", legXZ);
  legXZ->Draw();

  TCanvas* trcanv = new TCanvas("trcanv","trcanv",2400,800);
  trcanv->Divide(3,1);
  TLegend* legXY_tr = new TLegend(0.75,0.45,0.95,0.6);
  TLegend* legYZ_tr = new TLegend(0.75,0.45,0.95,0.6);
  TLegend* legXZ_tr = new TLegend(0.75,0.45,0.95,0.6);
  TH2D* trackXY = (TH2D*)nue->Get("trackXY"); 
  TH2D* trackYZ = (TH2D*)nue->Get("trackYZ");
  TH2D* trackXZ = (TH2D*)nue->Get("trackXZ");
  trcanv->cd(1);
  trackXY->Draw();
  //TPC 0
  DrawRectangle(-199.15,-200,-2.65,200,kBlue);//act vol
  DrawRectangle(-199.15+25,-200+25,-2.65-25,200-25, kRed);//fid vol
  //TPC 1
  DrawNamedRect(2.65,-200,199.15,200,kBlue, "active vol", legXY_tr);
  DrawNamedRect(2.65+25,-200+25,199.15-25,200-25, kRed, "fiducial vol", legXY_tr);
  DrawNamedRect(-260.1,-271.15,260.1,271.15,kYellow, "bounding box", legXY_tr);
  legXY_tr->Draw();

  trcanv->cd(2);
  trackYZ->Draw();
  DrawNamedRect(-200,0,200,500,kBlue, "active vol", legYZ_tr);
  DrawNamedRect(-200+25,25,200-25,500-25, kRed, "fiducial vol", legYZ_tr);
  DrawNamedRect(-271.15,-143.1,271.15,559.6,kYellow, "bounding box", legYZ_tr);
  legYZ_tr->Draw();
  
  trcanv->cd(3);
  trackXZ->Draw();
  //TPC 0
  DrawNamedRect(-199.15,0,-2.65,500,kBlue, "active vol", legXZ_tr);
  DrawNamedRect(-199.15+25,25,-2.65-25,500-25, kRed, "fiducial vol", legXZ_tr);
  //TPC 1
  DrawRectangle(2.65,0,199.15,500,kBlue);
  DrawRectangle(2.65+25,25,199.15-25,500-25, kRed);
  DrawNamedRect(-260.1,-143.1,260.1,559.6,kYellow, "bounding box", legXZ_tr);
  legXZ_tr->Draw();
        
  TCanvas* shcanv = new TCanvas("shcanv","shcanv",2400,800);
  shcanv->Divide(3,1);
  TLegend* legXY_sh = new TLegend(0.75,0.45,0.95,0.6);
  TLegend* legYZ_sh = new TLegend(0.75,0.45,0.95,0.6);
  TLegend* legXZ_sh = new TLegend(0.75,0.45,0.95,0.6);
  TH2D* showerXY = (TH2D*)nue->Get("showerXY"); 
  TH2D* showerYZ = (TH2D*)nue->Get("showerYZ");
  TH2D* showerXZ = (TH2D*)nue->Get("showerXZ");
  shcanv->cd(1);
  showerXY->Draw();
  //TPC 0
  DrawRectangle(-199.15,-200,-2.65,200,kBlue);//act vol
  DrawRectangle(-199.15+25,-200+25,-2.65-25,200-25, kRed);//fid vol
  //TPC 1
  DrawNamedRect(2.65,-200,199.15,200,kBlue, "active vol", legXY_sh);
  DrawNamedRect(2.65+25,-200+25,199.15-25,200-25, kRed, "fiducial vol", legXY_sh);
  DrawNamedRect(-260.1,-271.15,260.1,271.15,kYellow, "bounding box", legXY_sh);
  legXY_sh->Draw();

  shcanv->cd(2);
  showerYZ->Draw();
  DrawNamedRect(-200,0,200,500,kBlue, "active vol", legYZ_sh);
  DrawNamedRect(-200+25,25,200-25,500-25, kRed, "fiducial vol", legYZ_sh);
  DrawNamedRect(-271.15,-143.1,271.15,559.6,kYellow, "bounding box", legYZ_sh);
  legYZ_sh->Draw();
  
  shcanv->cd(3);
  showerXZ->Draw();
  //TPC 0
  DrawNamedRect(-199.15,0,-2.65,500,kBlue, "active vol", legXZ_sh);
  DrawNamedRect(-199.15+25,25,-2.65-25,500-25, kRed, "fiducial vol", legXZ_sh);
  //TPC 1
  DrawRectangle(2.65,0,199.15,500,kBlue);
  DrawRectangle(2.65+25,25,199.15-25,500-25, kRed);
  DrawNamedRect(-260.1,-143.1,260.1,559.6,kYellow, "bounding box", legXZ_sh);
  legXZ_sh->Draw();
}
