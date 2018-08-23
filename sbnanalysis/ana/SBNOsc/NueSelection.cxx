#include <iostream>
#include <vector>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include <string>
#include <cstring>
#include <TLorentzVector.h>
#include "TVector3.h"
#include <random>
#include <fstream>

namespace ana {
  namespace SBNOsc {

NueSelection::NueSelection() : SelectionBase(), fEventCounter(0), fNuCount(0) {}

std::vector<TH1D*> ColorFill(std::vector<TH1D*> hists, std::mt19937 rng, bool fill){
  //random number distribution in range [0,11]
  std::uniform_int_distribution<std::mt19937::result_type> color12(0,11);
  for(int i=0;i<hists.size();i++){
    bool color = true; //no repeated colors!
    do{
       auto num = color12(rng);
       if(fill == true){
         if(num==0) hists[i]->SetFillColor(kOrange);
         else if(num==1) hists[i]->SetFillColor(kRed+2);
         else if(num==2) hists[i]->SetFillColor(kPink+1);
         else if(num==3) hists[i]->SetFillColor(kMagenta-4);
         else if(num==4) hists[i]->SetFillColor(kViolet-5);
         else if(num==5) hists[i]->SetFillColor(kBlue);
         else if(num==6) hists[i]->SetFillColor(kAzure+1);
         else if(num==7) hists[i]->SetFillColor(kCyan+2);
         else if(num==8) hists[i]->SetFillColor(kTeal-2);
         else if(num==9) hists[i]->SetFillColor(kGreen+2);
         else if(num==10) hists[i]->SetFillColor(kSpring+8);
         else hists[i]->SetFillColor(kYellow+1);
         if(i!=0){
           for(int k=i-1; k>=0; k--){
             if(hists[i]->GetFillColor()==hists[k]->GetFillColor()){
               color=false;
               break;
             }
             else color=true;           
           }
         }
       }
       else{
         if(num==0) hists[i]->SetLineColor(kOrange);
         else if(num==1) hists[i]->SetLineColor(kRed+2);
         else if(num==2) hists[i]->SetLineColor(kPink+1);
         else if(num==3) hists[i]->SetLineColor(kMagenta-4);
         else if(num==4) hists[i]->SetLineColor(kViolet-5);
         else if(num==5) hists[i]->SetLineColor(kBlue);
         else if(num==6) hists[i]->SetLineColor(kAzure+1);
         else if(num==7) hists[i]->SetLineColor(kCyan+2);
         else if(num==8) hists[i]->SetLineColor(kTeal-2);
         else if(num==9) hists[i]->SetLineColor(kGreen+2);
         else if(num==10) hists[i]->SetLineColor(kSpring+8);
         else hists[i]->SetLineColor(kYellow+1);
         if(i!=0){
           for(int k=i-1; k>=0; k--){
             if(hists[i]->GetLineColor()==hists[k]->GetLineColor()){
               color=false;
               break;
             }
             else color=true;           
           }
         } 
       }
    } while(color==false);
  }
   return hists;
}

std::vector<TH1D*> InitializeHists(int nbins, double lowX, double highX, int size, std::string baseTitle, std::mt19937 rng, bool fill){
   std::vector<TH1D*> hists;
   char buffer[20];
   for(int i=0;i<size;i++){
      char* cstr = new char [baseTitle.length()+1];
      std::strcpy (cstr, baseTitle.c_str());
      std::sprintf(buffer, "%s%d", cstr, i); //in this case, using i to denote 1 (true = events that passed the cut) vs 0 (false)
      hists.push_back(new TH1D(buffer,"",nbins,lowX,highX));
   }
   hists = ColorFill(hists, rng, fill);
   return hists;
}

void WriteHists(std::vector<TH1D*> vec, THStack* stack){
  for(size_t i=0; i< vec.size(); i++){
    stack->Add(vec[i]);
  }
  stack->Write();
}

std::vector<TH2D*> Initialize2DHists(int nbinsx, double lowX, double highX, int nbinsy, double lowY, double highY, int size, std::string baseTitle, std::mt19937 rng){
   std::vector<TH2D*> hists;
   char buffer[20];
   for(int i=0;i<size;i++){
      char* cstr = new char [baseTitle.length()+1];
      std::strcpy (cstr, baseTitle.c_str());
      std::sprintf(buffer, "%s%d", cstr, i); //in this case, using i to denote 1 (true = events that passed the cut) vs 0 (false)
      hists.push_back(new TH2D(buffer,"",nbinsx,lowX,highX, nbinsy, lowY, highY));
   }
   return hists;
}

void Write2DHists(std::vector<TH2D*> vec, THStack* stack){
  for(size_t i=0; i< vec.size(); i++){
    stack->Add(vec[i]);
  }
  stack->Write();
}

void NueSelection::Initialize(Json::Value* config) {
  rng.seed(std::random_device()());

  // Load configuration parameters
  fTruthTag = { "generator" };
  fTrackTag = {"mcreco"};
  fShowerTag = {"mcreco"};

  if (config) {
    fTruthTag = { (*config)["SBNOsc"].get("MCTruthTag", "generator").asString() };
    fTrackTag = { (*config)["ExampleAnalysis"].get("MCTrackTag", "mcreco").asString() };
    fShowerTag = { (*config)["ExampleAnalysis"].get("MCShowerTag", "mcreco").asString() };
  }

  // Make a histogram
  prelimCuts = InitializeHists(30,0,5,5,"nu interactions",rng,false);
  prelim_stack_nue = new THStack("prelim_stack_nue",";E_#nu (GeV);count");

  dist_from_vertex = new TH1D("dist_from_vertex","Shower/track dist. from #nu vertex;Distance (cm);count",60,0,10);
  vertexDist_truth = InitializeHists(25,0,20,3,"nuegc",rng,true);
  truthVD_stack = new THStack("truthVD_stack","dist. from #nu vertex;dist(cm);count");
  vertexDist_reco = InitializeHists(60,0,10,2,"nuegc_reco",rng,false);
  recoVD_stack = new THStack("recoVD_stack","dist.from #nu vertex;dist(cm);count");

  nuE_vs_reco = new TH2D("nuE_vs_reco","truth v. reco E_#nu;E_#nu (GeV);Reconstructed E_#nu (GeV)", 30, 0,5,30,0,5); 
  nuereco_type = Initialize2DHists(30,0,5,60,0,10,10, "interaction type", rng);
  nuereco_stack = new THStack("nuereco_stack","truth v. reco E_#nu;E_#nu (GeV);Reconstructed E_#nu (GeV)");

  showerE = InitializeHists(60,0,1000,2,"assn",rng,false);
  showerE_stack = new THStack("showerE_stack","Shower energies + #nu assns.;E_shower (MeV);count");

  hello();
}

/** Returns whether the particle is within the fiducial volume of the detector */
bool AnybodyHome_SBND(TLorentzVector pos){
   if(((-199.15+25 < pos.X() && pos.X() < -2.65-25) || (2.65+25 < pos.X() && pos.X() < 199.15-25)) && (-200+25 < pos.Y() && pos.Y() < 200-25) && (0+25 < pos.Z() && pos.Z() < 500-25)) return true;
   else return false;
}

double FindDistance(TLorentzVector a, TLorentzVector b){
  auto distance = TMath::Sqrt(TMath::Power(a.X()-b.X(),2) + TMath::Power(a.Y()-b.Y(),2) + TMath::Power(a.Z()-b.Z(),2));
  return distance;
}

/**
 * * reconstructed neutrino energy based on cross-section for CCQE interactions w/ heavy nuclei
 * */
double FindRecoEnergy_nue(sim::MCShower mcshower){
  double m_n = 939.565; //neutron mass, in MeV
  double m_p = 938.272; //proton mass, in MeV
  double E_b = 340; //binding energy to bind with an Argon nucleus, in MeV
  double m_l = 0.510; //electron mass, in MeV
  auto E_l = mcshower.Start().E(); //lepton energy, in MeV
  auto p_l = mcshower.Start().Momentum().Vect().Mag(); //lep momentum, MeV/c
  auto theta_l = mcshower.Start().Position().Vect().Theta(); //azimuthal angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_b-E_l+p_l*TMath::Cos(theta_l));
  reco_energy /= 1000; //convert from MeV to GeV
  return reco_energy;
}

double FindRecoEnergy_nue(sim::MCTrack mctrack){
  double m_n = 939.565; //neutron mass, in MeV
  double m_p = 938.272; //proton mass, in MeV
  double E_b = 340; //binding energy to bind with an Argon nucleus, in MeV
  double m_l = 0.510; //electron mass, in MeV
  auto E_l = mctrack.Start().E(); //lepton energy, in MeV
  auto p_l = mctrack.Start().Momentum().Vect().Mag(); //lep momentum, MeV
  auto theta_l = mctrack.Start().Position().Vect().Theta(); //azimuthal angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_b-E_l+p_l*TMath::Cos(theta_l));
  reco_energy /= 1000; //convert from MeV to GeV
  return reco_energy;
}

double FindRecoEnergy_numu(sim::MCTrack mctrack){
  double m_n = 939.565; //neutron mass, in MeV
  double m_p = 938.272; //proton mass, in MeV
  double E_b = 340; //binding energy to bind with an Argon nucleus, in MeV
  auto m_l = 105.658; //muon mass, in MeV
  auto E_l = mctrack.Start().E(); //lepton energy, in MeV
  auto p_l = mctrack.Start().Momentum().Vect().Mag(); //lep momentum, MeV
  auto theta_l = mctrack.Start().Position().Vect().Theta(); //azimuthal angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_b-E_l+p_l*TMath::Cos(theta_l));
  reco_energy /= 1000;
  return reco_energy;
}

std::vector<sim::MCTrack> FindRelevantTracks(TLorentzVector nuVertex, std::vector<sim::MCTrack> mctracks){
  std::vector<sim::MCTrack> relTracks;
  //iterate through only each neutrino's "associated" tracks + showers (within 5 cm of neutrino interaction vertex)
  //(we're cheating - it's okay)
  for(size_t r=0;r<mctracks.size();r++){
    auto const& mctrack = mctracks.at(r);
    if(FindDistance(mctrack.Start().Position(),nuVertex)<=5){
      relTracks.push_back(mctrack);
    }
  }
  return relTracks;
}

std::vector<sim::MCShower> FindRelevantShowers(TLorentzVector nuVertex, std::vector<sim::MCShower> mcshowers){
  std::vector<sim::MCShower> relShowers;
  //iterate through only each neutrino's "associated" tracks + showers (within 5 cm of neutrino interaction vertex)
    for(size_t s=0;s<mcshowers.size();s++){
      auto const& mcshower = mcshowers.at(s);
      if(FindDistance(mcshower.Start().Position(),nuVertex)<=5){
        //shower associated with neutrino interaction
        showerE[0]->Fill(mcshower.Start().E(), 1);
        relShowers.push_back(mcshower);
      }
      else showerE[1]->Fill(mcshower.Start().E(),1);
    }
    return relShowers;
  }

void PrelimCuts_nue(simb::MCNeutrino nu, std::vector<sim::MCShower> mcshowers, std::vector<sim::MCTrack> mctracks, std::vector<TH1D*> prelimCuts){
  auto nuenergy = nu.Nu().E();
  auto nuVertex = nu.Nu().EndPosition();
  auto nuPdg = nu.Nu().PdgCode();
  auto ccnc = nu.CCNC();
  //make a table!
  prelimCuts[0]->Fill(nuenergy,1); //total interactions generated
  //# interactions in fiducial volume
  if(AnybodyHome_SBND(nuVertex)==true){
    prelimCuts[1]->Fill(nuenergy,1); 
    //# true nu_e CC 
    if(nuPdg==12 && ccnc==0){
      prelimCuts[2]->Fill(nuenergy,1);
      //# nu_e CC interactions with a matched shower, E_shower > 200 MeV
      //if the nu_e interaction has at least one matched shower
      int nMatched = 0;
      int primE = 0;
      for(size_t s=0; s<mcshowers.size();s++){
        auto mcshower = mcshowers.at(s); 
        if(FindDistance(nuVertex,mcshower.Start().Position()) <= 5){
          if(mcshower.Start().E() >= 200){  
            nMatched++;
            //efficiency: what fraction of nu_e CC events are accurately defined by my definition of a "signal event" (have a shower w/in 5 cm, energy above 200 MeV)?    
            if(nMatched == 1){ 
              prelimCuts[3]->Fill(nuenergy,1);
            }

            //Purity: what fraction of my "signals" (nu_e CC with E_shower > 200 MeV) are "true" signal events, i.e. primary electron showers? 
            if(mcshower.PdgCode()==11 && mcshower.Process()=="primary"){
              if(primE==0) prelimCuts[4]->Fill(nuenergy,1);
              primE++;
            }
          }
        }  
      }
    }
  }
}

void DistFromNuVertex(simb::MCNeutrino nu, std::vector<sim::MCShower> mcshowers, std::vector<sim::MCTrack> mctracks, TH1D* dist_from_vertex, std::vector<TH1D*> vertexDist_truth){
  
  for(size_t a=0; a<mctracks.size();a++){
    auto mctrack = mctracks.at(a);
    auto dist = FindDistance(mctrack.Start().Position(), nu.Nu().EndPosition());
    if(dist <= 20){
      dist_from_vertex->Fill(dist, 1);
      if(mctrack.Origin()==simb::kCosmicRay){
        vertexDist_truth[2]->Fill(dist,1);
      }
    }
  }

  for(size_t b=0; b<mcshowers.size();b++){
    auto mcshower = mcshowers.at(b);
    auto dist = FindDistance(mcshower.Start().Position(), nu.Nu().EndPosition());
    if(dist <= 20){
      dist_from_vertex->Fill(dist, 1);
      if(nu.Nu().PdgCode()==12 && nu.CCNC()==0 && mcshower.PdgCode()==11 && mcshower.Process()=="primary"){
        vertexDist_truth[0]->Fill(dist,1); //true nu_e CCNC
      }
      if(mcshower.PdgCode()==22){
        vertexDist_truth[1]->Fill(dist,1); //photons travel some distance away from neutrino vertex; distance has exponential distribution whose characteristic length is radiation length of the medium. For LAr, this characteristic length = 14 cm.
      }
      if(mcshower.Origin()==simb::kCosmicRay){
        vertexDist_truth[2]->Fill(dist,1);
      }
    }
  }
}

/**
 * Return vector containing unique codes of all types of neutrino interactions
 */
std::vector<int> UniqueNuTypes(std::vector<int> codes, simb::MCNeutrino nu){
  auto type = nu.Mode();
  //create vector of unique interaction type codes
  if(codes.empty()){
    codes.push_back(type);
  }
  else{
    if(find(codes.begin(),codes.end(),type)==codes.end()){
      codes.push_back(type);
    }
  }
  return codes;
}

void NuE_vs_RecoE(simb::MCNeutrino nu, std::vector<sim::MCShower> fRelShowers, std::vector<sim::MCTrack> fRelTracks){
  for(size_t a=0; a<fRelTracks.size();a++){
    auto mctrack = fRelTracks.at(a);
    auto reco = FindRecoEnergy_nue(mctrack);
    nuE_vs_reco->Fill(nu.Nu().E(), reco, 1); 
    //CCQE
    if(nu.CCNC()==0 && nu.Mode()==0) nuereco_type[0]->Fill(nu.Nu().E(), reco, 1);
    //NCQE
    else if(nu.CCNC()==1 && nu.Mode()==0) nuereco_type[1]->Fill(nu.Nu().E(), reco, 1);
    //resonance
    else if(nu.CCNC()==0 && nu.Mode()==1) nuereco_type[2]->Fill(nu.Nu().E(), reco, 1);
    else if(nu.CCNC()==1 && nu.Mode()==1) nuereco_type[3]->Fill(nu.Nu().E(), reco, 1);
    //deep inelastic scattering
    else if(nu.CCNC()==0 && nu.Mode()==2) nuereco_type[4]->Fill(nu.Nu().E(), reco, 1);
    else if(nu.CCNC()==1 && nu.Mode()==2) nuereco_type[5]->Fill(nu.Nu().E(), reco, 1);
    //coherent scattering
    else if(nu.CCNC()==0 && nu.Mode()==3) nuereco_type[6]->Fill(nu.Nu().E(), reco, 1);
    else if(nu.CCNC()==1 && nu.Mode()==3) nuereco_type[7]->Fill(nu.Nu().E(), reco, 1);
    //neutrino electron elastic scattering 
    else if(nu.Mode()==5) nuereco_type[8]->Fill(nu.Nu().E(), reco, 1);
    //Meson Exchange Current
    else if(nu.Mode()==10) nuereco_type[9]->Fill(nu.Nu().E(), reco, 1);
  }

  for(size_t b=0; b<fRelShowers.size();b++){
    auto mcshower = fRelShowers.at(b);
    auto reco = FindRecoEnergy_nue(mcshower);
    nuE_vs_reco->Fill(nu.Nu().E(), reco, 1);
    //QE
    if(nu.CCNC()==0 && nu.Mode()==0) nuereco_type[0]->Fill(nu.Nu().E(), reco, 1);
    else if(nu.CCNC()==1 && nu.Mode()==0) nuereco_type[1]->Fill(nu.Nu().E(), reco, 1);
    //resonance
    else if(nu.CCNC()==0 && nu.Mode()==1) nuereco_type[2]->Fill(nu.Nu().E(), reco, 1);
    else if(nu.CCNC()==1 && nu.Mode()==1) nuereco_type[3]->Fill(nu.Nu().E(), reco, 1);
    //deep inelastic scattering
    else if(nu.CCNC()==0 && nu.Mode()==2) nuereco_type[4]->Fill(nu.Nu().E(), reco, 1);
    else if(nu.CCNC()==1 && nu.Mode()==2) nuereco_type[5]->Fill(nu.Nu().E(), reco, 1);
    //coherent scattering
    else if(nu.CCNC()==0 && nu.Mode()==3) nuereco_type[6]->Fill(nu.Nu().E(), reco, 1);
    else if(nu.CCNC()==1 && nu.Mode()==3) nuereco_type[7]->Fill(nu.Nu().E(), reco, 1);
    //neutrino electron elastic scattering 
    else if(nu.Mode()==5) nuereco_type[8]->Fill(nu.Nu().E(), reco, 1);
    //Meson Exchange Current
    else if(nu.Mode()==10) nuereco_type[9]->Fill(nu.Nu().E(), reco, 1);
  }
}

void NueSelection::Finalize() {
  fOutputFile->cd();
  WriteHists(prelimCuts, prelim_stack_nue);
  dist_from_vertex->Write();
  WriteHists(vertexDist_truth, truthVD_stack);
  nuE_vs_reco->Write();
  Write2DHists(nuereco_type, nuereco_stack);
  WriteHists(showerE, showerE_stack); 
  std::cout << "Total num tracks+showers: " << fTotTracksShowers << std::endl;
  std::cout << "Neutrino interaction types: " << std::endl;
  for(int i=0; i<fCodes.size(); i++){
    std::cout << fCodes[i] << std::endl;
  }
}

bool NueSelection::ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco) {
  if (fEventCounter % 1000 == 0) {
    std::cout << "NueSelection: Processing event " << fEventCounter << " "
              << "(" << fNuCount << " neutrinos selected)"
              << std::endl;
  }
  fEventCounter++;

  // Grab a data product from the event
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
  auto const& mcfluxes = *ev.getValidHandle<std::vector<simb::MCFlux>>(fTruthTag);
  auto const& mctracks = \
    *ev.getValidHandle<std::vector<sim::MCTrack>>(fTrackTag);
  auto const& mcshowers = \
    *ev.getValidHandle<std::vector<sim::MCShower>>(fShowerTag);
  fTotTracksShowers += mctracks.size()+mcshowers.size();

  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();
    auto nuenergy = nu.Nu().E();
    auto nuVertex = nu.Nu().EndPosition();
    auto nuPdg = nu.Nu().PdgCode();
    auto ccnc = nu.CCNC();
    if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
      Event::Interaction interaction = TruthReco(mctruth);
      reco.push_back(interaction);
    }

    fRelShowers = FindRelevantShowers(nuVertex, mcshowers);
    fRelTracks = FindRelevantTracks(nuVertex, mctracks);
    PrelimCuts_nue(nu, mcshowers, mctracks, prelimCuts);    
    DistFromNuVertex(nu,mcshowers, mctracks, dist_from_vertex, vertexDist_truth);
    NuE_vs_RecoE(nu, fRelShowers, fRelTracks);  
    //fCodes = UniqueNuTypes(fCodes, nu);
  }

  bool selected = !reco.empty();
  // if neutrino is nue CCQE
  if (selected) {
    fNuCount++;
  }

  return selected;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelection)

