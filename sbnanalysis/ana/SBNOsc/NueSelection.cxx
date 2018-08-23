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
  fCut1 = InitializeHists(30, 0, 5, 3, "#nu_e CC", rng, true);  
  cut1stack = new THStack("cut1stack", "Intrinsic/signal #nu_e CC;E_#nu (GeV);count");
  fCut1_reco = InitializeHists(30,0,5,3,"nue1",rng, true);
  cut1stack_reco = new THStack("cut1stack_reco", "Cut 1 with reco E_#nu;Reconstructed E_#nu (GeV);count");
  fCut2 = InitializeHists(30,0,5,2, "events", rng, true);
  cut2stack = new THStack("cut2stack", "NC #gamma production;E_#nu (GeV);count");
  cut2a = InitializeHists(30,0,5,2,"2nd photon", rng, true);
  cut2a_stack = new THStack("cut2a_stack", "2nd photon cut;E_#nu (GeV);count");
  cut2b = InitializeHists(30,0,5,3,"conversion gap", rng, true);
  cut2b_stack = new THStack("cut2b_stack","conversion gap cut;E_#nu (GeV);count");
  pi2gamma = InitializeHists(30,0,500,6,"#gamma", rng, false);
  gamma2Stack = new THStack("2ndPhotonCut_stack", "Photon energies;E_#gamma (MeV);count");
  visiblevertex = InitializeHists(30,0,500,3,"vertex", rng, true);
  vertexstack = new THStack("visiblevertexstack","Visible vertex/conv. gap photons;#gamma shower energy (MeV);count"); 
  dEdx_gamma = InitializeHists(30,0,500,3,"gamma dEdx", rng, true);
  dEdx_gammastack = new THStack("gammastack","#gamma dE/dx cut (cut 2);#gamma shower energy (MeV);count");

  fCut3 = InitializeHists(30,0,5,2,"#nu_mu CC", rng, true);
  cut3stack = new THStack("cut3stack","CC #nu_#mu;E_#nu (GeV);count");


  fig11 = InitializeHists(30,0,5, 5,"thing", rng, true);
  fig11stack = new THStack("fig11stack",";Reconstructed Energy (GeV);Events/GeV");


  dEdx = new TH1D("dEdx","Shower dE/dx;dE/dx (MeV/cm);particle count", 30, 0, 5);
  dEdx_2 = InitializeHists(30,0,5,2,"particle",rng, true);
  dEdx_2_stack = new THStack("showerStack","Shower dE/dx;dE/dx(MeV/cm);count");

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

/*
 * reconstructed neutrino energy based on cross-section for CCQE interactions w/ heavy nuclei
 */
double FindRecoEnergy_nue(sim::MCShower mcshower){
  double m_n = 939.565; //neutron mass, in MeV
  double m_p = 938.272; //proton mass, in MeV
  double E_b = 340; //binding energy to bind with an Argon nucleus, in MeV
  double m_l = 0.510; //electron mass, in MeV
  auto E_l = mcshower.Start().E(); //lepton energy, in MeV
  auto p_l = mcshower.Start().Momentum().Vect().Mag(); //lep momentum, MeV/c
  auto theta_l = mcshower.Start().Position().Vect().Theta(); //angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_l+p_l*TMath::Cos(theta_l));
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
  auto theta_l = mctrack.Start().Position().Vect().Theta(); //angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_l+p_l*TMath::Cos(theta_l));
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
  auto theta_l = mctrack.Start().Position().Vect().Theta(); //angle b/w prod lepton + neutrino (z-axis)
  auto reco_energy = 0.5*(TMath::Power(m_p,2)-TMath::Power(m_n-E_b,2)-TMath::Power(m_l,2)+2*(m_n-E_b)*E_l)/(m_n-E_l+p_l*TMath::Cos(theta_l));
  reco_energy /= 1000;
  return reco_energy;
}

std::vector<sim::MCTrack> FindRelevantTracks(TLorentzVector nuVertex, std::vector<sim::MCTrack> mctracks){
  std::vector<sim::MCTrack> relTracks;
  //iterate through only each neutrino's "associated" tracks (within 5 cm of neutrino interaction vertex)
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
  //iterate through only each neutrino's "associated" showers (within 5 cm of neutrino interaction vertex)
  for(size_t s=0;s<mcshowers.size();s++){
      auto const& mcshower = mcshowers.at(s);
      if(FindDistance(mcshower.Start().Position(),nuVertex)<=5){
        //shower associated with neutrino interaction
        relShowers.push_back(mcshower);
      }
    }
    return relShowers;
  }

int LinkShower_HadronParent(simb::MCFlux mcflux, simb::MCTruth mctruth, sim::MCShower mcshower){
  int parent = 0;
  if(mctruth.GetNeutrino().Lepton().PdgCode()==mcshower.PdgCode()){
    if(mctruth.GetNeutrino().Lepton().Position().Vect()==mcshower.Start().Position().Vect()){
      parent = mcflux.fptype;
    }
  }
  return parent;
}

/*
 * Find active track length (track length w/in fid vol) of the track
 */
double ActiveTrackLength(sim::MCTrack mctrack){
  double tracklength = 0.0;
  for(size_t i=0; i<mctrack.size();i++){
    sim::MCStep& step = mctrack[i];
    if(i>0){
      if(AnybodyHome_SBND(step.Position())==true && AnybodyHome_SBND(mctrack[i-1].Position())==true){
        tracklength += FindDistance(step.Position(),mctrack[i-1].Position());
      }
    }
  }
  return tracklength;
}

/*
 * nu_e CC CUT 1 IMPLEMENTATION: nu_e CC signal, 80% efficiency
 */
void Cut1(simb::MCTruth mctruth, std::vector<sim::MCTrack> fRelTracks,std::vector<sim::MCShower> fRelShowers, std::vector<TH1D*> fCut1, std::mt19937 rng){
  std::uniform_int_distribution<std::mt19937::result_type> dist100(0,99); //distribution in range [0, 99]
  auto CCNC = mctruth.GetNeutrino().CCNC();
  auto energy = mctruth.GetNeutrino().Nu().E();
  if(mctruth.GetNeutrino().Nu().PdgCode()==12 && CCNC==0){
    fCut1[0]->Fill(energy,1); //nue CC interactions
    int shcount = 0; //counts num showers abv threshold = 200 MeV
    int ecount = 0; //counts num primary e- showers "
    for(size_t t=0; t<fRelShowers.size();t++){
      auto mcshower = fRelShowers.at(t);
      //nu_e CC with E_shower > 200 MeV (my defn of signal)
      if(mcshower.Start().E()>200){
        shcount++;
        if(shcount==1){
          fCut1[1]->Fill(energy,1); //only want to fill once per event
        }
        //# interactions producing prim e- sh w/ E_e>200  
        if(mcshower.PdgCode()==11 && mcshower.Process()=="primary"){
          ecount++;
          if(ecount==1){
            fCut1[2]->Fill(energy,1); 
            //implementing 80% identification efficiency
            //if(dist100(rng)<80) fCut1[3]->Fill(energy,1); 
            break; //done gathering info from this event!
          }
        }   
      }
    }
  }
}

void Cut1_reco(simb::MCFlux mcflux, simb::MCTruth mctruth, std::vector<sim::MCShower> mcshowers, std::vector<TH1D*> fCut1, std::vector<TH1D*> fig11, std::mt19937 rng){
  std::uniform_int_distribution<std::mt19937::result_type> dist100(0,99); //distribution in range [0, 99]
  int shcount = 0; 
  int ecount = 0;
  for(size_t i=0;i<mcshowers.size();i++){
    auto mcshower = mcshowers.at(i);
    auto energy = FindRecoEnergy_nue(mcshower);
    //all events w/ at least 1 matched shower above 200 MeV
    if(mcshower.Start().E() > 200){
      shcount++;
      if(shcount==1) fCut1[0]->Fill(energy,1);
      //# events w/ e- sh abv threshold
      if(mcshower.PdgCode()==11 && mcshower.Process()=="primary"){
        ecount++;
        if(ecount==1){
          fCut1[1]->Fill(energy,1);
          if(dist100(rng)<80){
            //implementing ID efficiency
            //fCut1[3]->Fill(energy, 1);
          }
          //muon parent
          if(LinkShower_HadronParent(mcflux,mctruth,mcshower)==13) fig11[0]->Fill(energy, 1);
          //K0 parent
          if(LinkShower_HadronParent(mcflux,mctruth,mcshower)==130 || LinkShower_HadronParent(mcflux,mctruth,mcshower)==310 || LinkShower_HadronParent(mcflux,mctruth,mcshower)==311) fig11[1]->Fill(energy, 1);
          //K+ parent 
          if(LinkShower_HadronParent(mcflux,mctruth,mcshower)==321){
            fig11[2]->Fill(energy, 1);
          }
        }
      } 
    }
  }
} 

/*
 * nu_e CC CUT 2 IMPLEMENTATION: NC photon production
 */
void Cut2(simb::MCNeutrino nu, std::vector<sim::MCTrack> fRelTracks, std::vector<sim::MCShower> fRelShowers, std::vector<TH1D*> hists, std::vector<TH1D*> fig11, std::mt19937 rng){
  std::uniform_int_distribution<std::mt19937::result_type> dist100(0,99); //distribution in range [0, 99]
  auto nuE = nu.Nu().E();
  auto nuVertex = nu.Nu().EndPosition(); 
  int count = 0; //count the number of photon showers produced by daughter pi0 decay (that pass the energy cut)
  double hadronE = 0.0; //total charged hadronic activity produced at nu vertex
  bool visible = false; //neutrino has "visible" interaction vertex if hadronE > 50 MeV at vertex
  bool reject = false; //does the neutrino interaction event pass the cuts?
  //keep events with photon showers above 200 MeV

  cut2a[0]->Fill(nuE,1); //all neutrino flux
  cut2b[0]->Fill(nuE,1);

  int nKids=0; //number of kids
  int momID;
  for(size_t n=0;n<fRelShowers.size();n++){
    auto const& mcshower = fRelShowers.at(n);
    dEdx->Fill(mcshower.dEdx(),1);
    
    if(mcshower.PdgCode()==22){
      pi2gamma[0]->Fill(mcshower.Start().E(),1); //all photon showers w/in 5 cm
      visiblevertex[0]->Fill(mcshower.Start().E(),1);
      dEdx_gamma[0]->Fill(mcshower.Start().E(),1);
      dEdx_2[0]->Fill(mcshower.dEdx(),1);
    }
    if(mcshower.PdgCode()==11) dEdx_2[1]->Fill(mcshower.dEdx(),1);

    //leading photon shower produced by pi0 must have E_gamma > 200 MeV; 2nd photon shower must have E_gamma > 100 MeV
    if(mcshower.PdgCode()==22 && mcshower.MotherPdgCode()==111 && mcshower.MotherProcess()=="primary"){
      if(nKids==0){
        pi2gamma[1]->Fill(mcshower.Start().E(),1); //# leading photon showers
        momID = mcshower.MotherTrackID(); //parent track ID of the first photon child
        if(mcshower.Start().E() > 200){
          pi2gamma[2]->Fill(mcshower.Start().E(), 1); //# primary g-showers that pass the energy cut
          count++;
        }
        nKids++;
      }      
      else if(nKids==1){
        if(mcshower.MotherTrackID()==momID){
          pi2gamma[3]->Fill(mcshower.Start().E(),1); //# 2nd photon showers
          if(count==1 && mcshower.Start().E() > 100){//if leading photon shower pass the energy cut
            pi2gamma[4]->Fill(mcshower.Start().E(), 1); //#2nd g-kids that pass the energy cut
            count++;
            //if the second photon converts (showers) within the active TPC, reject the event!
            if(AnybodyHome_SBND(mcshower.Start().Position())==true){
              pi2gamma[5]->Fill(mcshower.Start().E(),1);
              reject=true;
              cut2a[1]->Fill(nuE,1);//rejected events from 2nd photon cut
              break; 
            } 
          }
        }
      }
    }
  }

  for(size_t r=0;r<fRelTracks.size();r++){
    auto const& mctrack = fRelTracks.at(r);
    //nu_e CC SELECTION CUT 2B
    //Look for charged hadron primary particles, i.e. NOT leptons, protons, pi0, or K0. Find energy total (>50 MeV to be visible)
    if(mctrack.PdgCode()!=11 && mctrack.PdgCode()!=13 && mctrack.PdgCode()!=15 && mctrack.PdgCode()!=111 && mctrack.PdgCode()!=130 && mctrack.PdgCode()!= 310 && mctrack.PdgCode()!=311 && mctrack.Process()=="primary"){
      double anarG = mctrack.Start().E(); //in MeV
      hadronE += anarG;
    }
  }
  if(hadronE > 50){
    visible = true;
    cut2b[1]->Fill(nuE, 1);//num nu events w/ visible vertex
  }
  int nPhotonSh = 0; //total num of photon showers in the nu interaction
  int nfailed = 0; //num of photon showers that converted >3cm from vertex
  if(visible==true){
    for(size_t sh=0; sh<fRelShowers.size();sh++){
      auto const& mcshower = fRelShowers.at(sh);
      auto dist_from_vertex = FindDistance(mcshower.Start().Position(), nuVertex);
      if(mcshower.PdgCode()==22){
        visiblevertex[1]->Fill(mcshower.Start().E(),1); //how many photon showers from visible vertex
        nPhotonSh++;
        if(dist_from_vertex > 3){
          visiblevertex[2]->Fill(mcshower.Start().E(),1); //how many photon showers converting more than 3 cm away from vertex
          nfailed++;        
        } 
      } 
    }
    if(nPhotonSh!=0 && nfailed == nPhotonSh){      
      reject=true; 
      cut2b[2]->Fill(nuE,1);
    }
  } 
  //96% photon rejection rate from dE/dx cut - have to loop back thru all relevant showers
  if(reject==false){
    for(size_t sh=0; sh<fRelShowers.size();sh++){
      auto const& mcshower = fRelShowers.at(sh);
      if(mcshower.PdgCode()==22){
      dEdx_gamma[1]->Fill(mcshower.Start().E(),1); //total photon showers not rejected
        if(mcshower.dEdx()<1.6){
          hists[1]->Fill(nuE,1);
          fig11[3]->Fill(FindRecoEnergy_nue(mcshower),1);
          dEdx_gamma[2]->Fill(mcshower.Start().E(),1); //total # photon showers passing dEdx cut
        }
        else hists[0]->Fill(nuE,1);
      }
    }  
  }
  //if rejecting event, add all photons to "failed" hist
  else{
    for(size_t sh=0; sh<fRelShowers.size();sh++){
      auto const& mcshower = fRelShowers.at(sh);
      if(mcshower.PdgCode()==22) hists[0]->Fill(nuE, 1);  
    }
  } 
} 

/*
 * nu_e CC CUT 3 IMPLEMENTATION
 */
void Cut3(simb::MCTruth mctruth, std::vector<sim::MCTrack> fRelTracks, std::vector<sim::MCShower> fRelShowers, std::vector<TH1D*> hists, std::vector<TH1D*> fig11, std::mt19937 rng){
  if(mctruth.GetNeutrino().Nu().PdgCode()==14 && mctruth.GetNeutrino().CCNC()==0){
    for(size_t t=0; t<fRelTracks.size();t++){
      auto const& mctrack = fRelTracks.at(t);
      if(mctrack.PdgCode()==13){
        //if visible track length > 1 m or track exits (cannot tell what total track length is), then cut
        if(ActiveTrackLength(mctrack)>=100 || AnybodyHome_SBND(mctrack.End().Position())==false){
          hists[0]->Fill(mctruth.GetNeutrino().Nu().E(),1);  
        }
        else{
          int ct=0; //counts how many primary showers are associated w/ vertex
          for(size_t s=0; s<fRelShowers.size();s++){
            auto const& mcshower = fRelShowers.at(s);
            ct++;
          }
          //if only one shower associated with the CC event vertex 
          if(ct==1){
            //run Cut 2 photon selection criteria, events that are not rejected are retained as BG for nu_e CC sample.
            //I guess only cuts 2+3 (dEdx and produced charged hadron activity) are relevant?
            bool reject = false; 
            bool visible = false;
            double hadronE = 0.0;
            for(size_t r=0;r<fRelTracks.size();r++){
              auto const& track = fRelTracks.at(r);
              if(track.PdgCode()!=11 && track.PdgCode()!=13 && track.PdgCode()!=15 && track.PdgCode()!=111 && track.PdgCode()!=130 && track.PdgCode()!= 310 && track.PdgCode()!=311 && track.Process()=="primary"){
                double anarG = track.Start().E(); //in MeV
                hadronE += anarG;
              }
            }
            if(hadronE > 50){
              visible = true;
            }
            int nPhotonSh = 0; //total num of photon showers in the nu interaction
            int nfailed = 0; //num of photon showers that converted >3cm from vertex
            if(visible==true){
              for(size_t sh=0; sh<fRelShowers.size();sh++){
              auto const& mcshower = fRelShowers.at(sh);
              auto dist_from_vertex = FindDistance(mcshower.Start().Position(), mctruth.GetNeutrino().Nu().EndPosition());
                if(mcshower.PdgCode()==22){
                  nPhotonSh++;
                  if(dist_from_vertex > 3){
                    nfailed++; 
                  }  
                } 
              }
              if(nfailed == nPhotonSh){
              reject=true; 
              }
            }        
            if(reject==false){
              for(size_t sh=0; sh<fRelShowers.size();sh++){
                auto const& mcshower = fRelShowers.at(sh);
                if(mcshower.PdgCode()==22){
                  if(mcshower.dEdx()<1.6){
                    hists[1]->Fill(mctruth.GetNeutrino().Nu().E(),1);
                    fig11[4]->Fill(FindRecoEnergy_nue(mcshower),1);
                  }
                  else hists[0]->Fill(mctruth.GetNeutrino().Nu().E(),1);
                }
              }  
            }
            //if rejecting event, add all photons to "failed" hist
            else{
              for(size_t sh=0; sh<fRelShowers.size();sh++){
              auto const& mcshower = fRelShowers.at(sh);
              if(mcshower.PdgCode()==22) hists[0]->Fill(mctruth.GetNeutrino().Nu().E(), 1);  
              }
            }    
          }
        }
      }
    } 
  }
}

void NueSelection::Finalize() {
  fOutputFile->cd();
  WriteHists(fCut1, cut1stack);
  WriteHists(fCut1_reco, cut1stack_reco);

  WriteHists(fCut2, cut2stack); 
  WriteHists(cut2a, cut2a_stack);
  WriteHists(cut2b, cut2b_stack);
  WriteHists(pi2gamma, gamma2Stack);
  WriteHists(visiblevertex, vertexstack);
  WriteHists(dEdx_gamma, dEdx_gammastack);

  WriteHists(fCut3, cut3stack);
  WriteHists(fig11, fig11stack);
  dEdx->Write();
  WriteHists(dEdx_2, dEdx_2_stack);
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

  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;
    auto const& mctruth = mctruths.at(i);
    auto const& mcflux = mcfluxes.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();
    auto nuenergy = nu.Nu().E();
    auto nuVertex = nu.Nu().EndPosition();
    auto nuPdg = nu.Nu().PdgCode();
    auto ccnc = nu.CCNC();
    if (ccnc == simb::kCC && nu.Mode() == 0 && nuPdg == 12) {
      Event::Interaction interaction = TruthReco(mctruth);
      reco.push_back(interaction);
    }

    fRelTracks = FindRelevantTracks(nuVertex, mctracks);
    fRelShowers = FindRelevantShowers(nuVertex, mcshowers);
  
    Cut1(mctruth, fRelTracks, fRelShowers, fCut1, rng);
    Cut1_reco(mcflux, mctruth, mcshowers, fCut1_reco, fig11, rng);
    Cut2(nu, fRelTracks, fRelShowers, fCut2, fig11, rng);
    Cut3(mctruth, fRelTracks, fRelShowers, fCut3, fig11, rng);
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
