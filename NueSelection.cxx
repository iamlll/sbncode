#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"

namespace ana {
  namespace SBNOsc {

NueSelection::NueSelection() : SelectionBase(), fEventCounter(0), fNuCount(0) {}


void NueSelection::Initialize(Json::Value* config) {
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
  fNuVertexXZHist = new TH2D("nu_vtx_XZ", ";X pos (cm);Z pos (cm)",100, -1000, 1000, 100, -1000, 1000);
  fNuVertexXYHist = new TH2D("nu_vtx_XY", ";X pos (cm);Y pos (cm)",100, -1000, 1000, 100, -1000, 1000);
  fNuVertexYZHist = new TH2D("nu_vtx_YZ", ";Y pos (cm);Z pos (cm)",100, -1000, 1000, 100, -1000, 1000);
  fTrackXY = new TH2D("trackXY", ";X pos (cm);Z pos (cm)",100, -1000, 1000, 100, -1000, 1000);
  fTrackYZ = new TH2D("trackYZ", ";X pos (cm);Y pos (cm)",100, -1000, 1000, 100, -1000, 1000);
  fTrackXZ = new TH2D("trackXZ", ";Y pos (cm);Z pos (cm)",100, -1000, 1000, 100, -1000, 1000);
  fShowerXY = new TH2D("showerXY", ";X pos (cm);Z pos (cm)",100, -1000, 1000, 100, -1000, 1000);
  fShowerYZ = new TH2D("showerYZ", ";X pos (cm);Y pos (cm)",100, -1000, 1000, 100, -1000, 1000);
  fShowerXZ = new TH2D("showerXZ", ";Y pos (cm);Z pos (cm)",100, -1000, 1000, 100, -1000, 1000);

  hello();
}


void NueSelection::Finalize() {
  fOutputFile->cd(); 
  fNuVertexXZHist->Write();
  fNuVertexXYHist->Write();
  fNuVertexYZHist->Write();
  fTrackXY->Write();
  fTrackYZ->Write();
  fTrackXZ->Write();
  fShowerXY->Write();
  fShowerYZ->Write();
  fShowerXZ->Write();
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
  auto const& mctracks = \
    *ev.getValidHandle<std::vector<sim::MCTrack>>(fTrackTag);
  auto const& mcshowers = \
    *ev.getValidHandle<std::vector<sim::MCShower>>(fShowerTag);

  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    fNuVertexXZHist->Fill(mctruth.GetNeutrino().Nu().Vx(), mctruth.GetNeutrino().Nu().Vz());
    fNuVertexXYHist->Fill(mctruth.GetNeutrino().Nu().Vx(), mctruth.GetNeutrino().Nu().Vy());
    fNuVertexYZHist->Fill(mctruth.GetNeutrino().Nu().Vy(), mctruth.GetNeutrino().Nu().Vz());

    if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
      Event::Interaction interaction = TruthReco(mctruth);
      reco.push_back(interaction);
    }
  }

  for(size_t t=0; t<mctracks.size();t++){
    auto mctrack = mctracks.at(t);
    fTrackXY->Fill(mctrack.Start().X(), mctrack.Start().Y(),1);
    fTrackYZ->Fill(mctrack.Start().Y(), mctrack.Start().Z(),1);
    fTrackXZ->Fill(mctrack.Start().X(), mctrack.Start().Z(),1);
  }
  for(size_t s=0; s<mcshowers.size();s++){
    auto mcshower = mcshowers.at(s);
    fShowerXY->Fill(mcshower.Start().X(), mcshower.Start().Y(),1);
    fShowerYZ->Fill(mcshower.Start().Y(), mcshower.Start().Z(),1);
    fShowerXZ->Fill(mcshower.Start().X(), mcshower.Start().Z(),1);    
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

