#ifndef __sbnanalysis_ana_SBNOsc_NueSelection__
#define __sbnanalysis_ana_SBNOsc_NueSelection__

/**
 * \file NueSelection.h
 *
 * SBN nue selection.
 *
 * Author: 
 */

#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "lardataobj/MCBase/MCShower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include <random>

//Forward declarations
class TH2D;
class TH1D;
class THStack;

namespace ana {
  namespace SBNOsc {

/**
 * \class NueSelection
 * \brief Electron neutrino event selection
 */
class NueSelection : public core::SelectionBase {
public:
  /** Constructor. */
  NueSelection();

  /**
   * Initialization.
   *
   * \param config A configuration, as a JSON object
   */
  void Initialize(Json::Value* config=NULL);

  /** Finalize and write objects to the output file. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param reco Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco);

protected:
  unsigned fEventCounter;  //!< Count processed events
  unsigned fNuCount;  //!< Count num nu_e CCQE events
  unsigned fTotTracksShowers; //!<counts total num tracks and showers

  std::mt19937 rng; //!<random number generator for initializing vectors of histograms

  /** Configuration parameters */
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
};
  art::InputTag fTrackTag; //!< art tag for MCTrack information
  art::InputTag fShowerTag;//!< art tag for MCShower information

  /** vectors */
  std::vector<sim::MCTrack> fRelTracks; //!< MCTracks within 5 cm of neutrino interaction vertex
  std::vector<sim::MCShower> fRelShowers; //!< MCShowers within 5 cm of neutrino interaction vertex
  std::vector<int>  fCodes; //!< unique neutrino interaction codes

  /** Min + max distances from neutrino vertex */
  double fMinDist_e; //!< minimum distance of electron showers from vertex
  double fMaxDist_e; //!< maximum distance of electron showers from vertex
  double fMinDist_g; //!< minimum distance of photon showers from vertex
  double fMaxDist_g; //!< maximum distance of photon showers from vertex

  /**histograms!*/
  std::vector<TH1D*> prelimCuts; //!< preliminary cut hist
  THStack* prelim_stack_nue; //!< stacked histogram for preliminary cuts

  TH1D* dist_from_vertex; //!< 1D histogram of "matched" track + shower distances from neutrino interaction vertex
  std::vector<TH1D*> vertexDist_truth; //!< truth distances of electron and photon showers from neutrino interaction vertex
  THStack* truthVD_stack; //!< stacked histogram for differentiating e- & gamma shower distances from vertex

  TH2D* nuE_vs_reco; //!< 2D histogram plotting neutrino energy vs reconstructed energy using CCQE interaction equation
  std::vector<TH2D*> nuereco_type; //!< differentiating between different neutrino interaction types to tell accuracy of CCQE eqn for reconstructed neutrino energies
  THStack* nuereco_stack; //!< stacked histogram for neutrino energy vs reconstructed energy plots, separated by neutrino interaction types
  std::vector<TH1D*> showerE; //!< histogram of how many showers are "associated" with neutrino interactions (w/in 5 cm) 
  THStack* showerE_stack; //!< stacked histogram for shower associations

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelection__

