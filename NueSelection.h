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
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco);

protected:
  unsigned fEventCounter;  //!< Count processed events
  unsigned fNuCount;  //!< Count num nu_e CCQE events

  std::mt19937 rng;
  /** Configuration parameters */
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
};
  art::InputTag fTrackTag; //MCTrack
  art::InputTag fShowerTag;//MCShower

  /** vectors **/
  std::vector<sim::MCTrack> fRelTracks; //MCTracks within 5 cm of neutrino interaction vertex
  std::vector<sim::MCShower> fRelShowers; //MCShowers within 5 cm of neutrino interaction vertex
  std::vector<int>  fCodes;

  /**histograms!*/
  std::vector<TH1D*> prelimCuts;
  THStack* prelim_stack_nue;

  TH1D* dist_from_vertex;
  std::vector<TH1D*> vertexDist_truth;
  THStack* truthVD_stack;
  std::vector<TH1D*> vertexDist_reco;
  THStack* recoVD_stack;

  TH2D* nuE_vs_reco;
  std::vector<TH2D*> nuereco_type;
  THStack* nuereco_stack;
  std::vector<TH1D*> showerE;
  THStack* showerE_stack;

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelection__

