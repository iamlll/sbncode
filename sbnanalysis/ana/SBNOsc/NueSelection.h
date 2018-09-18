#ifndef __sbnanalysis_ana_SBNOsc_NueSelection__
#define __sbnanalysis_ana_SBNOsc_NueSelection__

/**
 * \file NueSelection.h
 *
 * SBN nue selection.
 *
 * Author: Lisa Lin
 */

#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "lardataobj/MCBase/MCShower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include <random>
#include <fstream>

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

  std::mt19937 rng; //!< Mersenne Twister random number generator
  std::ofstream myfile; //!< general file of gamma processes

  /** Configuration parameters */
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
};
  art::InputTag fTrackTag; //MCTrack
  art::InputTag fShowerTag;//MCShower
  double fdEdx;

  /** vectors */
  std::vector<std::string> gammaprocess; //!< vector of photon shower processes
  std::vector<int> parentcodes; //!< vector of neutrino parent PDG codes
  std::vector<sim::MCTrack> fRelTracks; //MCTracks within 5 cm of neutrino interaction vertex
  std::vector<sim::MCShower> fRelShowers; //MCShowers within 5 cm of neutrino interaction vertex

  /** Histograms of cuts */
  std::vector<TH1D*> fCut1; //cuts for selection criteria nu_e CC electron with E_e > 200 MeV
  THStack* cut1stack;
  std::vector<TH1D*> fCut2;
  THStack* cut2stack;
  std::vector<TH1D*> fCut3;
  THStack* cut3stack;

  std::vector<TH1D*> fCut1_reco; //Cut 1 with reconstructed neutrino energy
  THStack* cut1stack_reco;
  std::vector<TH1D*> cut2a;
  THStack* cut2a_stack;
  std::vector<TH1D*> cut2b;
  THStack* cut2b_stack;
  std::vector<TH1D*> pi2gamma; //ct vs photon shower energy
  THStack* gamma2Stack; 
  TH1D* distgamma; //!< distances of primary photon showers and photon showers from primary pi0s from the neutrino interaction vertex
  std::vector<TH1D*> visiblevertex;
  THStack* vertexstack;
  std::vector<TH1D*> dEdx_gamma;
  THStack* dEdx_gammastack; 

  std::vector<TH1D*> fig11;
  THStack* fig11stack;
  TH1D* dEdx;
  std::vector<TH1D*> dEdx_2; //differentiates b/w gamma and e- showers
  THStack* dEdx_2_stack;

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelection__

