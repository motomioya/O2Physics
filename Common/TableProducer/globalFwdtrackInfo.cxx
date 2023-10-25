// Copyright 2019-2020 CERN and copyright holds of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include <math.h>
#include <TLorentzVector.h>
#include <string>
#include <regex>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;

struct globalFwdtrackInfo {

  float etalow = -4;
  float etaup = -2.5;
  float pDCAcutrAtBsorberEndlow1 = 17.6;
  float pDCAcutrAtBsorberEndup1 = 26.5;
  float pDCAcutrAtBsorberEndlow2 = 26.5;
  float pDCAcutrAtBsorberEndup2 = 89.5;
  float pDCAcutdcaup1 = 594;
  float pDCAcutdcaup2 = 324;
  float chi2up = 1000000;
  float chi2MatchMCHMIDup = 1000000;
  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  HistogramRegistry registry{
    "registry", //
    {
      {"Pt_vs_chi2mftmch", "Pt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"SameCollisionPt_vs_chi2mftmch", "SameCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"DifferentCollisionPt_vs_chi2mftmch", "DifferentCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"chi2mft_vs_chi2mftmch", "chi2mft_vs_chi2mftmch", {HistType::kTH2F, {{400,0,200},{400,0,200}}}},
      {"chi2global_vs_chi2mftmch", "chi2mft_vs_chi2mftmch", {HistType::kTH2F, {{400,0,200},{400,0,200}}}},
      {"chi2mchmid_vs_chi2mftmch", "chi2mft_vs_chi2mftmch", {HistType::kTH2F, {{400,0,200},{400,0,200}}}},
      {"Pt_vs_phi_vs_chi2mftmch", "Pt_vs_phi_vs_chi2mftmch", {HistType::kTH3F, {{200,0,20},{100,-5,5},{400,0,200}}}},
      {"Pt_vs_Zvtx_vs_chi2mftmch", "Pt_vs_Zvtx_vs_chi2mftmch", {HistType::kTH3F, {{200,0,20},{200,-10,10},{400,0,200}}}},
      {"Pt_vs_chi2mft_vs_chi2mftmch", "Pt_vs_chi2mft_vs_chi2mftmch", {HistType::kTH3F, {{200,0,20},{400,0,200},{400,0,200}}}},
      {"Pt_vs_chi2mchmid_vs_chi2mftmch", "Pt_vs_chi2mchmid_vs_chi2mftmch", {HistType::kTH3F, {{200,0,20},{400,0,200},{400,0,200}}}},
      {"Pt_vs_chi2global_vs_chi2mftmch", "Pt_vs_chi2global_vs_chi2mftmch", {HistType::kTH3F, {{200,0,20},{400,0,200},{400,0,200}}}},

      {"McPt_vs_chi2mftmch", "McPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McTruePt_vs_chi2mftmch", "McTruePt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McTrueSameCollisionPt_vs_chi2mftmch", "McTrueSameCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McTrueDifferentCollisionPt_vs_chi2mftmch", "McTrueDifferentCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McTrueDifferentCollisiontime_vs_chi2mftmch", "McTrueDifferentCollisiontime_vs_chi2mftmch", {HistType::kTH2F, {{20000, -100, 100},{400,0,200}}}},
      {"McFalsePt_vs_chi2mftmch", "McFalsePt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McFalseSameCollisionPt_vs_chi2mftmch", "McFalseSameCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McFalseDifferentCollisionPt_vs_chi2mftmch", "McFalseDifferentCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McFalseDifferentCollisiontime_vs_chi2mftmch", "McFalseDifferentCollisiontime_vs_chi2mftmch", {HistType::kTH2F, {{20000, -100, 100},{400,0,200}}}},

      {"McMbPt_vs_chi2mftmch", "McMbPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McMbTruePt_vs_chi2mftmch", "McMbTruePt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McMbTrueSameCollisionPt_vs_chi2mftmch", "McMbTrueSameCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McMbTrueDifferentCollisionPt_vs_chi2mftmch", "McMbTrueDifferentCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McMbTrueDifferentCollisiontime_vs_chi2mftmch", "McMbTrueDifferentCollisiontime_vs_chi2mftmch", {HistType::kTH2F, {{20000, -100, 100},{400,0,200}}}},
      {"McMbFalsePt_vs_chi2mftmch", "McMbFalsePt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McMbFalseSameCollisionPt_vs_chi2mftmch", "McMbFalseSameCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McMbFalseDifferentCollisionPt_vs_chi2mftmch", "McMbFalseDifferentCollisionPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McMbFalseDifferentCollisiontime_vs_chi2mftmch", "McMbFalseDifferentCollisiontime_vs_chi2mftmch", {HistType::kTH2F, {{20000, -100, 100},{400,0,200}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions const& collisions, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto& fwdtrack: fwdtracks) {
      if (fwdtrack.has_collision()){
        if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack){
          registry.fill(HIST("Pt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
          registry.fill(HIST("chi2global_vs_chi2mftmch"),fwdtrack.chi2(),fwdtrack.chi2MatchMCHMFT());
          registry.fill(HIST("chi2mchmid_vs_chi2mftmch"),fwdtrack.chi2MatchMCHMID(),fwdtrack.chi2MatchMCHMFT());
          registry.fill(HIST("Pt_vs_chi2global_vs_chi2mftmch"), fwdtrack.pt(), fwdtrack.chi2(), fwdtrack.chi2MatchMCHMFT());
          registry.fill(HIST("Pt_vs_chi2mchmid_vs_chi2mftmch"), fwdtrack.pt(), fwdtrack.chi2MatchMCHMID(), fwdtrack.chi2MatchMCHMFT());
          registry.fill(HIST("Pt_vs_phi_vs_chi2mftmch"), fwdtrack.pt(), fwdtrack.phi(), fwdtrack.chi2MatchMCHMFT());
          registry.fill(HIST("Pt_vs_Zvtx_vs_chi2mftmch"), fwdtrack.pt(), fwdtrack.collision().posZ(), fwdtrack.chi2MatchMCHMFT());
          for (auto& mfttrack: mfttracks){
            if (mfttrack.has_collision()){
              if (fwdtrack.matchMFTTrackId() == mfttrack.globalIndex()){
                registry.fill(HIST("chi2mft_vs_chi2mftmch"), mfttrack.chi2(),fwdtrack.chi2MatchMCHMFT());
                registry.fill(HIST("Pt_vs_chi2mft_vs_chi2mftmch"), fwdtrack.pt(), mfttrack.chi2(), fwdtrack.chi2MatchMCHMFT());
              }
            }
          }
        }
      }
    }

  }

  void processGen(aod::Collisions const& collisions, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    for (auto& fwdtrack: fwdtracks) {
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack){
        for (auto& mfttrack: mfttracks){
          if (fwdtrack.has_collision() && mfttrack.has_collision() && fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
            auto fwdparticle = fwdtrack.mcParticle();
            auto mftparticle = mfttrack.mcParticle();
//            auto fwdbc = fwdtrack.collision().bc_as<aod::BCsWithTimestamps>();
 //           auto mftbc = mfttrack.collision().bc_as<aod::BCsWithTimestamps>();
            if (fwdparticle.fromBackgroundEvent() == 1 && mftparticle.fromBackgroundEvent() == 1) {
              if (fwdtrack.matchMFTTrackId() == mfttrack.globalIndex()){
                registry.fill(HIST("McPt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
                if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                  registry.fill(HIST("McTruePt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
                  if (fwdtrack.collisionId() == mfttrack.collisionId() ) {
                    registry.fill(HIST("McTrueSameCollisionPt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
                  } else {
                    registry.fill(HIST("McTrueDifferentCollisionPt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
                    registry.fill(HIST("McTrueDifferentCollisiontime_vs_chi2mftmch"), fwdtrack.collision().collisionTime() - mfttrack.collision().collisionTime(),fwdtrack.chi2MatchMCHMFT());
                  }
                } else {
                  registry.fill(HIST("McFalsePt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
                  if (fwdtrack.collisionId() == mfttrack.collisionId() ) {
                    registry.fill(HIST("McFalseSameCollisionPt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
                  } else {
                    registry.fill(HIST("McFalseDifferentCollisionPt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
                    registry.fill(HIST("McFalseDifferentCollisiontime_vs_chi2mftmch"), fwdtrack.collision().collisionTime() - mfttrack.collision().collisionTime(),fwdtrack.chi2MatchMCHMFT());
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(globalFwdtrackInfo, processGen, "Process generator-level info", false);
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<globalFwdtrackInfo>(cfgc)
  };
}
