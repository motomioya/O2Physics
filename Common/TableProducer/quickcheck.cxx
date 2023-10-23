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
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/mftmchMatchingML.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "GlobalTracking/MatchGlobalFwd.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct quickcheck {
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

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  HistogramRegistry registry{
    "registry",
    {
      {"TrueSameCollisionBcID", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"FalseSameCollisionBcID", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"TrueDifferentCollisionBcID", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"FalseDifferentCollisionBcID", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"TrueSameCollisionTimestamp", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"FalseSameCollisionTimestamp", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"TrueDifferentCollisionTimestamp", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"FalseDifferentCollisionTimestamp", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"TrueSameCollisionTime", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"FalseSameCollisionTime", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"TrueDifferentCollisionTime", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
      {"FalseDifferentCollisionTime", "TType", {HistType::kTH1F, {{20000, -10000, 10000}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions const& collisions, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks, aod::BCs const&)
  {
    for (auto& fwdtrack: fwdtracks) {
      if (fwdtrack.has_collision()){
        if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack){
          for (auto& mfttrack: mfttracks){
            if (mfttrack.has_collision()){
              if (fwdtrack.matchMFTTrackId() == mfttrack.globalIndex()){
              }
            }
          }
        }
      }
    }
  }

  void processGen(aod::Collisions const& collisions, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    for (auto& [fwdtrack, mfttrack] : combinations(o2::soa::CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        if (fwdtrack.has_collision() && mfttrack.has_collision() && fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
          auto fwdparticle = fwdtrack.mcParticle();
          auto mftparticle = mfttrack.mcParticle();
          auto fwdbc = fwdtrack.collision().bc_as<aod::BCsWithTimestamps>();
          auto mftbc = mfttrack.collision().bc_as<aod::BCsWithTimestamps>();
          if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
            if (fwdtrack.collisionId() == mfttrack.collisionId() ) {
              registry.fill(HIST("TrueSameCollisionBcID"), fwdbc.globalBC() - mftbc.globalBC());
              registry.fill(HIST("TrueSameCollisionTimestamp"), fwdbc.timestamp() - mftbc.timestamp());
              registry.fill(HIST("TrueSameCollisionTime"), fwdbc.timestamp() - mftbc.timestamp() + fwdtrack.collision().collisionTime() - mfttrack.collision().collisionTime());
            } else {
              registry.fill(HIST("TrueDifferentCollisionBcID"), fwdbc.globalBC() - mftbc.globalBC());
              registry.fill(HIST("TrueDifferentCollisionTimestamp"), fwdbc.timestamp() - mftbc.timestamp());
              registry.fill(HIST("TrueDifferentCollisionTime"), fwdbc.timestamp() - mftbc.timestamp() + fwdtrack.collision().collisionTime() - mfttrack.collision().collisionTime());
            }
          } else {
            if (fwdtrack.collisionId() == mfttrack.collisionId() ) {
              registry.fill(HIST("FalseSameCollisionBcID"), fwdbc.globalBC() - mftbc.globalBC());
              registry.fill(HIST("FalseSameCollisionTimestamp"), fwdbc.timestamp() - mftbc.timestamp());
              registry.fill(HIST("FalseSameCollisionTime"), fwdbc.timestamp() - mftbc.timestamp() + fwdtrack.collision().collisionTime() - mfttrack.collision().collisionTime());
            } else {
              registry.fill(HIST("FalseDifferentCollisionBcID"), fwdbc.globalBC() - mftbc.globalBC());
              registry.fill(HIST("FalseDifferentCollisionTimestamp"), fwdbc.timestamp() - mftbc.timestamp());
              registry.fill(HIST("FalseDifferentCollisionTime"), fwdbc.timestamp() - mftbc.timestamp() + fwdtrack.collision().collisionTime() - mfttrack.collision().collisionTime());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(quickcheck, processGen, "Process generator-level info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<quickcheck>(cfgc)
  };
}
