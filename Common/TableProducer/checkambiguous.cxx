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

struct checkambiguous {
  float etalow = -3.6;
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
      {"TruePositionZ", "TType", {HistType::kTH1F, {{100000, -50, 50}}}},
      {"FalsePositionZ", "TType", {HistType::kTH1F, {{100000, -50, 50}}}},
      {"TruePositionX", "TType", {HistType::kTH1F, {{10000, -5, 5}}}},
      {"FalsePositionX", "TType", {HistType::kTH1F, {{10000, -5, 5}}}},
      {"TruePositionY", "TType", {HistType::kTH1F, {{10000, -5, 5}}}},
      {"FalsePositionY", "TType", {HistType::kTH1F, {{10000, -5, 5}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions const& collisions, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    for (auto& [fwdtrack, mfttrack] : combinations(o2::soa::CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        if (fwdtrack.has_collision() && mfttrack.has_collision() && fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
          auto fwdparticle = fwdtrack.mcParticle();
          auto mftparticle = mfttrack.mcParticle();
          auto fwdbc = fwdtrack.collision().bc_as<aod::BCsWithTimestamps>();
          auto mftbc = mfttrack.collision().bc_as<aod::BCsWithTimestamps>();
          if (fwdtrack.collisionId() != mfttrack.collisionId()){
            if (fwdbc.timestamp() == mftbc.timestamp()){
              if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                registry.fill(HIST("TruePositionZ"), fwdtrack.collision().posZ()- mfttrack.collision().posZ());
                registry.fill(HIST("TruePositionX"), fwdtrack.collision().posX()- mfttrack.collision().posX());
                registry.fill(HIST("TruePositionY"), fwdtrack.collision().posY()- mfttrack.collision().posY());
              } else {
                registry.fill(HIST("FalsePositionZ"), fwdtrack.collision().posZ()- mfttrack.collision().posZ());
                registry.fill(HIST("FalsePositionX"), fwdtrack.collision().posX()- mfttrack.collision().posX());
                registry.fill(HIST("FalsePositionY"), fwdtrack.collision().posY()- mfttrack.collision().posY());
              }
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<checkambiguous>(cfgc)
  };
}
