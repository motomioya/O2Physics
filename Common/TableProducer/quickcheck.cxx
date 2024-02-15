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
#include "Common/DataModel/CollisionAssociationTables.h"
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
using namespace o2::soa;
using o2::globaltracking::MatchingFunc_t;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using MyMuonsColl = soa::Join<aod::FwdTracks, aod::FwdTrkCompColls>;
using MyMFTsColl = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;

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
      {"AllTrueDeltaCollId", "AllTrueDeltaCollId", {HistType::kTH1F, {{4001, -2000.5, 2000.5}}}},
      {"AfterAssociationTrueDeltaCollId", "AfterAssociationTrueDeltaCollId", {HistType::kTH1F, {{4001, -2000.5, 2000.5}}}},
      {"OnlyAssociatedTrueDeltaCollId", "OnlyAssociatedTrueDeltaCollId", {HistType::kTH1F, {{4001, -2000.5, 2000.5}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions const& collisions, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::MFTTrackAssoc const& mfttrackIndices)
  {
    for (auto& [fwdtrack, mfttrack] : combinations(o2::soa::CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        if (fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
          auto fwdparticle = fwdtrack.mcParticle();
          auto mftparticle = mfttrack.mcParticle();

          if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
            registry.fill(HIST("AllTrueDeltaCollId"), fwdtrack.collisionId() - mfttrack.collisionId());
          }
        }
      }
    }

    for (auto& [fwdtrackId, mfttrackId] : combinations(o2::soa::CombinationsFullIndexPolicy(fwdtrackIndices, mfttrackIndices))) {
      auto fwdtrack = fwdtrackId.fwdtrack_as<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>>();
      auto mfttrack = mfttrackId.mfttrack_as<soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>>();
      if (-4.0 < fwdtrack.eta() && fwdtrack.eta() < -2.5) {
        if ((17.6 < fwdtrack.rAtAbsorberEnd() && fwdtrack.rAtAbsorberEnd() < 26.5 && fwdtrack.pDca() < 594) || (26.5 < fwdtrack.rAtAbsorberEnd() && fwdtrack.rAtAbsorberEnd() < 89.5 && fwdtrack.pDca() < 324)) {
          if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
            if (fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
              auto fwdparticle = fwdtrack.mcParticle();
              auto mftparticle = mfttrack.mcParticle();

              if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                registry.fill(HIST("AfterAssociationTrueDeltaCollId"), fwdtrackId.collisionId() - mfttrackId.collisionId());
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
    adaptAnalysisTask<quickcheck>(cfgc)
  };
}
