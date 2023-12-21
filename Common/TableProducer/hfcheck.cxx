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
using MyMuonsColl = soa::Join<aod::FwdTracks, aod::FwdTrkCompColls>;
using MyMFTsColl = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;

struct hfcheck {
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
      {"hMuInclusivePt", "hMuInclusivePt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuFromHfPt", "hMuFromHfPt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuFromHadronPt", "hMuFromHadronPt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuFromHadronZpos", "hMuFromHadronZpos", {HistType::kTH1F, {{1000,-1000,0}}}},
      {"hMuFromMother", "hMuFromMother", {HistType::kTH1F, {{17, -8.5, 8.5}}}},
      {"hConter", "hConter", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions::iterator const& collision, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, aod::McParticles const& mcparticles)
  {
    for (auto& fwdtrack : fwdtracks) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        if (fwdtrack.has_mcParticle()){
          registry.fill(HIST("hMuInclusivePt"), fwdtrack.pt());
          registry.fill(HIST("hMuFromMother"), 0);
          auto fwdparticle = fwdtrack.mcParticle();
          const auto mcfwdtrackmothers = fwdtrack.mcParticle().mothersIds();
          registry.fill(HIST("hConter"), 0);
          if (mcfwdtrackmothers.size() > 0) {
            registry.fill(HIST("hConter"), 1);
            int mcfwdfirstmotherid = mcfwdtrackmothers[0];
            auto mcfwdfirstmother = mcparticles.iteratorAt(mcfwdfirstmotherid);
            //int mcfwdfirstmotherpdg = mcfwdfirstmother.pdgCode();
            int mcfwdlastmotherid = mcfwdtrackmothers[mcfwdtrackmothers.size() - 1];
            auto mcfwdlastmother = mcparticles.iteratorAt(mcfwdlastmotherid);
            int mcfwdlastmotherpdg = mcfwdlastmother.pdgCode();
            if (mcfwdfirstmother.fromBackgroundEvent() == false) {
              registry.fill(HIST("hConter"), 2);
              if (mcfwdlastmotherpdg >= 400 && mcfwdlastmotherpdg <= 439) {
                registry.fill(HIST("hMuFromHfPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromMother"), 1);
                registry.fill(HIST("hConter"), 3);
              } else if (mcfwdlastmotherpdg <= -400 && mcfwdlastmotherpdg >= -439) {
                registry.fill(HIST("hMuFromHfPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromMother"), -1);
                registry.fill(HIST("hConter"), 3);
              }
            } else {
              registry.fill(HIST("hConter"), 4);
              if (mcfwdlastmotherpdg >= 100 && mcfwdlastmotherpdg <= 119) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), 2);
                registry.fill(HIST("hConter"), 5);
              } else if (mcfwdlastmotherpdg <= -100 && mcfwdlastmotherpdg >= -119) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), -2);
                registry.fill(HIST("hConter"), 5);
              }
              if (mcfwdlastmotherpdg >= 1000 && mcfwdlastmotherpdg <= 1999) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), 3);
                registry.fill(HIST("hConter"), 5);
              } else if (mcfwdlastmotherpdg <= -1000 && mcfwdlastmotherpdg >= -1999) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), -3);
                registry.fill(HIST("hConter"), 5);
              }
              if (mcfwdlastmotherpdg >= 200 && mcfwdlastmotherpdg <= 299) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), 5);
                registry.fill(HIST("hConter"), 5);
              } else if (mcfwdlastmotherpdg <= -200 && mcfwdlastmotherpdg >= -299) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), -5);
                registry.fill(HIST("hConter"), 5);
              }
              if (mcfwdlastmotherpdg >= 2000 && mcfwdlastmotherpdg <= 2999) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), 6);
                registry.fill(HIST("hConter"), 5);
              } else if (mcfwdlastmotherpdg <= -2000 && mcfwdlastmotherpdg >= -2999) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), -6);
                registry.fill(HIST("hConter"), 5);
              }
              if (mcfwdlastmotherpdg >= 300 && mcfwdlastmotherpdg <= 399) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), 7);
                registry.fill(HIST("hConter"), 5);
              } else if (mcfwdlastmotherpdg <= -300 && mcfwdlastmotherpdg >= -399) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), -7);
                registry.fill(HIST("hConter"), 5);
              }
              if (mcfwdlastmotherpdg >= 3000 && mcfwdlastmotherpdg <= 3999) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), 8);
                registry.fill(HIST("hConter"), 5);
              } else if (mcfwdlastmotherpdg <= -3000 && mcfwdlastmotherpdg >= -3999) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hMuFromMother"), -8);
                registry.fill(HIST("hConter"), 5);
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
    adaptAnalysisTask<hfcheck>(cfgc)
  };
}
