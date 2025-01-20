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
#include "TDatabasePDG.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include <math.h>
#include <TLorentzVector.h>
#include <string>
#include <regex>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace std;
using MyEvents = soa::Join<aod::Collisions, aod::McCollisionLabels>;
using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;

struct mcCheck_cc {

  HistogramRegistry registry{
    "registry",
    {
      {"hMuonRecPtFromD0", "hMuonRecPtFromD0;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hMuonGenPtFromD0", "hMuonGenPtFromD0;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hMuonResPtFromD0", "hMuonResPtFromD0;VtxCol;entries", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hMuonResPhiFromD0", "hMuonResPhiFromD0;VtxCol;entries", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hMuonRecVtxColFromD0", "hMuonRecVtxColFromD0;VtxCol;entries", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hMuonGenVtxColFromD0", "hMuonGenVtxColFromD0;VtxCol;entries", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hMuonColPosZFromD0", "hMuonColPosZFromD0;VtxCol;entries", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hMuonMCColPosZFromD0", "hMuonMCColPosZFromD0;VtxCol;entries", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hMuonparproPosZFromD0", "hMuonparproPosZFromD0;VtxCol;entries", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hMuonColMCColPosZFromD0", "hMuonColMCColPosZFromD0;VtxCol;entries", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
    }};

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(MyEvents const& collisions, MyMuons const& muons, aod::FwdTrackAssoc const& muonAssocs, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    MCProng prongD0(2, {13, 421}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalD0;
    signalD0 = new MCSignal("muon from DO", "Muons from D0 decays", {prongD0}, {-1});

    for (auto& muonAssoc : muonAssocs) {
      auto muon = muonAssoc.template fwdtrack_as<MyMuons>();
      if (muon.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        if (muon.eta() >= -3.6 && muon.eta() <= -2.5 && (((17.6 < muon.rAtAbsorberEnd()) && (muon.rAtAbsorberEnd() < 26.5) && (muon.pDca() < 594)) || ((26.5 < muon.rAtAbsorberEnd()) && (muon.rAtAbsorberEnd() < 89.5) && (muon.pDca() < 324)))){
          if (muon.chi2MatchMCHMFT() < 40) {
            auto collisionassoc = muonAssoc.template collision_as<MyEvents>();
            //if (collisionassoc.posZ() > 10 || collisionassoc.posZ() < -10) continue;
            if (!muon.has_mcParticle()) continue;
            if (!muon.has_collision()) continue;
            auto collision = muon.template collision_as<MyEvents>();
            auto particle = muon.mcParticle();
            auto particlecollision = particle.mcCollision();
            if (!collision.has_mcCollision()) continue;
            if (!collisionassoc.has_mcCollision()) continue;
            auto mccollision = collision.mcCollision();
            auto mccollisionassoc = collisionassoc.mcCollision();
            if (particle.mcCollisionId() == mccollisionassoc.globalIndex()) {
              if((*signalD0).CheckSignal(true, particle)) {
                registry.fill(HIST("hMuonRecPtFromD0"), muon.pt());
                registry.fill(HIST("hMuonGenPtFromD0"), particle.pt());
                registry.fill(HIST("hMuonResPtFromD0"), muon.pt() - particle.pt());
                registry.fill(HIST("hMuonResPhiFromD0"), muon.phi() - particle.phi());
                registry.fill(HIST("hMuonColPosZFromD0"), collisionassoc.posZ());
                registry.fill(HIST("hMuonMCColPosZFromD0"), mccollisionassoc.posZ());
                registry.fill(HIST("hMuonparproPosZFromD0"), particle.vz());
                registry.fill(HIST("hMuonRecVtxColFromD0"), collisionassoc.posZ() - particle.vz());
                registry.fill(HIST("hMuonGenVtxColFromD0"), mccollisionassoc.posZ() - particle.vz());
                registry.fill(HIST("hMuonColMCColPosZFromD0"), collisionassoc.posZ() - mccollisionassoc.posZ());
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
    adaptAnalysisTask<mcCheck_cc>(cfgc)
  };
}
