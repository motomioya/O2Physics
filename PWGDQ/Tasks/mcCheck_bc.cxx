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

struct mccheckambiguity {

  HistogramRegistry registry{
    "registry",
    {
      {"hDeltaMcBcIdMuonFromPiplusReassoc", ";Delta McBcId;entries", {HistType::kTH1F, {{51, -25.5, 25.5}}}},
      {"hDeltaMcBcIdMuonFromLFReassoc", ";Delta McBcId;entries", {HistType::kTH1F, {{51, -25.5, 25.5}}}},
      {"hDeltaMcBcIdMuonFromCharmReassoc", ";Delta McBcId;entries", {HistType::kTH1F, {{51, -25.5, 25.5}}}},
      {"hDeltaMcBcIdMuonFromBottomReassoc", ";Delta McBcId;entries", {HistType::kTH1F, {{51, -25.5, 25.5}}}},
    }};

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(MyEvents const& collisions, MyMuons const& muons, aod::FwdTrackAssoc const& muonAssocs, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles, aod::BCs)
  {
    MCProng prongpiplus(2, {13, 211}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCSignal* signalpiplus;
    signalpiplus = new MCSignal("muFromPiPlus", "Electrons from piplus decays", {prongpiplus}, {-1});

    MCProng prongeta(2, {13, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCSignal* signaleta;
    signaleta = new MCSignal("muFromEta", "Electrons from eta decays", {prongeta}, {-1});

    MCProng prongetaprime(2, {13, 331}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCSignal* signaletaprime;
    signaletaprime = new MCSignal("muFromEtaprime", "Electrons from etaprime decays", {prongetaprime}, {-1});

    MCProng prongrho(2, {13, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCSignal* signalrho;
    signalrho = new MCSignal("muFromRho", "Electrons from rho decays", {prongrho}, {-1});

    MCProng prongomega(2, {13, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCSignal* signalomega;
    signalomega = new MCSignal("muFromOmega", "Electrons from omega decays", {prongomega}, {-1});

    MCProng prongphi(2, {13, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCSignal* signalphi;
    signalphi = new MCSignal("muFromPhi", "Electrons from phi decays", {prongphi}, {-1});

    MCProng prongHc(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongHc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalHc;
    signalHc = new MCSignal("muFromHc", "Electrons from open charmed hadron decays", {prongHc}, {-1});

    MCProng prongHb(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongHb.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalHb;
    signalHb = new MCSignal("muFromHb", "Electrons from open beauty hadron decays", {prongHb}, {-1});

    for (auto& muonAssoc : muonAssocs) {
      auto muon = muonAssoc.template fwdtrack_as<MyMuons>();
      if (muon.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        if (muon.eta() >= -4.0 && muon.eta() <= -2.5 && (((17.6 < muon.rAtAbsorberEnd()) && (muon.rAtAbsorberEnd() < 26.5) && (muon.pDca() < 594)) || ((26.5 < muon.rAtAbsorberEnd()) && (muon.rAtAbsorberEnd() < 89.5) && (muon.pDca() < 324)))){
          auto collisionassoc = muonAssoc.template collision_as<MyEvents>();
          if (!muon.has_mcParticle()) continue;
          if (!muon.has_bcs()) continue;
          auto bc = muon.bc();
          auto particle = muon.mcParticle();
          if (!particle.has_bcs()) continue;
          auto particlebc = particle.bc();
          if((*signalpiplus).CheckSignal(true, particle)) {
            registry.fill(HIST("hDeltaMcCollIdMuonFromPiplusReassoc"), bc.globalBC() - particlebc.globalBC());
          }
          if((*signaleta).CheckSignal(true, particle) || (*signaletaprime).CheckSignal(true, particle) || (*signalrho).CheckSignal(true, particle) || (*signalomega).CheckSignal(true, particle) || (*signalphi).CheckSignal(true, particle)) {
            registry.fill(HIST("hDeltaMcCollIdMuonFromLFReassoc"), bc.globalBC() - particlebc.globalBC());
          }
          if((*signalHc).CheckSignal(true, particle)) {
            registry.fill(HIST("hDeltaMcCollIdMuonFromCharmReassoc"), bc.globalBC() - particlebc.globalBC());
          }
          if((*signalHb).CheckSignal(true, particle)) {
            registry.fill(HIST("hDeltaMcCollIdMuonFromBottomReassoc"), bc.globalBC() - particlebc.globalBC());
          }
        }
      }
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mccheckambiguity>(cfgc)
  };
}
