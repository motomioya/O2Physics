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
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace std;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;

struct mccheck {

  Preslice<aod::McParticles> particlePerCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::FwdTrack> fwdtrackIndicesPerCollision = aod::fwdtrack::collisionId;

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(aod::McCollisions const& mcCollisions, aod::McParticles const& mcTracks, aod::Collisions const& collisions, MyMuons const& fwdtracks)
  {
    for (auto& collision : mcCollisions) {
      LOGF(info, "--- mccollision ID = %d ---", collision.globalIndex());
      auto particlethiscollision = mcTracks.sliceBy(particlePerCollision, collision.globalIndex());
      for (auto& particle : particlethiscollision) {
        if (particle.pdgCode() == 511 || particle.pdgCode() == 521 || particle.pdgCode() == 531 || particle.pdgCode() == 541 || particle.pdgCode() == 5112 || particle.pdgCode() == 5122 || particle.pdgCode() == 5232 || particle.pdgCode() == 5132 || particle.pdgCode() == 5332) {
          LOGF(info, "particle loop HF, pdg code = %d", particle.pdgCode());
          if (particle.eta() < -4.3 || particle.eta() > -2.2) continue;
          for (auto& particledaughter : particlethiscollision) {
            if (particledaughter.pdgCode() == 13 || particledaughter.pdgCode() == -13) {
              const auto mcmothersidlist = particledaughter.mothersIds();
              if (mcmothersidlist.size() > 0) {
                for (auto i = 0; i <= mcmothersidlist.size() - 1; i++) {
                  auto mcmother = mcTracks.iteratorAt(mcmothersidlist[i]);
                  if (mcmother.globalIndex() == particle.globalIndex()) {
                    if (particledaughter.pdgCode() == 13) {
                      LOGF(info, "Has mu+");
                    }
                    if (particledaughter.pdgCode() == -13) {
                      LOGF(info, "Has mu-");
                    }
                  }
                }
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
    adaptAnalysisTask<mccheck>(cfgc)
  };
}
