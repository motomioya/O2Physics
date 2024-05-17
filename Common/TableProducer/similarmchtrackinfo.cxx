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
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include <math.h>
#include <TLorentzVector.h>
#include <string>
#include <regex>
#include <iostream>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"

#include "DCAFitter/FwdDCAFitterN.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "Math/Vector3D.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace std;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTrkCompColls>;

struct simliarmchtrackinfo {
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = -5.0;
  float collisionZcut = 10.0f;
  Filter collisionFilter = nabs(aod::collision::posZ) < collisionZcut;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }


  //using BCsWithBcSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
  //Preslice<ColEvSels> perFoundBC = aod::evsel::foundBCId;
  void process(soa::Filtered<aod::Collisions> const& collisions, MyMuons const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, o2::aod::MFTTracks const& mfttracks)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      int numberofmuon = 0;
      for (auto& fwdtrackId : muonIdsThisCollision) {
        for (auto& fwdtrack : fwdtracks) {
          if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId.fwdtrackId() == fwdtrack.globalIndex() && fwdtrack.chi2MatchMCHMFT() < 50 && fwdtrack.compatibleCollIds().size() == 1) {
            if (fwdtrack.eta() > -4 && fwdtrack.eta() < -2.5 && (((17.6 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 26.5) && (fwdtrack.pDca() < 594)) || ((26.5 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 89.5) && (fwdtrack.pDca() < 324)))){
              numberofmuon++;
            }
          }
        }
      }
      if (numberofmuon > 1) {
        LOGF(info,"-----------collision-----------");
        for (auto& fwdtrackId : muonIdsThisCollision) {
          for (auto& fwdtrack : fwdtracks) {
            if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId.fwdtrackId() == fwdtrack.globalIndex() && fwdtrack.chi2MatchMCHMFT() < 50 && fwdtrack.compatibleCollIds().size() == 1) {
              if (fwdtrack.eta() > -4 && fwdtrack.eta() < -2.5 && (((17.6 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 26.5) && (fwdtrack.pDca() < 594)) || ((26.5 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 89.5) && (fwdtrack.pDca() < 324)))){
                LOGF(info,"---fwdtrack---");
                LOGF(info,"ID = %d",fwdtrack.globalIndex());
                LOGF(info,"Eta = %f",fwdtrack.eta());
                LOGF(info,"Phi = %f",fwdtrack.phi());
                LOGF(info,"Pt = %f",fwdtrack.pt());
                LOGF(info,"x = %f",fwdtrack.x());
                LOGF(info,"y = %f",fwdtrack.y());
                LOGF(info,"z = %f",fwdtrack.z());
              }
            }
          }
        }
        LOGF(info,"numberofmuon = %d",numberofmuon);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<simliarmchtrackinfo>(cfgc)
  };
}
