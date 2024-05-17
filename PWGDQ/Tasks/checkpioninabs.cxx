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
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace std;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTrkCompColls>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;

struct checkpioninabs{

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
  float collisionZcut = 10.0f;

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);
  Filter collisionFilter = nabs(aod::collision::posZ) < collisionZcut;
  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::MFTTracks> mfttracksPerCollision = aod::fwdtrack::collisionId;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = -5.0;
  //o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;


  HistogramRegistry registry{
    "registry", 
    {
      {"hcollision", "Collision PosZ;Z (cm)", {HistType::kTH1F, {{3000, -15, 15}}}},
      {"hfwdmftratio", "MCH-MID/MFT Standalone;N", {HistType::kTH1F, {{1000, 0.0, 1}}}},
      {"hnfwd", "Number of forward tracks;N", {HistType::kTH1F, {{10, 0.0, 10}}}},
      {"hnmft", "Number of mft tracks;N", {HistType::kTH1F, {{300, 0.0, 300}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(soa::Filtered<aod::Collisions> const& collisions, MyMuons const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, aod::MFTTracks const& mfttracks)
  {
    for (auto& collision : collisions) {
       registry.fill(HIST("hcollision"), collision.posZ());
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      auto mfttracksThisCollision = mfttracks.sliceBy(mfttracksPerCollision, collision.globalIndex());
      int numberofmft = 0;
      int numberoffwd = 0;
      double fwdmftratio = 0.0000;
      for (auto& mfttrack : mfttracksThisCollision) {
        if (mfttrack.eta() > -3.6 && mfttrack.eta() < -2.5) {
          numberofmft++;
        }
      }
      for (auto& fwdtrackId : muonIdsThisCollision) {
        for (auto& fwdtrack : fwdtracks) {
          if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack && fwdtrackId.fwdtrackId() == fwdtrack.globalIndex() && fwdtrack.chi2MatchMCHMFT() < 50 && fwdtrack.compatibleCollIds().size() == 1) {
            if (fwdtrack.eta() > -3.6 && fwdtrack.eta() < -2.5 && (((17.6 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 26.5) && (fwdtrack.pDca() < 594)) || ((26.5 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 89.5) && (fwdtrack.pDca() < 324)))){
              numberoffwd++;
            }
          }
        }
      }
     if (numberofmft == 0) continue;
     registry.fill(HIST("hnmft"), numberofmft);
     registry.fill(HIST("hnfwd"), numberoffwd);
     fwdmftratio = (double)numberoffwd/(double)numberofmft;
     registry.fill(HIST("hfwdmftratio"), fwdmftratio);
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<checkpioninabs>(cfgc)
  };
}
