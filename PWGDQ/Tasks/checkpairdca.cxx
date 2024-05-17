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
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTrkCompColls, aod::McFwdTrackLabels>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;

struct checkpairdca{

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

  Preslice<aod::McParticles> particlesIndicesPerCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = -5.0;
  //o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;


  HistogramRegistry registry{
    "registry", 
    {
      {"hDCATrueMatch", "DCA (mumu);DCA (cm)", {HistType::kTH1F, {{1000, 0.0, 100}}}},
      {"hDCAFakeMatch", "DCA (mumu);DCA (cm)", {HistType::kTH1F, {{1000, 0.0, 100}}}},
      {"hMFTDCATrueMatch", "DCA (mumu);DCA (cm)", {HistType::kTH1F, {{1000, 0.0, 100}}}},
      {"hMFTDCAFakeMatch", "DCA (mumu);DCA (cm)", {HistType::kTH1F, {{1000, 0.0, 100}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(soa::Filtered<aod::Collisions> const& collisions, MyMuons const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McCollisions const&, aod::McParticles const& particles)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& fwdtrackId : fwdtrackIndices) {
        for (auto& fwdtrack : fwdtracks) {
          if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId.fwdtrackId() == fwdtrack.globalIndex() && fwdtrack.chi2MatchMCHMFT() < 50 && fwdtrack.compatibleCollIds().size() == 1) {
            if (fwdtrack.eta() > -4 && fwdtrack.eta() < -2.5 && (((17.6 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 26.5) && (fwdtrack.pDca() < 594)) || ((26.5 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 89.5) && (fwdtrack.pDca() < 324)))){

              //set TrackParCovFwd for fwdtrack pairs
              double chi21 = fwdtrack.chi2();
              SMatrix5 tpars1(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
              std::vector<double> v11{fwdtrack.cXX(), fwdtrack.cXY(), fwdtrack.cYY(), fwdtrack.cPhiX(), fwdtrack.cPhiY(), fwdtrack.cPhiPhi(), fwdtrack.cTglX(), fwdtrack.cTglY(), fwdtrack.cTglPhi(), fwdtrack.cTglTgl(), fwdtrack.c1PtX(), fwdtrack.c1PtY(), fwdtrack.c1PtPhi(), fwdtrack.c1PtTgl(), fwdtrack.c1Pt21Pt2()};
              SMatrix55 tcovs1(v11.begin(), v11.end());

              o2::track::TrackParCovFwd pars1{fwdtrack.z(), tpars1, tcovs1, chi21};

              //fwdtrack propagation to collision
              pars1.propagateToZlinear(collision.posZ());

              //calculate pair DCA
              float fwd1dcaX = (pars1.getX() - collision.posX());
              float fwd1dcaY = (pars1.getY() - collision.posY());
              float DCA1 = std::sqrt(fwd1dcaX * fwd1dcaX + fwd1dcaY * fwd1dcaY);

              //MFT
              float mft1dcaX = -999;
              float mft1dcaY = -999;
              for (auto& mfttrack : mfttracks) {
                if (fwdtrack.matchMFTTrackId() == mfttrack.globalIndex()) {
                  //propagate mfttrack to matching position
                  double mftchi2 = mfttrack.chi2();
                  SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                  std::vector<double> mftv1;
                  SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                  o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                  //mftpars1.propagateToZlinear(MatchingPlaneZ);
                  mftpars1.propagateToZlinear(collision.posZ());
                  mft1dcaX = (mftpars1.getX() - collision.posX());
                  mft1dcaY = (mftpars1.getY() - collision.posY());
                  float MFTDCA1 = std::sqrt(mft1dcaX * mft1dcaX + mft1dcaY * mft1dcaY);
                  if (fwdtrack.has_mcParticle() == 1 && mfttrack.has_mcParticle() == 1) {
                    auto fwdparticle = fwdtrack.mcParticle();
                    auto mftparticle = mfttrack.mcParticle();
                    if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                      registry.fill(HIST("hDCATrueMatch"), DCA1);
                      registry.fill(HIST("hMFTDCATrueMatch"), MFTDCA1);
                    } else {
                      registry.fill(HIST("hDCAFakeMatch"), DCA1);
                      registry.fill(HIST("hMFTDCAFakeMatch"), MFTDCA1);
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
    adaptAnalysisTask<checkpairdca>(cfgc)
  };
}
