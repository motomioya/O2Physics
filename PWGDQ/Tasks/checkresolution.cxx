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

struct checkresolution {

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
      {"RelPtResArrCocktail", "p_{T} resolution;p_{T}^{gen} (Gev/c);p_{T}^{gen} - p_{T}^{rec} / p_{T}^{gen} (Gev/c)", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, -4, 4}}}},
      {"EtaResArr", "#eta resolution;p_{T}^{gen} (Gev/c);#eta^{gen} - #eta^{rec}", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, -2, 2}}}},
      {"EtaPosResArr", "#eta resolution for positive track;p_{T}^{gen} (Gev/c);#eta^{gen} - #eta^{rec}", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, -2, 2}}}},
      {"EtaNegResArr", "#eta resolution for negative track;p_{T}^{gen} (Gev/c);#eta^{gen} - #eta^{rec}", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, -2, 2}}}},
      {"PhiResArr", "#phi resolution;p_{T}^{gen} (Gev/c);#phi^{gen} - #phi^{rec}", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, -2, 2}}}},
      {"PhiPosResArr", "#phi resolution for positive track;p_{T}^{gen} (Gev/c);#phi^{gen} - #phi^{rec}", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, -2, 2}}}},
      {"PhiNegResArr", "#phi resolution for negative track;p_{T}^{gen} (Gev/c);#phi^{gen} - #phi^{rec}", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, -2, 2}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, soa::Filtered<MyMuons> const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const& particles)
  {
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& fwdtrackId : muonIdsThisCollision) {
        for (auto& fwdtrack : fwdtracks) {
          if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId.fwdtrackId() == fwdtrack.globalIndex() && fwdtrack.compatibleCollIds().size() == 1) {
            if (fwdtrack.has_mcParticle() == 1) {
              auto fwdparticle = fwdtrack.mcParticle();
              double resphi = fwdparticle.phi() - fwdtrack.phi();
              if (resphi > M_PI) resphi = resphi - 2 * M_PI;
              registry.fill(HIST("RelPtResArrCocktail"), fwdparticle.pt(), (fwdparticle.pt() - fwdtrack.pt())/fwdparticle.pt());
              registry.fill(HIST("EtaResArr"), fwdparticle.pt(), fwdparticle.eta() - fwdtrack.eta());
              registry.fill(HIST("PhiResArr"), fwdparticle.pt(), resphi);
              if (fwdtrack.sign() == 1) {
                registry.fill(HIST("EtaPosResArr"), fwdparticle.pt(), fwdparticle.eta() - fwdtrack.eta());
                registry.fill(HIST("PhiPosResArr"), fwdparticle.pt(), resphi);
              } else if (fwdtrack.sign() == -1) {
                registry.fill(HIST("EtaNegResArr"), fwdparticle.pt(), fwdparticle.eta() - fwdtrack.eta());
                registry.fill(HIST("PhiNegResArr"), fwdparticle.pt(), resphi);
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
    adaptAnalysisTask<checkresolution>(cfgc)
  };
}
