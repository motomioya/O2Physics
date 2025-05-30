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
// p

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
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
#include <THashList.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"

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
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;
namespace o2::aod
{
namespace dqanalysisflags
{
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);                       //! Muon track decisions (joinable to ReducedMuonsAssoc)
DECLARE_SOA_COLUMN(MuonAmbiguityInBunch, muonAmbiguityInBunch, int8_t);              //! Muon track in-bunch ambiguity
DECLARE_SOA_COLUMN(MuonAmbiguityOutOfBunch, muonAmbiguityOutOfBunch, int8_t);        //! Muon track out of bunch ambiguity
                                                                                                                                  }
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);                                                       //!  joinable to ReducedMuonsAssoc
DECLARE_SOA_TABLE(MuonAmbiguities, "AOD", "DQMUONAMB", dqanalysisflags::MuonAmbiguityInBunch, dqanalysisflags::MuonAmbiguityOutOfBunch);         //!  joinable to ReducedMuonTracks
                                                                                                                                  }

//using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::ReducedMuonsLabels>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
//using MyMuonTracksWithCovWithAmbiguities = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonAmbiguities, aod::ReducedMuonsLabels>;
using MyMuonTracksWithCovWithAmbiguities = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonAmbiguities>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;

struct myQAsingleData {

  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  float mMagField = -5.0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{
    "registry", 
    {
      {"hDCANoRefitlinear", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinear", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelix", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCANoRefitlinearLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCANoRefitlinearMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCANoRefitlinearHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}}
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(MyEventsVtxCov const&, MyMuonTracksWithCovWithAmbiguities const&, aod::ReducedMuonsAssoc const& assocs)
  {

    if (!o2::base::GeometryManager::isGeometryLoaded())
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");

    for (auto& fwdtrackAssocId : assocs) {
      auto fwdtrack = fwdtrackAssocId.template reducedmuon_as<MyMuonTracksWithCovWithAmbiguities>();
      //Check Ambiguity
      if (fwdtrack.muonAmbiguityInBunch() > 1) continue;
      if (fwdtrack.muonAmbiguityOutOfBunch() > 1) continue;
      auto collision = fwdtrackAssocId.template reducedevent_as<MyEventsVtxCov>();

      //calculate DCA
      float fwd0dcaX = fwdtrack.fwdDcaX();
      float fwd0dcaY = fwdtrack.fwdDcaY();
      float DCA0 = std::sqrt(fwd0dcaX * fwd0dcaX + fwd0dcaY * fwd0dcaY);
      //set TrackParCovFwd for fwdtrack pairs
      double chi21 = fwdtrack.chi2();
      SMatrix5 tpars1(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
      std::vector<double> v11{fwdtrack.cXX(), fwdtrack.cXY(), fwdtrack.cYY(), fwdtrack.cPhiX(), fwdtrack.cPhiY(), fwdtrack.cPhiPhi(), fwdtrack.cTglX(), fwdtrack.cTglY(), fwdtrack.cTglPhi(), fwdtrack.cTglTgl(), fwdtrack.c1PtX(), fwdtrack.c1PtY(), fwdtrack.c1PtPhi(), fwdtrack.c1PtTgl(), fwdtrack.c1Pt21Pt2()};
      SMatrix55 tcovs1(v11.begin(), v11.end());
      o2::track::TrackParCovFwd pars1{fwdtrack.z(), tpars1, tcovs1, chi21};
      pars1.propagateToZlinear(collision.posZ());
      float fwd1dcaX = (pars1.getX() - collision.posX());
      float fwd1dcaY = (pars1.getY() - collision.posY());
      float DCA1 = std::sqrt(fwd1dcaX * fwd1dcaX + fwd1dcaY * fwd1dcaY);

      //fwdtrack lenear propagation
      o2::track::TrackParCovFwd pars2{fwdtrack.z(), tpars1, tcovs1, chi21};
      o2::dataformats::GlobalFwdTrack propmuon;
      auto geoMan = o2::base::GeometryManager::meanMaterialBudget(fwdtrack.x(), fwdtrack.y(), fwdtrack.z(), collision.posX(), collision.posY(), collision.posZ());
      auto x2x0 = static_cast<float>(geoMan.meanX2X0);
      pars2.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, mMagField, x2x0);
      propmuon.setParameters(pars2.getParameters());
      propmuon.setZ(pars2.getZ());
      propmuon.setCovariances(pars2.getCovariances());
      //calculate DCA
      float fwd2dcaX = (pars2.getX() - collision.posX());
      float fwd2dcaY = (pars2.getY() - collision.posY());
      float DCA2 = std::sqrt(fwd2dcaX * fwd2dcaX + fwd2dcaY * fwd2dcaY);

      registry.fill(HIST("hDCANoRefitlinear"), DCA0);
      registry.fill(HIST("hDCARefitlinear"), DCA1);
      registry.fill(HIST("hDCARefithelix"), DCA2);
      if (fwdtrack.pt() < 0.5) {
        registry.fill(HIST("hDCANoRefitlinearLowPt"), DCA0);
        registry.fill(HIST("hDCARefitlinearLowPt"), DCA1);
        registry.fill(HIST("hDCARefithelixLowPt"), DCA2);
      } else if (fwdtrack.pt() >= 0.5 && fwdtrack.pt() < 0.8) {
        registry.fill(HIST("hDCANoRefitlinearMidPt"), DCA0);
        registry.fill(HIST("hDCARefitlinearMidPt"), DCA1);
        registry.fill(HIST("hDCARefithelixMidPt"), DCA2);
      } else if (fwdtrack.pt() >= 0.8) {
        registry.fill(HIST("hDCANoRefitlinearHighPt"), DCA0);
        registry.fill(HIST("hDCARefitlinearHighPt"), DCA1);
        registry.fill(HIST("hDCARefithelixHighPt"), DCA2);
      }
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myQAsingleData>(cfgc)
  };
}
