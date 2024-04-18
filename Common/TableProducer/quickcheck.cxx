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
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoNode.h>
#include <TGeoShape.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using o2::globaltracking::MatchingFunc_t;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
//using MyMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA, aod::FwdTrkCompColls, aod::McFwdTrackLabels>;
using MyMFTs = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::McMFTTrackLabels>;
using MyCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;

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

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  HistogramRegistry registry{
    "registry",
    {
      {"deltaCollisionId", "deltaCollisionId", {HistType::kTH1F, {{401, -200.5, 200.5}}}},
      {"MFTdeltaCollisionId", "MFTdeltaCollisionId", {HistType::kTH1F, {{401, -200.5, 200.5}}}},
      {"AssociateddeltaCollisionId", "AssociateddeltaCollisionId", {HistType::kTH1F, {{401, -200.5, 200.5}}}},
      {"AssociatedMFTdeltaCollisionId", "AssociatedMFTdeltaCollisionId", {HistType::kTH1F, {{401, -200.5, 200.5}}}},
      {"AllTrueDeltaCollId", "AllTrueDeltaCollId", {HistType::kTH1F, {{401, -200.5, 200.5}}}},
      {"AfterAssociationTrueDeltaCollId", "AfterAssociationTrueDeltaCollId", {HistType::kTH1F, {{401, -200.5, 200.5}}}},
      {"TrueDeltaXY", "TrueDeltaXY", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"TrueDeltaPhiTanl", "TrueDeltaPhiTanl", {HistType::kTH1F, {{1000, -25, 25}}}},
      {"AssocTrueDeltaXY", "TrueDeltaXY", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"AssocTrueDeltaPhiTanl", "TrueDeltaPhiTanl", {HistType::kTH1F, {{1000, -25, 25}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
    fCCDB->setURL("http://alice-ccdb.cern.ch");
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDBApi.init("http://alice-ccdb.cern.ch");
  }

  void process(o2::aod::BCs const& bcs, MyCollisions const& collisions, soa::Filtered<MyMuons> const& fwdtracks, MyMFTs const& mfttracks, aod::McParticles const&, aod::FwdTrackAssoc const& fwdtrackIndices, aod::MFTTrackAssoc const& mfttrackIndices, aod::McCollisions const&)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;
    float Bz = 0;                                         // Magnetic field for MFT
    static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
    static o2::globaltracking::MatchGlobalFwd mMatching;
    int run = bcs.begin().runNumber();

    std::map<string, string> metadata, headers;
    headers = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", run), metadata, -1);
    int64_t ts = std::atol(headers["SOR"].c_str());
    auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded())
      fCCDB->get<TGeoManager>("GLO/Config/GeometryAligned");
    o2::mch::TrackExtrap::setField();
    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    Bz = field->getBz(centerMFT);
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz;
    /*
    for (auto& fwdtrack : fwdtracks) {
      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        if (fwdtrack.has_mcParticle() & fwdtrack.has_collision()){
          auto fwdCollision = fwdtrack.collision_as<MyCollisions>();
          if (fwdCollision.has_mcCollision()){
            auto fwdCollisionId = fwdCollision.mcCollision_as<aod::McCollisions>().globalIndex();
            auto parCollisionId = fwdtrack.mcParticle().mcCollisionId();
            auto deltaCollisionId = fwdCollisionId - parCollisionId;
            registry.fill(HIST("deltaCollisionId"), deltaCollisionId);
          }
        }
      }
    }
    for (auto& mfttrack : mfttracks) {
      if (mfttrack.has_mcParticle() & mfttrack.has_collision()){
        auto mftCollision = mfttrack.collision_as<MyCollisions>();
        if (mftCollision.has_mcCollision()){
          auto mftCollisionId = mftCollision.mcCollision_as<aod::McCollisions>().globalIndex();
          auto parCollisionId = mfttrack.mcParticle().mcCollisionId();
          auto deltaCollisionId = mftCollisionId - parCollisionId;
          registry.fill(HIST("MFTdeltaCollisionId"), deltaCollisionId);
        }
      }
    }

    for (auto& fwdtrackId : fwdtrackIndices) {
      for (auto& fwdtrack : fwdtracks) {
        if (fwdtrackId.fwdtrackId() == fwdtrack.globalIndex()) {
          if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
            int num = 0;
            num = fwdtrack.compatibleCollIds().size();
            if (num == 1) {
              if (fwdtrack.has_mcParticle() && fwdtrack.has_collision()){
                for (auto& col : collisions) {
                  if (col.globalIndex() == fwdtrackId.collisionId()){
                    if (col.has_mcCollision()){
                      auto fwdCollisionId = col.mcCollision_as<aod::McCollisions>().globalIndex();
                      auto parCollisionId = fwdtrack.mcParticle().mcCollisionId();
                      auto deltaCollisionId = fwdCollisionId - parCollisionId;
                      registry.fill(HIST("AssociateddeltaCollisionId"), deltaCollisionId);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    for (auto& mfttrackId : mfttrackIndices) {
      for (auto& mfttrack : mfttracks) {
        if (mfttrackId.mfttrackId() == mfttrack.globalIndex()) {
          int num = 0;
          num = mfttrack.compatibleCollIds().size();
          if (num == 1) {
            if (mfttrack.has_mcParticle() && mfttrack.has_collision()){
              for (auto& col : collisions) {
                if (col.globalIndex() == mfttrackId.collisionId()){
                  if (col.has_mcCollision()){
                    auto mftCollisionId = col.mcCollision_as<aod::McCollisions>().globalIndex();
                    auto parCollisionId = mfttrack.mcParticle().mcCollisionId();
                    auto deltaCollisionId = mftCollisionId - parCollisionId;
                    registry.fill(HIST("AssociatedMFTdeltaCollisionId"), deltaCollisionId);
                  }
                }
              }
            }
          }
        }
      }
    }
    */

    for (auto& [fwdtrack, mfttrack] : combinations(o2::soa::CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {
      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        if (fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
          auto fwdparticle = fwdtrack.mcParticle();
          auto mftparticle = mfttrack.mcParticle();

          //Propagate MCH-MID track
          double muonchi2 = fwdtrack.chi2();
          SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
          std::vector<double> fwdtrackv1{fwdtrack.cXX(), fwdtrack.cXY(), fwdtrack.cYY(), fwdtrack.cPhiX(), fwdtrack.cPhiY(), fwdtrack.cPhiPhi(), fwdtrack.cTglX(), fwdtrack.cTglY(), fwdtrack.cTglPhi(), fwdtrack.cTglTgl(), fwdtrack.c1PtX(), fwdtrack.c1PtY(), fwdtrack.c1PtPhi(), fwdtrack.c1PtTgl(), fwdtrack.c1Pt21Pt2()};
          SMatrix55 muoncovs(fwdtrackv1.begin(), fwdtrackv1.end());
          o2::track::TrackParCovFwd muonpartrack{fwdtrack.z(), muonpars, muoncovs, muonchi2};
          muonpartrack.propagateToZ(MatchingPlaneZ, 0);
          /*
          o2::dataformats::GlobalFwdTrack propmuon;
          o2::dataformats::GlobalFwdTrack gtrack;
          gtrack.setParameters(muonpars);
          gtrack.setZ(muonpartrack.getZ());
          gtrack.setCovariances(muoncovs);
          auto mchTrack = mMatching.FwdtoMCH(gtrack);
          o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, MatchingPlaneZ);
          auto proptrack = mMatching.MCHtoFwd(mchTrack);
          propmuon.setParameters(proptrack.getParameters());
          propmuon.setZ(proptrack.getZ());
          propmuon.setCovariances(proptrack.getCovariances());
          */

          //Propagate MFT track
          double mftchi2 = mfttrack.chi2();
          SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
          std::vector<double> mftv1;
          SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
          o2::track::TrackParCovFwd mftpartrack{mfttrack.z(), mftpars, mftcovs, mftchi2};
          mftpartrack.propagateToZ(MatchingPlaneZ, Bz);

          //Get parameters
          Float_t MFT_X = mftpartrack.getX();
          Float_t MFT_Y = mftpartrack.getY();
          Float_t MFT_Phi = mftpartrack.getPhi();
          Float_t MFT_Tanl = mftpartrack.getTanl();

          /*
          Float_t MCH_X = propmuon.getX();
          Float_t MCH_Y = propmuon.getY();
          Float_t MCH_Phi = propmuon.getPhi();
          Float_t MCH_Tanl = propmuon.getTanl();
          */
          Float_t MCH_X = muonpartrack.getX();
          Float_t MCH_Y = muonpartrack.getY();
          Float_t MCH_Phi = muonpartrack.getPhi();
          Float_t MCH_Tanl = muonpartrack.getTanl();

          Float_t Delta_X = MFT_X - MCH_X;
          Float_t Delta_Y = MFT_Y - MCH_Y;
          Float_t Delta_Phi = MFT_Phi - MCH_Phi;
          Float_t Delta_Tanl = MFT_Tanl - MCH_Tanl;

          Float_t Delta_XY = sqrt(Delta_X * Delta_X + Delta_Y * Delta_Y);
          Float_t Delta_PhiTanl = sqrt(Delta_Phi * Delta_Phi + Delta_Tanl * Delta_Tanl);

          if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
            registry.fill(HIST("AllTrueDeltaCollId"), fwdtrack.collisionId() - mfttrack.collisionId());
            if (fwdtrack.collisionId() - mfttrack.collisionId() < 6) {
              registry.fill(HIST("TrueDeltaXY"), Delta_XY);
              registry.fill(HIST("TrueDeltaPhiTanl"), Delta_PhiTanl);
            }
          }
        }
      }
    }

    for (auto& [fwdtrackId, mfttrack] : combinations(o2::soa::CombinationsFullIndexPolicy(fwdtrackIndices, mfttracks))) {
      for (auto& fwdtrack : fwdtracks) {
        if (fwdtrackId.fwdtrackId() == fwdtrack.globalIndex()) {
          if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
            int num = 0;
            num = fwdtrack.compatibleCollIds().size();
            if (num == 1) {
              if (fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
                auto fwdparticle = fwdtrack.mcParticle();
                auto mftparticle = mfttrack.mcParticle();

                //Propagate MCH-MID track
                double muonchi2 = fwdtrack.chi2();
                SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
                std::vector<double> fwdtrackv1{fwdtrack.cXX(), fwdtrack.cXY(), fwdtrack.cYY(), fwdtrack.cPhiX(), fwdtrack.cPhiY(), fwdtrack.cPhiPhi(), fwdtrack.cTglX(), fwdtrack.cTglY(), fwdtrack.cTglPhi(), fwdtrack.cTglTgl(), fwdtrack.c1PtX(), fwdtrack.c1PtY(), fwdtrack.c1PtPhi(), fwdtrack.c1PtTgl(), fwdtrack.c1Pt21Pt2()};
                SMatrix55 muoncovs(fwdtrackv1.begin(), fwdtrackv1.end());
                o2::track::TrackParCovFwd muonpartrack{fwdtrack.z(), muonpars, muoncovs, muonchi2};
                o2::dataformats::GlobalFwdTrack propmuon;
                o2::dataformats::GlobalFwdTrack gtrack;
                gtrack.setParameters(muonpars);
                gtrack.setZ(muonpartrack.getZ());
                gtrack.setCovariances(muoncovs);
                auto mchTrack = mMatching.FwdtoMCH(gtrack);
                o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, MatchingPlaneZ);
                auto proptrack = mMatching.MCHtoFwd(mchTrack);
                propmuon.setParameters(proptrack.getParameters());
                propmuon.setZ(proptrack.getZ());
                propmuon.setCovariances(proptrack.getCovariances());

                //Propagate MFT track
                double mftchi2 = mfttrack.chi2();
                SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                std::vector<double> mftv1;
                SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                o2::track::TrackParCovFwd mftpartrack{mfttrack.z(), mftpars, mftcovs, mftchi2};
                mftpartrack.propagateToZ(MatchingPlaneZ, Bz);

                //Get parameters
                Float_t MFT_X = mftpartrack.getX();
                Float_t MFT_Y = mftpartrack.getY();
                Float_t MFT_Phi = mftpartrack.getPhi();
                Float_t MFT_Tanl = mftpartrack.getTanl();

                Float_t MCH_X = propmuon.getX();
                Float_t MCH_Y = propmuon.getY();
                Float_t MCH_Phi = propmuon.getPhi();
                Float_t MCH_Tanl = propmuon.getTanl();

                Float_t Delta_X = MFT_X - MCH_X;
                Float_t Delta_Y = MFT_Y - MCH_Y;
                Float_t Delta_Phi = MFT_Phi - MCH_Phi;
                Float_t Delta_Tanl = MFT_Tanl - MCH_Tanl;

                Float_t Delta_XY = sqrt(Delta_X * Delta_X + Delta_Y * Delta_Y);
                Float_t Delta_PhiTanl = sqrt(Delta_Phi * Delta_Phi + Delta_Tanl * Delta_Tanl);
                if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                  registry.fill(HIST("AfterAssociationTrueDeltaCollId"), fwdtrackId.collisionId() - mfttrack.collisionId());
                  if (fwdtrack.collisionId() - mfttrack.collisionId() < 6) {
                    registry.fill(HIST("AssocTrueDeltaXY"), Delta_XY);
                    registry.fill(HIST("AssocTrueDeltaPhiTanl"), Delta_PhiTanl);
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
    adaptAnalysisTask<quickcheck>(cfgc)
  };
}
