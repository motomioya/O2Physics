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
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
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
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTrkCompColls, aod::McFwdTrackLabels>;

struct checkmftdca{

  float collisionZcut = 10.0f;
  Filter collisionFilter = nabs(aod::collision::posZ) < collisionZcut;
  Preslice<aod::McParticles> particlesIndicesPerCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  float mMagField = -5.0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{
    "registry", 
    {
      {"hDCAMFTlinear", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAMFTlinearPrimary", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAMFTlinearFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAMFTlinearPion", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAMFTlinearKaon", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hRxyz", "Rxyz;DCA (cm)", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"hDCAGloballinear", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGloballinearFromPion", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGloballinearFromKaon", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGloballinearFromDtoPion", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGloballinearFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGloballinearFromB", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGloballinearFromnJpsi", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGloballinearFrompJpsi", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGloballinearFromLF", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelix", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelixFromPion", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelixFromKaon", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelixFromDtoPion", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelixFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelixFromB", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelixFromnJpsi", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelixFrompJpsi", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
      {"hDCAGlobalhelixFromLF", "DCA;DCA (cm)", {HistType::kTH1F, {{2000, 0, 1}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(soa::Filtered<aod::Collisions> const& collisions, soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McCollisions const&, aod::McParticles const& particles, MyMuons const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    if (!o2::base::GeometryManager::isGeometryLoaded())
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
    //Everything From D
    MCSignal* signalEfromD;
    MCProng prongEfromD(2, {0, 403}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signalEfromD = new MCSignal("signalEfromD", "Everything from charm", {prongEfromD}, {-1});
    //Pion not from D
    MCSignal* signalpion;
    MCProng pionprong(2, {211, 403}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false}); 
    signalpion = new MCSignal("signalpion", "Primary Muons", {pionprong}, {-1});
    //Kaon not from D
    MCSignal* signalkaon;
    MCProng kaonprong(2, {321, 403}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false}); 
    signalkaon = new MCSignal("signalkaon", "Primary Muons", {kaonprong}, {-1});
    //muon from pion
    MCSignal* signalmufrompion;
    MCProng mufrompionprong(2, {13, 211}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufrompion = new MCSignal("signalmufrompion", "Primary Muons", {mufrompionprong}, {-1});
    //muon from kaon
    MCSignal* signalmufromkaon;
    MCProng mufromkaonprong(2, {13, 321}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromkaon = new MCSignal("signalmufromkaon", "Primary Muons", {mufromkaonprong}, {-1});
    //muon from DtoPion
    MCSignal* signalmufromdtopion;
    MCProng mufromdtopionprong(3, {13, 211, 403}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}); 
    signalmufromdtopion = new MCSignal("signaldtopion", "Primary Muons", {mufromdtopionprong}, {-1});
    //muon from D
    MCSignal* signalmufromD;
    MCProng mufromDprong(2, {13, 403}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromD = new MCSignal("signalmufromD", "Primary Muons", {mufromDprong}, {-1});
    //muon from B
    MCSignal* signalmufromB;
    MCProng mufromBprong(2, {13, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromB = new MCSignal("signalmufromB", "Primary Muons", {mufromBprong}, {-1});
    //muon from omega, eta, phi
    MCSignal* signalmufromphi;
    MCProng mufromphiprong(2, {13, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromphi = new MCSignal("signalmufromphi", "Primary Muons", {mufromphiprong}, {-1});
    MCSignal* signalmufromomega;
    MCProng mufromomegaprong(2, {13, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromomega = new MCSignal("signalmufromomega", "Primary Muons", {mufromomegaprong}, {-1});
    MCSignal* signalmufrometa;
    MCProng mufrometaprong(2, {13, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufrometa = new MCSignal("signalmufrometa", "Primary Muons", {mufrometaprong}, {-1});
    MCSignal* signalmufromrho;
    MCProng mufromrhoprong(2, {13, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromrho = new MCSignal("signalmufromrho", "Primary Muons", {mufromrhoprong}, {-1});
    //muon from J/psi
    MCSignal* signalmufromjpsi;
    MCProng mufromjpsiprong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromjpsi = new MCSignal("signalmufromjpsi", "Primary Muons", {mufromjpsiprong}, {-1});
    //muon from nonprompt J/psi
    MCSignal* signalmufromnjpsi;
    MCProng mufromnjpsiprong(3, {13, 443, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}); 
    signalmufromnjpsi = new MCSignal("signalmufromnjpsi", "Primary Muons", {mufromnjpsiprong}, {-1});

    for (auto& mfttrack : mfttracks) {
      if (mfttrack.eta() > -3.6 && mfttrack.eta() < -2.5) {
        if (mfttrack.has_mcParticle()) {
          auto mftparticle = mfttrack.mcParticle();
          if (mftparticle.has_mcCollision()) {
            auto mccollision = mftparticle.mcCollision();

            float mftdcaX = -999;
            float mftdcaY = -999;

            //calculate DCA
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            std::vector<double> mftv1;
            SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
            o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
            //mftpars1.propagateToZlinear(MatchingPlaneZ);
            mftpars1.propagateToZlinear(mccollision.posZ());
            mftdcaX = (mftpars1.getX() - mccollision.posX());
            mftdcaY = (mftpars1.getY() - mccollision.posY());
            float MFTDCA = std::sqrt(mftdcaX * mftdcaX + mftdcaY * mftdcaY);
            registry.fill(HIST("hDCAMFTlinear"), MFTDCA);
            float rxyz = std::sqrt((mccollision.posX() - mftparticle.vx()) * (mccollision.posX() - mftparticle.vx()) + (mccollision.posX() - mftparticle.vx()) * (mccollision.posX() - mftparticle.vx()) + (mccollision.posX() - mftparticle.vx()) * (mccollision.posX() - mftparticle.vx()));


            registry.fill(HIST("hRxyz"), rxyz);
            if (rxyz < 0.01) {
              registry.fill(HIST("hDCAMFTlinearPrimary"), MFTDCA);
            }
            if (signalEfromD->CheckSignal(true, mftparticle)) {
              registry.fill(HIST("hDCAMFTlinearFromD"), MFTDCA);
            } else if (signalpion->CheckSignal(true, mftparticle)) {
              registry.fill(HIST("hDCAMFTlinearPion"), MFTDCA);
            } else if (signalkaon->CheckSignal(true, mftparticle)) {
              registry.fill(HIST("hDCAMFTlinearKaon"), MFTDCA);
            }
          }
        }
      }
    }
    for (auto& collision : collisions) {
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& fwdtrackId : fwdtrackIndices) {
        for (auto& fwdtrack : fwdtracks) {
          if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId.fwdtrackId() == fwdtrack.globalIndex() && fwdtrack.chi2MatchMCHMFT() < 50 && fwdtrack.compatibleCollIds().size() == 1) {
            if (fwdtrack.eta() > -4 && fwdtrack.eta() < -2.5 && (((17.6 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 26.5) && (fwdtrack.pDca() < 594)) || ((26.5 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 89.5) && (fwdtrack.pDca() < 324)))){
              //
              //set TrackParCovFwd for fwdtrack pairs
              double chi21 = fwdtrack.chi2();
              SMatrix5 tpars1(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
              std::vector<double> v11{fwdtrack.cXX(), fwdtrack.cXY(), fwdtrack.cYY(), fwdtrack.cPhiX(), fwdtrack.cPhiY(), fwdtrack.cPhiPhi(), fwdtrack.cTglX(), fwdtrack.cTglY(), fwdtrack.cTglPhi(), fwdtrack.cTglTgl(), fwdtrack.c1PtX(), fwdtrack.c1PtY(), fwdtrack.c1PtPhi(), fwdtrack.c1PtTgl(), fwdtrack.c1Pt21Pt2()};
              SMatrix55 tcovs1(v11.begin(), v11.end());

              o2::track::TrackParCovFwd pars1{fwdtrack.z(), tpars1, tcovs1, chi21};

              //fwdtrack linear propagation
              pars1.propagateToZlinear(collision.posZ());
              //calculate DCA
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
              if (fwdtrack.has_mcParticle()) {
                auto fwdparticle = fwdtrack.mcParticle();
                //Fill histograms
                registry.fill(HIST("hDCAGloballinear"), DCA1);
                registry.fill(HIST("hDCAGlobalhelix"), DCA2);
                if (signalmufrompion->CheckSignal(true, fwdparticle)) {
                  registry.fill(HIST("hDCAGloballinearFromPion"), DCA1);
                  registry.fill(HIST("hDCAGlobalhelixFromPion"), DCA2);
                }
                if (signalmufromkaon->CheckSignal(true, fwdparticle)) {
                  registry.fill(HIST("hDCAGloballinearFromKaon"), DCA1);
                  registry.fill(HIST("hDCAGlobalhelixFromKaon"), DCA2);
                }
                if (signalmufromdtopion->CheckSignal(true, fwdparticle)) {
                  registry.fill(HIST("hDCAGloballinearFromDtoPion"), DCA1);
                  registry.fill(HIST("hDCAGlobalhelixFromDtoPion"), DCA2);
                }
                if (signalmufromD->CheckSignal(true, fwdparticle)) {
                  registry.fill(HIST("hDCAGloballinearFromD"), DCA1);
                  registry.fill(HIST("hDCAGlobalhelixFromD"), DCA2);
                }
                if (signalmufromB->CheckSignal(true, fwdparticle)) {
                  registry.fill(HIST("hDCAGloballinearFromB"), DCA1);
                  registry.fill(HIST("hDCAGlobalhelixFromB"), DCA2);
                }
                if (signalmufromphi->CheckSignal(true, fwdparticle) || signalmufromomega->CheckSignal(true, fwdparticle) || signalmufrometa->CheckSignal(true, fwdparticle) || signalmufromrho->CheckSignal(true, fwdparticle)) {
                  registry.fill(HIST("hDCAGloballinearFromLF"), DCA1);
                  registry.fill(HIST("hDCAGlobalhelixFromLF"), DCA2);
                }
                if (signalmufromnjpsi->CheckSignal(true, fwdparticle)) {
                  registry.fill(HIST("hDCAGloballinearFromnJpsi"), DCA1);
                  registry.fill(HIST("hDCAGlobalhelixFromnJpsi"), DCA2);
                } else if (signalmufromjpsi->CheckSignal(true, fwdparticle)) {
                  registry.fill(HIST("hDCAGloballinearFrompJpsi"), DCA1);
                  registry.fill(HIST("hDCAGlobalhelixFrompJpsi"), DCA2);
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
    adaptAnalysisTask<checkmftdca>(cfgc)
  };
}
