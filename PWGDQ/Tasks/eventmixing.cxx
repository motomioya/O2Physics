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
#include "Framework/ASoAHelpers.h"
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
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
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
using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;

struct eventmixing {

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
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = -5.0;
  //o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  Configurable<int> ndepth{"ndepth", 5, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -7.5f, -5.0f, -2.5f, 0.0f, 2.5f, 5.0f, 7.5f, 10.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfVtxBinsOne{"ConfVtxBinsOne", {VARIABLE_WIDTH, -10.0f, 10.0f}, "Mixing bins - z-vertex"};

  BinningType colBinning{{ConfVtxBins}, true};
  BinningType colBinningOne{{ConfVtxBinsOne}, true};

  HistogramRegistry registry{
    "registry", 
    {
      {"hPMMassPtEM", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hPMMassPtEMVtx", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hPMSV", "SV;SV(z);SV(z)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hPMMassPtEMSVcut", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hPMMassPtEMSVcutVtx", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hPMSVVtx", "SV;SV(z);SV(z)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hPPMassPtEM", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hPPMassPtEMVtx", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hPPSV", "SV;SV(z);SV(z)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hPPMassPtEMSVcut", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hPPMassPtEMSVcutVtx", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hPPSVVtx", "SV;SV(z);SV(z)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hMMMassPtEM", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hMMMassPtEMVtx", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hMMSV", "SV;SV(z);SV(z)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hMMMassPtEMSVcut", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hMMMassPtEMSVcutVtx", "Invariant;Invariant Mass (GeV/c^{2});mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"hMMSVVtx", "SV;SV(z);SV(z)", {HistType::kTH1F, {{1100, -1000, 100}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
    }
  }

  void process(soa::Filtered<aod::Collisions> const& collisions, soa::Filtered<MyMuons> const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;

    bool propagateToPCA = true;
    float maxR = 200;
    float minParamChange = 1.0e-3;
    float minRelChi2Change = 0.9;
    bool useAbsDCA = false;

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    fgFitterTwoProngFwd.setTGeoMat(false);
    fgFitterTwoProngFwd.setMatLUT(lut);
    fgFitterTwoProngFwd.setBz(mMagField);
    fgFitterTwoProngFwd.setPropagateToPCA(propagateToPCA);
    fgFitterTwoProngFwd.setMaxR(maxR);
    fgFitterTwoProngFwd.setMinParamChange(minParamChange);
    fgFitterTwoProngFwd.setMinRelChi2Change(minRelChi2Change);
    fgFitterTwoProngFwd.setUseAbsDCA(useAbsDCA);


    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) {
      auto muonIdsThisCollision1 = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision1.globalIndex());
      auto muonIdsThisCollision2 = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision2.globalIndex());
      for (auto& [fwdtrackId1, fwdtrackId2] : combinations(o2::soa::CombinationsFullIndexPolicy(muonIdsThisCollision1, muonIdsThisCollision2))) {
        for (auto& fwdtrack1 : fwdtracks) {
          if (fwdtrack1.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId1.fwdtrackId() == fwdtrack1.globalIndex() && fwdtrack1.chi2MatchMCHMFT() < 50 && fwdtrack1.compatibleCollIds().size() == 1) {
            for (auto& fwdtrack2 : fwdtracks) {
              if (fwdtrack2.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId2.fwdtrackId() == fwdtrack2.globalIndex() && fwdtrack2.chi2MatchMCHMFT() < 50 && fwdtrack2.compatibleCollIds().size() == 1) {
                //calculate mass
                TLorentzVector lv1, lv2, lv;
                lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
                lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
                lv = lv1 + lv2;

                //set TrackParCovFwd for fwdtrack pairs
                double chi21 = fwdtrack1.chi2();
                SMatrix5 tpars1(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.phi(), fwdtrack1.tgl(), fwdtrack1.signed1Pt());
                std::vector<double> v11{fwdtrack1.cXX(), fwdtrack1.cXY(), fwdtrack1.cYY(), fwdtrack1.cPhiX(), fwdtrack1.cPhiY(), fwdtrack1.cPhiPhi(), fwdtrack1.cTglX(), fwdtrack1.cTglY(), fwdtrack1.cTglPhi(), fwdtrack1.cTglTgl(), fwdtrack1.c1PtX(), fwdtrack1.c1PtY(), fwdtrack1.c1PtPhi(), fwdtrack1.c1PtTgl(), fwdtrack1.c1Pt21Pt2()};
                SMatrix55 tcovs1(v11.begin(), v11.end());
                o2::track::TrackParCovFwd pars1{fwdtrack1.z(), tpars1, tcovs1, chi21};

                double chi22 = fwdtrack2.chi2();
                SMatrix5 tpars2(fwdtrack2.x(), fwdtrack2.y(), fwdtrack2.phi(), fwdtrack2.tgl(), fwdtrack2.signed1Pt());
                std::vector<double> v12{fwdtrack2.cXX(), fwdtrack2.cXY(), fwdtrack2.cYY(), fwdtrack2.cPhiX(), fwdtrack2.cPhiY(), fwdtrack2.cPhiPhi(), fwdtrack2.cTglX(), fwdtrack2.cTglY(), fwdtrack2.cTglPhi(), fwdtrack2.cTglTgl(), fwdtrack2.c1PtX(), fwdtrack2.c1PtY(), fwdtrack2.c1PtPhi(), fwdtrack2.c1PtTgl(), fwdtrack2.c1Pt21Pt2()};
                SMatrix55 tcovs2(v12.begin(), v12.end());
                o2::track::TrackParCovFwd pars2{fwdtrack2.z(), tpars2, tcovs2, chi22};

                //Get secondary vertex using DCAFitterN
                int procCode = fgFitterTwoProngFwd.process(pars1, pars2);
                double VertexingSV = -999;
                if (procCode == 0 ) {
                  Vec3D secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
                  fgFitterTwoProngFwd.calcPCACovMatrixFlat();
                  VertexingSV = secondaryVertex[2];
                }

                if (fwdtrack1.signed1Pt()*fwdtrack2.signed1Pt() < 0) {
                  registry.fill(HIST("hPMMassPtEMVtx"), lv.M(), lv.Pt());
                  registry.fill(HIST("hPMSVVtx"), VertexingSV);
                  if (VertexingSV >= -20 || VertexingSV == -999) {
                    registry.fill(HIST("hPMMassPtEMSVcutVtx"), lv.M(), lv.Pt());
                  }
                } else if (fwdtrack1.signed1Pt() >= 0 && fwdtrack2.signed1Pt() >= 0) {
                  registry.fill(HIST("hPPMassPtEMVtx"), lv.M(), lv.Pt());
                  registry.fill(HIST("hPPSVVtx"), VertexingSV);
                  if (VertexingSV >= -20 || VertexingSV == -999) {
                    registry.fill(HIST("hPPMassPtEMSVcutVtx"), lv.M(), lv.Pt());
                  }
                } else if (fwdtrack1.signed1Pt() < 0 && fwdtrack2.signed1Pt() < 0) {
                  registry.fill(HIST("hMMMassPtEMVtx"), lv.M(), lv.Pt());
                  registry.fill(HIST("hMMSVVtx"), VertexingSV);
                  if (VertexingSV >= -20 || VertexingSV == -999) {
                    registry.fill(HIST("hMMMassPtEMSVcutVtx"), lv.M(), lv.Pt());
                  }
                }
              }
            }
          }
        }
      }
    }
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinningOne, ndepth, -1, collisions, collisions)) {
      auto muonIdsThisCollision1 = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision1.globalIndex());
      auto muonIdsThisCollision2 = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision2.globalIndex());
      for (auto& [fwdtrackId1, fwdtrackId2] : combinations(o2::soa::CombinationsFullIndexPolicy(muonIdsThisCollision1, muonIdsThisCollision2))) {
        for (auto& fwdtrack1 : fwdtracks) {
          if (fwdtrack1.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId1.fwdtrackId() == fwdtrack1.globalIndex() && fwdtrack1.chi2MatchMCHMFT() < 50 && fwdtrack1.compatibleCollIds().size() == 1) {
            for (auto& fwdtrack2 : fwdtracks) {
              if (fwdtrack2.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId2.fwdtrackId() == fwdtrack2.globalIndex() && fwdtrack2.chi2MatchMCHMFT() < 50 && fwdtrack2.compatibleCollIds().size() == 1) {
                //calculate mass
                TLorentzVector lv1, lv2, lv;
                lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
                lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
                lv = lv1 + lv2;

                //set TrackParCovFwd for fwdtrack pairs
                double chi21 = fwdtrack1.chi2();
                SMatrix5 tpars1(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.phi(), fwdtrack1.tgl(), fwdtrack1.signed1Pt());
                std::vector<double> v11{fwdtrack1.cXX(), fwdtrack1.cXY(), fwdtrack1.cYY(), fwdtrack1.cPhiX(), fwdtrack1.cPhiY(), fwdtrack1.cPhiPhi(), fwdtrack1.cTglX(), fwdtrack1.cTglY(), fwdtrack1.cTglPhi(), fwdtrack1.cTglTgl(), fwdtrack1.c1PtX(), fwdtrack1.c1PtY(), fwdtrack1.c1PtPhi(), fwdtrack1.c1PtTgl(), fwdtrack1.c1Pt21Pt2()};
                SMatrix55 tcovs1(v11.begin(), v11.end());
                o2::track::TrackParCovFwd pars1{fwdtrack1.z(), tpars1, tcovs1, chi21};

                double chi22 = fwdtrack2.chi2();
                SMatrix5 tpars2(fwdtrack2.x(), fwdtrack2.y(), fwdtrack2.phi(), fwdtrack2.tgl(), fwdtrack2.signed1Pt());
                std::vector<double> v12{fwdtrack2.cXX(), fwdtrack2.cXY(), fwdtrack2.cYY(), fwdtrack2.cPhiX(), fwdtrack2.cPhiY(), fwdtrack2.cPhiPhi(), fwdtrack2.cTglX(), fwdtrack2.cTglY(), fwdtrack2.cTglPhi(), fwdtrack2.cTglTgl(), fwdtrack2.c1PtX(), fwdtrack2.c1PtY(), fwdtrack2.c1PtPhi(), fwdtrack2.c1PtTgl(), fwdtrack2.c1Pt21Pt2()};
                SMatrix55 tcovs2(v12.begin(), v12.end());
                o2::track::TrackParCovFwd pars2{fwdtrack2.z(), tpars2, tcovs2, chi22};

                //Get secondary vertex using DCAFitterN
                int procCode = fgFitterTwoProngFwd.process(pars1, pars2);
                double VertexingSV = -999;
                if (procCode == 0 ) {
                  Vec3D secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
                  fgFitterTwoProngFwd.calcPCACovMatrixFlat();
                  VertexingSV = secondaryVertex[2];
                }

                if (fwdtrack1.signed1Pt()*fwdtrack2.signed1Pt() < 0) {
                  registry.fill(HIST("hPMMassPtEM"), lv.M(), lv.Pt());
                  registry.fill(HIST("hPMSV"), VertexingSV);
                  if (VertexingSV >= -20 || VertexingSV == -999) {
                    registry.fill(HIST("hPMMassPtEMSVcut"), lv.M(), lv.Pt());
                  }
                } else if (fwdtrack1.signed1Pt() >= 0 && fwdtrack2.signed1Pt() >= 0) {
                  registry.fill(HIST("hPPMassPtEM"), lv.M(), lv.Pt());
                  registry.fill(HIST("hPPSV"), VertexingSV);
                  if (VertexingSV >= -20 || VertexingSV == -999) {
                    registry.fill(HIST("hPPMassPtEMSVcut"), lv.M(), lv.Pt());
                  }
                } else if (fwdtrack1.signed1Pt() < 0 && fwdtrack2.signed1Pt() < 0) {
                  registry.fill(HIST("hMMMassPtEM"), lv.M(), lv.Pt());
                  registry.fill(HIST("hMMSV"), VertexingSV);
                  if (VertexingSV >= -20 || VertexingSV == -999) {
                    registry.fill(HIST("hMMMassPtEMSVcut"), lv.M(), lv.Pt());
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
    adaptAnalysisTask<eventmixing>(cfgc)
  };
}
