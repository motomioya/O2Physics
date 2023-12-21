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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using o2::globaltracking::MatchingFunc_t;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using MyMuonsColl = soa::Join<aod::FwdTracks, aod::FwdTrkCompColls>;
using MyMFTsColl = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;

struct quickcheck {
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

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  HistogramRegistry registry{
    "registry",
    {
      {"TrueDeltaCollId", "TrueDeltaCollId", {HistType::kTH1F, {{4001, -2000.5, 2000.5}}}},
      {"FalseDeltaCollId", "FalseDeltaCollId", {HistType::kTH1F, {{4001, -2000.5, 2000.5}}}},
      {"TrueDeltaXY", "TrueDeltaXY", {HistType::kTH1F, {{6000, 0, 60}}}},
      {"FalseDeltaXY", "FalseDeltaXY", {HistType::kTH1F, {{6000, 0, 60}}}},
      {"TrueDeltaPhiTanl", "TrueDeltaPhiTanl", {HistType::kTH1F, {{4000, 0, 40}}}},
      {"FalseDeltaPhiTanl", "FalseDeltaPhiTanl", {HistType::kTH1F, {{4000, 0, 40}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions const& collisions, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&)
  {
    for (auto& [fwdtrack, mfttrack] : combinations(o2::soa::CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        if (fwdtrack.has_collision() && mfttrack.has_collision() && fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
          auto fwdparticle = fwdtrack.mcParticle();
          auto mftparticle = mfttrack.mcParticle();
          static constexpr Double_t MatchingPlaneZ = -77.5;

          // propagate muontrack to matching position
          double muonchi2 = fwdtrack.chi2();
          SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
          std::vector<double> muonv1;
          SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
          o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
          muonpars1.propagateToZlinear(MatchingPlaneZ);

          // propagate mfttrack to matching position
          double mftchi2 = mfttrack.chi2();
          SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
          std::vector<double> mftv1;
          SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
          o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
          mftpars1.propagateToZlinear(MatchingPlaneZ);

          Float_t MFT_X = mftpars1.getX();
          Float_t MFT_Y = mftpars1.getY();
          Float_t MFT_Phi = mftpars1.getPhi();
          Float_t MFT_Tanl = mftpars1.getTanl();

          Float_t MCH_X = muonpars1.getX();
          Float_t MCH_Y = muonpars1.getY();
          Float_t MCH_Phi = muonpars1.getPhi();
          Float_t MCH_Tanl = muonpars1.getTanl();


          Float_t Delta_X = MFT_X - MCH_X;
          Float_t Delta_Y = MFT_Y - MCH_Y;
          Float_t Delta_Phi = MFT_Phi - MCH_Phi;
          Float_t Delta_Tanl = MFT_Tanl - MCH_Tanl;

          Float_t Delta_XY = sqrt(Delta_X * Delta_X + Delta_Y * Delta_Y);
          Float_t Delta_PhiTanl = sqrt(Delta_Phi * Delta_Phi + Delta_Tanl * Delta_Tanl);
          if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
            registry.fill(HIST("TrueDeltaCollId"), fwdtrack.collisionId() - mfttrack.collisionId());
            registry.fill(HIST("TrueDeltaXY"), Delta_XY);
            registry.fill(HIST("TrueDeltaPhiTanl"), Delta_PhiTanl);
          } else {
            registry.fill(HIST("FalseDeltaCollId"), fwdtrack.collisionId() - mfttrack.collisionId());
            registry.fill(HIST("FalseDeltaXY"), Delta_XY);
            registry.fill(HIST("FalseDeltaPhiTanl"), Delta_PhiTanl);
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
