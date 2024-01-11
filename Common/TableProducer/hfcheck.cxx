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
#include "DCAFitter/FwdDCAFitterN.h"
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
#include <TLorentzVector.h>
#include "TDatabasePDG.h"
#include "Math/Vector3D.h"
#include <vector>
#include <cmath>
#include <utility>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using MyMuonsColl = soa::Join<aod::FwdTracks, aod::FwdTrkCompColls>;
using MyMFTsColl = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;
using Vec3D = ROOT::Math::SVector<double, 3>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct hfcheck {
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

  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();
  int kaonPDGCode = 311;
  TParticlePDG* kaonParticle = TDatabasePDG::Instance()->GetParticle(kaonPDGCode);
  double kaonMass = kaonParticle->Mass();

  HistogramRegistry registry{
    "registry",
    {
      {"hMuInclusivePt", "hMuInclusivePt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuInclusiveP", "hMuInclusiveP", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuFromHfPt", "hMuFromHfPt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuFromHfP", "hMuFromHfP", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuMfromHfLxyz", "hMuMfromHfLxyz", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuMfromHfPCA", "hMuMfromHfPCA", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuFromHadronPt", "hMuFromHadronPt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuFromHadronP", "hMuFromHadronP", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuFromHadronZpos", "hMuFromHadronZpos", {HistType::kTH1F, {{1000,-1000,0}}}},
      {"hConter", "hConter", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
      {"hInvariantMassHF", "hInvariantMassHF", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hInvariantMassLF", "hInvariantMassLF", {HistType::kTH1F, {{5000, 0, 50}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  template <typename TFwdTrack, typename TMFTTrack, typename TCollision>
  std::pair<double, float> fitDCA(TFwdTrack const& fwdtrack, TMFTTrack const& mfttrack, TCollision const& collision)
  {
    std::pair<double, float> pair_fitDCA;

    o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;
    Vec3D secondaryVertex;

    fgFitterTwoProngFwd.setBz(5.0f);
    fgFitterTwoProngFwd.setPropagateToPCA(true);
    fgFitterTwoProngFwd.setMaxR(200.0f);
    fgFitterTwoProngFwd.setMinParamChange(1.0e-3f);
    fgFitterTwoProngFwd.setMinRelChi2Change(0.9f);
    fgFitterTwoProngFwd.setUseAbsDCA(true);

    double chi21 = fwdtrack.chi2();
    double chi22 = mfttrack.chi2();

    SMatrix5 fwdtrackpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
    std::vector<double> v1;
    SMatrix55 fwdtrackcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd pars1{fwdtrack.z(), fwdtrackpars, fwdtrackcovs, chi21};

    SMatrix5 mfttrackpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
    std::vector<double> v2;
    SMatrix55 mfttrackcovs(v2.begin(), v2.end());
    o2::track::TrackParCovFwd pars2{mfttrack.z(), mfttrackpars, mfttrackcovs, chi22};

    int procCode = fgFitterTwoProngFwd.process(pars1,pars2);

    if (procCode != 0) {
      secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
  double kVertexingLxy = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) + (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
      double kVertexingLz = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
      double kVertexingLxyz = kVertexingLxy + kVertexingLz;
      kVertexingLxy = std::sqrt(kVertexingLxy);
      kVertexingLz = std::sqrt(kVertexingLz);
      kVertexingLxyz = std::sqrt(kVertexingLxyz);
      float kPCAChi2 = fgFitterTwoProngFwd.getChi2AtPCACandidate();
      pair_fitDCA =  std::make_pair(kVertexingLxyz, kPCAChi2);
      return pair_fitDCA;
    }
  }

  void process(aod::Collisions::iterator const& collision, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, aod::MFTTracks const& mfttracks, aod::McParticles const& mcparticles)
  {
    for (auto& fwdtrack : fwdtracks) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        if (fwdtrack.has_mcParticle()){
          registry.fill(HIST("hMuInclusivePt"), fwdtrack.pt());
          registry.fill(HIST("hMuInclusiveP"), fwdtrack.p());
          auto fwdparticle = fwdtrack.mcParticle();
          const auto mcfwdtrackmothers = fwdtrack.mcParticle().mothersIds();
          registry.fill(HIST("hConter"), 0);
          if (mcfwdtrackmothers.size() > 0) {
            registry.fill(HIST("hConter"), 1);
            int mcfwdfirstmotherid = mcfwdtrackmothers[0];
            auto mcfwdfirstmother = mcparticles.iteratorAt(mcfwdfirstmotherid);
            //int mcfwdfirstmotherpdg = mcfwdfirstmother.pdgCode();
            int mcfwdlastmotherid = mcfwdtrackmothers[mcfwdtrackmothers.size() - 1];
            auto mcfwdlastmother = mcparticles.iteratorAt(mcfwdlastmotherid);
            int mcfwdlastmotherpdg = mcfwdlastmother.pdgCode();
            if (mcfwdfirstmother.fromBackgroundEvent() == false) {
              registry.fill(HIST("hConter"), 2);
              if ( (mcfwdlastmotherpdg >= 400 && mcfwdlastmotherpdg <= 439) || (mcfwdlastmotherpdg <= -400 && mcfwdlastmotherpdg >= -439) ) {
                registry.fill(HIST("hMuFromHfPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHfP"), fwdtrack.p());
                registry.fill(HIST("hConter"), 3);
                for (auto& mfttrack : mfttracks) {

                  std::pair<double, float> res_fitDCA = fitDCA(fwdtrack, mfttrack, collision);
                  registry.fill(HIST("hMuMfromHfLxyz"), res_fitDCA.first);
                  registry.fill(HIST("hMuMfromHfPCA"), res_fitDCA.second);

                  TLorentzVector lv1, lv2, lv;
                  lv1.SetPtEtaPhiM(fwdtrack.pt(), fwdtrack.eta(), fwdtrack.phi(), muonMass);
                  lv2.SetPtEtaPhiM(mfttrack.pt(), mfttrack.eta(), mfttrack.phi(), kaonMass);
                  lv = lv1 + lv2;
                  registry.fill(HIST("hInvariantMassHF"), lv.M());
                }
              }
            } else {
              registry.fill(HIST("hConter"), 4);
              if ( (mcfwdlastmotherpdg >= 100 && mcfwdlastmotherpdg <= 119) || (mcfwdlastmotherpdg <= -100 && mcfwdlastmotherpdg >= -119) || (mcfwdlastmotherpdg >= 1000 && mcfwdlastmotherpdg <= 1999) || (mcfwdlastmotherpdg <= -1000 && mcfwdlastmotherpdg >= -1999) || (mcfwdlastmotherpdg >= 200 && mcfwdlastmotherpdg <= 299) || (mcfwdlastmotherpdg <= -200 && mcfwdlastmotherpdg >= -299) || (mcfwdlastmotherpdg >= 2000 && mcfwdlastmotherpdg <= 2999) || (mcfwdlastmotherpdg <= -2000 && mcfwdlastmotherpdg >= -2999) || (mcfwdlastmotherpdg >= 300 && mcfwdlastmotherpdg <= 399) || (mcfwdlastmotherpdg <= -300 && mcfwdlastmotherpdg >= -399) || (mcfwdlastmotherpdg >= 3000 && mcfwdlastmotherpdg <= 3999) || (mcfwdlastmotherpdg <= -3000 && mcfwdlastmotherpdg >= -3999) ) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronP"), fwdtrack.p());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hConter"), 5);
                for (auto& mfttrack : mfttracks) {
                  TLorentzVector lv1, lv2, lv;
                  lv1.SetPtEtaPhiM(fwdtrack.pt(), fwdtrack.eta(), fwdtrack.phi(), muonMass);
                  lv2.SetPtEtaPhiM(mfttrack.pt(), mfttrack.eta(), mfttrack.phi(), kaonMass);
                  lv = lv1 + lv2;
                  registry.fill(HIST("hInvariantMassLF"), lv.M());
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
    adaptAnalysisTask<hfcheck>(cfgc)
  };
}
