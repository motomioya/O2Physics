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

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct masscheck {
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

  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  HistogramRegistry registry{
    "registry",
    {
      {"hInvariantMass", "Invariant Mass of MCH standalone track;Invariant Mass (GeV/c^{2});Counts", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hCounter", "hCounter", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"Pt_vs_chi2mftmch", "Pt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"Eta_vs_chi2mftmch", "Eta_vs_chi2mftmch", {HistType::kTH2F, {{300, -5, -2},{400,0,200}}}},
      {"Pdca_vs_chi2mftmch", "Pdca_vs_chi2mftmch", {HistType::kTH2F, {{5000, 0, 5000},{400,0,200}}}},
      {"RAbs_vs_chi2mftmch", "RAbs_vs_chi2mftmch", {HistType::kTH2F, {{1500, 0, 150},{400,0,200}}}},
      {"McPt_vs_chi2mftmch", "McPt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McEta_vs_chi2mftmch", "McEta_vs_chi2mftmch", {HistType::kTH2F, {{300, -5, -2},{400,0,200}}}},
      {"McPdca_vs_chi2mftmch", "McPdca_vs_chi2mftmch", {HistType::kTH2F, {{5000, 0, 5000},{400,0,200}}}},
      {"McRAbs_vs_chi2mftmch", "McRAbs_vs_chi2mftmch", {HistType::kTH2F, {{1500, 0, 150},{400,0,200}}}},
      {"McTruePt_vs_chi2mftmch", "McTruePt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McTrueEta_vs_chi2mftmch", "McTrueEta_vs_chi2mftmch", {HistType::kTH2F, {{300, -5, -2},{400,0,200}}}},
      {"McTruePdca_vs_chi2mftmch", "McTruePdca_vs_chi2mftmch", {HistType::kTH2F, {{5000, 0, 5000},{400,0,200}}}},
      {"McTrueRAbs_vs_chi2mftmch", "McTrueRAbs_vs_chi2mftmch", {HistType::kTH2F, {{1500, 0, 150},{400,0,200}}}},
      {"McFalsePt_vs_chi2mftmch", "McFalsePt_vs_chi2mftmch", {HistType::kTH2F, {{200, 0, 20},{400,0,200}}}},
      {"McFalseEta_vs_chi2mftmch", "McFalseEta_vs_chi2mftmch", {HistType::kTH2F, {{300, -5, -2},{400,0,200}}}},
      {"McFalsePdca_vs_chi2mftmch", "McFalsePdca_vs_chi2mftmch", {HistType::kTH2F, {{5000, 0, 5000},{400,0,200}}}},
      {"McFalseRAbs_vs_chi2mftmch", "McFalseRAbs_vs_chi2mftmch", {HistType::kTH2F, {{1500, 0, 150},{400,0,200}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions const& collisions, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto& [fwdtrack1,fwdtrack2] : combinations(CombinationsStrictlyUpperIndexPolicy(fwdtracks, fwdtracks))) {
      if (fwdtrack1.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack2.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        if (fwdtrack1.collisionId() == fwdtrack2.collisionId()) {
          TLorentzVector lv1, lv2, lv;
          lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
          lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
          lv = lv1 + lv2;
          registry.fill(HIST("hInvariantMass"), lv.M());
          if (lv.M() < 0.4 && 0.2 < lv.M()) {
            registry.fill(HIST("Pt_vs_chi2mftmch"), fwdtrack1.pt(),fwdtrack1.chi2MatchMCHMFT());
            registry.fill(HIST("Pt_vs_chi2mftmch"), fwdtrack2.pt(),fwdtrack2.chi2MatchMCHMFT());
            registry.fill(HIST("Eta_vs_chi2mftmch"), fwdtrack1.eta(),fwdtrack1.chi2MatchMCHMFT());
            registry.fill(HIST("Eta_vs_chi2mftmch"), fwdtrack2.eta(),fwdtrack2.chi2MatchMCHMFT());
            registry.fill(HIST("Pdca_vs_chi2mftmch"), fwdtrack1.pDca(),fwdtrack1.chi2MatchMCHMFT());
            registry.fill(HIST("Pdca_vs_chi2mftmch"), fwdtrack2.pDca(),fwdtrack2.chi2MatchMCHMFT());
            registry.fill(HIST("RAbs_vs_chi2mftmch"), fwdtrack1.rAtAbsorberEnd(),fwdtrack1.chi2MatchMCHMFT());
            registry.fill(HIST("RAbs_vs_chi2mftmch"), fwdtrack2.rAtAbsorberEnd(),fwdtrack2.chi2MatchMCHMFT());
          }
        }
      }
    }
  }

  void processGen(aod::Collisions const& collisions, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    for (auto& [fwdtrack1,fwdtrack2] : combinations(CombinationsStrictlyUpperIndexPolicy(fwdtracks, fwdtracks))) {
      if (fwdtrack1.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack2.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        if (fwdtrack1.collisionId() == fwdtrack2.collisionId()) {
          TLorentzVector lv1, lv2, lv;
          lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
          lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
          lv = lv1 + lv2;
          registry.fill(HIST("hInvariantMass"), lv.M());
          if (lv.M() < 0.4 && 0.2 < lv.M()) {
            for (auto& mfttrack: mfttracks){
              if (fwdtrack1.has_mcParticle() && fwdtrack2.has_mcParticle() && mfttrack.has_mcParticle()){
                auto fwdparticle1 = fwdtrack1.mcParticle();
                auto fwdparticle2 = fwdtrack2.mcParticle();
                auto mftparticle = mfttrack.mcParticle();
                if (fwdparticle1.fromBackgroundEvent() == 1 && fwdparticle2.fromBackgroundEvent() == 1 && mftparticle.fromBackgroundEvent() == 1) {
                  if (fwdtrack1.matchMFTTrackId() == mfttrack.globalIndex()){
                    registry.fill(HIST("McPt_vs_chi2mftmch"), fwdtrack1.pt(),fwdtrack1.chi2MatchMCHMFT());
                    registry.fill(HIST("McEta_vs_chi2mftmch"), fwdtrack1.eta(),fwdtrack1.chi2MatchMCHMFT());
                    registry.fill(HIST("McPdca_vs_chi2mftmch"), fwdtrack1.pDca(),fwdtrack1.chi2MatchMCHMFT());
                    registry.fill(HIST("McRAbs_vs_chi2mftmch"), fwdtrack1.rAtAbsorberEnd(),fwdtrack1.chi2MatchMCHMFT());
                    if (fwdparticle1.globalIndex() == mftparticle.globalIndex()){
                      registry.fill(HIST("McTruePt_vs_chi2mftmch"), fwdtrack1.pt(),fwdtrack1.chi2MatchMCHMFT());
                      registry.fill(HIST("McTrueEta_vs_chi2mftmch"), fwdtrack1.eta(),fwdtrack1.chi2MatchMCHMFT());
                      registry.fill(HIST("McTruePdca_vs_chi2mftmch"), fwdtrack1.pDca(),fwdtrack1.chi2MatchMCHMFT());
                      registry.fill(HIST("McTrueRAbs_vs_chi2mftmch"), fwdtrack1.rAtAbsorberEnd(),fwdtrack1.chi2MatchMCHMFT());
                    } else {
                      registry.fill(HIST("McFalsePt_vs_chi2mftmch"), fwdtrack1.pt(),fwdtrack1.chi2MatchMCHMFT());
                      registry.fill(HIST("McFalseEta_vs_chi2mftmch"), fwdtrack1.eta(),fwdtrack1.chi2MatchMCHMFT());
                      registry.fill(HIST("McFalsePdca_vs_chi2mftmch"), fwdtrack1.pDca(),fwdtrack1.chi2MatchMCHMFT());
                      registry.fill(HIST("McFalseRAbs_vs_chi2mftmch"), fwdtrack1.rAtAbsorberEnd(),fwdtrack1.chi2MatchMCHMFT());
                    }
                  }
                  if (fwdtrack2.matchMFTTrackId() == mfttrack.globalIndex()){
                    registry.fill(HIST("McPt_vs_chi2mftmch"), fwdtrack2.pt(),fwdtrack2.chi2MatchMCHMFT());
                    registry.fill(HIST("McEta_vs_chi2mftmch"), fwdtrack2.eta(),fwdtrack2.chi2MatchMCHMFT());
                    registry.fill(HIST("McPdca_vs_chi2mftmch"), fwdtrack2.pDca(),fwdtrack2.chi2MatchMCHMFT());
                    registry.fill(HIST("McRAbs_vs_chi2mftmch"), fwdtrack2.rAtAbsorberEnd(),fwdtrack2.chi2MatchMCHMFT());
                    if (fwdparticle2.globalIndex() == mftparticle.globalIndex()){
                      registry.fill(HIST("McTruePt_vs_chi2mftmch"), fwdtrack2.pt(),fwdtrack2.chi2MatchMCHMFT());
                      registry.fill(HIST("McTrueEta_vs_chi2mftmch"), fwdtrack2.eta(),fwdtrack2.chi2MatchMCHMFT());
                      registry.fill(HIST("McTruePdca_vs_chi2mftmch"), fwdtrack2.pDca(),fwdtrack2.chi2MatchMCHMFT());
                      registry.fill(HIST("McTrueRAbs_vs_chi2mftmch"), fwdtrack2.rAtAbsorberEnd(),fwdtrack2.chi2MatchMCHMFT());
                    } else {
                      registry.fill(HIST("McFalsePt_vs_chi2mftmch"), fwdtrack2.pt(),fwdtrack2.chi2MatchMCHMFT());
                      registry.fill(HIST("McFalseEta_vs_chi2mftmch"), fwdtrack2.eta(),fwdtrack2.chi2MatchMCHMFT());
                      registry.fill(HIST("McFalsePdca_vs_chi2mftmch"), fwdtrack2.pDca(),fwdtrack2.chi2MatchMCHMFT());
                      registry.fill(HIST("McFalseRAbs_vs_chi2mftmch"), fwdtrack2.rAtAbsorberEnd(),fwdtrack2.chi2MatchMCHMFT());
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
  PROCESS_SWITCH(masscheck, processGen, "Process generator-level info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<masscheck>(cfgc)
  };

}
