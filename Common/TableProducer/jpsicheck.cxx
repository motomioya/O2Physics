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

struct jpsicheck {
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
      {"hInvariantMass", "Invariant Mass of mch-mid track;Invariant Mass (GeV/c^{2});Counts", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hInvariantMassGlobalIndependent", "Invariant Mass of mft-mch-mid track;Invariant Mass (GeV/c^{2});Counts", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hCounter", "hCounter", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"hPtGlobalIndependent", "hPtGlobalIndependent", {HistType::kTH1F, {{100, 0, 10}}}},
      {"hPtGlobal", "hPtGlobal", {HistType::kTH1F, {{100, 0, 10}}}},
      {"hPtMuon", "hPtMuon", {HistType::kTH1F, {{100, 0, 10}}}},
      {"hRapidityMuon", "hRapidityMuon", {HistType::kTH1F, {{50, -5, 0}}}},
      {"hRapidityGlobal", "hRapidityGlobal", {HistType::kTH1F, {{50, -5, 0}}}},
      {"hRapidityGlobalIndependent", "hRapidityGlobalIndependent", {HistType::kTH1F, {{50, -5, 0}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions::iterator const& collision, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto& [fwdtrack1,fwdtrack2] : combinations(CombinationsStrictlyUpperIndexPolicy(fwdtracks, fwdtracks))) {
      if (fwdtrack1.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack && fwdtrack2.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        if (fwdtrack1.sign() != fwdtrack2.sign()) {
          TLorentzVector lv1, lv2, lv;
          lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
          lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
          lv = lv1 + lv2;
          registry.fill(HIST("hInvariantMass"), lv.M());
          if (lv.M() < 3.093 + 0.072*3 && 3.093 - 0.072*3 < lv.M()) {
            registry.fill(HIST("hCounter"), 0);
            registry.fill(HIST("hCounter"), 0);
            registry.fill(HIST("hPtMuon"), fwdtrack1.pt());
            registry.fill(HIST("hPtMuon"), fwdtrack2.pt());
            registry.fill(HIST("hRapidityMuon"),fwdtrack1.eta());
            registry.fill(HIST("hRapidityMuon"),fwdtrack2.eta());
            for (auto& fwdtrackjpsi : fwdtracks) {
              if (fwdtrackjpsi.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
                if (fwdtrackjpsi.matchMCHTrackId() == fwdtrack1.globalIndex()) {
                  registry.fill(HIST("hCounter"), 1);
                  registry.fill(HIST("hPtGlobal"), fwdtrack1.pt());
                  registry.fill(HIST("hRapidityGlobal"),fwdtrack1.eta());
                }
                if (fwdtrackjpsi.matchMCHTrackId() == fwdtrack2.globalIndex()) {
                  registry.fill(HIST("hCounter"), 1);
                  registry.fill(HIST("hPtGlobal"), fwdtrack2.pt());
                  registry.fill(HIST("hRapidityGlobal"),fwdtrack2.eta());
                }
              }
            }
          }
        }
      }
      if (fwdtrack1.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack2.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        if (fwdtrack1.sign() != fwdtrack2.sign()) {
          TLorentzVector lv1, lv2, lv;
          lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
          lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
          lv = lv1 + lv2;
          registry.fill(HIST("hInvariantMassGlobalIndependent"), lv.M());
          if (lv.M() < 3.093 + 0.072*3 && 3.093 - 0.072*3 < lv.M()) {
            registry.fill(HIST("hPtGlobalIndependent"), fwdtrack1.pt());
            registry.fill(HIST("hPtGlobalIndependent"), fwdtrack2.pt());
            registry.fill(HIST("hRapidityGlobalIndependent"), fwdtrack1.eta());
            registry.fill(HIST("hRapidityGlobalIndependent"), fwdtrack2.eta());
          }
        }
      }
    }
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<jpsicheck>(cfgc)
  };

}
