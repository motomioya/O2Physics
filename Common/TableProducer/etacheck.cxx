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

struct etacheck {
  float etalow = -4.0;
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
      {"hEtaGlobal", "hEtaGlobal", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hEtaGlobalTrue", "hEtaGlobalTrue", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hEtaGlobalFalse", "hEtaGlobalFalse", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hEtaMatchedMCH", "hEtaMatchedMCH", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hEtaMatchedMCHTrue", "hEtaMatchedMCHTrue", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hEtaMatchedMCHFalse", "hEtaMatchedMCHFalse", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hEtaMatchedMFT", "hEtaMatchedMFT", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hEtaMatchedMFTTrue", "hEtaMatchedMFTTrue", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hEtaMatchedMFTFalse", "hEtaMatchedMFTFalse", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hDeltaEtaMCH", "hDeltaEtaMCH", {HistType::kTH1F, {{500, -2.5, 2.5}}}},
      {"hDeltaEtaMCHTrue", "hDeltaEtaMCHTrue", {HistType::kTH1F, {{500, -2.5, 2.5}}}},
      {"hDeltaEtaMCHFalse", "hDeltaEtaMCHFalse", {HistType::kTH1F, {{500, -2.5, 2.5}}}},
      {"hDeltaEtaMFT", "hDeltaEtaMFT", {HistType::kTH1F, {{500, -2.5, 2.5}}}},
      {"hDeltaEtaMFTTrue", "hDeltaEtaMFTTrue", {HistType::kTH1F, {{500, -2.5, 2.5}}}},
      {"hDeltaEtaMFTFalse", "hDeltaEtaMFTFalse", {HistType::kTH1F, {{500, -2.5, 2.5}}}},
      {"hPtGlobal", "hPtGlobal", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hPtGlobalTrue", "hPtGlobalTrue", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hPtGlobalFalse", "hPtGlobalFalse", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hPtMatchedMCH", "hPtMatchedMCH", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hPtMatchedMCHTrue", "hPtMatchedMCHTrue", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hPtMatchedMCHFalse", "hPtMatchedMCHFalse", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hPtMatchedMFT", "hPtMatchedMFT", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hPtMatchedMFTTrue", "hPtMatchedMFTTrue", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hPtMatchedMFTFalse", "hPtMatchedMFTFalse", {HistType::kTH1F, {{500, 0, 10}}}},
      {"hDeltaPtMCH", "hDeltaPtMCH", {HistType::kTH1F, {{500, -5, 5}}}},
      {"hDeltaPtMCHTrue", "hDeltaPtMCHTrue", {HistType::kTH1F, {{500, -5, 5}}}},
      {"hDeltaPtMCHFalse", "hDeltaPtMCHFalse", {HistType::kTH1F, {{500, -5, 5}}}},
      {"hDeltaPtMFT", "hDeltaPtMFT", {HistType::kTH1F, {{500, -5, 5}}}},
      {"hDeltaPtMFTTrue", "hDeltaPtMFTTrue", {HistType::kTH1F, {{500, -5, 5}}}},
      {"hDeltaPtMFTFalse", "hDeltaPtMFTFalse", {HistType::kTH1F, {{500, -5, 5}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }


  void process(soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&)
  {
    for (auto& fwdtrack : fwdtracks) {
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack.chi2MatchMCHMFT() < 40) {
        if (fwdtrack.eta() > -4 && fwdtrack.eta() < -2.5 && (((17.6 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 26.5) && (fwdtrack.pDca() < 594)) || ((26.5 < fwdtrack.rAtAbsorberEnd()) && (fwdtrack.rAtAbsorberEnd() < 89.5) && (fwdtrack.pDca() < 324)))){
          for (auto& mchtrack : fwdtracks) {
            if (fwdtrack.matchMCHTrackId() == mchtrack.globalIndex()) {
              registry.fill(HIST("hEtaGlobal"), fwdtrack.eta());
              registry.fill(HIST("hEtaMatchedMCH"), mchtrack.eta());
              registry.fill(HIST("hDeltaEtaMCH"), fwdtrack.eta() - mchtrack.eta());
              registry.fill(HIST("hPtGlobal"), fwdtrack.pt());
              registry.fill(HIST("hPtMatchedMCH"), mchtrack.pt());
              registry.fill(HIST("hDeltaPtMCH"), fwdtrack.pt() - mchtrack.pt());
              for (auto& mfttrack : mfttracks) {
                if (fwdtrack.matchMFTTrackId() == mfttrack.globalIndex()) {
                  registry.fill(HIST("hEtaMatchedMFT"), mfttrack.eta());
                  registry.fill(HIST("hDeltaEtaMFT"), fwdtrack.eta() - mfttrack.eta());
                  registry.fill(HIST("hPtMatchedMFT"), mfttrack.pt());
                  registry.fill(HIST("hDeltaPtMFT"), fwdtrack.pt() - mfttrack.pt());
                  if (mchtrack.mcParticleId() == mfttrack.mcParticleId()) {
                    registry.fill(HIST("hEtaGlobalTrue"), fwdtrack.eta());
                    registry.fill(HIST("hEtaMatchedMCHTrue"), mchtrack.eta());
                    registry.fill(HIST("hDeltaEtaMCHTrue"), fwdtrack.eta() - mchtrack.eta());
                    registry.fill(HIST("hEtaMatchedMFTTrue"), mfttrack.eta());
                    registry.fill(HIST("hDeltaEtaMFTTrue"), fwdtrack.eta() - mfttrack.eta());
                    registry.fill(HIST("hPtGlobalTrue"), fwdtrack.pt());
                    registry.fill(HIST("hPtMatchedMCHTrue"), mchtrack.pt());
                    registry.fill(HIST("hDeltaPtMCHTrue"), fwdtrack.pt() - mchtrack.pt());
                    registry.fill(HIST("hPtMatchedMFTTrue"), mfttrack.pt());
                    registry.fill(HIST("hDeltaPtMFTTrue"), fwdtrack.pt() - mfttrack.pt());
                  } else {
                    registry.fill(HIST("hEtaGlobalFalse"), fwdtrack.eta());
                    registry.fill(HIST("hEtaMatchedMCHFalse"), mchtrack.eta());
                    registry.fill(HIST("hDeltaEtaMCHFalse"), fwdtrack.eta() - mchtrack.eta());
                    registry.fill(HIST("hEtaMatchedMFTFalse"), mfttrack.eta());
                    registry.fill(HIST("hDeltaEtaMFTFalse"), fwdtrack.eta() - mfttrack.eta());
                    registry.fill(HIST("hPtGlobalFalse"), fwdtrack.pt());
                    registry.fill(HIST("hPtMatchedMCHFalse"), mchtrack.pt());
                    registry.fill(HIST("hDeltaPtMCHFalse"), fwdtrack.pt() - mchtrack.pt());
                    registry.fill(HIST("hPtMatchedMFTFalse"), mfttrack.pt());
                    registry.fill(HIST("hDeltaPtMFTFalse"), fwdtrack.pt() - mfttrack.pt());
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
    adaptAnalysisTask<etacheck>(cfgc)
  };

}
