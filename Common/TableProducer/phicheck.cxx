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

struct phicheck {
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
      {"PphiGlobal", "PphiGlobal", {HistType::kTH1F, {{2000, -10, 10}}}},
      {"MphiGlobal", "MphiGlobal", {HistType::kTH1F, {{2000, -10, 10}}}},
      {"PphiMuon", "PphiMuon", {HistType::kTH1F, {{2000, -10, 10}}}},
      {"MphiMuon", "MphiMuon", {HistType::kTH1F, {{2000, -10, 10}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions const& collisions, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto& [fwdtrack1,fwdtrack2] : combinations(CombinationsStrictlyUpperIndexPolicy(fwdtracks, fwdtracks))) {
      if (fwdtrack1.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack2.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        if (fwdtrack1.sign() > 0) registry.fill(HIST("PphiGlobal"), fwdtrack1.phi());
        if (fwdtrack2.sign() < 0) registry.fill(HIST("MphiGlobal"), fwdtrack1.phi());

      } else if (fwdtrack1.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack && fwdtrack2.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        if (fwdtrack2.sign() > 0) registry.fill(HIST("PphiMuon"), fwdtrack2.phi());
        if (fwdtrack2.sign() < 0) registry.fill(HIST("MphiMuon"), fwdtrack2.phi());
      }
    }
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phicheck>(cfgc)
  };
}
