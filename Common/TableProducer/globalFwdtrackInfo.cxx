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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace evsel;

struct globalFwdtrackInfo {

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
  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (etaup < aod::fwdtrack::eta ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) || (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) || (aod::fwdtrack::pDca < pDCAcutdcaup1)) && ((aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndlow2) || (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) || (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  HistogramRegistry registry{
    "registry", //
    {
      {"Pt_vs_chi2mftmch", "Pt_vs_chi2mftmch", {HistType::kTH2F, {{300, 0, 30},{300,0,300}}}},
      {"chi2mft_vs_chi2mftmch", "chi2mft_vs_chi2mftmch", {HistType::kTH2F, {{300,0,300},{300,0,300}}}},
      {"chi2global_vs_chi2mftmch", "chi2mft_vs_chi2mftmch", {HistType::kTH2F, {{300,0,300},{300,0,300}}}},
      {"chi2mchmid_vs_chi2mftmch", "chi2mft_vs_chi2mftmch", {HistType::kTH2F, {{300,0,300},{300,0,300}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions::iterator const& collision, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto& fwdtrack: fwdtracks) {
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack){
        registry.fill(HIST("Pt_vs_chi2mftmch"), fwdtrack.pt(),fwdtrack.chi2MatchMCHMFT());
        registry.fill(HIST("chi2global_vs_chi2mftmch"), fwdtrack.chi2(),fwdtrack.chi2MatchMCHMFT());
        registry.fill(HIST("chi2mchmid_vs_chi2mftmch"), fwdtrack.chi2MatchMCHMID(),fwdtrack.chi2MatchMCHMFT());
        for (auto& mfttrack: mfttracks){
          if (fwdtrack.matchMFTTrackId() == mfttrack.globalIndex()){
            registry.fill(HIST("chi2mft_vs_chi2mftmch"), mfttrack.chi2(),fwdtrack.chi2MatchMCHMFT());
          }
        }
      }
    }

  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<globalFwdtrackInfo>(cfgc)
  };
}
