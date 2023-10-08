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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

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

  HistogramRegistry registry{
    "registry",
    {
      {"TrackType", "TType", {HistType::kTH1F, {{5, 0, 5}}}},
      {"GlobalMuonTrackPt", "GMTPT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"MuonStandaloneTrackPt", "MSPT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"PDca", "PDca", {HistType::kTH1F, {{10000, 0, 1000}}}},
      {"NFTrackType", "TType", {HistType::kTH1F, {{5, 0, 5}}}},
      {"NFGlobalMuonTrackPt", "GMTPT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"NFMuonStandaloneTrackPt", "MSPT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"NFPDca", "NFPDca", {HistType::kTH1F, {{10000, 0, 1000}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions::iterator const& collision, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto const& fwdtrack : fwdtracks){
      registry.fill(HIST("TrackType"), fwdtrack.trackType());
      registry.fill(HIST("PDca"), fwdtrack.pDca());
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack){
        registry.fill(HIST("GlobalMuonTrackPt"), fwdtrack.pt());
      }
      else if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack){
        registry.fill(HIST("MuonStandaloneTrackPt"), fwdtrack.pt());
      }
    }
  }

  void processNoFilter(aod::Collisions::iterator const& collision, aod::FwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto const& fwdtrack : fwdtracks){
      registry.fill(HIST("NFTrackType"), fwdtrack.trackType());
      registry.fill(HIST("NFPDca"), fwdtrack.pDca());
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack){
        registry.fill(HIST("NFGlobalMuonTrackPt"), fwdtrack.pt());
      }
      else if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack){
        registry.fill(HIST("NFMuonStandaloneTrackPt"), fwdtrack.pt());
      }
    }
  }
  PROCESS_SWITCH(quickcheck, processNoFilter, "Process generator-level info", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<quickcheck>(cfgc)
  };
}
