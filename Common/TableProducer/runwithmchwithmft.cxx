// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <math.h>
#include <string>
#include <regex>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/mftmchMatchingML.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "CCDB/CcdbApi.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct runwithmchwithmft {
  Produces<aod::FwdTracksML> fwdtrackml;

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
      {"hCounter", "hCounter", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions::iterator const& collision, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto& fwdtrack : fwdtracks) {
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        registry.fill(HIST("hCounter"), 0);
        for (auto& fwdtrackglobal : fwdtracks) {
          if (fwdtrackglobal.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
            if (fwdtrackglobal.matchMCHTrackId() == fwdtrack.globalIndex()) {
              registry.fill(HIST("hCounter"), 1);
              fwdtrackml(fwdtrack.collisionId(),2,fwdtrack.x(),fwdtrack.y(),fwdtrack.z(),fwdtrack.phi(),fwdtrack.tgl(),fwdtrack.signed1Pt(),fwdtrack.nClusters(),fwdtrack.pDca(),fwdtrack.rAtAbsorberEnd(),fwdtrack.chi2(),fwdtrack.chi2MatchMCHMID(),fwdtrack.chi2MatchMCHMID(),fwdtrack.matchScoreMCHMFT(),fwdtrack.matchMFTTrackId(),fwdtrack.matchMCHTrackId(),fwdtrack.mchBitMap(),fwdtrack.midBitMap(),fwdtrack.midBoards(),fwdtrack.trackTime(),fwdtrack.trackTimeRes(), fwdtrack.eta(), fwdtrack.pt(),fwdtrack.p(),0,0);
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
    adaptAnalysisTask<runwithmchwithmft>(cfgc)
  };
}
