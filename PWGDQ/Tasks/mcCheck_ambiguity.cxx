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
#include "TDatabasePDG.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"

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
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace std;
using MyEventsVtxCovSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov, aod::ReducedMCEventLabels>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::ReducedMuonsLabels>;

struct mccheckambiguity {

  Preslice<aod::ReducedMuonsAssoc> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(MyEventsVtxCovSelected const& events, aod::ReducedMuonsAssoc const& muonAssocs, MyMuonTracksWithCov const& muons, ReducedMCEvents const& mcEvents, ReducedMCTracks const& mcTracks)
  {
    MCProng prongHc(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongHc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalHc;
    signalHc = new MCSignal("muFromHc", "Electrons from open charmed hadron decays", {prongHc}, {-1});

    MCProng prongHb(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongHb.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalHb;
    signalHb = new MCSignal("muFromHb", "Electrons from open beauty hadron decays", {prongHb}, {-1});

    auto groupedMuonAssocs = assocs.sliceBy(muonAssocsPerCollision, event.globalIndex());

    for (auto& muonAssoc : muonAssocs) {
      auto muon = muonAssoc.template reducedtrack_as<MyMuonTracksWithCov>();
      LOGF(info, "particle loop HF, pdg code = %d", particle.pdgCode());
      checked = sig.CheckSignal(true, mctrack);
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mccheckambiguity>(cfgc)
  };
}
