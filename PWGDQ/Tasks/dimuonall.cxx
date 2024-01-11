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
#include "PWGDQ/DataModel/ReducedInfoTables.h"

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

struct dimuonall {

  HistogramRegistry registry{
    /*
    "registry", 
    {
      {"massPM1", "massPM1", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"massMMPP1", "massMMPP1", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"massPM24", "massPM24", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"massMMPP24", "massMMPP24", {HistType::kTH1F, {{2000, 0, 20}}}}
    },
    */
    "registry", 
    {
      {"massPM0", "massPM0", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"massMMPP0", "massMMPP0", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"massPM1", "massPM1", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"oaPM1", "oaPM1", {HistType::kTH1F, {{2000, -10, 10}}}},
      {"oaMMPP1", "oaMMPP1", {HistType::kTH1F, {{2000, -10, 10}}}},
      {"massMMPP1", "massMMPP1", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"massPM3", "massPM3", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"massMMPP3", "massMMPP3", {HistType::kTH1F, {{2000, 0, 20}}}}
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  /*
  void process(aod::DimuonsAll const& dimuonsall)
  {
    for (auto& dimuon : dimuonsall) {
      if (dimuon.sign() == 0) {
        if (dimuon.mcDecision() == 1) {
          registry.fill(HIST("massPM1"), dimuon.mass());
        } else if (dimuon.mcDecision() == 2 || dimuon.mcDecision() == 4 ) {
          registry.fill(HIST("massPM24"), dimuon.mass());
        }
      } else {
        if (dimuon.mcDecision() == 1) {
          registry.fill(HIST("massMMPP1"), dimuon.mass());
        } else if (dimuon.mcDecision() == 2 || dimuon.mcDecision() == 4 ) {
          registry.fill(HIST("massMMPP24"), dimuon.mass());
        }
      }
    }
  }
  */

  void process(aod::DimuonsAll const& dimuonsall)
  {
    for (auto& dimuon : dimuonsall) {
      if (dimuon.sign() == 0) {
        if (dimuon.mcDecision() == 1) {
          registry.fill(HIST("massPM1"), dimuon.mass());
          registry.fill(HIST("oaPM1"), dimuon.phi1() - dimuon.phi2());
        } else if (dimuon.mcDecision() == 3) {
          registry.fill(HIST("massPM3"), dimuon.mass());
        } else if (dimuon.mcDecision() == 0) {
          registry.fill(HIST("massPM0"), dimuon.mass());
        }
      } else {
        if (dimuon.mcDecision() == 1) {
          registry.fill(HIST("massMMPP1"), dimuon.mass());
          registry.fill(HIST("oaMMPP1"), dimuon.phi1() - dimuon.phi2());
        } else if (dimuon.mcDecision() == 3) {
          registry.fill(HIST("massMMPP3"), dimuon.mass());
        } else if (dimuon.mcDecision() == 0) {
          registry.fill(HIST("massMMPP0"), dimuon.mass());
        }
      }
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dimuonall>(cfgc)
  };
}
