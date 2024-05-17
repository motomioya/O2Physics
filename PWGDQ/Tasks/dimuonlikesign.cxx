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

struct dimuonlikesign {

  HistogramRegistry registry{
    "registry", 
    {
      {"SV_LS", "SV_LS", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"SV_LS_lowmass", "SV_LS_lowmass", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"SV_LS_uppermass", "SV_LS_uppermass", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"mass_LS", "mass_LS", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"mass_LS_lowSV", "mass_LS_lowSV", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"mass_LS_upperSV", "mass_LS_upperSV", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"mass_LS_failedSV", "mass_LS_failedSV", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"SV_US", "SV_US", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"SV_US_lowmass", "SV_US_lowmass", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"SV_US_uppermass", "SV_US_uppermass", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"mass_US", "mass_US", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"mass_US_lowSV", "mass_US_lowSV", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"mass_US_upperSV", "mass_US_upperSV", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"mass_US_failedSV", "mass_US_failedSV", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::DimuonsAll const& dimuons)
  {
    for (auto& dimuon : dimuons) {
      if (dimuon.sign() == 0) {
        registry.fill(HIST("mass_US"), dimuon.mass());
        registry.fill(HIST("SV_US"), dimuon.sVertex());
        if (dimuon.mass() < 0.26) {
          registry.fill(HIST("SV_US_lowmass"), dimuon.sVertex());
        } else {
          registry.fill(HIST("SV_US_uppermass"), dimuon.sVertex());
        }
        if (dimuon.sVertex() == -999) {
          registry.fill(HIST("mass_US_failedSV"), dimuon.mass());
        } else if (dimuon.sVertex() < -40) {
          registry.fill(HIST("mass_US_lowSV"), dimuon.mass());
        } else {
          registry.fill(HIST("mass_US_upperSV"), dimuon.mass());
        }
      } else if (dimuon.sign() == -2 || dimuon.sign() == 2) {
        registry.fill(HIST("mass_LS"), dimuon.mass());
        registry.fill(HIST("SV_LS"), dimuon.sVertex());
        if (dimuon.mass() < 0.26) {
          registry.fill(HIST("SV_LS_lowmass"), dimuon.sVertex());
        } else {
          registry.fill(HIST("SV_LS_uppermass"), dimuon.sVertex());
        }
        if (dimuon.sVertex() == -999) {
          registry.fill(HIST("mass_LS_failedSV"), dimuon.mass());
        } else if (dimuon.sVertex() < -40) {
          registry.fill(HIST("mass_LS_lowSV"), dimuon.mass());
        } else {
          registry.fill(HIST("mass_LS_upperSV"), dimuon.mass());
        }
      }
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dimuonlikesign>(cfgc)
  };
}
