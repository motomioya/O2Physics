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
    "registry", 
    {
      {"mass_DCAmumuPM0", "mass_DCAmumuPM0", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM0", "mass_DCAmumuMM0", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP0", "mass_DCAmumuPP0", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM1", "mass_DCAmumuPM1", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM1", "mass_DCAmumuMM1", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP1", "mass_DCAmumuPP1", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM2", "mass_DCAmumuPM2", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM2", "mass_DCAmumuMM2", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP2", "mass_DCAmumuPP2", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM4", "mass_DCAmumuPM4", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM4", "mass_DCAmumuMM4", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP4", "mass_DCAmumuPP4", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM8", "mass_DCAmumuPM8", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM8", "mass_DCAmumuMM8", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP8", "mass_DCAmumuPP8", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM10", "mass_DCAmumuPM10", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM10", "mass_DCAmumuMM10", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP10", "mass_DCAmumuPP10", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM12", "mass_DCAmumuPM12", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM12", "mass_DCAmumuMM12", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP12", "mass_DCAmumuPP12", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM15", "mass_DCAmumuPM15", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM15", "mass_DCAmumuMM15", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP15", "mass_DCAmumuPP15", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM18", "mass_DCAmumuPM18", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM18", "mass_DCAmumuMM18", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP18", "mass_DCAmumuPP18", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM32", "mass_DCAmumuPM32", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM32", "mass_DCAmumuMM32", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP32", "mass_DCAmumuPP32", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM36", "mass_DCAmumuPM36", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM36", "mass_DCAmumuMM36", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP36", "mass_DCAmumuPP36", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM40", "mass_DCAmumuPM40", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM40", "mass_DCAmumuMM40", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP40", "mass_DCAmumuPP40", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM42", "mass_DCAmumuPM42", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM42", "mass_DCAmumuMM42", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP42", "mass_DCAmumuPP42", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM48", "mass_DCAmumuPM48", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM48", "mass_DCAmumuMM48", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP48", "mass_DCAmumuPP48", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM52", "mass_DCAmumuPM52", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM52", "mass_DCAmumuMM52", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP52", "mass_DCAmumuPP52", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPM57", "mass_DCAmumuPM57", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuMM57", "mass_DCAmumuMM57", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
      {"mass_DCAmumuPP57", "mass_DCAmumuPP57", {HistType::kTH2F, {{1000, 0, 20}, {1000, 0, 10}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::DimuonsAll const& dimuons)
  {
    for (auto& dimuon : dimuons) {
      float DCA1 = std::sqrt(dimuon.fwdDcaX1() * dimuon.fwdDcaX1() + dimuon.fwdDcaY1() * dimuon.fwdDcaY1());
      float DCA2 = std::sqrt(dimuon.fwdDcaX2() * dimuon.fwdDcaX2() + dimuon.fwdDcaY2() * dimuon.fwdDcaY2());
      float DCAmumu = std::sqrt((DCA1 * DCA1 + DCA2 * DCA2)/2);
      if (dimuon.sign() == 0) {
        if (dimuon.mcDecision() == 0) {
          registry.fill(HIST("mass_DCAmumuPM0"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 1) {
          registry.fill(HIST("mass_DCAmumuPM1"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 2) {
          registry.fill(HIST("mass_DCAmumuPM2"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 4) {
          registry.fill(HIST("mass_DCAmumuPM4"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 8) {
          registry.fill(HIST("mass_DCAmumuPM8"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 10) {
          registry.fill(HIST("mass_DCAmumuPM10"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 12) {
          registry.fill(HIST("mass_DCAmumuPM12"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 15) {
          registry.fill(HIST("mass_DCAmumuPM15"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 18) {
          registry.fill(HIST("mass_DCAmumuPM18"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 32) {
          registry.fill(HIST("mass_DCAmumuPM32"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 36) {
          registry.fill(HIST("mass_DCAmumuPM36"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 40) {
          registry.fill(HIST("mass_DCAmumuPM40"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 42) {
          registry.fill(HIST("mass_DCAmumuPM42"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 48) {
          registry.fill(HIST("mass_DCAmumuPM48"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 52) {
          registry.fill(HIST("mass_DCAmumuPM52"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 57) {
          registry.fill(HIST("mass_DCAmumuPM57"), dimuon.mass(), DCAmumu);
        }
      } else if (dimuon.sign() == -2) {
        if (dimuon.mcDecision() == 0) {
          registry.fill(HIST("mass_DCAmumuMM0"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 1) {
          registry.fill(HIST("mass_DCAmumuMM1"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 2) {
          registry.fill(HIST("mass_DCAmumuMM2"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 4) {
          registry.fill(HIST("mass_DCAmumuMM4"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 8) {
          registry.fill(HIST("mass_DCAmumuMM8"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 10) {
          registry.fill(HIST("mass_DCAmumuMM10"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 12) {
          registry.fill(HIST("mass_DCAmumuMM12"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 15) {
          registry.fill(HIST("mass_DCAmumuMM15"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 18) {
          registry.fill(HIST("mass_DCAmumuMM18"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 32) {
          registry.fill(HIST("mass_DCAmumuMM32"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 36) {
          registry.fill(HIST("mass_DCAmumuMM36"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 40) {
          registry.fill(HIST("mass_DCAmumuMM40"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 42) {
          registry.fill(HIST("mass_DCAmumuMM42"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 48) {
          registry.fill(HIST("mass_DCAmumuMM48"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 52) {
          registry.fill(HIST("mass_DCAmumuMM52"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 57) {
          registry.fill(HIST("mass_DCAmumuMM57"), dimuon.mass(), DCAmumu);
        }
      } else if (dimuon.sign() == 2) {
        if (dimuon.mcDecision() == 0) {
          registry.fill(HIST("mass_DCAmumuPP0"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 1) {
          registry.fill(HIST("mass_DCAmumuPP1"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 2) {
          registry.fill(HIST("mass_DCAmumuPP2"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 4) {
          registry.fill(HIST("mass_DCAmumuPP4"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 8) {
          registry.fill(HIST("mass_DCAmumuPP8"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 10) {
          registry.fill(HIST("mass_DCAmumuPP10"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 12) {
          registry.fill(HIST("mass_DCAmumuPP12"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 15) {
          registry.fill(HIST("mass_DCAmumuPP15"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 18) {
          registry.fill(HIST("mass_DCAmumuPP18"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 32) {
          registry.fill(HIST("mass_DCAmumuPP32"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 36) {
          registry.fill(HIST("mass_DCAmumuPP36"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 40) {
          registry.fill(HIST("mass_DCAmumuPP40"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 42) {
          registry.fill(HIST("mass_DCAmumuPP42"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 48) {
          registry.fill(HIST("mass_DCAmumuPP48"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 52) {
          registry.fill(HIST("mass_DCAmumuPP52"), dimuon.mass(), DCAmumu);
        } else if (dimuon.mcDecision() == 57) {
          registry.fill(HIST("mass_DCAmumuPP57"), dimuon.mass(), DCAmumu);
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
