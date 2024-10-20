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

struct myanalysis_dimuonAll_mc_cc {

  HistogramRegistry registry{
    "registry", 
    {
      {"mass_DCAmumuPM0", "mass_DCAmumuPM0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM0", "mass_DCAmumuMM0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP0", "mass_DCAmumuPP0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM0", "mass_DCAmumunoambiPM0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM0", "mass_DCAmumunoambiMM0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP0", "mass_DCAmumunoambiPP0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM0", "mass_LxyzmumuPM0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM0", "mass_LxyzmumuMM0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP0", "mass_LxyzmumuPP0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM0", "mass_LxyzmumunoambiPM0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM0", "mass_LxyzmumunoambiMM0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP0", "mass_LxyzmumunoambiPP0", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM0", "mass_PtmumuPM0", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM0", "mass_PtmumuMM0", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP0", "mass_PtmumuPP0", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM0", "mass_PtmumunoambiPM0", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM0", "mass_PtmumunoambiMM0", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP0", "mass_PtmumunoambiPP0", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM1", "mass_DCAmumuPM1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM1", "mass_DCAmumuMM1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP1", "mass_DCAmumuPP1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM1", "mass_DCAmumunoambiPM1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM1", "mass_DCAmumunoambiMM1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP1", "mass_DCAmumunoambiPP1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM1", "mass_LxyzmumuPM1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM1", "mass_LxyzmumuMM1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP1", "mass_LxyzmumuPP1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM1", "mass_LxyzmumunoambiPM1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM1", "mass_LxyzmumunoambiMM1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP1", "mass_LxyzmumunoambiPP1", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM1", "mass_PtmumuPM1", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM1", "mass_PtmumuMM1", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP1", "mass_PtmumuPP1", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM1", "mass_PtmumunoambiPM1", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM1", "mass_PtmumunoambiMM1", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP1", "mass_PtmumunoambiPP1", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM3", "mass_DCAmumuPM3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM3", "mass_DCAmumuMM3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP3", "mass_DCAmumuPP3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM3", "mass_DCAmumunoambiPM3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM3", "mass_DCAmumunoambiMM3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP3", "mass_DCAmumunoambiPP3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM3", "mass_LxyzmumuPM3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM3", "mass_LxyzmumuMM3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP3", "mass_LxyzmumuPP3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM3", "mass_LxyzmumunoambiPM3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM3", "mass_LxyzmumunoambiMM3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP3", "mass_LxyzmumunoambiPP3", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM3", "mass_PtmumuPM3", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM3", "mass_PtmumuMM3", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP3", "mass_PtmumuPP3", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM3", "mass_PtmumunoambiPM3", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM3", "mass_PtmumunoambiMM3", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP3", "mass_PtmumunoambiPP3", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM5", "mass_DCAmumuPM5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM5", "mass_DCAmumuMM5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP5", "mass_DCAmumuPP5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM5", "mass_DCAmumunoambiPM5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM5", "mass_DCAmumunoambiMM5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP5", "mass_DCAmumunoambiPP5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM5", "mass_LxyzmumuPM5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM5", "mass_LxyzmumuMM5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP5", "mass_LxyzmumuPP5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM5", "mass_LxyzmumunoambiPM5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM5", "mass_LxyzmumunoambiMM5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP5", "mass_LxyzmumunoambiPP5", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM5", "mass_PtmumuPM5", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM5", "mass_PtmumuMM5", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP5", "mass_PtmumuPP5", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM5", "mass_PtmumunoambiPM5", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM5", "mass_PtmumunoambiMM5", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP5", "mass_PtmumunoambiPP5", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM9", "mass_DCAmumuPM9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM9", "mass_DCAmumuMM9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP9", "mass_DCAmumuPP9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM9", "mass_DCAmumunoambiPM9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM9", "mass_DCAmumunoambiMM9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP9", "mass_DCAmumunoambiPP9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM9", "mass_LxyzmumuPM9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM9", "mass_LxyzmumuMM9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP9", "mass_LxyzmumuPP9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM9", "mass_LxyzmumunoambiPM9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM9", "mass_LxyzmumunoambiMM9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP9", "mass_LxyzmumunoambiPP9", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM9", "mass_PtmumuPM9", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM9", "mass_PtmumuMM9", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP9", "mass_PtmumuPP9", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM9", "mass_PtmumunoambiPM9", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM9", "mass_PtmumunoambiMM9", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP9", "mass_PtmumunoambiPP9", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM17", "mass_DCAmumuPM17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM17", "mass_DCAmumuMM17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP17", "mass_DCAmumuPP17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM17", "mass_DCAmumunoambiPM17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM17", "mass_DCAmumunoambiMM17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP17", "mass_DCAmumunoambiPP17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM17", "mass_LxyzmumuPM17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM17", "mass_LxyzmumuMM17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP17", "mass_LxyzmumuPP17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM17", "mass_LxyzmumunoambiPM17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM17", "mass_LxyzmumunoambiMM17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP17", "mass_LxyzmumunoambiPP17", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM17", "mass_PtmumuPM17", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM17", "mass_PtmumuMM17", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP17", "mass_PtmumuPP17", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM17", "mass_PtmumunoambiPM17", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM17", "mass_PtmumunoambiMM17", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP17", "mass_PtmumunoambiPP17", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM33", "mass_DCAmumuPM33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM33", "mass_DCAmumuMM33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP33", "mass_DCAmumuPP33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM33", "mass_DCAmumunoambiPM33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM33", "mass_DCAmumunoambiMM33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP33", "mass_DCAmumunoambiPP33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM33", "mass_LxyzmumuPM33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM33", "mass_LxyzmumuMM33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP33", "mass_LxyzmumuPP33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM33", "mass_LxyzmumunoambiPM33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM33", "mass_LxyzmumunoambiMM33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP33", "mass_LxyzmumunoambiPP33", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM33", "mass_PtmumuPM33", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM33", "mass_PtmumuMM33", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP33", "mass_PtmumuPP33", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM33", "mass_PtmumunoambiPM33", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM33", "mass_PtmumunoambiMM33", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP33", "mass_PtmumunoambiPP33", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM65", "mass_DCAmumuPM65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM65", "mass_DCAmumuMM65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP65", "mass_DCAmumuPP65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM65", "mass_DCAmumunoambiPM65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM65", "mass_DCAmumunoambiMM65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP65", "mass_DCAmumunoambiPP65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM65", "mass_LxyzmumuPM65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM65", "mass_LxyzmumuMM65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP65", "mass_LxyzmumuPP65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM65", "mass_LxyzmumunoambiPM65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM65", "mass_LxyzmumunoambiMM65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP65", "mass_LxyzmumunoambiPP65", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM65", "mass_PtmumuPM65", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM65", "mass_PtmumuMM65", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP65", "mass_PtmumuPP65", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM65", "mass_PtmumunoambiPM65", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM65", "mass_PtmumunoambiMM65", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP65", "mass_PtmumunoambiPP65", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM129", "mass_DCAmumuPM129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM129", "mass_DCAmumuMM129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP129", "mass_DCAmumuPP129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM129", "mass_DCAmumunoambiPM129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM129", "mass_DCAmumunoambiMM129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP129", "mass_DCAmumunoambiPP129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM129", "mass_LxyzmumuPM129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM129", "mass_LxyzmumuMM129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP129", "mass_LxyzmumuPP129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM129", "mass_LxyzmumunoambiPM129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM129", "mass_LxyzmumunoambiMM129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP129", "mass_LxyzmumunoambiPP129", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM129", "mass_PtmumuPM129", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM129", "mass_PtmumuMM129", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP129", "mass_PtmumuPP129", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM129", "mass_PtmumunoambiPM129", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM129", "mass_PtmumunoambiMM129", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP129", "mass_PtmumunoambiPP129", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM257", "mass_DCAmumuPM257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM257", "mass_DCAmumuMM257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP257", "mass_DCAmumuPP257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM257", "mass_DCAmumunoambiPM257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM257", "mass_DCAmumunoambiMM257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP257", "mass_DCAmumunoambiPP257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM257", "mass_LxyzmumuPM257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM257", "mass_LxyzmumuMM257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP257", "mass_LxyzmumuPP257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM257", "mass_LxyzmumunoambiPM257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM257", "mass_LxyzmumunoambiMM257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP257", "mass_LxyzmumunoambiPP257", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM257", "mass_PtmumuPM257", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM257", "mass_PtmumuMM257", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP257", "mass_PtmumuPP257", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM257", "mass_PtmumunoambiPM257", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM257", "mass_PtmumunoambiMM257", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP257", "mass_PtmumunoambiPP257", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM513", "mass_DCAmumuPM513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM513", "mass_DCAmumuMM513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP513", "mass_DCAmumuPP513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM513", "mass_DCAmumunoambiPM513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM513", "mass_DCAmumunoambiMM513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP513", "mass_DCAmumunoambiPP513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM513", "mass_LxyzmumuPM513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM513", "mass_LxyzmumuMM513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP513", "mass_LxyzmumuPP513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM513", "mass_LxyzmumunoambiPM513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM513", "mass_LxyzmumunoambiMM513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP513", "mass_LxyzmumunoambiPP513", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM513", "mass_PtmumuPM513", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM513", "mass_PtmumuMM513", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP513", "mass_PtmumuPP513", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM513", "mass_PtmumunoambiPM513", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM513", "mass_PtmumunoambiMM513", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP513", "mass_PtmumunoambiPP513", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM1025", "mass_DCAmumuPM1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM1025", "mass_DCAmumuMM1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP1025", "mass_DCAmumuPP1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM1025", "mass_DCAmumunoambiPM1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM1025", "mass_DCAmumunoambiMM1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP1025", "mass_DCAmumunoambiPP1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM1025", "mass_LxyzmumuPM1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM1025", "mass_LxyzmumuMM1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP1025", "mass_LxyzmumuPP1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM1025", "mass_LxyzmumunoambiPM1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM1025", "mass_LxyzmumunoambiMM1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP1025", "mass_LxyzmumunoambiPP1025", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM1025", "mass_PtmumuPM1025", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM1025", "mass_PtmumuMM1025", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP1025", "mass_PtmumuPP1025", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM1025", "mass_PtmumunoambiPM1025", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM1025", "mass_PtmumunoambiMM1025", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP1025", "mass_PtmumunoambiPP1025", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"SVFailurePM", "SVFailurePM", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"SVFailureMM", "SVFailureMM", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"SVFailurePP", "SVFailurePP", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"SVFailurenoambiPM", "SVFailurenoambiPM", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"SVFailurenoambiMM", "SVFailurenoambiMM", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"SVFailurenoambiPP", "SVFailurenoambiPP", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
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
      float lxy = dimuon.tauxy() * dimuon.pt() * o2::constants::physics::LightSpeedCm2NS/dimuon.mass();
      float lz = std::sqrt((dimuon.posZ() - dimuon.sVertex()) * (dimuon.posZ() - dimuon.sVertex()));
      float lxyz = std::sqrt(lxy * lxy + lz * lz);
      if (dimuon.eta1() >= -3.6 && dimuon.eta1() <= -2.5 && dimuon.eta2() >= -3.6 && dimuon.eta2() <= -2.5) {
        if (dimuon.sign() == 0) {
          if (dimuon.mcDecision() == 0) {
            registry.fill(HIST("mass_DCAmumuPM0"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM0"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM0"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1) {
            registry.fill(HIST("mass_DCAmumuPM1"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM1"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM1"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 3) {
            registry.fill(HIST("mass_DCAmumuPM3"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM3"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM3"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 5) {
            registry.fill(HIST("mass_DCAmumuPM5"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM5"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM5"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 9) {
            registry.fill(HIST("mass_DCAmumuPM9"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM9"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM9"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 17) {
            registry.fill(HIST("mass_DCAmumuPM17"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM17"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM17"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 33) {
            registry.fill(HIST("mass_DCAmumuPM33"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM33"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM33"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 65) {
            registry.fill(HIST("mass_DCAmumuPM65"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM65"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM65"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 129) {
            registry.fill(HIST("mass_DCAmumuPM129"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM129"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM129"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 257) {
            registry.fill(HIST("mass_DCAmumuPM257"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM257"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM257"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 513) {
            registry.fill(HIST("mass_DCAmumuPM513"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM513"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM513"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1025) {
            registry.fill(HIST("mass_DCAmumuPM1025"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM1025"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM1025"), dimuon.mass(), dimuon.pt());
          }
          if (lxyz == -999) {
            registry.fill(HIST("SVFailurePM"), 1);
          } else {
            registry.fill(HIST("SVFailurePM"), 0);
          }
        } else if (dimuon.sign() == -2) {
          if (dimuon.mcDecision() == 0) {
            registry.fill(HIST("mass_DCAmumuMM0"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM0"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM0"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1) {
            registry.fill(HIST("mass_DCAmumuMM1"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM1"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM1"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 3) {
            registry.fill(HIST("mass_DCAmumuMM3"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM3"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM3"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 5) {
            registry.fill(HIST("mass_DCAmumuMM5"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM5"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM5"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 9) {
            registry.fill(HIST("mass_DCAmumuMM9"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM9"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM9"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 17) {
            registry.fill(HIST("mass_DCAmumuMM17"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM17"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM17"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 33) {
            registry.fill(HIST("mass_DCAmumuMM33"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM33"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM33"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 65) {
            registry.fill(HIST("mass_DCAmumuMM65"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM65"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM65"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 129) {
            registry.fill(HIST("mass_DCAmumuMM129"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM129"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM129"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 257) {
            registry.fill(HIST("mass_DCAmumuMM257"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM257"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM257"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 513) {
            registry.fill(HIST("mass_DCAmumuMM513"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM513"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM513"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1025) {
            registry.fill(HIST("mass_DCAmumuMM1025"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM1025"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM1025"), dimuon.mass(), dimuon.pt());
          }
          if (lxyz == -999) {
            registry.fill(HIST("SVFailureMM"), 1);
          } else {
            registry.fill(HIST("SVFailureMM"), 0);
          }
        } else if (dimuon.sign() == 2) {
          if (dimuon.mcDecision() == 0) {
            registry.fill(HIST("mass_DCAmumuPP0"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP0"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP0"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1) {
            registry.fill(HIST("mass_DCAmumuPP1"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP1"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP1"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 3) {
            registry.fill(HIST("mass_DCAmumuPP3"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP3"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP3"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 5) {
            registry.fill(HIST("mass_DCAmumuPP5"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP5"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP5"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 9) {
            registry.fill(HIST("mass_DCAmumuPP9"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP9"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP9"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 17) {
            registry.fill(HIST("mass_DCAmumuPP17"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP17"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP17"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 33) {
            registry.fill(HIST("mass_DCAmumuPP33"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP33"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP33"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 65) {
            registry.fill(HIST("mass_DCAmumuPP65"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP65"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP65"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 129) {
            registry.fill(HIST("mass_DCAmumuPP129"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP129"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP129"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 257) {
            registry.fill(HIST("mass_DCAmumuPP257"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP257"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP257"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 513) {
            registry.fill(HIST("mass_DCAmumuPP513"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP513"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP513"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1025) {
            registry.fill(HIST("mass_DCAmumuPP1025"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP1025"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP1025"), dimuon.mass(), dimuon.pt());
          }
          if (lxyz == -999) {
            registry.fill(HIST("SVFailurePP"), 1);
          } else {
            registry.fill(HIST("SVFailurePP"), 0);
          }
        }
        //ambiguous tracks
        if (dimuon.isAmbig1() == 0 && dimuon.isAmbig2() == 0) {
          if (dimuon.sign() == 0) {
          if (dimuon.mcDecision() == 0) {
            registry.fill(HIST("mass_DCAmumunoambiPM0"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM0"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM0"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1) {
            registry.fill(HIST("mass_DCAmumunoambiPM1"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM1"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM1"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 3) {
            registry.fill(HIST("mass_DCAmumunoambiPM3"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM3"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM3"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 5) {
            registry.fill(HIST("mass_DCAmumunoambiPM5"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM5"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM5"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 9) {
            registry.fill(HIST("mass_DCAmumunoambiPM9"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM9"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM9"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 17) {
            registry.fill(HIST("mass_DCAmumunoambiPM17"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM17"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM17"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 33) {
            registry.fill(HIST("mass_DCAmumunoambiPM33"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM33"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM33"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 65) {
            registry.fill(HIST("mass_DCAmumunoambiPM65"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM65"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM65"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 129) {
            registry.fill(HIST("mass_DCAmumunoambiPM129"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM129"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM129"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 257) {
            registry.fill(HIST("mass_DCAmumunoambiPM257"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM257"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM257"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 513) {
            registry.fill(HIST("mass_DCAmumunoambiPM513"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM513"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM513"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1025) {
            registry.fill(HIST("mass_DCAmumunoambiPM1025"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM1025"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM1025"), dimuon.mass(), dimuon.pt());
          }
            if (lxyz == -999) {
              registry.fill(HIST("SVFailurenoambiPM"), 1);
            } else {
              registry.fill(HIST("SVFailurenoambiPM"), 0);
            }
          } else if (dimuon.sign() == -2) {
          if (dimuon.mcDecision() == 0) {
            registry.fill(HIST("mass_DCAmumunoambiMM0"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM0"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM0"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1) {
            registry.fill(HIST("mass_DCAmumunoambiMM1"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM1"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM1"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 3) {
            registry.fill(HIST("mass_DCAmumunoambiMM3"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM3"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM3"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 5) {
            registry.fill(HIST("mass_DCAmumunoambiMM5"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM5"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM5"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 9) {
            registry.fill(HIST("mass_DCAmumunoambiMM9"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM9"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM9"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 17) {
            registry.fill(HIST("mass_DCAmumunoambiMM17"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM17"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM17"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 33) {
            registry.fill(HIST("mass_DCAmumunoambiMM33"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM33"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM33"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 65) {
            registry.fill(HIST("mass_DCAmumunoambiMM65"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM65"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM65"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 129) {
            registry.fill(HIST("mass_DCAmumunoambiMM129"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM129"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM129"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 257) {
            registry.fill(HIST("mass_DCAmumunoambiMM257"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM257"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM257"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 513) {
            registry.fill(HIST("mass_DCAmumunoambiMM513"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM513"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM513"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1025) {
            registry.fill(HIST("mass_DCAmumunoambiMM1025"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM1025"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM1025"), dimuon.mass(), dimuon.pt());
          }
            if (lxyz == -999) {
              registry.fill(HIST("SVFailurenoambiMM"), 1);
            } else {
              registry.fill(HIST("SVFailurenoambiMM"), 0);
            }
          } else if (dimuon.sign() == 2) {
          if (dimuon.mcDecision() == 0) {
            registry.fill(HIST("mass_DCAmumunoambiPP0"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP0"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP0"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1) {
            registry.fill(HIST("mass_DCAmumunoambiPP1"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP1"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP1"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 3) {
            registry.fill(HIST("mass_DCAmumunoambiPP3"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP3"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP3"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 5) {
            registry.fill(HIST("mass_DCAmumunoambiPP5"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP5"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP5"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 9) {
            registry.fill(HIST("mass_DCAmumunoambiPP9"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP9"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP9"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 17) {
            registry.fill(HIST("mass_DCAmumunoambiPP17"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP17"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP17"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 33) {
            registry.fill(HIST("mass_DCAmumunoambiPP33"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP33"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP33"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 65) {
            registry.fill(HIST("mass_DCAmumunoambiPP65"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP65"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP65"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 129) {
            registry.fill(HIST("mass_DCAmumunoambiPP129"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP129"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP129"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 257) {
            registry.fill(HIST("mass_DCAmumunoambiPP257"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP257"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP257"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 513) {
            registry.fill(HIST("mass_DCAmumunoambiPP513"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP513"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP513"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 1025) {
            registry.fill(HIST("mass_DCAmumunoambiPP1025"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP1025"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP1025"), dimuon.mass(), dimuon.pt());
          }
            if (lxyz == -999) {
              registry.fill(HIST("SVFailurenoambiPP"), 1);
            } else {
              registry.fill(HIST("SVFailurenoambiPP"), 0);
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
    adaptAnalysisTask<myanalysis_dimuonAll_mc_cc>(cfgc)
  };
}
