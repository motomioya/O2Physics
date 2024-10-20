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

struct myanalysis_dimuonAll_mc {

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
      {"mass_DCAmumuPM2", "mass_DCAmumuPM2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM2", "mass_DCAmumuMM2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP2", "mass_DCAmumuPP2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM2", "mass_DCAmumunoambiPM2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM2", "mass_DCAmumunoambiMM2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP2", "mass_DCAmumunoambiPP2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM2", "mass_LxyzmumuPM2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM2", "mass_LxyzmumuMM2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP2", "mass_LxyzmumuPP2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM2", "mass_LxyzmumunoambiPM2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM2", "mass_LxyzmumunoambiMM2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP2", "mass_LxyzmumunoambiPP2", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM2", "mass_PtmumuPM2", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM2", "mass_PtmumuMM2", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP2", "mass_PtmumuPP2", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM2", "mass_PtmumunoambiPM2", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM2", "mass_PtmumunoambiMM2", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP2", "mass_PtmumunoambiPP2", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM4", "mass_DCAmumuPM4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM4", "mass_DCAmumuMM4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP4", "mass_DCAmumuPP4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM4", "mass_DCAmumunoambiPM4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM4", "mass_DCAmumunoambiMM4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP4", "mass_DCAmumunoambiPP4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM4", "mass_LxyzmumuPM4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM4", "mass_LxyzmumuMM4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP4", "mass_LxyzmumuPP4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM4", "mass_LxyzmumunoambiPM4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM4", "mass_LxyzmumunoambiMM4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP4", "mass_LxyzmumunoambiPP4", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM4", "mass_PtmumuPM4", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM4", "mass_PtmumuMM4", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP4", "mass_PtmumuPP4", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM4", "mass_PtmumunoambiPM4", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM4", "mass_PtmumunoambiMM4", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP4", "mass_PtmumunoambiPP4", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM8", "mass_DCAmumuPM8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM8", "mass_DCAmumuMM8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP8", "mass_DCAmumuPP8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM8", "mass_DCAmumunoambiPM8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM8", "mass_DCAmumunoambiMM8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP8", "mass_DCAmumunoambiPP8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM8", "mass_LxyzmumuPM8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM8", "mass_LxyzmumuMM8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP8", "mass_LxyzmumuPP8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM8", "mass_LxyzmumunoambiPM8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM8", "mass_LxyzmumunoambiMM8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP8", "mass_LxyzmumunoambiPP8", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM8", "mass_PtmumuPM8", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM8", "mass_PtmumuMM8", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP8", "mass_PtmumuPP8", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM8", "mass_PtmumunoambiPM8", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM8", "mass_PtmumunoambiMM8", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP8", "mass_PtmumunoambiPP8", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM20", "mass_DCAmumuPM20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM20", "mass_DCAmumuMM20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP20", "mass_DCAmumuPP20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM20", "mass_DCAmumunoambiPM20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM20", "mass_DCAmumunoambiMM20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP20", "mass_DCAmumunoambiPP20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM20", "mass_LxyzmumuPM20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM20", "mass_LxyzmumuMM20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP20", "mass_LxyzmumuPP20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM20", "mass_LxyzmumunoambiPM20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM20", "mass_LxyzmumunoambiMM20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP20", "mass_LxyzmumunoambiPP20", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM20", "mass_PtmumuPM20", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM20", "mass_PtmumuMM20", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP20", "mass_PtmumuPP20", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM20", "mass_PtmumunoambiPM20", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM20", "mass_PtmumunoambiMM20", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP20", "mass_PtmumunoambiPP20", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM40", "mass_DCAmumuPM40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM40", "mass_DCAmumuMM40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP40", "mass_DCAmumuPP40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM40", "mass_DCAmumunoambiPM40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM40", "mass_DCAmumunoambiMM40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP40", "mass_DCAmumunoambiPP40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM40", "mass_LxyzmumuPM40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM40", "mass_LxyzmumuMM40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP40", "mass_LxyzmumuPP40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM40", "mass_LxyzmumunoambiPM40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM40", "mass_LxyzmumunoambiMM40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP40", "mass_LxyzmumunoambiPP40", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM40", "mass_PtmumuPM40", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM40", "mass_PtmumuMM40", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP40", "mass_PtmumuPP40", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM40", "mass_PtmumunoambiPM40", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM40", "mass_PtmumunoambiMM40", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP40", "mass_PtmumunoambiPP40", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM68", "mass_DCAmumuPM68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM68", "mass_DCAmumuMM68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP68", "mass_DCAmumuPP68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM68", "mass_DCAmumunoambiPM68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM68", "mass_DCAmumunoambiMM68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP68", "mass_DCAmumunoambiPP68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM68", "mass_LxyzmumuPM68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM68", "mass_LxyzmumuMM68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP68", "mass_LxyzmumuPP68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM68", "mass_LxyzmumunoambiPM68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM68", "mass_LxyzmumunoambiMM68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP68", "mass_LxyzmumunoambiPP68", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM68", "mass_PtmumuPM68", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM68", "mass_PtmumuMM68", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP68", "mass_PtmumuPP68", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM68", "mass_PtmumunoambiPM68", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM68", "mass_PtmumunoambiMM68", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP68", "mass_PtmumunoambiPP68", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_DCAmumuPM136", "mass_DCAmumuPM136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuMM136", "mass_DCAmumuMM136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumuPP136", "mass_DCAmumuPP136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPM136", "mass_DCAmumunoambiPM136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiMM136", "mass_DCAmumunoambiMM136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_DCAmumunoambiPP136", "mass_DCAmumunoambiPP136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 3}}}},
      {"mass_LxyzmumuPM136", "mass_LxyzmumuPM136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuMM136", "mass_LxyzmumuMM136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumuPP136", "mass_LxyzmumuPP136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPM136", "mass_LxyzmumunoambiPM136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiMM136", "mass_LxyzmumunoambiMM136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_LxyzmumunoambiPP136", "mass_LxyzmumunoambiPP136", {HistType::kTH2F, {{750, 0, 15}, {1000, 0, 16}}}},
      {"mass_PtmumuPM136", "mass_PtmumuPM136", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuMM136", "mass_PtmumuMM136", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumuPP136", "mass_PtmumuPP136", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPM136", "mass_PtmumunoambiPM136", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiMM136", "mass_PtmumunoambiMM136", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
      {"mass_PtmumunoambiPP136", "mass_PtmumunoambiPP136", {HistType::kTH2F, {{750, 0, 15}, {750, 0, 15}}}},
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
          } else if (dimuon.mcDecision() == 2) {
            registry.fill(HIST("mass_DCAmumuPM2"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM2"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM2"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 4) {
            registry.fill(HIST("mass_DCAmumuPM4"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM4"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM4"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 8) {
            registry.fill(HIST("mass_DCAmumuPM8"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM8"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM8"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 20) {
            registry.fill(HIST("mass_DCAmumuPM20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM20"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM20"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 40) {
            registry.fill(HIST("mass_DCAmumuPM40"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM40"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM40"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 68) {
            registry.fill(HIST("mass_DCAmumuPM68"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM68"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM68"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 136) {
            registry.fill(HIST("mass_DCAmumuPM136"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPM136"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPM136"), dimuon.mass(), dimuon.pt());
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
          } else if (dimuon.mcDecision() == 2) {
            registry.fill(HIST("mass_DCAmumuMM2"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM2"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM2"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 4) {
            registry.fill(HIST("mass_DCAmumuMM4"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM4"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM4"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 8) {
            registry.fill(HIST("mass_DCAmumuMM8"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM8"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM8"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 20) {
            registry.fill(HIST("mass_DCAmumuMM20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM20"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM20"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 40) {
            registry.fill(HIST("mass_DCAmumuMM40"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM40"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM40"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 68) {
            registry.fill(HIST("mass_DCAmumuMM68"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM68"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM68"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 136) {
            registry.fill(HIST("mass_DCAmumuMM136"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuMM136"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuMM136"), dimuon.mass(), dimuon.pt());
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
          } else if (dimuon.mcDecision() == 2) {
            registry.fill(HIST("mass_DCAmumuPP2"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP2"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP2"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 4) {
            registry.fill(HIST("mass_DCAmumuPP4"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP4"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP4"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 8) {
            registry.fill(HIST("mass_DCAmumuPP8"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP8"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP8"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 20) {
            registry.fill(HIST("mass_DCAmumuPP20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP20"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP20"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 40) {
            registry.fill(HIST("mass_DCAmumuPP40"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP40"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP40"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 68) {
            registry.fill(HIST("mass_DCAmumuPP68"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP68"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP68"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 136) {
            registry.fill(HIST("mass_DCAmumuPP136"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumuPP136"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumuPP136"), dimuon.mass(), dimuon.pt());
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
          } else if (dimuon.mcDecision() == 2) {
            registry.fill(HIST("mass_DCAmumunoambiPM2"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM2"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM2"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 4) {
            registry.fill(HIST("mass_DCAmumunoambiPM4"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM4"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM4"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 8) {
            registry.fill(HIST("mass_DCAmumunoambiPM8"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM8"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM8"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 20) {
            registry.fill(HIST("mass_DCAmumunoambiPM20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM20"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM20"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 40) {
            registry.fill(HIST("mass_DCAmumunoambiPM40"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM40"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM40"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 68) {
            registry.fill(HIST("mass_DCAmumunoambiPM68"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM68"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM68"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 136) {
            registry.fill(HIST("mass_DCAmumunoambiPM136"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPM136"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPM136"), dimuon.mass(), dimuon.pt());
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
          } else if (dimuon.mcDecision() == 2) {
            registry.fill(HIST("mass_DCAmumunoambiMM2"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM2"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM2"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 4) {
            registry.fill(HIST("mass_DCAmumunoambiMM4"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM4"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM4"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 8) {
            registry.fill(HIST("mass_DCAmumunoambiMM8"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM8"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM8"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 20) {
            registry.fill(HIST("mass_DCAmumunoambiMM20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM20"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM20"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 40) {
            registry.fill(HIST("mass_DCAmumunoambiMM40"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM40"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM40"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 68) {
            registry.fill(HIST("mass_DCAmumunoambiMM68"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM68"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM68"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 136) {
            registry.fill(HIST("mass_DCAmumunoambiMM136"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiMM136"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiMM136"), dimuon.mass(), dimuon.pt());
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
          } else if (dimuon.mcDecision() == 2) {
            registry.fill(HIST("mass_DCAmumunoambiPP2"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP2"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP2"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 4) {
            registry.fill(HIST("mass_DCAmumunoambiPP4"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP4"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP4"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 8) {
            registry.fill(HIST("mass_DCAmumunoambiPP8"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP8"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP8"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 20) {
            registry.fill(HIST("mass_DCAmumunoambiPP20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP20"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP20"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 40) {
            registry.fill(HIST("mass_DCAmumunoambiPP40"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP40"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP40"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 68) {
            registry.fill(HIST("mass_DCAmumunoambiPP68"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP68"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP68"), dimuon.mass(), dimuon.pt());
          } else if (dimuon.mcDecision() == 136) {
            registry.fill(HIST("mass_DCAmumunoambiPP136"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzmumunoambiPP136"), dimuon.mass(), lxyz);
            registry.fill(HIST("mass_PtmumunoambiPP136"), dimuon.mass(), dimuon.pt());
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
    adaptAnalysisTask<myanalysis_dimuonAll_mc>(cfgc)
  };
}
