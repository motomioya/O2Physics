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

struct dimuonall_data {

  HistogramRegistry registry{
    "registry", 
    {
      {"dca1PM", "dca1PM", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"dca2PM", "dca2PM", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"pairdcaPM", "pairdcaPM", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"pairdcaambiPM", "pairdcaambiPM", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"svPM", "svPM", {HistType::kTH1F, {{1000, 0., 20}}}},
      {"pcachi2PM", "pcachi2PM", {HistType::kTH1F, {{1000, 0., 6}}}},
      {"lxyzPM", "lxyzPM", {HistType::kTH1F, {{1000, 0, 8}}}},

      {"dca1PP", "dca1PP", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"dca2PP", "dca2PP", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"pairdcaPP", "pairdcaPP", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"svPP", "svPP", {HistType::kTH1F, {{1000, 0., 20}}}},
      {"pcachi2PP", "pcachi2PP", {HistType::kTH1F, {{1000, 0., 6}}}},
      {"lxyzPP", "lxyzPP", {HistType::kTH1F, {{1000, 0, 8}}}},

      {"dca1MM", "dca1MM", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"dca2MM", "dca2MM", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"pairdcaMM", "pairdcaMM", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"svMM", "svMM", {HistType::kTH1F, {{1000, 0., 20}}}},
      {"pcachi2MM", "pcachi2MM", {HistType::kTH1F, {{1000, 0., 6}}}},
      {"lxyzMM", "lxyzMM", {HistType::kTH1F, {{1000, 0, 8}}}},

      {"mass_ptPM", "mass_ptPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM", "mass_ptMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP", "mass_ptPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_dcaxyPM", "mass_dcaxyPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 2}}}},
      {"mass_dcaxyMM", "mass_dcaxyMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 2}}}},
      {"mass_dcaxyPP", "mass_dcaxyPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 2}}}},
      {"mass_pcachiPM", "mass_pcachiPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"mass_pcachiMM", "mass_pcachiMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"mass_pcachiPP", "mass_pcachiPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"dca1_dca2PM", "dca1_dca2PM", {HistType::kTH2F, {{1000, 0.0, 1.0}, {1000, 0.0, 1.0}}}},
      {"dca1_dca2MM", "dca1_dca2MM", {HistType::kTH2F, {{1000, 0.0, 1.0}, {1000, 0.0, 1.0}}}},
      {"dca1_dca2PP", "dca1_dca2PP", {HistType::kTH2F, {{1000, 0.0, 1.0}, {1000, 0.0, 1.0}}}},

      {"mass_ptPM_dcaxylow", "mass_ptPM_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_dcaxylow", "mass_ptMM_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_dcaxylow", "mass_ptPP_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_dcaxyhigh", "mass_ptPM_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_dcaxyhigh", "mass_ptMM_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_dcaxyhigh", "mass_ptPP_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_singledcalow", "mass_ptPM_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_singledcalow", "mass_ptMM_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_singledcalow", "mass_ptPP_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_singledcalowhigh", "mass_ptPM_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_singledcalowhigh", "mass_ptMM_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_singledcalowhigh", "mass_ptPP_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_singledcahigh", "mass_ptPM_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_singledcahigh", "mass_ptMM_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_singledcahigh", "mass_ptPP_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_pcachilow", "mass_ptPM_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_pcachilow", "mass_ptMM_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_pcachilow", "mass_ptPP_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_pcachihigh", "mass_ptPM_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_pcachihigh", "mass_ptMM_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_pcachihigh", "mass_ptPP_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_pcachifail", "mass_ptPM_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_pcachifail", "mass_ptMM_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_pcachifail", "mass_ptPP_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_vectormesoncut", "mass_ptPM_vectormesoncut", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_vectormesoncut", "mass_ptMM_vectormesoncut", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_vectormesoncut", "mass_ptPP_vectormesoncut", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_hfcut", "mass_ptPM_hfcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_hfcut", "mass_ptMM_hfcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_hfcut", "mass_ptPP_hfcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_lxyzlow", "mass_ptPM_lxyzlow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_lxyzlow", "mass_ptMM_lxyzlow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_lxyzlow", "mass_ptPP_lxyzlow", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_lxyzhigh", "mass_ptPM_lxyzhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_lxyzhigh", "mass_ptMM_lxyzhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_lxyzhigh", "mass_ptPP_lxyzhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPM_lxyzfail", "mass_ptPM_lxyzfail", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptMM_lxyzfail", "mass_ptMM_lxyzfail", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_ptPP_lxyzfail", "mass_ptPP_lxyzfail", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},

      {"pcaPM_singledcalowhigh", "pcaPM_singledcalowhigh", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"pcaPM_singledcabothsame", "pcaPM_singledcabothsame", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"svPM_singledcalowhigh", "svPM_singledcalowhigh", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"svPM_singledcabothsame", "svPM_singledcabothsame", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"diffetaPM_singledcalowhigh", "diffetaPM_singledcalowhigh", {HistType::kTH1F, {{1000, -10, 10}}}},
      {"diffetaPM_singledcabothsame", "diffetaPM_singledcabothsame", {HistType::kTH1F, {{1000, -10, 10}}}},
      {"diffphiPM_singledcalowhigh", "diffphiPM_singledcalowhigh", {HistType::kTH1F, {{1000, -10, 10}}}},
      {"diffphiPM_singledcabothsame", "diffphiPM_singledcabothsame", {HistType::kTH1F, {{1000, -10, 10}}}},

      /*
      {"dca1PM_singleptcut", "dca1PM_singleptcut", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"dca2PM_singleptcut", "dca2PM_singleptcut", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"pairdcaPM_singleptcut", "pairdcaPM_singleptcut", {HistType::kTH1F, {{2250, 0.0, 3.0}}}},
      {"svPM_singleptcut", "svPM_singleptcut", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"pcachi2PM_singleptcut", "pcachi2PM_singleptcut", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"lxyzPM_singleptcut", "lxyzPM_singleptcut", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"dca1PP_singleptcut", "dca1PP_singleptcut", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"dca2PP_singleptcut", "dca2PP_singleptcut", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"pairdcaPP_singleptcut", "pairdcaPP_singleptcut", {HistType::kTH1F, {{2250, 0.0, 3.0}}}},
      {"svPP_singleptcut", "svPP_singleptcut", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"pcachi2PP_singleptcut", "pcachi2PP_singleptcut", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"lxyzPP_singleptcut", "lxyzPP_singleptcut", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"dca1MM_singleptcut", "dca1MM_singleptcut", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"dca2MM_singleptcut", "dca2MM_singleptcut", {HistType::kTH1F, {{1000, 0.0, 1.0}}}},
      {"pairdcaMM_singleptcut", "pairdcaMM_singleptcut", {HistType::kTH1F, {{2250, 0.0, 3.0}}}},
      {"svMM_singleptcut", "svMM_singleptcut", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"pcachi2MM_singleptcut", "pcachi2MM_singleptcut", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"lxyzMM_singleptcut", "lxyzMM_singleptcut", {HistType::kTH1F, {{3000, 0, 30}}}},

      {"mass_ptPM_singleptcut", "mass_ptPM_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut", "mass_ptMM_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut", "mass_ptPP_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_dcaxyPM_singleptcut", "mass_dcaxyPM_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {2250, 0.0, 3.0}}}},
      {"mass_dcaxyMM_singleptcut", "mass_dcaxyMM_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {2250, 0.0, 3.0}}}},
      {"mass_dcaxyPP_singleptcut", "mass_dcaxyPP_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {2250, 0.0, 3.0}}}},
      {"mass_pcachiPM_singleptcut", "mass_pcachiPM_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {10000, 0, 100}}}},
      {"mass_pcachiMM_singleptcut", "mass_pcachiMM_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {10000, 0, 100}}}},
      {"mass_pcachiPP_singleptcut", "mass_pcachiPP_singleptcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {10000, 0, 100}}}},
      {"dca1_dca2PM_singleptcut", "dca1_dca2PM_singleptcut", {HistType::kTH2F, {{1000, 0.0, 1.0}, {1000, 0.0, 1.0}}}},
      {"dca1_dca2MM_singleptcut", "dca1_dca2MM_singleptcut", {HistType::kTH2F, {{1000, 0.0, 1.0}, {1000, 0.0, 1.0}}}},
      {"dca1_dca2PP_singleptcut", "dca1_dca2PP_singleptcut", {HistType::kTH2F, {{1000, 0.0, 1.0}, {1000, 0.0, 1.0}}}},

      {"mass_ptPM_singleptcut_dcaxylow", "mass_ptPM_singleptcut_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_dcaxylow", "mass_ptMM_singleptcut_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_dcaxylow", "mass_ptPP_singleptcut_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_dcaxyhigh", "mass_ptPM_singleptcut_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_dcaxyhigh", "mass_ptMM_singleptcut_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_dcaxyhigh", "mass_ptPP_singleptcut_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_singledcalow", "mass_ptPM_singleptcut_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_singledcalow", "mass_ptMM_singleptcut_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_singledcalow", "mass_ptPP_singleptcut_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_singledcalowhigh", "mass_ptPM_singleptcut_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_singledcalowhigh", "mass_ptMM_singleptcut_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_singledcalowhigh", "mass_ptPP_singleptcut_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_singledcahigh", "mass_ptPM_singleptcut_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_singledcahigh", "mass_ptMM_singleptcut_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_singledcahigh", "mass_ptPP_singleptcut_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_pcachilow", "mass_ptPM_singleptcut_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_pcachilow", "mass_ptMM_singleptcut_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_pcachilow", "mass_ptPP_singleptcut_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_pcachihigh", "mass_ptPM_singleptcut_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_pcachihigh", "mass_ptMM_singleptcut_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_pcachihigh", "mass_ptPP_singleptcut_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_pcachifail", "mass_ptPM_singleptcut_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_pcachifail", "mass_ptMM_singleptcut_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_pcachifail", "mass_ptPP_singleptcut_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_vectormesoncut", "mass_ptPM_singleptcut_vectormesoncut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_vectormesoncut", "mass_ptMM_singleptcut_vectormesoncut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_vectormesoncut", "mass_ptPP_singleptcut_vectormesoncut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_hfcut", "mass_ptPM_singleptcut_hfcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptMM_singleptcut_hfcut", "mass_ptMM_singleptcut_hfcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPP_singleptcut_hfcut", "mass_ptPP_singleptcut_hfcut", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30}}}},
      {"mass_ptPM_singleptcut_lxyzlow", "mass_ptPM_singleptcut_lxyzlow", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},
      {"mass_ptMM_singleptcut_lxyzlow", "mass_ptMM_singleptcut_lxyzlow", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},
      {"mass_ptPP_singleptcut_lxyzlow", "mass_ptPP_singleptcut_lxyzlow", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},
      {"mass_ptPM_singleptcut_lxyzhigh", "mass_ptPM_singleptcut_lxyzhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},
      {"mass_ptMM_singleptcut_lxyzhigh", "mass_ptMM_singleptcut_lxyzhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},
      {"mass_ptPP_singleptcut_lxyzhigh", "mass_ptPP_singleptcut_lxyzhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},
      {"mass_ptPM_singleptcut_lxyzfail", "mass_ptPM_singleptcut_lxyzfail", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},
      {"mass_ptMM_singleptcut_lxyzfail", "mass_ptMM_singleptcut_lxyzfail", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},
      {"mass_ptPP_singleptcut_lxyzfail", "mass_ptPP_singleptcut_lxyzfail", {HistType::kTH2F, {{750, 0.0, 15.0}, {730, 0.0, 15.0}}}},

      {"pcaPM_singleptcut_singledcalowhigh", "pcaPM_singleptcut_singledcalowhigh", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"pcaPM_singleptcut_singledcabothsame", "pcaPM_singleptcut_singledcabothsame", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"svPM_singleptcut_singledcalowhigh", "svPM_singleptcut_singledcalowhigh", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"svPM_singleptcut_singledcabothsame", "svPM_singleptcut_singledcabothsame", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"diffetaPM_singleptcut_singledcalowhigh", "diffetaPM_singleptcut_singledcalowhigh", {HistType::kTH1F, {{1000, -10, 10}}}},
      {"diffetaPM_singleptcut_singledcabothsame", "diffetaPM_singleptcut_singledcabothsame", {HistType::kTH1F, {{1000, -10, 10}}}},
      {"diffphiPM_singleptcut_singledcalowhigh", "diffphiPM_singleptcut_singledcalowhigh", {HistType::kTH1F, {{1000, -10, 10}}}},
      {"diffphiPM_singleptcut_singledcabothsame", "diffphiPM_singleptcut_singledcabothsame", {HistType::kTH1F, {{1000, -10, 10}}}},
      */
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
      if (dimuon.chi2pca() == -999) lxyz = -999;
      double pcacut = 0.05;
      double lxyzcut = 0.5;
      double dcacut = 0.2;
      //double singleptcut = 0.5;
      if (dimuon.sVertex() >= -20 || dimuon.sVertex() == -999) {
        if (dimuon.isAmbig1() == 0 && dimuon.isAmbig2() == 0) {
          if (dimuon.chi2MatchMCHMFT1() < 50 && dimuon.chi2MatchMCHMFT2() < 50) {
            if (dimuon.sign() == 0) {
              registry.fill(HIST("dca1PM"), DCA1);
              registry.fill(HIST("dca2PM"), DCA2);
              registry.fill(HIST("pairdcaPM"), DCAmumu);
              registry.fill(HIST("svPM"), dimuon.sVertex());
              registry.fill(HIST("pcachi2PM"), dimuon.chi2pca());
              registry.fill(HIST("lxyzPM"), lxyz);
              registry.fill(HIST("dca1_dca2PM"), DCA1, DCA2);
              registry.fill(HIST("mass_ptPM"), dimuon.mass(), dimuon.pt());
              registry.fill(HIST("mass_dcaxyPM"), dimuon.mass(), DCAmumu);
              registry.fill(HIST("mass_pcachiPM"), dimuon.mass(), dimuon.chi2pca());
              if (DCAmumu <= dcacut) {
                registry.fill(HIST("mass_ptPM_dcaxylow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPM_dcaxyhigh"), dimuon.mass(), dimuon.pt());
              }
              if (dimuon.chi2pca() == -999) {
                registry.fill(HIST("mass_ptPM_pcachifail"), dimuon.mass(), dimuon.pt());
              } else if (dimuon.chi2pca() <= pcacut){
                registry.fill(HIST("mass_ptPM_pcachilow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPM_pcachihigh"), dimuon.mass(), dimuon.pt());
              }
              if (DCA1 <= dcacut && DCA2 <= dcacut) {
                registry.fill(HIST("mass_ptPM_singledcalow"), dimuon.mass(), dimuon.pt());
                registry.fill(HIST("pcaPM_singledcabothsame"), dimuon.chi2pca());
                registry.fill(HIST("svPM_singledcabothsame"), dimuon.sVertex());
                registry.fill(HIST("diffetaPM_singledcabothsame"), dimuon.eta1() - dimuon.eta2());
                registry.fill(HIST("diffphiPM_singledcabothsame"), dimuon.phi1() - dimuon.phi2());
              } else if (DCA1 > dcacut && DCA2 > dcacut) {
                registry.fill(HIST("mass_ptPM_singledcahigh"), dimuon.mass(), dimuon.pt());
                registry.fill(HIST("pcaPM_singledcabothsame"), dimuon.chi2pca());
                registry.fill(HIST("svPM_singledcabothsame"), dimuon.sVertex());
                registry.fill(HIST("diffetaPM_singledcabothsame"), dimuon.eta1() - dimuon.eta2());
                registry.fill(HIST("diffphiPM_singledcabothsame"), dimuon.phi1() - dimuon.phi2());
              } else {
                registry.fill(HIST("mass_ptPM_singledcalowhigh"), dimuon.mass(), dimuon.pt());
                registry.fill(HIST("pcaPM_singledcalowhigh"), dimuon.chi2pca());
                registry.fill(HIST("svPM_singledcalowhigh"), dimuon.sVertex());
                registry.fill(HIST("diffetaPM_singledcalowhigh"), dimuon.eta1() - dimuon.eta2());
                registry.fill(HIST("diffphiPM_singledcalowhigh"), dimuon.phi1() - dimuon.phi2());
              }
              if (lxyz == -999) {
                registry.fill(HIST("mass_ptPM_lxyzfail"), dimuon.mass(), dimuon.pt());
              } else if (lxyz <= lxyzcut) {
                registry.fill(HIST("mass_ptPM_lxyzlow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPM_lxyzhigh"), dimuon.mass(), dimuon.pt());
              }
              if ((dimuon.chi2pca() <= pcacut) && (dimuon.chi2pca() != -999) && (DCAmumu <= dcacut) && (lxyz <= lxyzcut)) {
                registry.fill(HIST("mass_ptPM_vectormesoncut"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPM_hfcut"), dimuon.mass(), dimuon.pt());
              }
            } else if (dimuon.sign() == -2) {
              registry.fill(HIST("dca1MM"), DCA1);
              registry.fill(HIST("dca2MM"), DCA2);
              registry.fill(HIST("pairdcaMM"), DCAmumu);
              registry.fill(HIST("svMM"), dimuon.sVertex());
              registry.fill(HIST("pcachi2MM"), dimuon.chi2pca());
              registry.fill(HIST("lxyzMM"), lxyz);
              registry.fill(HIST("dca1_dca2MM"), DCA1, DCA2);
              registry.fill(HIST("mass_ptMM"), dimuon.mass(), dimuon.pt());
              registry.fill(HIST("mass_dcaxyMM"), dimuon.mass(), DCAmumu);
              registry.fill(HIST("mass_pcachiMM"), dimuon.mass(), dimuon.chi2pca());
              if (DCAmumu <= dcacut) {
                registry.fill(HIST("mass_ptMM_dcaxylow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptMM_dcaxyhigh"), dimuon.mass(), dimuon.pt());
              }
              if (dimuon.chi2pca() == -999) {
                registry.fill(HIST("mass_ptMM_pcachifail"), dimuon.mass(), dimuon.pt());
              } else if (dimuon.chi2pca() <= pcacut){
                registry.fill(HIST("mass_ptMM_pcachilow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptMM_pcachihigh"), dimuon.mass(), dimuon.pt());
              }
              if (DCA1 <= dcacut && DCA2 <= dcacut) {
                registry.fill(HIST("mass_ptMM_singledcalow"), dimuon.mass(), dimuon.pt());
              } else if (DCA1 > dcacut && DCA2 > dcacut) {
                registry.fill(HIST("mass_ptMM_singledcahigh"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptMM_singledcalowhigh"), dimuon.mass(), dimuon.pt());
              }
              if ((dimuon.chi2pca() <= pcacut) && (dimuon.chi2pca() != -999) && (DCAmumu <= dcacut) && (lxyz <= lxyzcut)) {
                registry.fill(HIST("mass_ptMM_vectormesoncut"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptMM_hfcut"), dimuon.mass(), dimuon.pt());
              }
              if (lxyz == -999) {
                registry.fill(HIST("mass_ptMM_lxyzfail"), dimuon.mass(), dimuon.pt());
              } else if (lxyz <= lxyzcut) {
                registry.fill(HIST("mass_ptMM_lxyzlow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptMM_lxyzhigh"), dimuon.mass(), dimuon.pt());
              }
            } else if (dimuon.sign() == 2) {
              registry.fill(HIST("dca1PP"), DCA1);
              registry.fill(HIST("dca2PP"), DCA2);
              registry.fill(HIST("pairdcaPP"), DCAmumu);
              registry.fill(HIST("svPP"), dimuon.sVertex());
              registry.fill(HIST("pcachi2PP"), dimuon.chi2pca());
              registry.fill(HIST("lxyzPP"), lxyz);
              registry.fill(HIST("dca1_dca2PP"), DCA1, DCA2);
              registry.fill(HIST("mass_ptPP"), dimuon.mass(), dimuon.pt());
              registry.fill(HIST("mass_dcaxyPP"), dimuon.mass(), DCAmumu);
              registry.fill(HIST("mass_pcachiPP"), dimuon.mass(), dimuon.chi2pca());
              if (DCAmumu <= dcacut) {
                registry.fill(HIST("mass_ptPP_dcaxylow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPP_dcaxyhigh"), dimuon.mass(), dimuon.pt());
              }
              if (dimuon.chi2pca() == -999) {
                registry.fill(HIST("mass_ptPP_pcachifail"), dimuon.mass(), dimuon.pt());
              } else if (dimuon.chi2pca() <= pcacut){
                registry.fill(HIST("mass_ptPP_pcachilow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPP_pcachihigh"), dimuon.mass(), dimuon.pt());
              }
              if (DCA1 <= dcacut && DCA2 <= dcacut) {
                registry.fill(HIST("mass_ptPP_singledcalow"), dimuon.mass(), dimuon.pt());
              } else if (DCA1 > dcacut && DCA2 > dcacut) {
                registry.fill(HIST("mass_ptPP_singledcahigh"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPP_singledcalowhigh"), dimuon.mass(), dimuon.pt());
              }
              if ((dimuon.chi2pca() <= pcacut) && (dimuon.chi2pca() != -999) && (DCAmumu <= dcacut) && (lxyz <= lxyzcut)) {
                registry.fill(HIST("mass_ptPP_vectormesoncut"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPP_hfcut"), dimuon.mass(), dimuon.pt());
              }
              if (lxyz == -999) {
                registry.fill(HIST("mass_ptPP_lxyzfail"), dimuon.mass(), dimuon.pt());
              } else if (lxyz <= lxyzcut) {
                registry.fill(HIST("mass_ptPP_lxyzlow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPP_lxyzhigh"), dimuon.mass(), dimuon.pt());
              }
            }
          }
        } else {
          if (dimuon.chi2MatchMCHMFT1() < 50 && dimuon.chi2MatchMCHMFT2() < 50) {
            if (dimuon.sign() == 0) {
              registry.fill(HIST("pairdcaambiPM"), DCAmumu);
            }
          }
        }
      }
      /*
      if (dimuon.pt1() >= singleptcut && dimuon.pt2() >= singleptcut) {
        if (dimuon.sVertex() >= -20 || dimuon.sVertex() == -999) {
          if (dimuon.isAmbig1() == 0 && dimuon.isAmbig2() == 0) {
            if (dimuon.chi2MatchMCHMFT1() < 50 && dimuon.chi2MatchMCHMFT2() < 50) {
              if (dimuon.sign() == 0) {
                registry.fill(HIST("dca1PM_singleptcut"), DCA1);
                registry.fill(HIST("dca2PM_singleptcut"), DCA2);
                registry.fill(HIST("pairdcaPM_singleptcut"), DCAmumu);
                registry.fill(HIST("svPM_singleptcut"), dimuon.sVertex());
                registry.fill(HIST("lxyzPM_singleptcut"), lxyz);
                registry.fill(HIST("pcachi2PM_singleptcut"), dimuon.chi2pca());
                registry.fill(HIST("dca1_dca2PM_singleptcut"), DCA1, DCA2);
                registry.fill(HIST("mass_ptPM_singleptcut"), dimuon.mass(), dimuon.pt());
                registry.fill(HIST("mass_dcaxyPM_singleptcut"), dimuon.mass(), DCAmumu);
                registry.fill(HIST("mass_pcachiPM_singleptcut"), dimuon.mass(), dimuon.chi2pca());
                if (DCAmumu <= dcacut) {
                  registry.fill(HIST("mass_ptPM_singleptcut_dcaxylow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPM_singleptcut_dcaxyhigh"), dimuon.mass(), dimuon.pt());
                }
                if (dimuon.chi2pca() == -999) {
                  registry.fill(HIST("mass_ptPM_singleptcut_pcachifail"), dimuon.mass(), dimuon.pt());
                } else if (dimuon.chi2pca() <= pcacut){
                  registry.fill(HIST("mass_ptPM_singleptcut_pcachilow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPM_singleptcut_pcachihigh"), dimuon.mass(), dimuon.pt());
                }
                if (DCA1 <= dcacut && DCA2 <= dcacut) {
                  registry.fill(HIST("mass_ptPM_singleptcut_singledcalow"), dimuon.mass(), dimuon.pt());
                  registry.fill(HIST("pcaPM_singleptcut_singledcabothsame"), dimuon.chi2pca());
                  registry.fill(HIST("svPM_singleptcut_singledcabothsame"), dimuon.sVertex());
                  registry.fill(HIST("diffetaPM_singleptcut_singledcabothsame"), dimuon.eta1() - dimuon.eta2());
                  registry.fill(HIST("diffphiPM_singleptcut_singledcabothsame"), dimuon.phi1() - dimuon.phi2());
                } else if (DCA1 > dcacut && DCA2 > dcacut) {
                  registry.fill(HIST("mass_ptPM_singleptcut_singledcahigh"), dimuon.mass(), dimuon.pt());
                  registry.fill(HIST("pcaPM_singleptcut_singledcabothsame"), dimuon.chi2pca());
                  registry.fill(HIST("svPM_singleptcut_singledcabothsame"), dimuon.sVertex());
                  registry.fill(HIST("diffetaPM_singleptcut_singledcabothsame"), dimuon.eta1() - dimuon.eta2());
                  registry.fill(HIST("diffphiPM_singleptcut_singledcabothsame"), dimuon.phi1() - dimuon.phi2());
                } else {
                  registry.fill(HIST("mass_ptPM_singleptcut_singledcalowhigh"), dimuon.mass(), dimuon.pt());
                  registry.fill(HIST("pcaPM_singleptcut_singledcalowhigh"), dimuon.chi2pca());
                  registry.fill(HIST("svPM_singleptcut_singledcalowhigh"), dimuon.sVertex());
                  registry.fill(HIST("diffetaPM_singleptcut_singledcalowhigh"), dimuon.eta1() - dimuon.eta2());
                  registry.fill(HIST("diffphiPM_singleptcut_singledcalowhigh"), dimuon.phi1() - dimuon.phi2());
                }
                if (lxyz == -999) {
                  registry.fill(HIST("mass_ptPM_singleptcut_lxyzfail"), dimuon.mass(), dimuon.pt());
                } else if (lxyz <= lxyzcut) {
                  registry.fill(HIST("mass_ptPM_singleptcut_lxyzlow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPM_singleptcut_lxyzhigh"), dimuon.mass(), dimuon.pt());
                }
                if ((dimuon.chi2pca() <= pcacut && dimuon.chi2pca() != -999) && ((DCA1 > dcacut && DCA2 > dcacut) || (DCA1 <= dcacut && DCA2 <= dcacut)) && (lxyz <= lxyzcut)) {
                  registry.fill(HIST("mass_ptPM_singleptcut_vectormesoncut"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPM_singleptcut_hfcut"), dimuon.mass(), dimuon.pt());
                }
              } else if (dimuon.sign() == -2) {
                registry.fill(HIST("dca1MM_singleptcut"), DCA1);
                registry.fill(HIST("dca2MM_singleptcut"), DCA2);
                registry.fill(HIST("pairdcaMM_singleptcut"), DCAmumu);
                registry.fill(HIST("svMM_singleptcut"), dimuon.sVertex());
                registry.fill(HIST("pcachi2MM_singleptcut"), dimuon.chi2pca());
                registry.fill(HIST("lxyzMM_singleptcut"), lxyz);
                registry.fill(HIST("dca1_dca2MM_singleptcut"), DCA1, DCA2);
                registry.fill(HIST("mass_ptMM_singleptcut"), dimuon.mass(), dimuon.pt());
                registry.fill(HIST("mass_dcaxyMM_singleptcut"), dimuon.mass(), DCAmumu);
                registry.fill(HIST("mass_pcachiMM_singleptcut"), dimuon.mass(), dimuon.chi2pca());
                if (DCAmumu <= dcacut) {
                  registry.fill(HIST("mass_ptMM_singleptcut_dcaxylow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptMM_singleptcut_dcaxyhigh"), dimuon.mass(), dimuon.pt());
                }
                if (dimuon.chi2pca() == -999) {
                  registry.fill(HIST("mass_ptMM_singleptcut_pcachifail"), dimuon.mass(), dimuon.pt());
                } else if (dimuon.chi2pca() <= pcacut){
                  registry.fill(HIST("mass_ptMM_singleptcut_pcachilow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptMM_singleptcut_pcachihigh"), dimuon.mass(), dimuon.pt());
                }
                if (DCA1 <= dcacut && DCA2 <= dcacut) {
                  registry.fill(HIST("mass_ptMM_singleptcut_singledcalow"), dimuon.mass(), dimuon.pt());
                } else if (DCA1 > dcacut && DCA2 > dcacut) {
                  registry.fill(HIST("mass_ptMM_singleptcut_singledcahigh"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptMM_singleptcut_singledcalowhigh"), dimuon.mass(), dimuon.pt());
                }
                if (lxyz == -999) {
                  registry.fill(HIST("mass_ptMM_singleptcut_lxyzfail"), dimuon.mass(), dimuon.pt());
                } else if (lxyz <= lxyzcut) {
                  registry.fill(HIST("mass_ptMM_singleptcut_lxyzlow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptMM_singleptcut_lxyzhigh"), dimuon.mass(), dimuon.pt());
                }
                if ((dimuon.chi2pca() <= pcacut && dimuon.chi2pca() != -999) && ((DCA1 > dcacut && DCA2 > dcacut) || (DCA1 <= dcacut && DCA2 <= dcacut)) && (lxyz <= lxyzcut)) {
                  registry.fill(HIST("mass_ptMM_singleptcut_vectormesoncut"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptMM_singleptcut_hfcut"), dimuon.mass(), dimuon.pt());
                }
              } else if (dimuon.sign() == 2) {
                registry.fill(HIST("dca1PP_singleptcut"), DCA1);
                registry.fill(HIST("dca2PP_singleptcut"), DCA2);
                registry.fill(HIST("pairdcaPP_singleptcut"), DCAmumu);
                registry.fill(HIST("svPP_singleptcut"), dimuon.sVertex());
                registry.fill(HIST("pcachi2PP_singleptcut"), dimuon.chi2pca());
                registry.fill(HIST("lxyzPP_singleptcut"), lxyz);
                registry.fill(HIST("dca1_dca2PP_singleptcut"), DCA1, DCA2);
                registry.fill(HIST("mass_ptPP_singleptcut"), dimuon.mass(), dimuon.pt());
                registry.fill(HIST("mass_dcaxyPP_singleptcut"), dimuon.mass(), DCAmumu);
                registry.fill(HIST("mass_pcachiPP_singleptcut"), dimuon.mass(), dimuon.chi2pca());
                if (DCAmumu <= dcacut) {
                  registry.fill(HIST("mass_ptPP_singleptcut_dcaxylow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPP_singleptcut_dcaxyhigh"), dimuon.mass(), dimuon.pt());
                }
                if (dimuon.chi2pca() == -999) {
                  registry.fill(HIST("mass_ptPP_singleptcut_pcachifail"), dimuon.mass(), dimuon.pt());
                } else if (dimuon.chi2pca() <= pcacut){
                  registry.fill(HIST("mass_ptPP_singleptcut_pcachilow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPP_singleptcut_pcachihigh"), dimuon.mass(), dimuon.pt());
                }
                if (DCA1 <= dcacut && DCA2 <= dcacut) {
                  registry.fill(HIST("mass_ptPP_singleptcut_singledcalow"), dimuon.mass(), dimuon.pt());
                } else if (DCA1 > dcacut && DCA2 > dcacut) {
                  registry.fill(HIST("mass_ptPP_singleptcut_singledcahigh"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPP_singleptcut_singledcalowhigh"), dimuon.mass(), dimuon.pt());
                }
                if (lxyz == -999) {
                  registry.fill(HIST("mass_ptPP_singleptcut_lxyzfail"), dimuon.mass(), dimuon.pt());
                } else if (lxyz <= lxyzcut) {
                  registry.fill(HIST("mass_ptPP_singleptcut_lxyzlow"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPP_singleptcut_lxyzhigh"), dimuon.mass(), dimuon.pt());
                }
                if ((dimuon.chi2pca() <= pcacut && dimuon.chi2pca() != -999) && ((DCA1 > dcacut && DCA2 > dcacut) || (DCA1 <= dcacut && DCA2 <= dcacut)) && (lxyz <= lxyzcut)) {
                  registry.fill(HIST("mass_ptPP_singleptcut_vectormesoncut"), dimuon.mass(), dimuon.pt());
                } else {
                  registry.fill(HIST("mass_ptPP_singleptcut_hfcut"), dimuon.mass(), dimuon.pt());
                }
              }
            }
          }
        }
      }
    */
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dimuonall_data>(cfgc)
  };
}
