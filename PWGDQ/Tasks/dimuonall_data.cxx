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
      {"mass_ptPM", "mass_ptPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM", "mass_ptMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP", "mass_ptPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_dcaxyPM", "mass_dcaxyPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0, 15}}}},
      {"mass_dcaxyMM", "mass_dcaxyMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0, 15}}}},
      {"mass_dcaxyPP", "mass_dcaxyPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0, 15}}}},
      {"mass_pcachiPM", "mass_pcachiPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {11000, -1000, 100}}}},
      {"mass_pcachiMM", "mass_pcachiMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {11000, -1000, 100}}}},
      {"mass_pcachiPP", "mass_pcachiPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {11000, -1000, 100}}}},

      {"mass_ptPM_dcaxylow", "mass_ptPM_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM_dcaxylow", "mass_ptMM_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP_dcaxylow", "mass_ptPP_dcaxylow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPM_dcaxyhigh", "mass_ptPM_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM_dcaxyhigh", "mass_ptMM_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP_dcaxyhigh", "mass_ptPP_dcaxyhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPM_singledcalow", "mass_ptPM_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM_singledcalow", "mass_ptMM_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP_singledcalow", "mass_ptPP_singledcalow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPM_singledcalowhigh", "mass_ptPM_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM_singledcalowhigh", "mass_ptMM_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP_singledcalowhigh", "mass_ptPP_singledcalowhigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPM_singledcahigh", "mass_ptPM_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM_singledcahigh", "mass_ptMM_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP_singledcahigh", "mass_ptPP_singledcahigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPM_pcachilow", "mass_ptPM_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM_pcachilow", "mass_ptMM_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP_pcachilow", "mass_ptPP_pcachilow", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPM_pcachihigh", "mass_ptPM_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM_pcachihigh", "mass_ptMM_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP_pcachihigh", "mass_ptPP_pcachihigh", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPM_pcachifail", "mass_ptPM_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptMM_pcachifail", "mass_ptMM_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"mass_ptPP_pcachifail", "mass_ptPP_pcachifail", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    LOGF(info, "----init----");
  }

  void process(aod::DimuonsAll const& dimuons)
  {
    LOGF(info, "----process----");
    for (auto& dimuon : dimuons) {
      LOGF(info, "----dimuon loop----");
      float DCA1 = std::sqrt(dimuon.fwdDcaX1() * dimuon.fwdDcaX1() + dimuon.fwdDcaY1() * dimuon.fwdDcaY1());
      float DCA2 = std::sqrt(dimuon.fwdDcaX2() * dimuon.fwdDcaX2() + dimuon.fwdDcaY2() * dimuon.fwdDcaY2());
      float DCAmumu = std::sqrt((DCA1 * DCA1 + DCA2 * DCA2)/2);
      LOGF(info, "before svertex if statement");
      if (dimuon.sVertex() >= -20 || dimuon.sVertex() == -999) {
        LOGF(info, "before siAmbigu if statement");
        if (dimuon.isAmbig1() == 0 && dimuon.isAmbig2() == 0) {
          LOGF(info, "before chi2 if statement");
          if (dimuon.chi2MatchMCHMFT1() < 50 && dimuon.chi2MatchMCHMFT2() < 50) {
            LOGF(info, "before sign if statement");
            if (dimuon.sign() == 0) {
              LOGF(info, "inside if sign() == 0");
              registry.fill(HIST("mass_ptPM"), dimuon.mass(), dimuon.pt());
              registry.fill(HIST("mass_dcaxyPM"), dimuon.mass(), DCAmumu);
              registry.fill(HIST("mass_pcachiPM"), dimuon.mass(), dimuon.chi2pca());
              if (DCAmumu <= 0.02) {
                registry.fill(HIST("mass_ptPM_dcaxylow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPM_dcaxyhigh"), dimuon.mass(), dimuon.pt());
              }
              if (dimuon.chi2pca() == -999) {
                registry.fill(HIST("mass_ptPM_pcachifail"), dimuon.mass(), dimuon.pt());
              } else if (dimuon.chi2pca() <= 2){
                registry.fill(HIST("mass_ptPM_pcachilow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPM_pcachihigh"), dimuon.mass(), dimuon.pt());
              }
              if (DCA1 <= 0.02 && DCA2 <= 0.02) {
                registry.fill(HIST("mass_ptPM_singledcalow"), dimuon.mass(), dimuon.pt());
              } else if (DCA1 > 0.02 && DCA2 > 0.02) {
                registry.fill(HIST("mass_ptPM_singledcahigh"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPM_singledcalowhigh"), dimuon.mass(), dimuon.pt());
              }
            } else if (dimuon.sign() == -2) {
              registry.fill(HIST("mass_ptMM"), dimuon.mass(), dimuon.pt());
              registry.fill(HIST("mass_dcaxyMM"), dimuon.mass(), DCAmumu);
              registry.fill(HIST("mass_pcachiMM"), dimuon.mass(), dimuon.chi2pca());
              if (DCAmumu <= 0.02) {
                registry.fill(HIST("mass_ptMM_dcaxylow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptMM_dcaxyhigh"), dimuon.mass(), dimuon.pt());
              }
              if (dimuon.chi2pca() == -999) {
                registry.fill(HIST("mass_ptMM_pcachifail"), dimuon.mass(), dimuon.pt());
              } else if (dimuon.chi2pca() <= 2){
                registry.fill(HIST("mass_ptMM_pcachilow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptMM_pcachihigh"), dimuon.mass(), dimuon.pt());
              }
              if (DCA1 <= 0.02 && DCA2 <= 0.02) {
                registry.fill(HIST("mass_ptMM_singledcalow"), dimuon.mass(), dimuon.pt());
              } else if (DCA1 > 0.02 && DCA2 > 0.02) {
                registry.fill(HIST("mass_ptMM_singledcahigh"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptMM_singledcalowhigh"), dimuon.mass(), dimuon.pt());
              }
            } else if (dimuon.sign() == 2) {
              registry.fill(HIST("mass_ptPP"), dimuon.mass(), dimuon.pt());
              registry.fill(HIST("mass_dcaxyPP"), dimuon.mass(), DCAmumu);
              registry.fill(HIST("mass_pcachiPP"), dimuon.mass(), dimuon.chi2pca());
              if (DCAmumu <= 0.02) {
                registry.fill(HIST("mass_ptPP_dcaxylow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPP_dcaxyhigh"), dimuon.mass(), dimuon.pt());
              }
              if (dimuon.chi2pca() == -999) {
                registry.fill(HIST("mass_ptPP_pcachifail"), dimuon.mass(), dimuon.pt());
              } else if (dimuon.chi2pca() <= 2){
                registry.fill(HIST("mass_ptPP_pcachilow"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPP_pcachihigh"), dimuon.mass(), dimuon.pt());
              }
              if (DCA1 <= 0.02 && DCA2 <= 0.02) {
                registry.fill(HIST("mass_ptPP_singledcalow"), dimuon.mass(), dimuon.pt());
              } else if (DCA1 > 0.02 && DCA2 > 0.02) {
                registry.fill(HIST("mass_ptPP_singledcahigh"), dimuon.mass(), dimuon.pt());
              } else {
                registry.fill(HIST("mass_ptPP_singledcalowhigh"), dimuon.mass(), dimuon.pt());
              }
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
    adaptAnalysisTask<dimuonall_data>(cfgc)
  };
}
