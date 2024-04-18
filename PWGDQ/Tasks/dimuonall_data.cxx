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
      {"mass_svPM", "mass_svPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"mass_svMM", "mass_svMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"mass_svPP", "mass_svPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      /*
      {"mass_LxyzPM", "mass_LxyzPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_LxyzMM", "mass_LxyzMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_LxyzPP", "mass_LxyzPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      */
      {"mass_pcachiPM", "mass_pcachiPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 100}}}},
      {"mass_pcachiMM", "mass_pcachiMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 100}}}},
      {"mass_pcachiPP", "mass_pcachiPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 100}}}},

      {"chi2cutmass_ptPM", "chi2cutmass_ptPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"chi2cutmass_ptMM", "chi2cutmass_ptMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"chi2cutmass_ptPP", "chi2cutmass_ptPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {1200, 0.0, 30.0}}}},
      {"chi2cutmass_dcaxyPM", "chi2cutmass_dcaxyPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0, 15}}}},
      {"chi2cutmass_dcaxyMM", "chi2cutmass_dcaxyMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0, 15}}}},
      {"chi2cutmass_dcaxyPP", "chi2cutmass_dcaxyPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0, 15}}}},
      {"chi2cutmass_svPM", "chi2cutmass_svPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"chi2cutmass_svMM", "chi2cutmass_svMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"chi2cutmass_svPP", "chi2cutmass_svPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      /*
      {"chi2cutmass_LxyzPM", "chi2cutmass_LxyzPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_LxyzMM", "chi2cutmass_LxyzMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_LxyzPP", "chi2cutmass_LxyzPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      */
      {"chi2cutmass_pcachiPM", "chi2cutmass_pcachiPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 100}}}},
      {"chi2cutmass_pcachiMM", "chi2cutmass_pcachiMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 100}}}},
      {"chi2cutmass_pcachiPP", "chi2cutmass_pcachiPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 100}}}},
      {"chi2cutpcachiPM", "chi2cutpcachiPM", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"chi2cutpcachiMM", "chi2cutpcachiMM", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"chi2cutpcachiPP", "chi2cutpcachiPP", {HistType::kTH1F, {{1100, -1000, 100}}}},
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
      //float Tauxy = dimuon.tauxy();
      //float Tauz = dimuon.tauz();
      //float Lxy = (Tauz * dimuon.p() * o2::constants::physics::LightSpeedCm2NS) / dimuon.mass();;
      if (dimuon.sign() == 0) {
        registry.fill(HIST("mass_ptPM"), dimuon.mass(), dimuon.pt1());
        registry.fill(HIST("mass_dcaxyPM"), dimuon.mass(), DCAmumu);
        registry.fill(HIST("mass_pcachiPM"), dimuon.mass(), dimuon.chi2pca());
        registry.fill(HIST("mass_svPM"), dimuon.mass(), dimuon.sVertex());
      } else if (dimuon.sign() == -2) {
        registry.fill(HIST("mass_ptMM"), dimuon.mass(), dimuon.pt1());
        registry.fill(HIST("mass_dcaxyMM"), dimuon.mass(), DCAmumu);
        registry.fill(HIST("mass_pcachiMM"), dimuon.mass(), dimuon.chi2pca());
        registry.fill(HIST("mass_svMM"), dimuon.mass(), dimuon.sVertex());
      } else if (dimuon.sign() == 2) {
        registry.fill(HIST("mass_ptPP"), dimuon.mass(), dimuon.pt1());
        registry.fill(HIST("mass_dcaxyPP"), dimuon.mass(), DCAmumu);
        registry.fill(HIST("mass_pcachiPP"), dimuon.mass(), dimuon.chi2pca());
        registry.fill(HIST("mass_svPP"), dimuon.mass(), dimuon.sVertex());
      }
      if (dimuon.chi2MatchMCHMFT1() < 50 && dimuon.chi2MatchMCHMFT2() < 50) {
        if (dimuon.sign() == 0) {
          registry.fill(HIST("chi2cutmass_ptPM"), dimuon.mass(), dimuon.pt1());
          registry.fill(HIST("chi2cutmass_dcaxyPM"), dimuon.mass(), DCAmumu);
          registry.fill(HIST("chi2cutmass_pcachiPM"), dimuon.mass(), dimuon.chi2pca());
          registry.fill(HIST("chi2cutmass_svPM"), dimuon.mass(), dimuon.sVertex());
          registry.fill(HIST("chi2cutpcachiPM"), dimuon.chi2pca());
        } else if (dimuon.sign() == -2) {
          registry.fill(HIST("chi2cutmass_ptMM"), dimuon.mass(), dimuon.pt1());
          registry.fill(HIST("chi2cutmass_dcaxyMM"), dimuon.mass(), DCAmumu);
          registry.fill(HIST("chi2cutmass_pcachiMM"), dimuon.mass(), dimuon.chi2pca());
          registry.fill(HIST("chi2cutmass_svMM"), dimuon.mass(), dimuon.sVertex());
          registry.fill(HIST("chi2cutpcachiMM"), dimuon.chi2pca());
        } else if (dimuon.sign() == 2) {
          registry.fill(HIST("chi2cutmass_ptPP"), dimuon.mass(), dimuon.pt1());
          registry.fill(HIST("chi2cutmass_dcaxyPP"), dimuon.mass(), DCAmumu);
          registry.fill(HIST("chi2cutmass_pcachiPP"), dimuon.mass(), dimuon.chi2pca());
          registry.fill(HIST("chi2cutmass_svPP"), dimuon.mass(), dimuon.sVertex());
          registry.fill(HIST("chi2cutpcachiPP"), dimuon.chi2pca());
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
