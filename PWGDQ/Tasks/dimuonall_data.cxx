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
      {"pt1PM", "pt1PM", {HistType::kTH1F, {{120, 0.0, 30.0}}}},
      {"pt2PM", "pt2PM", {HistType::kTH1F, {{120, 0.0, 30.0}}}},
      {"massPM", "massPM", {HistType::kTH1F, {{120, 0.0, 30.0}}}},
      {"mass_ptPM", "mass_ptPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_ptMM", "mass_ptMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_ptPP", "mass_ptPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy0PM", "mass_pt_dcaxy0PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy0MM", "mass_pt_dcaxy0MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy0PP", "mass_pt_dcaxy0PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy1PM", "mass_pt_dcaxy1PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy1MM", "mass_pt_dcaxy1MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy1PP", "mass_pt_dcaxy1PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy2PM", "mass_pt_dcaxy2PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy2MM", "mass_pt_dcaxy2MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy2PP", "mass_pt_dcaxy2PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy3PM", "mass_pt_dcaxy3PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy3MM", "mass_pt_dcaxy3MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy3PP", "mass_pt_dcaxy3PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy4PM", "mass_pt_dcaxy4PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy4MM", "mass_pt_dcaxy4MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"mass_pt_dcaxy4PP", "mass_pt_dcaxy4PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"dcaxyPM", "dcaxyPM", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"dcaxyMM", "dcaxyMM", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"dcaxyPP", "dcaxyPP", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"chi2cutmass_ptPM", "chi2cutmass_ptPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_ptMM", "chi2cutmass_ptMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_ptPP", "chi2cutmass_ptPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy0PM", "chi2cutmass_pt_dcaxy0PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy0MM", "chi2cutmass_pt_dcaxy0MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy0PP", "chi2cutmass_pt_dcaxy0PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy1PM", "chi2cutmass_pt_dcaxy1PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy1MM", "chi2cutmass_pt_dcaxy1MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy1PP", "chi2cutmass_pt_dcaxy1PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy2PM", "chi2cutmass_pt_dcaxy2PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy2MM", "chi2cutmass_pt_dcaxy2MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy2PP", "chi2cutmass_pt_dcaxy2PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy3PM", "chi2cutmass_pt_dcaxy3PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy3MM", "chi2cutmass_pt_dcaxy3MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy3PP", "chi2cutmass_pt_dcaxy3PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy4PM", "chi2cutmass_pt_dcaxy4PM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy4MM", "chi2cutmass_pt_dcaxy4MM", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutmass_pt_dcaxy4PP", "chi2cutmass_pt_dcaxy4PP", {HistType::kTH2F, {{750, 0.0, 15.0}, {120, 0.0, 30.0}}}},
      {"chi2cutdcaxyPM", "chi2cutdcaxyPM", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"chi2cutdcaxyMM", "chi2cutdcaxyMM", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"chi2cutdcaxyPP", "chi2cutdcaxyPP", {HistType::kTH1F, {{1000, 0, 10}}}},
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
        registry.fill(HIST("pt1PM"), dimuon.pt1());
        registry.fill(HIST("pt2PM"), dimuon.pt2());
        registry.fill(HIST("massPM"), dimuon.mass());
        registry.fill(HIST("dcaxyPM"), DCAmumu);
        registry.fill(HIST("mass_ptPM"), dimuon.mass(), dimuon.pt());
        if (DCAmumu < 0.2) {
          registry.fill(HIST("mass_pt_dcaxy0PM"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 0.5) {
          registry.fill(HIST("mass_pt_dcaxy1PM"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 1) {
          registry.fill(HIST("mass_pt_dcaxy2PM"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 2) {
          registry.fill(HIST("mass_pt_dcaxy3PM"), dimuon.mass(), dimuon.pt());
        } else {
          registry.fill(HIST("mass_pt_dcaxy4PM"), dimuon.mass(), dimuon.pt());
        }
      } else if (dimuon.sign() == -2) {
        registry.fill(HIST("dcaxyMM"), DCAmumu);
        registry.fill(HIST("mass_ptMM"), dimuon.mass(), dimuon.pt());
        if (DCAmumu < 0.2) {
          registry.fill(HIST("mass_pt_dcaxy0MM"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 0.5) {
          registry.fill(HIST("mass_pt_dcaxy1MM"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 1) {
          registry.fill(HIST("mass_pt_dcaxy2MM"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 2) {
          registry.fill(HIST("mass_pt_dcaxy3MM"), dimuon.mass(), dimuon.pt());
        } else {
          registry.fill(HIST("mass_pt_dcaxy4MM"), dimuon.mass(), dimuon.pt());
        }
      } else if (dimuon.sign() == 2) {
        registry.fill(HIST("dcaxyPP"), DCAmumu);
        registry.fill(HIST("mass_ptPP"), dimuon.mass(), dimuon.pt());
        if (DCAmumu < 0.2) {
          registry.fill(HIST("mass_pt_dcaxy0PP"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 0.5) {
          registry.fill(HIST("mass_pt_dcaxy1PP"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 1) {
          registry.fill(HIST("mass_pt_dcaxy2PP"), dimuon.mass(), dimuon.pt());
        } else if (DCAmumu < 2) {
          registry.fill(HIST("mass_pt_dcaxy3PP"), dimuon.mass(), dimuon.pt());
        } else {
          registry.fill(HIST("mass_pt_dcaxy4PP"), dimuon.mass(), dimuon.pt());
        }
      }

      if (dimuon.chi2MatchMCHMFT1() < 40) {
        if (dimuon.sign() == 0) {
          registry.fill(HIST("chi2cutdcaxyPM"), DCAmumu);
          registry.fill(HIST("chi2cutmass_ptPM"), dimuon.mass(), dimuon.pt());
          if (DCAmumu < 0.2) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy0PM"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 0.5) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy1PM"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 1) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy2PM"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 2) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy3PM"), dimuon.mass(), dimuon.pt());
          } else {
            registry.fill(HIST("chi2cutmass_pt_dcaxy4PM"), dimuon.mass(), dimuon.pt());
          }
        } else if (dimuon.sign() == -2) {
          registry.fill(HIST("chi2cutdcaxyMM"), DCAmumu);
          registry.fill(HIST("chi2cutmass_ptMM"), dimuon.mass(), dimuon.pt());
          if (DCAmumu < 0.2) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy0MM"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 0.5) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy1MM"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 1) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy2MM"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 2) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy3MM"), dimuon.mass(), dimuon.pt());
          } else {
            registry.fill(HIST("chi2cutmass_pt_dcaxy4MM"), dimuon.mass(), dimuon.pt());
          }
        } else if (dimuon.sign() == 2) {
          registry.fill(HIST("chi2cutdcaxyPP"), DCAmumu);
          registry.fill(HIST("chi2cutmass_ptPP"), dimuon.mass(), dimuon.pt());
          if (DCAmumu < 0.2) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy0PP"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 0.5) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy1PP"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 1) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy2PP"), dimuon.mass(), dimuon.pt());
          } else if (DCAmumu < 2) {
            registry.fill(HIST("chi2cutmass_pt_dcaxy3PP"), dimuon.mass(), dimuon.pt());
          } else {
            registry.fill(HIST("chi2cutmass_pt_dcaxy4PP"), dimuon.mass(), dimuon.pt());
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
