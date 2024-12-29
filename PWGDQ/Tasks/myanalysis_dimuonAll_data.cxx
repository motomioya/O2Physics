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

struct myanalysis_dimuonAll_data {

  HistogramRegistry registry{
    "registry", 
    {
      {"mass_pTPM", "mass_pTPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_pTPMnoambi", "mass_pTPMnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_DCAPM", "mass_DCAPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 10.0}}}},
      {"mass_DCAPMnoambi", "mass_DCAPMnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 10.0}}}},
      {"mass_LxyzPM", "mass_LxyzPM", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 20.0}}}},
      {"mass_LxyzPMnoambi", "mass_LxyzPMnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 20.0}}}},

      {"mass_pTPP", "mass_pTPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_pTPPnoambi", "mass_pTPPnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_DCAPP", "mass_DCAPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 10.0}}}},
      {"mass_DCAPPnoambi", "mass_DCAPPnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 10.0}}}},
      {"mass_LxyzPP", "mass_LxyzPP", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 20.0}}}},
      {"mass_LxyzPPnoambi", "mass_LxyzPPnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 20.0}}}},

      {"mass_pTMM", "mass_pTMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_pTMMnoambi", "mass_pTMMnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"mass_DCAMM", "mass_DCAMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 10.0}}}},
      {"mass_DCAMMnoambi", "mass_DCAMMnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 10.0}}}},
      {"mass_LxyzMM", "mass_LxyzMM", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 20.0}}}},
      {"mass_LxyzMMnoambi", "mass_LxyzMMnoambi", {HistType::kTH2F, {{750, 0.0, 15.0}, {2000, 0.0, 20.0}}}},
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
      if (!(dimuon.eta1() >= -3.6 && dimuon.eta1() <= -2.5)) continue;
      if (!(dimuon.eta2() >= -3.6 && dimuon.eta2() <= -2.5)) continue;
      if (dimuon.sign() == 0) {
        registry.fill(HIST("mass_pTPM"), dimuon.mass(), dimuon.pt());
        registry.fill(HIST("mass_DCAPM"), dimuon.mass(), DCAmumu);
        registry.fill(HIST("mass_LxyzPM"), dimuon.mass(), lxyz);
      } else if (dimuon.sign() == -2) {
        registry.fill(HIST("mass_pTMM"), dimuon.mass(), dimuon.pt());
        registry.fill(HIST("mass_DCAMM"), dimuon.mass(), DCAmumu);
        registry.fill(HIST("mass_LxyzMM"), dimuon.mass(), lxyz);
      } else if (dimuon.sign() == 2) {
        registry.fill(HIST("mass_pTPP"), dimuon.mass(), dimuon.pt());
        registry.fill(HIST("mass_DCAPP"), dimuon.mass(), DCAmumu);
        registry.fill(HIST("mass_LxyzPP"), dimuon.mass(), lxyz);
      }
      //ambiguous tracks
      if (dimuon.isAmbig1() == 0 && dimuon.isAmbig2() == 0) {
        if (dimuon.sign() == 0) {
          registry.fill(HIST("mass_pTPMnoambi"), dimuon.mass(), dimuon.pt());
          registry.fill(HIST("mass_DCAPMnoambi"), dimuon.mass(), DCAmumu);
          registry.fill(HIST("mass_LxyzPMnoambi"), dimuon.mass(), lxyz);
        } else if (dimuon.sign() == -2) {
          registry.fill(HIST("mass_pTMMnoambi"), dimuon.mass(), dimuon.pt());
          registry.fill(HIST("mass_DCAMMnoambi"), dimuon.mass(), DCAmumu);
          registry.fill(HIST("mass_LxyzMMnoambi"), dimuon.mass(), lxyz);
        } else if (dimuon.sign() == 2) {
          registry.fill(HIST("mass_pTPPnoambi"), dimuon.mass(), dimuon.pt());
          registry.fill(HIST("mass_DCAPPnoambi"), dimuon.mass(), DCAmumu);
          registry.fill(HIST("mass_LxyzPPnoambi"), dimuon.mass(), lxyz);
        }
      }
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myanalysis_dimuonAll_data>(cfgc)
  };
}
