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
      {"mass_pTPMnoambi", "mass_pTPMnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPMnoambi", "mass_DCAPMnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPMnoambi", "mass_LxyzPMnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPPnoambi", "mass_pTPPnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPPnoambi", "mass_DCAPPnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPPnoambi", "mass_LxyzPPnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTMMnoambi", "mass_pTMMnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAMMnoambi", "mass_DCAMMnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzMMnoambi", "mass_LxyzMMnoambi", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPMnoambi_over05", "mass_pTPMnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPMnoambi_over05", "mass_DCAPMnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPMnoambi_over05", "mass_LxyzPMnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPPnoambi_over05", "mass_pTPPnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPPnoambi_over05", "mass_DCAPPnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPPnoambi_over05", "mass_LxyzPPnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTMMnoambi_over05", "mass_pTMMnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAMMnoambi_over05", "mass_DCAMMnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzMMnoambi_over05", "mass_LxyzMMnoambi_over05", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPMnoambi_over07", "mass_pTPMnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPMnoambi_over07", "mass_DCAPMnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPMnoambi_over07", "mass_LxyzPMnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPPnoambi_over07", "mass_pTPPnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPPnoambi_over07", "mass_DCAPPnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPPnoambi_over07", "mass_LxyzPPnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTMMnoambi_over07", "mass_pTMMnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAMMnoambi_over07", "mass_DCAMMnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzMMnoambi_over07", "mass_LxyzMMnoambi_over07", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPMnoambi_over10", "mass_pTPMnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPMnoambi_over10", "mass_DCAPMnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPMnoambi_over10", "mass_LxyzPMnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPPnoambi_over10", "mass_pTPPnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPPnoambi_over10", "mass_DCAPPnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPPnoambi_over10", "mass_LxyzPPnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTMMnoambi_over10", "mass_pTMMnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAMMnoambi_over10", "mass_DCAMMnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzMMnoambi_over10", "mass_LxyzMMnoambi_over10", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPMnoambi_over15", "mass_pTPMnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPMnoambi_over15", "mass_DCAPMnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPMnoambi_over15", "mass_LxyzPMnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPPnoambi_over15", "mass_pTPPnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPPnoambi_over15", "mass_DCAPPnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPPnoambi_over15", "mass_LxyzPPnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTMMnoambi_over15", "mass_pTMMnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAMMnoambi_over15", "mass_DCAMMnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzMMnoambi_over15", "mass_LxyzMMnoambi_over15", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPMnoambi_over20", "mass_pTPMnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPMnoambi_over20", "mass_DCAPMnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPMnoambi_over20", "mass_LxyzPMnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTPPnoambi_over20", "mass_pTPPnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAPPnoambi_over20", "mass_DCAPPnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzPPnoambi_over20", "mass_LxyzPPnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
      {"mass_pTMMnoambi_over20", "mass_pTMMnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {120, 0.0, 30.0}}}},
      {"mass_DCAMMnoambi_over20", "mass_DCAMMnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {900, 0.0, 3.0}}}},
      {"mass_LxyzMMnoambi_over20", "mass_LxyzMMnoambi_over20", {HistType::kTH2F, {{250, 0.0, 5.0}, {1000, 0.0, 5}}}},
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
        if (dimuon.pt1() > 0.5 && dimuon.pt2() > 0.5) {
          if (dimuon.sign() == 0) {
            registry.fill(HIST("mass_pTPMnoambi_over05"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPMnoambi_over05"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPMnoambi_over05"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == -2) {
            registry.fill(HIST("mass_pTMMnoambi_over05"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAMMnoambi_over05"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzMMnoambi_over05"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == 2) {
            registry.fill(HIST("mass_pTPPnoambi_over05"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPPnoambi_over05"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPPnoambi_over05"), dimuon.mass(), lxyz);
          }
        }
        if (dimuon.pt1() > 0.7 && dimuon.pt2() > 0.7) {
          if (dimuon.sign() == 0) {
            registry.fill(HIST("mass_pTPMnoambi_over07"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPMnoambi_over07"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPMnoambi_over07"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == -2) {
            registry.fill(HIST("mass_pTMMnoambi_over07"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAMMnoambi_over07"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzMMnoambi_over07"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == 2) {
            registry.fill(HIST("mass_pTPPnoambi_over07"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPPnoambi_over07"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPPnoambi_over07"), dimuon.mass(), lxyz);
          }
        }
        if (dimuon.pt1() > 1.0 && dimuon.pt2() > 1.0) {
          if (dimuon.sign() == 0) {
            registry.fill(HIST("mass_pTPMnoambi_over10"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPMnoambi_over10"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPMnoambi_over10"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == -2) {
            registry.fill(HIST("mass_pTMMnoambi_over10"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAMMnoambi_over10"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzMMnoambi_over10"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == 2) {
            registry.fill(HIST("mass_pTPPnoambi_over10"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPPnoambi_over10"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPPnoambi_over10"), dimuon.mass(), lxyz);
          }
        }
        if (dimuon.pt1() > 1.5 && dimuon.pt2() > 1.5) {
          if (dimuon.sign() == 0) {
            registry.fill(HIST("mass_pTPMnoambi_over15"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPMnoambi_over15"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPMnoambi_over15"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == -2) {
            registry.fill(HIST("mass_pTMMnoambi_over15"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAMMnoambi_over15"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzMMnoambi_over15"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == 2) {
            registry.fill(HIST("mass_pTPPnoambi_over15"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPPnoambi_over15"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPPnoambi_over15"), dimuon.mass(), lxyz);
          }
        }
        if (dimuon.pt1() > 2.0 && dimuon.pt2() > 2.0) {
          if (dimuon.sign() == 0) {
            registry.fill(HIST("mass_pTPMnoambi_over20"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPMnoambi_over20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPMnoambi_over20"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == -2) {
            registry.fill(HIST("mass_pTMMnoambi_over20"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAMMnoambi_over20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzMMnoambi_over20"), dimuon.mass(), lxyz);
          } else if (dimuon.sign() == 2) {
            registry.fill(HIST("mass_pTPPnoambi_over20"), dimuon.mass(), dimuon.pt());
            registry.fill(HIST("mass_DCAPPnoambi_over20"), dimuon.mass(), DCAmumu);
            registry.fill(HIST("mass_LxyzPPnoambi_over20"), dimuon.mass(), lxyz);
          }
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
