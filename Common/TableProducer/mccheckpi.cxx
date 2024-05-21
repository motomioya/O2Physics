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
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"

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

struct mccheckpi {

  //MCSignal* mySignal;
  HistogramRegistry registry{
    "registry", //
    {
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(aod::McParticles const& mcTracks)
  {

    for (auto& mctrack : mcTracks) {
      if (mctrack.eta() > -3.6 && mctrack.eta() < -2.5){
        if (mctrack.pdgCode() == 211 || mctrack.pdgCode() == -211) {
          LOGF(info, "------find forward pion(+/-)------");
          if (mctrack.has_daughters()){
            LOGF(info, "----has daughters----");
            const auto daughters = mctrack.daughters_as<aod::McParticles>();
            for (auto& daughter : daughters) {
              LOGF(info, "--daughter particle--");
              LOGF(info, "production point (z): %f", daughter.vz());
              LOGF(info, "pdg code: %d", daughter.pdgCode());
              LOGF(info, "VMC physics code: %d", daughter.getProcess());
              LOGF(info, "pion, isPhysicsPrimary: %i", mctrack.isPhysicalPrimary());
            }
          } else {
            LOGF(info, "----no daughter----");
          }
        }
      }
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mccheckpi>(cfgc)
  };
}
