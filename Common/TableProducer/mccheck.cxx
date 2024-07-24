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

struct mccheck {

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
  void process(aod::McCollisions const& mcCollisions, aod::McParticles const& mcTracks)
  {
    for (auto& collision : mcCollisions) {
      LOGF(info, "---------collision---------");
      for (auto& mctrack : mcTracks) {
        if (mctrack.mcCollisionId() == collision.globalIndex()){
          if (mctrack.pdgCode() == 221 || mctrack.pdgCode() == -221 || mctrack.pdgCode() == 331 || mctrack.pdgCode() == -331 || mctrack.pdgCode() == 223 || mctrack.pdgCode() == -223 || mctrack.pdgCode() == 333 || mctrack.pdgCode() == -333 || mctrack.pdgCode() == 113 || mctrack.pdgCode() == -113) {
            //LOGF(info, "---check from vector meson---", mctrack.globalIndex());
            //LOGF(info, "mctrack.pdgCode() = %d", mctrack.pdgCode());
            //LOGF(info, "mctrack.pdgCode() = %f", mctrack.pt());
            //LOGF(info, "mctrack.pdgCode() = %f", mctrack.eta());
            //if (mctrack.has_daughters()) {
              //for (auto& d : mctrack.daughters_as<aod::McParticles>()) {
                //LOGF(info, "daughter: pdg Code = %d", d.pdgCode());
              //}
            //}
          }
          if (mctrack.pdgCode() == 13 || mctrack.pdgCode() == -13) {
            LOGF(info, "---check from muon---");
            LOGF(info, "mctrack.pdgCode() = %d", mctrack.pdgCode());
            for (auto& d : mctrack.mothers_as<aod::McParticles>()) {
              LOGF(info, "mother: pdg Code = %d", d.pdgCode());
              LOGF(info, "mother: eta = %f", d.eta());
              LOGF(info, "mother: pt = %f", d.pt());
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
    adaptAnalysisTask<mccheck>(cfgc)
  };
}
