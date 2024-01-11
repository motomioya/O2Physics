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
      {"muoncounter", "muoncounter", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
      {"daughtercounter", "daughtercounter", {HistType::kTH1F, {{5, -0.5, 4.5}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    //MCProng prongElectronNonPromptJpsi(3,{11,443,502},{true,true,true},{false,false,false},{0,0,0},{0,0,0},{false,false,false});
    //mySignal = new MCSignal("jpsiBeautyElectron", "Electrons from beauty jpsi decays", {prongElectronNonPromptJpsi}, {-1});
    //MCProng prongAllFromBeauty(2,{0,503},{true,true},{false,false},{0,0},{0,0},{false,false});
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(aod::McParticles const& mcTracks)
  {
    for (auto& mctrack : mcTracks) {
      
      if (-3.6 < mctrack.eta() && mctrack.eta() < -2.5) {
        if (mctrack.pdgCode() == 13 || mctrack.pdgCode() == -13) {
          registry.fill(HIST("muoncounter"), 0);
          const auto mcmothersidlist = mctrack.mothersIds();
          if (mcmothersidlist.size() > 0) {
            for (int i = 0; i <= mcmothersidlist.size() - 1; i++) {
              auto mcmother = mcTracks.iteratorAt(mcmothersidlist[i]);
              int mcmotherpdg = mcmother.pdgCode();
              //int mcmotherid = mcmother.globalIndex();
              if ( (mcmotherpdg >= 500 && mcmotherpdg <= 549) || (mcmotherpdg <= -500 && mcmotherpdg >= -549) ) {
                registry.fill(HIST("muoncounter"), 1);
                const auto mcgmothersidlist1 = mcmother.mothersIds();
                if (mcgmothersidlist1.size() > 0) {
                  for (int k = 0; k <= mcgmothersidlist1.size() - 1; k++) {
                    auto mcgmother = mcTracks.iteratorAt(mcgmothersidlist1[k]);
                    int mcgmotherpdg = mcgmother.pdgCode();
                    if ( (mcgmotherpdg == mcmotherpdg * -1) ) {
                      registry.fill(HIST("muoncounter"), 2);
                    }
                  }
                }
              }
              if ( (mcmotherpdg >= 400 && mcmotherpdg <= 439) || (mcmotherpdg <= -400 && mcmotherpdg >= -439) ) {
                registry.fill(HIST("muoncounter"), 3);
                const auto mcgmothersidlist2 = mcmother.mothersIds();
                if (mcgmothersidlist2.size() > 0) {
                  for (int k = 0; k <= mcgmothersidlist2.size() - 1; k++) {
                    auto mcgmother = mcTracks.iteratorAt(mcgmothersidlist2[k]);
                    int mcgmotherpdg = mcgmother.pdgCode();
                    if ( (mcgmotherpdg >= 500 && mcgmotherpdg <= 549) || (mcgmotherpdg <= -500 && mcgmotherpdg >= -549) ) {
                      registry.fill(HIST("muoncounter"), 4);
                    }
                  }
                }
              }
            }
          }
        }
        if ( (mctrack.pdgCode() >= 500 && mctrack.pdgCode() <= 549) || (mctrack.pdgCode() <= -500 && mctrack.pdgCode() >= -549) ) {
          registry.fill(HIST("daughtercounter"), 0);
          const auto mcdaughtersidlist1 = mctrack.daughtersIds();
          if (mcdaughtersidlist1.size() > 0) {
            for (int l = 0; l <= mcdaughtersidlist1.size() - 1; l++) {
              auto mcdaughter1 = mcTracks.iteratorAt(mcdaughtersidlist1[l]);
              int mcdaughterpdg = mcdaughter1.pdgCode();
              if ( (mcdaughterpdg == 13 || mcdaughterpdg == -13) ) {
                registry.fill(HIST("daughtercounter"), 1);
              }
            }
          }
        }
        if ( (mctrack.pdgCode() >= 400 && mctrack.pdgCode() <= 439) || (mctrack.pdgCode() <= -400 && mctrack.pdgCode() >= -439) ) {
          registry.fill(HIST("daughtercounter"), 2);
          const auto mcdaughtersidlist2 = mctrack.daughtersIds();
          if (mcdaughtersidlist2.size() > 0) {
            for (int l = 0; l <= mcdaughtersidlist2.size() - 1; l++) {
              auto mcdaughter2 = mcTracks.iteratorAt(mcdaughtersidlist2[l]);
              int mcdaughterpdg = mcdaughter2.pdgCode();
              if ( (mcdaughterpdg == 13 || mcdaughterpdg == -13) ) {
                registry.fill(HIST("daughtercounter"), 3);
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
    adaptAnalysisTask<mccheck>(cfgc)
  };
}
