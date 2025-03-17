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
#include "Common/DataModel/CollisionAssociationTables.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
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

struct mcCheck_cc_single {

  HistogramRegistry registry{
    "registry",
    {
      {"hD0Pt", "hD0Pt;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hD0PtHasMuon", "hD0PtHasMuon;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hDplusPt", "hDplusPt;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hDplusPtHasMuon", "hDplusPtHasMuon;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hDsPt", "hDsPt;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hDsPtHasMuon", "hDsPtHasMuon;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hLambdacPt", "hLambdacPt;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hLambdacPtHasMuon", "hLambdacPtHasMuon;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hXiC0Pt", "hXiC0Pt;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hXiC0PtHasMuon", "hXiC0PtHasMuon;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hXiCplusPt", "hXiCplusPt;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hXiCplusPtHasMuon", "hXiCplusPtHasMuon;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hOmegacPt", "hOmegacPt;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
      {"hOmegacPtHasMuon", "hOmegacPtHasMuon;pT;entries", {HistType::kTH1F, {{1000, 0, 10.0}}}},
    }};

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    MCProng prongD0(2, {13, 421}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalD0;
    signalD0 = new MCSignal("muon from DO", "Muons from D0 decays", {prongD0}, {-1});

    MCProng prongDplus(2, {13, 411}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongDplus.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalDplus;
    signalDplus = new MCSignal("muon from Dplus", "Muons from Dplus decays", {prongDplus}, {-1});

    MCProng prongDs(2, {13, 431}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongDs.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalDs;
    signalDs = new MCSignal("muon from Ds", "Muons from Ds decays", {prongDs}, {-1});

    MCProng prongLambdac(2, {13, 4122}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLambdac.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalLambdac;
    signalLambdac = new MCSignal("muon from Lambdac", "Muons from Lambdac decays", {prongLambdac}, {-1});

    MCProng prongXiC0(2, {13, 4132}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongXiC0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalXiC0;
    signalXiC0 = new MCSignal("muon from XiC0", "Muons from XiC0 decays", {prongXiC0}, {-1});

    MCProng prongXiCplus(2, {13, 4232}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongXiCplus.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalXiCplus;
    signalXiCplus = new MCSignal("muon from XiCplus", "Muons from XiCplus decays", {prongXiCplus}, {-1});

    MCProng prongOmegac(2, {13, 4332}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongOmegac.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCSignal* signalOmegac;
    signalOmegac = new MCSignal("muon from Omegac", "Muons from Omegac decays", {prongOmegac}, {-1});

    for (auto& particle : mcParticles) {
      if (particle.eta() >= -3.6 && particle.eta() <= -2.5) {
        if (particle.pdgCode() == 421 || particle.pdgCode() == -421) {
          registry.fill(HIST("hD0Pt"), particle.pt());
        }
        if (particle.pdgCode() == 411 || particle.pdgCode() == -411) {
          registry.fill(HIST("hDplusPt"), particle.pt());
        }
        if (particle.pdgCode() == 431 || particle.pdgCode() == -431) {
          registry.fill(HIST("hDsPt"), particle.pt());
        }
        if (particle.pdgCode() == 4122 || particle.pdgCode() == -4122) {
          registry.fill(HIST("hLambdacPt"), particle.pt());
        }
        if (particle.pdgCode() == 4132 || particle.pdgCode() == -4132) {
          registry.fill(HIST("hXiC0Pt"), particle.pt());
        }
        if (particle.pdgCode() == 4232 || particle.pdgCode() == -4232) {
          registry.fill(HIST("hXiCplusPt"), particle.pt());
        }
        if (particle.pdgCode() == 4332 || particle.pdgCode() == -4332) {
          registry.fill(HIST("hOmegacPt"), particle.pt());
        }
        if (particle.pdgCode() == 13 || particle.pdgCode() == -13) {
          if((*signalD0).CheckSignal(true, particle)) {
            registry.fill(HIST("hD0PtHasMuon"), particle.pt());
          }
          if((*signalDplus).CheckSignal(true, particle)) {
            registry.fill(HIST("hDplusPtHasMuon"), particle.pt());
          }
          if((*signalDs).CheckSignal(true, particle)) {
            registry.fill(HIST("hDsPtHasMuon"), particle.pt());
          }
          if((*signalLambdac).CheckSignal(true, particle)) {
            registry.fill(HIST("hLambdacPtHasMuon"), particle.pt());
          }
          if((*signalXiC0).CheckSignal(true, particle)) {
            registry.fill(HIST("hXiC0PtHasMuon"), particle.pt());
          }
          if((*signalXiCplus).CheckSignal(true, particle)) {
            registry.fill(HIST("hXiCplusPtHasMuon"), particle.pt());
          }
          if((*signalOmegac).CheckSignal(true, particle)) {
            registry.fill(HIST("hOmegacPtHasMuon"), particle.pt());
          }
        }
      }
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mcCheck_cc_single>(cfgc)
  };
}
