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

  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  //MCSignal* mySignal;
  HistogramRegistry registry{
    "registry", //
    {
      {"hetaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hetaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hetaprimeMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hetaprimeMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hrhoMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hrhoMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"homegaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"homegaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hphiMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hphiMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(aod::McCollisions const& mcCollision, aod::McParticles const& mcTracks)
  {
    for (auto& mctrack : mcTracks) {
      if (mctrack.pdgCode() == 221 || mctrack.pdgCode() == 331 || mctrack.pdgCode() == 113 || mctrack.pdgCode() == 223 || mctrack.pdgCode() == 333) {
        TLorentzVector lv1, lv2, lv;

        const auto mcdaughtersidlist = mctrack.daughtersIds();
        bool hasmuon1 = 0;
        bool hasmuon2 = 0;
        if (mcdaughtersidlist.size() > 0 && mcdaughtersidlist.size() > 0) {
          for (auto i = 0; i <= mcdaughtersidlist.size() - 1; i++) {
            auto mcdaughter = mcTracks.iteratorAt(mcdaughtersidlist[i]);
            if (mcdaughter.pdgCode() == 13) {
              lv1.SetPtEtaPhiM(mcdaughter.pt(), mcdaughter.eta(), mcdaughter.phi(), muonMass);
              hasmuon1 = 1;
            }
            if (mcdaughter.pdgCode() == -13) {
              lv2.SetPtEtaPhiM(mcdaughter.pt(), mcdaughter.eta(), mcdaughter.phi(), muonMass);
              hasmuon2 = 1;
            }
          }
          if (hasmuon1 == 0 || hasmuon2 == 0) continue;
          lv = lv1 + lv2;
          if (mctrack.pdgCode() == 221) {
            registry.fill(HIST("hetaMcMassPM"), lv.M());
            registry.fill(HIST("hetaMcPtPM"), mctrack.pt());
          }
          if (mctrack.pdgCode() == 331) {
            registry.fill(HIST("hetaprimeMcMassPM"), lv.M());
            registry.fill(HIST("hetaprimeMcPtPM"), mctrack.pt());
          }
          if (mctrack.pdgCode() == 113) {
            registry.fill(HIST("hrhoMcMassPM"), lv.M());
            registry.fill(HIST("hrhoMcPtPM"), mctrack.pt());
          }
          if (mctrack.pdgCode() == 223) {
            registry.fill(HIST("homegaMcMassPM"), lv.M());
            registry.fill(HIST("homegaMcPtPM"), mctrack.pt());
          }
          if (mctrack.pdgCode() == 333) {
            registry.fill(HIST("hphiMcMassPM"), lv.M());
            registry.fill(HIST("hphiMcPtPM"), mctrack.pt());
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
