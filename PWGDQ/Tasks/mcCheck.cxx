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
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;

struct mccheck {

  Preslice<aod::McParticles> particlePerCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::FwdTrack> fwdtrackIndicesPerCollision = aod::fwdtrack::collisionId;
  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  //MCSignal* mySignal;
  HistogramRegistry registry{
    "registry", //
    {
      {"hetaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hetaHasGammaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hetaHasPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hetaHasGammaPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hetaprimeMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hetaprimeHasGammaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hetaprimeHasPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hetaprimeHasGammaPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hrhoMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hrhoHasGammaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hrhoHasPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hrhoHasGammaPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"homegaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"homegaHasGammaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"homegaHasPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"homegaHasGammaPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hphiMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hphiHasGammaMcPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hphiHasPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},
      {"hphiHasGammaPi0McPtPM", "Invariant;Invariant McPt (GeV/c^{2});McPt", {HistType::kTH1F, {{1000, 0.0, 10.0}}}},

      {"hetaMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hetaHasGammaMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hetaHasPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hetaHasGammaPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hetaprimeMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hetaprimeHasGammaMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hetaprimeHasPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hetaprimeHasGammaPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hrhoMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hrhoHasGammaMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hrhoHasPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hrhoHasGammaPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"homegaMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"homegaHasGammaMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"homegaHasPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"homegaHasGammaPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hphiMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hphiHasGammaMcEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hphiHasPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},
      {"hphiHasGammaPi0McEtaPM", "Invariant;Invariant McEta (GeV/c^{2});McEta", {HistType::kTH1F, {{1000, -5.0, 0.0}}}},

      {"hetaMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hetaHasGammaMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hetaHasPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hetaHasGammaPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hetaprimeMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hetaprimeHasGammaMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hetaprimeHasPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hetaprimeHasGammaPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hrhoMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hrhoHasGammaMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hrhoHasPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hrhoHasGammaPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"homegaMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"homegaHasGammaMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"homegaHasPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"homegaHasGammaPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hphiMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hphiHasGammaMcPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hphiHasPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},
      {"hphiHasGammaPi0McPhiPM", "Invariant;Invariant McPhi (GeV/c^{2});McPhi", {HistType::kTH1F, {{1000, -10.0, 10.0}}}},

      {"hetaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hetaHasGammaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hetaHasPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hetaHasGammaPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hetaprimeMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hetaprimeHasGammaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hetaprimeHasPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hetaprimeHasGammaPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hrhoMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hrhoHasGammaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hrhoHasPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hrhoHasGammaPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"homegaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"homegaHasGammaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"homegaHasPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"homegaHasGammaPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hphiMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hphiHasGammaMcMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hphiHasPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
      {"hphiHasGammaPi0McMassPM", "Invariant;Invariant McMass (GeV/c^{2});McMass", {HistType::kTH1F, {{1000, 0.0, 5.0}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(aod::McCollisions const& mcCollisions, aod::McParticles const& mcTracks, aod::Collisions const& collisions, MyMuons const& fwdtracks)
  {
    for (auto& collision : mcCollisions) {
      auto particlethiscollision = mcTracks.sliceBy(particlePerCollision, collision.globalIndex());
      for (auto& particle : particlethiscollision) {
        if (particle.pdgCode() == 221 || particle.pdgCode() == 331 || particle.pdgCode() == 113 || particle.pdgCode() == 223 || particle.pdgCode() == 333) {
          if (particle.eta() < -4.0 || particle.eta() > -2.5) continue;
          bool hasmuplus = 0;
          bool hasmuminus = 0;
          bool haspi0 = 0;
          bool hasgamma = 0;
          int processCodemuplus = 0;
          int processCodemuminus = 0;
          int processCodepion = 0;
          int processCodegamma = 0;
          TLorentzVector lv1, lv2, lv;
          for (auto& particledaughter : particlethiscollision) {
            const auto mcmothersidlist = particledaughter.mothersIds();
            if (mcmothersidlist.size() > 0) {
              for (auto i = 0; i <= mcmothersidlist.size() - 1; i++) {
                auto mcmother = mcTracks.iteratorAt(mcmothersidlist[i]);
                if (mcmother.globalIndex() == particle.globalIndex()) {
                  if (particledaughter.pdgCode() == 13 && particledaughter.getProcess() == 0) {
                    hasmuplus = 1;
                    lv1.SetPtEtaPhiM(particledaughter.pt(), particledaughter.eta(), particledaughter.phi(), muonMass);
                    processCodemuplus = particledaughter.getProcess();
                  }
                  if (particledaughter.pdgCode() == -13 && particledaughter.getProcess() == 0) {
                    hasmuminus = 1;
                    lv2.SetPtEtaPhiM(particledaughter.pt(), particledaughter.eta(), particledaughter.phi(), muonMass);
                    processCodemuminus = particledaughter.getProcess();
                  }
                  if (particledaughter.pdgCode() == 111 && particledaughter.getProcess() == 0) {
                    haspi0 = 1;
                    processCodepion = particledaughter.getProcess();
                  }
                  if (particledaughter.pdgCode() == 22 && particledaughter.getProcess() == 0) {
                    hasgamma = 1;
                    processCodegamma = particledaughter.getProcess();
                  }
                }
              }
            }
          }
          lv = lv1 + lv2;
          if (hasmuplus == 1 && hasmuminus == 1 && haspi0 == 0 && hasgamma == 0) {
            /*
            LOGF(info, "--------");
            LOGF(info, "particle.pdgCode() == %d", particle.pdgCode());
            LOGF(info, "hasmuplus == 1 && hasmuminus == 1 && haspi0 == 0 && hasgamma == 0");
            LOGF(info, "processCodemuplus == %d", processCodemuplus);
            LOGF(info, "processCodemuminus == %d", processCodemuminus);
            */
            if (particle.pdgCode() == 221) {
              registry.fill(HIST("hetaMcPtPM"), particle.pt());
              registry.fill(HIST("hetaMcEtaPM"), particle.eta());
              registry.fill(HIST("hetaMcPhiPM"), particle.phi());
              registry.fill(HIST("hetaMcMassPM"), lv.M());
            } else if (particle.pdgCode() == 331) {
              registry.fill(HIST("hetaprimeMcPtPM"), particle.pt());
              registry.fill(HIST("hetaprimeMcEtaPM"), particle.eta());
              registry.fill(HIST("hetaprimeMcPhiPM"), particle.phi());
              registry.fill(HIST("hetaprimeMcMassPM"), lv.M());
            } else if (particle.pdgCode() == 113) {
              registry.fill(HIST("hrhoMcPtPM"), particle.pt());
              registry.fill(HIST("hrhoMcEtaPM"), particle.eta());
              registry.fill(HIST("hrhoMcPhiPM"), particle.phi());
              registry.fill(HIST("hrhoMcMassPM"), lv.M());
            } else if (particle.pdgCode() == 223) {
              registry.fill(HIST("homegaMcPtPM"), particle.pt());
              registry.fill(HIST("homegaMcEtaPM"), particle.eta());
              registry.fill(HIST("homegaMcPhiPM"), particle.phi());
              registry.fill(HIST("homegaMcMassPM"), lv.M());
            } else if (particle.pdgCode() == 333) {
              registry.fill(HIST("hphiMcPtPM"), particle.pt());
              registry.fill(HIST("hphiMcEtaPM"), particle.eta());
              registry.fill(HIST("hphiMcPhiPM"), particle.phi());
              registry.fill(HIST("hphiMcMassPM"), lv.M());
            }
          }
          if (hasmuplus == 1 && hasmuminus == 1 && haspi0 == 1 && hasgamma == 0) {
            /*
            LOGF(info, "--------");
            LOGF(info, "particle.pdgCode() == %d", particle.pdgCode());
            LOGF(info, "hasmuplus == 1 && hasmuminus == 1 && haspi0 == 1 && hasgamma == 0");
            LOGF(info, "processCodemuplus == %d", processCodemuplus);
            LOGF(info, "processCodemuminus == %d", processCodemuminus);
            LOGF(info, "processCodepion == %d", processCodepion);
            */
            if (processCodepion != 0) continue;
            if (particle.pdgCode() == 221) {
              registry.fill(HIST("hetaHasPi0McPtPM"), particle.pt());
              registry.fill(HIST("hetaHasPi0McEtaPM"), particle.eta());
              registry.fill(HIST("hetaHasPi0McPhiPM"), particle.phi());
              registry.fill(HIST("hetaHasPi0McMassPM"), lv.M());
            } else if (particle.pdgCode() == 331) {
              registry.fill(HIST("hetaprimeHasPi0McPtPM"), particle.pt());
              registry.fill(HIST("hetaprimeHasPi0McEtaPM"), particle.eta());
              registry.fill(HIST("hetaprimeHasPi0McPhiPM"), particle.phi());
              registry.fill(HIST("hetaprimeHasPi0McMassPM"), lv.M());
            } else if (particle.pdgCode() == 113) {
              registry.fill(HIST("hrhoHasPi0McPtPM"), particle.pt());
              registry.fill(HIST("hrhoHasPi0McEtaPM"), particle.eta());
              registry.fill(HIST("hrhoHasPi0McPhiPM"), particle.phi());
              registry.fill(HIST("hrhoHasPi0McMassPM"), lv.M());
            } else if (particle.pdgCode() == 223) {
              registry.fill(HIST("homegaHasPi0McPtPM"), particle.pt());
              registry.fill(HIST("homegaHasPi0McEtaPM"), particle.eta());
              registry.fill(HIST("homegaHasPi0McPhiPM"), particle.phi());
              registry.fill(HIST("homegaHasPi0McMassPM"), lv.M());
            } else if (particle.pdgCode() == 333) {
              registry.fill(HIST("hphiHasPi0McPtPM"), particle.pt());
              registry.fill(HIST("hphiHasPi0McEtaPM"), particle.eta());
              registry.fill(HIST("hphiHasPi0McPhiPM"), particle.phi());
              registry.fill(HIST("hphiHasPi0McMassPM"), lv.M());
            }
          }
          if (hasmuplus == 1 && hasmuminus == 1 && haspi0 == 0 && hasgamma == 1) {
            /*
            LOGF(info, "--------");
            LOGF(info, "particle.pdgCode() == %d", particle.pdgCode());
            LOGF(info, "hasmuplus == 1 && hasmuminus == 1 && haspi0 == 0 && hasgamma == 1");
            LOGF(info, "processCodemuplus == %d", processCodemuplus);
            LOGF(info, "processCodemuminus == %d", processCodemuminus);
            LOGF(info, "processCodegamma == %d", processCodegamma);
            */
            if (particle.pdgCode() == 221) {
              registry.fill(HIST("hetaHasGammaMcPtPM"), particle.pt());
              registry.fill(HIST("hetaHasGammaMcEtaPM"), particle.eta());
              registry.fill(HIST("hetaHasGammaMcPhiPM"), particle.phi());
              registry.fill(HIST("hetaHasGammaMcMassPM"), lv.M());
            } else if (particle.pdgCode() == 331) {
              registry.fill(HIST("hetaprimeHasGammaMcPtPM"), particle.pt());
              registry.fill(HIST("hetaprimeHasGammaMcEtaPM"), particle.eta());
              registry.fill(HIST("hetaprimeHasGammaMcPhiPM"), particle.phi());
              registry.fill(HIST("hetaprimeHasGammaMcMassPM"), lv.M());
            } else if (particle.pdgCode() == 113) {
              registry.fill(HIST("hrhoHasGammaMcPtPM"), particle.pt());
              registry.fill(HIST("hrhoHasGammaMcEtaPM"), particle.eta());
              registry.fill(HIST("hrhoHasGammaMcPhiPM"), particle.phi());
              registry.fill(HIST("hrhoHasGammaMcMassPM"), lv.M());
            } else if (particle.pdgCode() == 223) {
              registry.fill(HIST("homegaHasGammaMcPtPM"), particle.pt());
              registry.fill(HIST("homegaHasGammaMcEtaPM"), particle.eta());
              registry.fill(HIST("homegaHasGammaMcPhiPM"), particle.phi());
              registry.fill(HIST("homegaHasGammaMcMassPM"), lv.M());
            } else if (particle.pdgCode() == 333) {
              registry.fill(HIST("hphiHasGammaMcPtPM"), particle.pt());
              registry.fill(HIST("hphiHasGammaMcEtaPM"), particle.eta());
              registry.fill(HIST("hphiHasGammaMcPhiPM"), particle.phi());
              registry.fill(HIST("hphiHasGammaMcMassPM"), lv.M());
            }
          }
          if (hasmuplus == 1 && hasmuminus == 1 && haspi0 == 1 && hasgamma == 1) {
            /*
            LOGF(info, "--------");
            LOGF(info, "particle.pdgCode() == %d", particle.pdgCode());
            LOGF(info, "hasmuplus == 1 && hasmuminus == 1 && haspi0 == 1 && hasgamma == 1");
            LOGF(info, "processCodemuplus == %d", processCodemuplus);
            LOGF(info, "processCodemuminus == %d", processCodemuminus);
            LOGF(info, "processCodepion == %d", processCodepion);
            LOGF(info, "processCodegamma == %d", processCodegamma);
            */
            if (particle.pdgCode() == 221) {
              registry.fill(HIST("hetaHasGammaPi0McPtPM"), particle.pt());
              registry.fill(HIST("hetaHasGammaPi0McEtaPM"), particle.eta());
              registry.fill(HIST("hetaHasGammaPi0McPhiPM"), particle.phi());
              registry.fill(HIST("hetaHasGammaPi0McMassPM"), lv.M());
            } else if (particle.pdgCode() == 331) {
              registry.fill(HIST("hetaprimeHasGammaPi0McPtPM"), particle.pt());
              registry.fill(HIST("hetaprimeHasGammaPi0McEtaPM"), particle.eta());
              registry.fill(HIST("hetaprimeHasGammaPi0McPhiPM"), particle.phi());
              registry.fill(HIST("hetaprimeHasGammaPi0McMassPM"), lv.M());
            } else if (particle.pdgCode() == 113) {
              registry.fill(HIST("hrhoHasGammaPi0McPtPM"), particle.pt());
              registry.fill(HIST("hrhoHasGammaPi0McEtaPM"), particle.eta());
              registry.fill(HIST("hrhoHasGammaPi0McPhiPM"), particle.phi());
              registry.fill(HIST("hrhoHasGammaPi0McMassPM"), lv.M());
            } else if (particle.pdgCode() == 223) {
              registry.fill(HIST("homegaHasGammaPi0McPtPM"), particle.pt());
              registry.fill(HIST("homegaHasGammaPi0McEtaPM"), particle.eta());
              registry.fill(HIST("homegaHasGammaPi0McPhiPM"), particle.phi());
              registry.fill(HIST("homegaHasGammaPi0McMassPM"), lv.M());
            } else if (particle.pdgCode() == 333) {
              registry.fill(HIST("hphiHasGammaPi0McPtPM"), particle.pt());
              registry.fill(HIST("hphiHasGammaPi0McEtaPM"), particle.eta());
              registry.fill(HIST("hphiHasGammaPi0McPhiPM"), particle.phi());
              registry.fill(HIST("hphiHasGammaPi0McMassPM"), lv.M());
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
