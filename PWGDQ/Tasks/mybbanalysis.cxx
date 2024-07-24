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
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include <math.h>
#include <TLorentzVector.h>
#include <string>
#include <regex>
#include <iostream>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"

#include "DCAFitter/FwdDCAFitterN.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "Math/Vector3D.h"
#include "CommonConstants/PhysicsConstants.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace std;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTrkCompColls, aod::McFwdTrackLabels>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;

struct mybbanalysis {

  float etalow = -4;
  float etaup = -2.5;
  float pDCAcutrAtBsorberEndlow1 = 17.6;
  float pDCAcutrAtBsorberEndup1 = 26.5;
  float pDCAcutrAtBsorberEndlow2 = 26.5;
  float pDCAcutrAtBsorberEndup2 = 89.5;
  float pDCAcutdcaup1 = 594;
  float pDCAcutdcaup2 = 324;
  float chi2up = 1000000;
  float chi2MatchMCHMIDup = 1000000;
  float collisionZcut = 10.0f;

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);
  Filter collisionFilter = nabs(aod::collision::posZ) < collisionZcut;
  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  Preslice<aod::McParticles> particlesIndicesPerCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = -5.0;
  //o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;


  HistogramRegistry registry{
    "registry", 
    {
      {"hMassvsDCAPM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAMM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAPP", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAbbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAvecPM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAvecMM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAvecPP", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},

      {"hMassvsPtPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtMM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtPP", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtbbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtvecPM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtvecMM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtvecPP", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},

      {"hMassvsLxyzPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzMM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzPP", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzbbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzvecPM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzvecMM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzvecPP", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},

      {"hMassvsChi2PM", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2MM", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2PP", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2bbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2vecPM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2vecMM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2vecPP", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},

    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, soa::Filtered<MyMuons> const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const& particles)
  {
    o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;

    bool propagateToPCA = true;
    float maxR = 200;
    float minParamChange = 1.0e-3;
    float minRelChi2Change = 0.9;
    bool useAbsDCA = false;
    //double singleptcut = 0.5;

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    fgFitterTwoProngFwd.setTGeoMat(false);
    fgFitterTwoProngFwd.setMatLUT(lut);
    fgFitterTwoProngFwd.setBz(mMagField);
    fgFitterTwoProngFwd.setPropagateToPCA(propagateToPCA);
    fgFitterTwoProngFwd.setMaxR(maxR);
    fgFitterTwoProngFwd.setMinParamChange(minParamChange);
    fgFitterTwoProngFwd.setMinRelChi2Change(minRelChi2Change);
    fgFitterTwoProngFwd.setUseAbsDCA(useAbsDCA);

    //Configuring signals bb -> mumu
    MCSignal* signalbb;
    MCProng prongbmu(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongbmu.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb = new MCSignal("bbTomumu", "mumu pairs from b->mu and b->mu", {prongbmu, prongbmu}, {-1, -1});
    //Configuring signals b->cmu->mumu
    MCSignal* signalbb1;
    MCProng prongbmu1A(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongbmu1B(3, {13, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prongbmu1A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu1B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb1 = new MCSignal("bbcTomumusingle", "mumu pairs from b->e and b->c->e (single b)", {prongbmu1A, prongbmu1B}, {1, 2}); // signal at pair level
    //Configuring signals b->cmu->mumu, inversed
    MCSignal* signalbb2;
    MCProng prongbmu2A(3, {13, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongbmu2B(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongbmu2A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu2B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb2 = new MCSignal("bcbTomumusingle", "mumu pairs from b->e and b->c->e (single b), inversed", {prongbmu2A, prongbmu2B}, {2, 1}); // signal at pair level
    //Configuring signals b->cmu->mumu
    MCSignal* signalbb3;
    MCProng prongbmu3A(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongbmu3B(4, {13, 402, 402, 502}, {true, true, true, true}, {false, false, false, false}, {0, 0, 0, 0}, {0, 0, 0, 0}, {false, false, false, false});
    prongbmu3A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu3B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb3 = new MCSignal("bbcTomumusingle", "mumu pairs from b->e and b->c->e (single b)", {prongbmu3A, prongbmu3B}, {1, 3}); // signal at pair level
    //Configuring signals b->cmu->mumu, inversed
    MCSignal* signalbb4;
    MCProng prongbmu4A(4, {13, 402, 402, 502}, {true, true, true, true}, {false, false, false, false}, {0, 0, 0, 0}, {0, 0, 0, 0}, {false, false, false, false});
    MCProng prongbmu4B(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongbmu4A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu4B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb4 = new MCSignal("bcbTomumusingle", "mumu pairs from b->e and b->c->e (single b), inversed", {prongbmu4A, prongbmu4B}, {3, 1}); // signal at pair level
                                                                                                                                            //
    MCSignal* signalbb5;
    MCProng prongbmu5(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prongbmu5.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb5 = new MCSignal("bcbcTomumu", "mumu pairs from b->c->mu and b->c->mu", {prongbmu5, prongbmu5}, {-1, -1});

    MCSignal* signalbb6;
    MCProng prongbmu6A(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); // check if mother pdg code is in history
    MCProng prongbmu6B(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prongbmu6A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu6B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb6 = new MCSignal("bbcTomumu", "mumu pairs from b->mu and b->c->mu", {prongbmu6A, prongbmu6B}, {-1, -1}); // signal at pair level
                                                                                                                    
    MCSignal* signalbb7;
    MCProng prongbmu7A(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    MCProng prongbmu7B(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); // check if mother pdg code is in history
    prongbmu7A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu7B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb7 = new MCSignal("bcbTomumu", "mumu pairs from b->c->mu and b->mu", {prongbmu7A, prongbmu7B}, {-1, -1}); // signal at pair level
                                                                                                                     //
    MCSignal* signalbb8;
    MCProng prongbmu8(1, {13}, {true}, {false}, {0}, {0}, {false}, false, {502}, {false}); // check if mother pdg code is in history
    prongbmu8.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb8 = new MCSignal("otherbbTomumu", "mumu pairs from bb", {prongbmu8, prongbmu8}, {-1, -1}); // signal at pair level


    MCSignal* signalmufromphi;
    MCProng mufromphiprong(2, {13, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromphi = new MCSignal("signalmufromphi", "Primary Muons", {mufromphiprong, mufromphiprong}, {-1, -1});

    MCSignal* signalmufromomega;
    MCProng mufromomegaprong(2, {13, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromomega = new MCSignal("signalmufromomega", "Primary Muons", {mufromomegaprong, mufromomegaprong}, {-1, -1});

    MCSignal* signalmufrometa;
    MCProng mufrometaprong(2, {13, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufrometa = new MCSignal("signalmufrometa", "Primary Muons", {mufrometaprong, mufrometaprong}, {-1, -1});

    MCSignal* signalmufromrho;
    MCProng mufromrhoprong(2, {13, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromrho = new MCSignal("signalmufromrho", "Primary Muons", {mufromrhoprong, mufromrhoprong}, {-1, -1});

    MCSignal* signalmufromjpsi;
    MCProng mufromjpsiprong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromjpsi = new MCSignal("signalmufromjpsi", "Primary Muons", {mufromjpsiprong, mufromjpsiprong}, {-1, -1});

    for (auto& collision : collisions) {
      //auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      //grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      //mMagField = grpmag->getNominalL3Field();
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& [fwdtrackId1, fwdtrackId2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(muonIdsThisCollision, muonIdsThisCollision))) {
        for (auto& fwdtrack1 : fwdtracks) {
          //if (fwdtrack1.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId1.fwdtrackId() == fwdtrack1.globalIndex() && fwdtrack1.chi2MatchMCHMFT() < 50 && fwdtrack1.compatibleCollIds().size() == 1) {
          if (fwdtrack1.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId1.fwdtrackId() == fwdtrack1.globalIndex()) {
            for (auto& fwdtrack2 : fwdtracks) {
              //if (fwdtrack2.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId2.fwdtrackId() == fwdtrack2.globalIndex() && fwdtrack2.chi2MatchMCHMFT() < 50 && fwdtrack2.compatibleCollIds().size() == 1) {
              if (fwdtrack2.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId2.fwdtrackId() == fwdtrack2.globalIndex()) {
                //calculate mass
                TLorentzVector lv1, lv2, lv;
                lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
                lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
                lv = lv1 + lv2;

                //set TrackParCovFwd for fwdtrack pairs
                double chi21 = fwdtrack1.chi2();
                SMatrix5 tpars1(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.phi(), fwdtrack1.tgl(), fwdtrack1.signed1Pt());
                std::vector<double> v11{fwdtrack1.cXX(), fwdtrack1.cXY(), fwdtrack1.cYY(), fwdtrack1.cPhiX(), fwdtrack1.cPhiY(), fwdtrack1.cPhiPhi(), fwdtrack1.cTglX(), fwdtrack1.cTglY(), fwdtrack1.cTglPhi(), fwdtrack1.cTglTgl(), fwdtrack1.c1PtX(), fwdtrack1.c1PtY(), fwdtrack1.c1PtPhi(), fwdtrack1.c1PtTgl(), fwdtrack1.c1Pt21Pt2()};
                SMatrix55 tcovs1(v11.begin(), v11.end());
                /*
                float detXY1 = fwdtrack1.cXX() * fwdtrack1.cYY() - fwdtrack1.cXY() * fwdtrack1.cXY();
                if (detXY1 <= 0) {
                  continue;
                }
                */
                o2::track::TrackParCovFwd pars1{fwdtrack1.z(), tpars1, tcovs1, chi21};

                double chi22 = fwdtrack2.chi2();
                SMatrix5 tpars2(fwdtrack2.x(), fwdtrack2.y(), fwdtrack2.phi(), fwdtrack2.tgl(), fwdtrack2.signed1Pt());
                std::vector<double> v12{fwdtrack2.cXX(), fwdtrack2.cXY(), fwdtrack2.cYY(), fwdtrack2.cPhiX(), fwdtrack2.cPhiY(), fwdtrack2.cPhiPhi(), fwdtrack2.cTglX(), fwdtrack2.cTglY(), fwdtrack2.cTglPhi(), fwdtrack2.cTglTgl(), fwdtrack2.c1PtX(), fwdtrack2.c1PtY(), fwdtrack2.c1PtPhi(), fwdtrack2.c1PtTgl(), fwdtrack2.c1Pt21Pt2()};
                SMatrix55 tcovs2(v12.begin(), v12.end());
                /*
                float detXY2 = fwdtrack2.cXX() * fwdtrack2.cYY() - fwdtrack2.cXY() * fwdtrack2.cXY();
                if (detXY2 <= 0) {
                  continue;
                }
                */
                o2::track::TrackParCovFwd pars2{fwdtrack2.z(), tpars2, tcovs2, chi22};

                //Get secondary vertex using DCAFitterN
                int procCode = fgFitterTwoProngFwd.process(pars1, pars2);
                double chi2PCA = -999;
                //double VertexingSV = -999;
                double VertexingLxyz = -999;
                if (procCode == 0 ) {
                  continue;
                }
                Vec3D secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
                fgFitterTwoProngFwd.calcPCACovMatrixFlat();
                chi2PCA = fgFitterTwoProngFwd.getChi2AtPCACandidate();
                auto VertexingLxy = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) + (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
                auto VertexingLz = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
                VertexingLxyz = VertexingLxy + VertexingLz;
                VertexingLxy = std::sqrt(VertexingLxy);
                VertexingLz = std::sqrt(VertexingLz);
                VertexingLxyz = std::sqrt(VertexingLxyz);
                //VertexingSV = secondaryVertex[2];

                //Get secondary vertex using KFparticles
                /*
                KFPTrack kfpTrack0 = createKFPFwdTrackFromFwdTrack(fwdtrack1);
                trk0KF = KFParticle(kfpTrack0, -13 * fwdtrack1.sign());
                KFPTrack kfpTrack1 = createKFPFwdTrackFromFwdTrack(fwdtrack2);
                trk1KF = KFParticle(kfpTrack1, -13 * fwdtrack2.sign());
                KFGeoTwoProng.SetConstructMethod(2);
                KFGeoTwoProng.AddDaughter(trk0KF);
                KFGeoTwoProng.AddDaughter(trk1KF);
                
                KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
                KFParticle KFPV(kfpVertex);
                double dxPair2PV = KFGeoTwoProng.GetX() - KFPV.GetX();
                double dyPair2PV = KFGeoTwoProng.GetY() - KFPV.GetY();
                double dzPair2PV = KFGeoTwoProng.GetZ() - KFPV.GetZ();
                auto KFVertexingLxy = std::sqrt(dxPair2PV * dxPair2PV + dyPair2PV * dyPair2PV);
                auto KFVertexingLz = std::sqrt(dzPair2PV * dzPair2PV);
                auto KFVertexingLxyz = std::sqrt(dxPair2PV * dxPair2PV + dyPair2PV * dyPair2PV + dzPair2PV * dzPair2PV);
                auto KFVertexingSV = KFGeoTwoProng.GetZ();
                auto kKFDCAxyzBetweenProngs = trk0KF.GetDistanceFromParticle(trk1KF);
                */

                //fwdtrack propagation to collision
                //pars1.propagateToZlinear(collision.posZ());
                //pars2.propagateToZlinear(collision.posZ());
                //double centerMFT[3] = {0, 0, -61.4};
                //o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
                //auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
                auto Bz = -5.;

                o2::dataformats::GlobalFwdTrack propmuon1;
                auto geoMan1 = o2::base::GeometryManager::meanMaterialBudget(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.z(), collision.posX(), collision.posY(), collision.posZ());
                auto x2x01 = static_cast<float>(geoMan1.meanX2X0);
                pars1.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x01);
                propmuon1.setParameters(pars1.getParameters());
                propmuon1.setZ(pars1.getZ());
                propmuon1.setCovariances(pars1.getCovariances());

                o2::dataformats::GlobalFwdTrack propmuon2;
                auto geoMan2 = o2::base::GeometryManager::meanMaterialBudget(fwdtrack2.x(), fwdtrack2.y(), fwdtrack2.z(), collision.posX(), collision.posY(), collision.posZ());
                auto x2x02 = static_cast<float>(geoMan2.meanX2X0);
                pars2.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x02);
                propmuon2.setParameters(pars2.getParameters());
                propmuon2.setZ(pars2.getZ());
                propmuon2.setCovariances(pars2.getCovariances());

                //calculate pair DCA
                //float fwd1dcaX = (pars1.getX() - collision.posX());
                //float fwd1dcaY = (pars1.getY() - collision.posY());
                //float fwd2dcaX = (pars2.getX() - collision.posX());
                //float fwd2dcaY = (pars2.getY() - collision.posY());
                float fwd1dcaX = (propmuon1.getX() - collision.posX());
                float fwd1dcaY = (propmuon1.getY() - collision.posY());
                float fwd2dcaX = (propmuon2.getX() - collision.posX());
                float fwd2dcaY = (propmuon2.getY() - collision.posY());
                float DCA1 = std::sqrt(fwd1dcaX * fwd1dcaX + fwd1dcaY * fwd1dcaY);
                float DCA2 = std::sqrt(fwd2dcaX * fwd2dcaX + fwd2dcaY * fwd2dcaY);
                float DCAmumu = std::sqrt((DCA1 * DCA1 + DCA2 * DCA2)/2);

                //get MC particles
                if (fwdtrack1.has_mcParticle() == 1 && fwdtrack2.has_mcParticle() == 1) {
                  auto fwdparticle1 = fwdtrack1.mcParticle();
                  auto fwdparticle2 = fwdtrack2.mcParticle();
                  const auto mcmothersidlist1 = fwdparticle1.mothersIds();
                  const auto mcmothersidlist2 = fwdparticle2.mothersIds();
                  if (fwdparticle1.globalIndex() != fwdparticle2.globalIndex()) {
                    //svcut
                    //if(VertexingSV > -20 || VertexingSV == -999) {
                      if (fwdtrack1.signed1Pt()*fwdtrack2.signed1Pt() < 0) {
                        //fill mass vs DCA
                        registry.fill(HIST("hMassvsPtPM"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAPM"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2PM"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsLxyzPM"), lv.M(), VertexingLxyz);
                        /*
                        if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                          registry.fill(HIST("hMassvsPtPtcutPM"), lv.M(), lv.Pt());
                        }
                        */
                        //select muons from HF semileptonic decays
                        if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM1"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM1"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM1"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM1"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM1"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM2"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM2"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM2"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM2"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM2"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM3"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM3"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM3"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM3"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM3"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM4"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM4"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM4"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM4"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM4"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM5"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM5"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM5"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM5"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM5"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM6"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM6"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM6"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM6"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM6"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM7"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM7"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM7"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM7"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM7"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPM8"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPM8"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPM8"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPM8"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPM8"), lv.M(), lv.Pt());
                          }
                          */
                        } 
                        if (signalmufromphi->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromomega->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufrometa->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromrho->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromjpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtvecPM"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAvecPM"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2vecPM"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzvecPM"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutvecPM"), lv.M(), lv.Pt());
                          }
                          */
                        }
                      } else if (fwdtrack1.signed1Pt() >= 0 && fwdtrack2.signed1Pt() >= 0) {
                        //fill mass vs DCA
                        registry.fill(HIST("hMassvsPtPP"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAPP"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2PP"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsLxyzPP"), lv.M(), VertexingLxyz);
                        /*
                        if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                          registry.fill(HIST("hMassvsPtPtcutPP"), lv.M(), lv.Pt());
                        }
                        */
                        //select muons from HF semileptonic decays
                        if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP1"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP1"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP1"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP1"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP1"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP2"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP2"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP2"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP2"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP2"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP3"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP3"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP3"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP3"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP3"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP4"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP4"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP4"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP4"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP4"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP5"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP5"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP5"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP5"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP5"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP6"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP6"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP6"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP6"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP6"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP7"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP7"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP7"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP7"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP7"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbPP8"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbPP8"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbPP8"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbPP8"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbPP8"), lv.M(), lv.Pt());
                          }
                          */
                        }
                        if (signalmufromphi->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromomega->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufrometa->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromrho->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromjpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtvecPP"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAvecPP"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2vecPP"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzvecPP"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutvecPP"), lv.M(), lv.Pt());
                          }
                          */
                        }
                      } else if (fwdtrack1.signed1Pt() < 0 && fwdtrack2.signed1Pt() < 0) {
                        //fill mass vs DCA
                        registry.fill(HIST("hMassvsPtMM"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAMM"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2MM"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsLxyzMM"), lv.M(), VertexingLxyz);
                        /*
                        if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                          registry.fill(HIST("hMassvsPtPtcutMM"), lv.M(), lv.Pt());
                        }
                        */
                        //select muons from HF semileptonic decays
                        if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM1"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM1"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM1"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM1"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM1"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM2"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM2"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM2"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM2"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM2"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM3"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM3"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM3"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM3"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM3"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM4"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM4"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM4"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM4"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM4"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM5"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM5"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM5"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM5"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM5"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM6"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM6"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM6"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM6"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM6"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM7"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM7"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM7"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM7"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM7"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalbb8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtbbMM8"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAbbMM8"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2bbMM8"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzbbMM8"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutbbMM8"), lv.M(), lv.Pt());
                          }
                          */
                        }
                        if (signalmufromphi->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromomega->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufrometa->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromrho->CheckSignal(true, fwdparticle1, fwdparticle2) || signalmufromjpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtvecMM"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAvecMM"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2vecMM"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzvecMM"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutvecMM"), lv.M(), lv.Pt());
                          }
                          */
                        }
                      }
                    //}
                  }
                }
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
    adaptAnalysisTask<mybbanalysis>(cfgc)
  };
}
