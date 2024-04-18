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

struct myhfanalysis {

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
      {"hMassvsPtPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},

      {"hMassvsDCAPM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAMM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPP", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},

      {"hMassvsDCAPt0PM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0MM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0PP", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},

      {"hMassvsDCAPt1PM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1MM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1PP", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},

      {"hMassvsDCAPt2PM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2MM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2PP", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},

      {"hMassvsLxyzPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzbbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsLxyzccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Lxyz", {HistType::kTH2F, {{750, 0.0, 15.0}, {300, 0, 30}}}},
      {"hMassvsChi2PM", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2MM", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2PP", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPP", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPP1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM1", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPP2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM2", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPP3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM3", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPP4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM4", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPP5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM5", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPP6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM6", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPP7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM7", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsSVZPM", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZMM", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZPP", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPM", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbMM", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPP", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPP1", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPM1", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbMM1", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPP2", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPM2", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbMM2", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPP3", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPM3", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbMM3", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPP4", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPM4", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbMM4", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPP5", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPM5", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbMM5", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPP6", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPM6", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbMM6", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPP7", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbPM7", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZbbMM7", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPM", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccMM", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPP", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPM1", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccMM1", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPP1", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPM2", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccMM2", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPP2", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPM3", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccMM3", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPP3", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPM4", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccMM4", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPP4", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPM5", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccMM5", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPP5", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPM6", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccMM6", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
      {"hMassvsSVZccPP6", "Secondary vertex", {HistType::kTH2F, {{750, 0.0, 15.0}, {1100, -55, 55}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, soa::Filtered<MyMuons> const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const& particles)
  {
    o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;

    bool propagateToPCA = true;
    float maxR = 200;
    float minParamChange = 1.0e-3;
    float minRelChi2Change = 0.9;
    bool useAbsDCA = false;

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

    //Configuring signals cc -> mumu
    MCSignal* signalcc;
    MCProng prongcmu(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc = new MCSignal("ccTomumu", "mumu pairs from c->mu and c->mu", {prongcmu, prongcmu}, {-1, -1});
    //Configuring signals D+D- -> mumu
    MCSignal* signalcc1;
    MCProng prongcmu1(2, {13, 411}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu1.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc1 = new MCSignal("D+D-Tomumu", "mumu pairs from D+->mu and D-->mu", {prongcmu1, prongcmu1}, {-1, -1});
    //Configuring signals D0D0 -> mumu
    MCSignal* signalcc2;
    MCProng prongcmu2(2, {13, 421}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu2.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc2 = new MCSignal("D0D0Tomumu", "mumu pairs from D0->mu and D0->mu", {prongcmu2, prongcmu2}, {-1, -1});

    //Configuring signals D+D0 -> mumu
    MCSignal* signalcc3;
    MCProng prongcmu3A(2, {13, 411}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongcmu3B(2, {13, 421}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu3A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu3B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc3 = new MCSignal("D+D0Tomumu", "mumu pairs from D+->mu and D0->mu", {prongcmu3A, prongcmu3B}, {-1, -1});
    //Configuring signals D+D0 -> mumu, inversed
    MCSignal* signalcc4;
    MCProng prongcmu4A(2, {13, 421}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongcmu4B(2, {13, 411}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu4A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu4B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc4 = new MCSignal("D0D+Tomumu", "mumu pairs from D+->mu and D0->mu, inversed", {prongcmu4A, prongcmu4B}, {-1, -1});
    //Configuring signals D-lambdac -> mumu
    MCSignal* signalcc5;
    MCProng prongcmu5A(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongcmu5B(2, {13, 4122}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu5A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu5B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc5 = new MCSignal("DLambdacTomumu", "mumu pairs from D->mu and Lambda+c->mu", {prongcmu5A, prongcmu5B}, {-1, -1});
    //Configuring signals D-lambdac -> mumu, inversed
    MCSignal* signalcc6;
    MCProng prongcmu6A(2, {13, 4122}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongcmu6B(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu6A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu6B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc6 = new MCSignal("LambdacDTomumu", "mumu pairs from D->mu and Lambda+c->mu, inversed", {prongcmu6A, prongcmu6B}, {-1, -1});


    for (auto& collision : collisions) {
      //auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      //grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      //mMagField = grpmag->getNominalL3Field();
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& [fwdtrackId1, fwdtrackId2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(muonIdsThisCollision, muonIdsThisCollision))) {
        for (auto& fwdtrack1 : fwdtracks) {
          if (fwdtrack1.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId1.fwdtrackId() == fwdtrack1.globalIndex() && fwdtrack1.chi2MatchMCHMFT() < 50 && fwdtrack1.compatibleCollIds().size() == 1) {
            for (auto& fwdtrack2 : fwdtracks) {
              if (fwdtrack2.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId2.fwdtrackId() == fwdtrack2.globalIndex() && fwdtrack2.chi2MatchMCHMFT() < 50 && fwdtrack2.compatibleCollIds().size() == 1) {
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
                double VertexingSV = -999;
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
                VertexingSV = secondaryVertex[2];

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
                pars1.propagateToZlinear(collision.posZ());
                pars2.propagateToZlinear(collision.posZ());

                //calculate pair DCA
                float fwd1dcaX = (pars1.getX() - collision.posX());
                float fwd1dcaY = (pars1.getY() - collision.posY());
                float fwd2dcaX = (pars2.getX() - collision.posX());
                float fwd2dcaY = (pars2.getY() - collision.posY());
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
                    if (fwdtrack1.signed1Pt()*fwdtrack2.signed1Pt() < 0) {
                      //fill mass vs DCA
                      registry.fill(HIST("hMassvsPtPM"), lv.M(), lv.Pt());
                      registry.fill(HIST("hMassvsDCAPM"), lv.M(), DCAmumu);
                      registry.fill(HIST("hMassvsLxyzPM"), lv.M(), VertexingLxyz);
                      registry.fill(HIST("hMassvsChi2PM"), lv.M(), chi2PCA);
                      registry.fill(HIST("hMassvsSVZPM"), lv.M(), VertexingSV);
                      if (lv.Pt() >= 0 && lv.Pt() < 1) {
                        registry.fill(HIST("hMassvsDCAPt0PM"), lv.M(), DCAmumu);
                      } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                        registry.fill(HIST("hMassvsDCAPt1PM"), lv.M(), DCAmumu);
                      } else if (lv.Pt() >= 2) {
                        registry.fill(HIST("hMassvsDCAPt2PM"), lv.M(), DCAmumu);
                      }
                      //select muons from HF semileptonic decays
                      if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPM"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPM"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPM"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM1"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM1"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPM1"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPM1"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPM1"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM1"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM1"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM1"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM2"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM2"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPM2"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPM2"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPM2"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM2"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM2"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM2"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM3"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM3"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPM3"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPM3"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPM3"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM3"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM3"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM3"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM4"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM4"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPM4"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPM4"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPM4"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM4"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM4"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM4"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM5"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM5"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPM5"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPM5"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPM5"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM5"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM5"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM5"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM6"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM6"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPM6"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPM6"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPM6"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM6"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM6"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM6"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM7"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM7"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPM7"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPM7"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPM7"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM7"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzccPM"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2ccPM"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZccPM"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM"), lv.M(), DCAmumu);
                        }
                        if (signalcc1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM1"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM1"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPM1"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPM1"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPM1"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPM1"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPM1"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPM1"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM2"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM2"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPM2"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPM2"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPM2"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPM2"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPM2"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPM2"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM3"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM3"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPM3"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPM3"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPM3"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPM3"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPM3"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPM3"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM4"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM4"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPM4"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPM4"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPM4"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPM4"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPM4"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPM4"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM5"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM5"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPM5"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPM5"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPM5"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPM5"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPM5"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPM5"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM6"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM6"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPM6"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPM6"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPM6"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPM6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPM6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPM6"), lv.M(), DCAmumu);
                          }
                        }
                      }
                    } else if (fwdtrack1.signed1Pt() >= 0 && fwdtrack2.signed1Pt() > 0) {
                      //fill mass vs DCA
                      registry.fill(HIST("hMassvsPtPP"), lv.M(), lv.Pt());
                      registry.fill(HIST("hMassvsDCAPP"), lv.M(), DCAmumu);
                      registry.fill(HIST("hMassvsLxyzPP"), lv.M(), VertexingLxyz);
                      registry.fill(HIST("hMassvsChi2PP"), lv.M(), chi2PCA);
                      registry.fill(HIST("hMassvsSVZPP"), lv.M(), VertexingSV);
                      if (lv.Pt() >= 0 && lv.Pt() < 1) {
                        registry.fill(HIST("hMassvsDCAPt0PP"), lv.M(), DCAmumu);
                      } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                        registry.fill(HIST("hMassvsDCAPt1PP"), lv.M(), DCAmumu);
                      } else if (lv.Pt() >= 2) {
                        registry.fill(HIST("hMassvsDCAPt2PP"), lv.M(), DCAmumu);
                      }
                      //select muons from HF semileptonic decays
                      if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPP"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPP"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPP"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP1"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP1"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPP1"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPP1"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPP1"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP1"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP1"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP1"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP2"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP2"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPP2"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPP2"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPP2"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP2"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP2"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP2"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP3"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP3"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPP3"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPP3"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPP3"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP3"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP3"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP3"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP4"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP4"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPP4"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPP4"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPP4"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP4"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP4"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP4"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP5"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP5"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPP5"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPP5"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPP5"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP5"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP5"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP5"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP6"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP6"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPP6"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPP6"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPP6"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP6"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP6"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP6"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP7"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP7"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbPP7"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbPP7"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbPP7"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP7"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzccPP"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2ccPP"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZccPP"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP"), lv.M(), DCAmumu);
                        }
                        if (signalcc1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP1"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP1"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPP1"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPP1"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPP1"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPP1"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPP1"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPP1"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP2"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP2"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPP2"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPP2"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPP2"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPP2"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPP2"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPP2"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP3"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP3"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPP3"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPP3"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPP3"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPP3"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPP3"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPP3"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP4"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP4"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPP4"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPP4"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPP4"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPP4"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPP4"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPP4"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP5"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP5"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPP5"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPP5"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPP5"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPP5"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPP5"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPP5"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP6"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP6"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccPP6"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccPP6"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccPP6"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPP6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPP6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPP6"), lv.M(), DCAmumu);
                          }
                        }
                      }
                    } else if (fwdtrack1.signed1Pt() < 0 && fwdtrack2.signed1Pt() < 0) {
                      //fill mass vs DCA
                      registry.fill(HIST("hMassvsPtMM"), lv.M(), lv.Pt());
                      registry.fill(HIST("hMassvsDCAMM"), lv.M(), DCAmumu);
                      registry.fill(HIST("hMassvsLxyzMM"), lv.M(), VertexingLxyz);
                      registry.fill(HIST("hMassvsChi2MM"), lv.M(), chi2PCA);
                      registry.fill(HIST("hMassvsSVZMM"), lv.M(), VertexingSV);
                      if (lv.Pt() >= 0 && lv.Pt() < 1) {
                        registry.fill(HIST("hMassvsDCAPt0MM"), lv.M(), DCAmumu);
                      } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                        registry.fill(HIST("hMassvsDCAPt1MM"), lv.M(), DCAmumu);
                      } else if (lv.Pt() >= 2) {
                        registry.fill(HIST("hMassvsDCAPt2MM"), lv.M(), DCAmumu);
                      }
                      //select muons from HF semileptonic decays
                      if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbMM"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbMM"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbMM"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM1"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM1"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbMM1"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbMM1"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbMM1"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM1"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM1"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM1"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM2"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM2"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbMM2"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbMM2"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbMM2"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM2"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM2"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM2"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM3"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM3"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbMM3"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbMM3"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbMM3"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM3"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM3"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM3"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM4"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM4"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbMM4"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbMM4"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbMM4"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM4"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM4"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM4"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM5"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM5"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbMM5"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbMM5"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbMM5"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM5"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM5"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM5"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM6"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM6"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbMM6"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbMM6"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbMM6"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM6"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM6"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM6"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM7"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM7"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzbbMM7"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2bbMM7"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZbbMM7"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM7"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsLxyzccMM"), lv.M(), VertexingLxyz);
                        registry.fill(HIST("hMassvsChi2ccMM"), lv.M(), chi2PCA);
                        registry.fill(HIST("hMassvsSVZccMM"), lv.M(), VertexingSV);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM"), lv.M(), DCAmumu);
                        }
                        if (signalcc1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM1"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM1"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccMM1"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccMM1"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccMM1"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccMM1"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccMM1"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccMM1"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM2"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM2"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccMM2"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccMM2"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccMM2"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccMM2"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccMM2"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccMM2"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM3"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM3"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccMM3"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccMM3"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccMM3"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccMM3"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccMM3"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccMM3"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM4"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM4"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccMM4"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccMM4"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccMM4"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccMM4"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccMM4"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccMM4"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM5"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM5"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccMM5"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccMM5"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccMM5"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccMM5"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccMM5"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccMM5"), lv.M(), DCAmumu);
                          }
                        } else if (signalcc6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM6"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM6"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsLxyzccMM6"), lv.M(), VertexingLxyz);
                          registry.fill(HIST("hMassvsChi2ccMM6"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsSVZccMM6"), lv.M(), VertexingSV);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccMM6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccMM6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccMM6"), lv.M(), DCAmumu);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      auto mccolision = collision.mcCollision();
      auto particlessThisCollision = particles.sliceBy(particlesIndicesPerCollision, collision.globalIndex());
      for (auto& [fwdparticle1, fwdparticle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(particlessThisCollision, particlessThisCollision))) {
        if (fwdparticle1.pdgCode() == 13 || fwdparticle1.pdgCode() == -13 ) {
          if (fwdparticle2.pdgCode() == 13 || fwdparticle2.pdgCode() == -13 ) {
            //calculate mass
            if (fwdparticle1.eta() > -4 && fwdparticle1.eta() < -2.5) {
              if (fwdparticle2.eta() > -4 && fwdparticle2.eta() < -2.5) {
                TLorentzVector lv1, lv2, lv;
                lv1.SetPtEtaPhiM(fwdparticle1.pt(), fwdparticle1.eta(), fwdparticle1.phi(), muonMass);
                lv2.SetPtEtaPhiM(fwdparticle2.pt(), fwdparticle2.eta(), fwdparticle2.phi(), muonMass);
                lv = lv1 + lv2;

                if (fwdparticle1.globalIndex() != fwdparticle2.globalIndex()) {
                  if (fwdparticle1.pdgCode()*fwdparticle2.pdgCode() < 0) {
                    //fill mass vs DCA
                    registry.fill(HIST("hGenMassvsPtPM"), lv.M(), lv.Pt());
                    //select muons from HF semileptonic decays
                    if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM"), lv.M(), lv.Pt());
                    } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM1"), lv.M(), lv.Pt());
                    } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM2"), lv.M(), lv.Pt());
                    } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM3"), lv.M(), lv.Pt());
                    } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM4"), lv.M(), lv.Pt());
                    } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM5"), lv.M(), lv.Pt());
                    } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM6"), lv.M(), lv.Pt());
                    } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM7"), lv.M(), lv.Pt());
                    } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtccPM"), lv.M(), lv.Pt());
                    }
                  } else if (fwdparticle1.pdgCode() > 0 && fwdparticle2.pdgCode() > 0) {
                    //fill mass vs DCA
                    registry.fill(HIST("hGenMassvsPtPP"), lv.M(), lv.Pt());
                    //select muons from HF semileptonic decays
                    if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP"), lv.M(), lv.Pt());
                    } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP1"), lv.M(), lv.Pt());
                    } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP2"), lv.M(), lv.Pt());
                    } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP3"), lv.M(), lv.Pt());
                    } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP4"), lv.M(), lv.Pt());
                    } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP5"), lv.M(), lv.Pt());
                    } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP6"), lv.M(), lv.Pt());
                    } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP7"), lv.M(), lv.Pt());
                    } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtccPP"), lv.M(), lv.Pt());
                    }
                  } else if (fwdparticle1.pdgCode() < 0 && fwdparticle2.pdgCode() < 0) {
                    //fill mass vs DCA
                    registry.fill(HIST("hGenMassvsPtMM"), lv.M(), lv.Pt());
                    //select muons from HF semileptonic decays
                    if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM"), lv.M(), lv.Pt());
                    } else if (signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM1"), lv.M(), lv.Pt());
                    } else if (signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM2"), lv.M(), lv.Pt());
                    } else if (signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM3"), lv.M(), lv.Pt());
                    } else if (signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM4"), lv.M(), lv.Pt());
                    } else if (signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM5"), lv.M(), lv.Pt());
                    } else if (signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM6"), lv.M(), lv.Pt());
                    } else if (signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM7"), lv.M(), lv.Pt());
                    } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtccMM"), lv.M(), lv.Pt());
                    }
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
    adaptAnalysisTask<myhfanalysis>(cfgc)
  };
}
