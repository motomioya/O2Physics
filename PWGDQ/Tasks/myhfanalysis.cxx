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
      {"hMassvsPtbbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtbbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
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
      {"hMassvsPtccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hMassvsPtccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},

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
      {"hMassvsDCAbbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAbbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
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
      {"hMassvsDCAccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},

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
      {"hMassvsDCAPt0bbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0bbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
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
      {"hMassvsDCAPt0ccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt0ccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},

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
      {"hMassvsDCAPt1bbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1bbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
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
      {"hMassvsDCAPt1ccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt1ccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},

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
      {"hMassvsDCAPt2bbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2bbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
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
      {"hMassvsDCAPt2ccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},
      {"hMassvsDCAPt2ccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {600, 0, 3}}}},

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
      {"hMassvsChi2bbPP8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbPM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2bbMM8", "Invariant;Invariant Mass (GeV/c^{2}) from bb;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
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
      {"hMassvsChi2ccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
      {"hMassvsChi2ccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 100}}}},
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
                                                                                                                     //
    MCSignal* signalbb8;
    MCProng prongbmu8(1, {13}, {true}, {false}, {0}, {0}, {false}, false, {502}, {false}); // check if mother pdg code is in history
    prongbmu8.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb8 = new MCSignal("otherbbTomumu", "mumu pairs from bb", {prongbmu8, prongbmu8}, {-1, -1}); // signal at pair level

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

    // -----
    //Configuring signals b->cmu->mumu
    MCSignal* signalcc7;
    MCProng prongcmu7A(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongcmu7B(3, {13, 300, 402}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prongcmu7A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu7B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc7 = new MCSignal("ccsTomumusingle", "mumu pairs from b->e and b->c->e (single b)", {prongcmu7A, prongcmu7B}, {1, 2}); // signal at pair level
    //Configuring signals b->cmu->mumu, inversed
    MCSignal* signalcc8;
    MCProng prongcmu8A(3, {13, 300, 402}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongcmu8B(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu8A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu8B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc8 = new MCSignal("cscTomumusingle", "mumu pairs from b->e and b->c->e (single b), inversed", {prongcmu8A, prongcmu8B}, {2, 1}); // signal at pair level
    //Configuring signals b->cmu->mumu
    MCSignal* signalcc9;
    MCProng prongcmu9A(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongcmu9B(4, {13, 300, 300, 402}, {true, true, true, true}, {false, false, false, false}, {0, 0, 0, 0}, {0, 0, 0, 0}, {false, false, false, false});
    prongcmu9A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu9B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc9 = new MCSignal("ccssTomumusingle", "mumu pairs from b->e and b->c->e (single b)", {prongcmu9A, prongcmu9B}, {1, 3}); // signal at pair level
    //Configuring signals b->cmu->mumu, inversed
    MCSignal* signalcc10;
    MCProng prongcmu10A(4, {13, 300, 300, 402}, {true, true, true, true}, {false, false, false, false}, {0, 0, 0, 0}, {0, 0, 0, 0}, {false, false, false, false});
    MCProng prongcmu10B(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu10A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu10B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc10 = new MCSignal("csscTomumusingle", "mumu pairs from b->e and b->c->e (single b), inversed", {prongcmu10A, prongcmu10B}, {3, 1}); // signal at pair level
                                                                                                                                            //
    MCSignal* signalcc11;
    MCProng prongcmu11(2, {13, 300}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {402}, {false}); // check if mother pdg code is in history
    prongcmu11.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc11 = new MCSignal("cscsTomumu", "mumu pairs from b->c->mu and b->c->mu", {prongcmu11, prongcmu11}, {-1, -1});

    MCSignal* signalcc12;
    MCProng prongcmu12A(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); // check if mother pdg code is in history
    MCProng prongcmu12B(2, {13, 300}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {402}, {false}); // check if mother pdg code is in history
    prongcmu12A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu12B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc12 = new MCSignal("ccsTomumu", "mumu pairs from b->mu and b->c->mu", {prongcmu12A, prongcmu12B}, {-1, -1}); // signal at pair level
                                                                                                                    
    MCSignal* signalcc13;
    MCProng prongcmu13A(2, {13, 300}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {402}, {false}); // check if mother pdg code is in history
    MCProng prongcmu13B(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); // check if mother pdg code is in history
    prongcmu13A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongcmu13B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc13 = new MCSignal("cscTomumu", "mumu pairs from b->c->mu and b->mu", {prongcmu13A, prongcmu13B}, {-1, -1}); // signal at pair level
                                                                                                                     //
    MCSignal* signalcc14;
    MCProng prongcmu14(1, {13}, {true}, {false}, {0}, {0}, {false}, false, {402}, {false}); // check if mother pdg code is in history
    prongcmu14.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc14 = new MCSignal("otherccTomumu", "mumu pairs from bb", {prongcmu14, prongcmu14}, {-1, -1}); // signal at pair level


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
                      registry.fill(HIST("hMassvsChi2PM"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPM"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPM1"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPM2"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPM3"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPM4"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPM5"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPM6"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPM7"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM7"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPM8"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPM8"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2bbPM8"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPM8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPM8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPM8"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPM1"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPM2"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPM3"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPM4"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPM5"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPM6"), lv.M(), chi2PCA);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPM6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPM6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPM6"), lv.M(), DCAmumu);
                          }
                        }
                      } else if (signalcc7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM7"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM7"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM7"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM7"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM8"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM8"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM8"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM8"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc9->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM9"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM9"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM9"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM9"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM9"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM9"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc10->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM10"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM10"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM10"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM10"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM10"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM10"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc11->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM11"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM11"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM11"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM11"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM11"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM11"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc12->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM12"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM12"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM12"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM12"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM12"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM12"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc13->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM13"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM13"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM13"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM13"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM13"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM13"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc14->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPM14"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPM14"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPM14"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPM14"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPM14"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPM14"), lv.M(), DCAmumu);
                        }
                      }
                    } else if (fwdtrack1.signed1Pt() >= 0 && fwdtrack2.signed1Pt() >= 0) {
                      //fill mass vs DCA
                      registry.fill(HIST("hMassvsPtPP"), lv.M(), lv.Pt());
                      registry.fill(HIST("hMassvsDCAPP"), lv.M(), DCAmumu);
                      registry.fill(HIST("hMassvsChi2PP"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPP"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPP1"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPP2"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPP3"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPP4"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPP5"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPP6"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbPP7"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP7"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbPP8"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbPP8"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2bbPP8"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbPP8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbPP8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbPP8"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPP1"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPP2"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPP3"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPP4"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPP5"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccPP6"), lv.M(), chi2PCA);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccPP6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccPP6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccPP6"), lv.M(), DCAmumu);
                          }
                        }
                      } else if (signalcc7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP7"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP7"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP7"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP7"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP8"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP8"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP8"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP8"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc9->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP9"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP9"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP9"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP9"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP9"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP9"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc10->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP10"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP10"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP10"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP10"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP10"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP10"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc11->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP11"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP11"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP11"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP11"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP11"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP11"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc12->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP12"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP12"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP12"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP12"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP12"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP12"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc13->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP13"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP13"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP13"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP13"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP13"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP13"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc14->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccPP14"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccPP14"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccPP14"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccPP14"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccPP14"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccPP14"), lv.M(), DCAmumu);
                        }
                      }
                    } else if (fwdtrack1.signed1Pt() < 0 && fwdtrack2.signed1Pt() < 0) {
                      //fill mass vs DCA
                      registry.fill(HIST("hMassvsPtMM"), lv.M(), lv.Pt());
                      registry.fill(HIST("hMassvsDCAMM"), lv.M(), DCAmumu);
                      registry.fill(HIST("hMassvsChi2MM"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbMM"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbMM1"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbMM2"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbMM3"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbMM4"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbMM5"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbMM6"), lv.M(), chi2PCA);
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
                        registry.fill(HIST("hMassvsChi2bbMM7"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM7"), lv.M(), DCAmumu);
                        }
                      } else if (signalbb8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtbbMM8"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAbbMM8"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2bbMM8"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0bbMM8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1bbMM8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2bbMM8"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccMM1"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccMM2"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccMM3"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccMM4"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccMM5"), lv.M(), chi2PCA);
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
                          registry.fill(HIST("hMassvsChi2ccMM6"), lv.M(), chi2PCA);
                          if (lv.Pt() >= 0 && lv.Pt() < 1) {
                            registry.fill(HIST("hMassvsDCAPt0ccMM6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                            registry.fill(HIST("hMassvsDCAPt1ccMM6"), lv.M(), DCAmumu);
                          } else if (lv.Pt() >= 2) {
                            registry.fill(HIST("hMassvsDCAPt2ccMM6"), lv.M(), DCAmumu);
                          }
                        }
                      } else if (signalcc7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM7"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM7"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM7"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM7"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM7"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM8"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM8"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM8"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM8"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM8"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc9->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM9"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM9"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM9"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM9"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM9"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM9"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc10->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM10"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM10"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM10"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM10"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM10"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM10"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc11->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM11"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM11"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM11"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM11"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM11"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM11"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc12->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM12"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM12"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM12"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM12"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM12"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM12"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc13->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM13"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM13"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM13"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM13"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM13"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM13"), lv.M(), DCAmumu);
                        }
                      } else if (signalcc14->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hMassvsPtccMM14"), lv.M(), lv.Pt());
                        registry.fill(HIST("hMassvsDCAccMM14"), lv.M(), DCAmumu);
                        registry.fill(HIST("hMassvsChi2ccMM14"), lv.M(), chi2PCA);
                        if (lv.Pt() >= 0 && lv.Pt() < 1) {
                          registry.fill(HIST("hMassvsDCAPt0ccMM14"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 1 && lv.Pt() < 2) {
                          registry.fill(HIST("hMassvsDCAPt1ccMM14"), lv.M(), DCAmumu);
                        } else if (lv.Pt() >= 2) {
                          registry.fill(HIST("hMassvsDCAPt2ccMM14"), lv.M(), DCAmumu);
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
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myhfanalysis>(cfgc)
  };
}
