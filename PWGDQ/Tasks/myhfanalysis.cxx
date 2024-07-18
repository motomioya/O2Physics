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
      {"hMassvsDCAPM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAMM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAPP", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAvecPM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAvecMM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},
      {"hMassvsDCAvecPP", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Pair DCA", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0.0, 0.1}}}},

      {"hMassvsPtPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtMM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtPP", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},

      {"hMassvsPtccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtvecPM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtvecMM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},
      {"hMassvsPtvecPP", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {750, 0.0, 15.0}}}},

      {"hMassvsLxyzPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzMM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzPP", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzvecPM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzvecMM", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},
      {"hMassvsLxyzvecPP", "Invariant;Invariant Mass (GeV/c^{2}) from vec;Mass", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0, 15}}}},

      {"hMassvsChi2PM", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2MM", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2PP", "Invariant;Invariant Mass (GeV/c^{2});Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP1", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP2", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP3", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP4", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP5", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP6", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP7", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP8", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP9", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP10", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP11", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP12", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP13", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccMM14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
      {"hMassvsChi2ccPP14", "Invariant;Invariant Mass (GeV/c^{2}) from cc;Pair Chi2", {HistType::kTH2F, {{750, 0.0, 15.0}, {1000, 0., 20}}}},
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
                                                                                                          // 

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
                    if(VertexingSV > -20 || VertexingSV == -999) {
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
                        if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM"), lv.M(), lv.Pt());
                          }
                          */
                          if (signalcc1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPM1"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPM1"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPM1"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPM1"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPM1"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPM2"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPM2"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPM2"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPM2"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPM2"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPM3"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPM3"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPM3"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPM3"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPM2"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPM4"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPM4"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPM4"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPM4"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPM4"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPM5"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPM5"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPM5"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPM5"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPM5"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPM6"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPM6"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPM6"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPM6"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPM6"), lv.M(), lv.Pt());
                            }
                            */
                          }
                        } else if (signalcc7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM7"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM7"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM7"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM7"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM7"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM8"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM8"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM8"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM8"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM8"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc9->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM9"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM9"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM9"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM9"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM9"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc10->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM10"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM10"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM10"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM10"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM10"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc11->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM11"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM11"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM11"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM11"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM11"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc12->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM12"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM12"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM12"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM12"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM12"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc13->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM13"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM13"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM13"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM13"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM13"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc14->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPM14"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPM14"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPM14"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPM14"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPM14"), lv.M(), lv.Pt());
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
                       if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP"), lv.M(), lv.Pt());
                          }
                          */
                          if (signalcc1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPP1"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPP1"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPP1"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPP1"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPP1"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPP2"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPP2"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPP2"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPP2"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPP2"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPP3"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPP3"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPP3"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPP3"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPP2"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPP4"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPP4"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPP4"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPP4"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPP4"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPP5"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPP5"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPP5"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPP5"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPP5"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccPP6"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccPP6"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccPP6"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccPP6"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccPP6"), lv.M(), lv.Pt());
                            }
                            */
                          }
                        } else if (signalcc7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP7"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP7"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP7"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP7"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP7"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP8"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP8"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP8"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP8"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP8"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc9->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP9"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP9"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP9"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP9"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP9"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc10->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP10"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP10"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP10"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP10"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP10"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc11->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP11"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP11"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP11"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP11"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP11"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc12->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP12"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP12"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP12"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP12"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP12"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc13->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP13"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP13"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP13"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP13"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP13"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc14->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccPP14"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccPP14"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccPP14"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccPP14"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccPP14"), lv.M(), lv.Pt());
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
                       if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM"), lv.M(), lv.Pt());
                          }
                          */
                          if (signalcc1->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccMM1"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccMM1"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccMM1"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccMM1"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccMM1"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc2->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccMM2"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccMM2"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccMM2"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccMM2"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccMM2"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc3->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccMM3"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccMM3"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccMM3"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccMM3"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccMM2"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc4->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccMM4"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccMM4"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccMM4"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccMM4"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccMM4"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc5->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccMM5"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccMM5"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccMM5"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccMM5"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccMM5"), lv.M(), lv.Pt());
                            }
                            */
                          } else if (signalcc6->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                            registry.fill(HIST("hMassvsPtccMM6"), lv.M(), lv.Pt());
                            registry.fill(HIST("hMassvsDCAccMM6"), lv.M(), DCAmumu);
                            registry.fill(HIST("hMassvsChi2ccMM6"), lv.M(), chi2PCA);
                            registry.fill(HIST("hMassvsLxyzccMM6"), lv.M(), VertexingLxyz);
                            /*
                            if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                              registry.fill(HIST("hMassvsPtPtcutccMM6"), lv.M(), lv.Pt());
                            }
                            */
                          }
                        } else if (signalcc7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM7"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM7"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM7"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM7"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM7"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc8->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM8"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM8"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM8"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM8"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM8"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc9->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM9"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM9"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM9"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM9"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM9"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc10->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM10"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM10"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM10"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM10"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM10"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc11->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM11"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM11"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM11"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM11"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM11"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc12->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM12"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM12"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM12"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM12"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM12"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc13->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM13"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM13"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM13"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM13"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM13"), lv.M(), lv.Pt());
                          }
                          */
                        } else if (signalcc14->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hMassvsPtccMM14"), lv.M(), lv.Pt());
                          registry.fill(HIST("hMassvsDCAccMM14"), lv.M(), DCAmumu);
                          registry.fill(HIST("hMassvsChi2ccMM14"), lv.M(), chi2PCA);
                          registry.fill(HIST("hMassvsLxyzccMM14"), lv.M(), VertexingLxyz);
                          /*
                          if (fwdtrack1.pt() >= singleptcut && fwdtrack2.pt() >= singleptcut) {
                            registry.fill(HIST("hMassvsPtPtcutccMM14"), lv.M(), lv.Pt());
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
