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
// p

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
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
#include <math.h>
#include <TLorentzVector.h>
#include <string>
#include <regex>
#include <iostream>
#include <THashList.h>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"

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
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;
namespace o2::aod
{
namespace dqanalysisflags
{
DECLARE_SOA_BITMAP_COLUMN(IsMuonSelected, isMuonSelected, 32);                       //! Muon track decisions (joinable to ReducedMuonsAssoc)
DECLARE_SOA_COLUMN(MuonAmbiguityInBunch, muonAmbiguityInBunch, int8_t);              //! Muon track in-bunch ambiguity
DECLARE_SOA_COLUMN(MuonAmbiguityOutOfBunch, muonAmbiguityOutOfBunch, int8_t);        //! Muon track out of bunch ambiguity
                                                                                                                                  }
DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "DQANAMUONCUTS", dqanalysisflags::IsMuonSelected);                                                       //!  joinable to ReducedMuonsAssoc
DECLARE_SOA_TABLE(MuonAmbiguities, "AOD", "DQMUONAMB", dqanalysisflags::MuonAmbiguityInBunch, dqanalysisflags::MuonAmbiguityOutOfBunch);         //!  joinable to ReducedMuonTracks
                                                                                                                                  }

using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::ReducedMuonsLabels>;
using MyMuonTracksWithCovWithAmbiguities = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::MuonAmbiguities, aod::ReducedMuonsLabels>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyMFTTrack = soa::Join<aod::ReducedMFTs, aod::ReducedMFTsExtra, aod::ReducedMFTLabels>;
o2::base::MatLayerCylSet* lut = nullptr;

struct myQApair {

  float collisionZcut = 10.0f;
  Preslice<aod::McParticles> particlesIndicesPerCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::ReducedMuonsAssoc> fwdtrackIndicesPerCollision = aod::reducedtrack_association::reducedeventId;
  float mMagField = -5.0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  HistogramRegistry registry{
    "registry", 
    {
      {"hDCANoRefitlinearNoHF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromB", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromB", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromB", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHF", "LXYZ;LXYZ (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromD", "LXYZ;LXYZ (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromB", "LXYZ;LXYZ (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFLowPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDLowPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBLowPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFLowPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDLowPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBLowPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFLowPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDLowPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBLowPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFMidPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDMidPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBMidPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFMidPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDMidPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBMidPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFMidPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDMidPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBMidPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFHighPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDHighPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBHighPt", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFHighPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDHighPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBHighPtCorrectMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},

      {"hDCANoRefitlinearNoHFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hLxyzRefithelixNoHFHighPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromDHighPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
      {"hLxyzRefithelixFromBHighPtWrongMatch", "Lxyz;Lxyz (cm)", {HistType::kTH1F, {{10000, 0, 2.0}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(MyEventsVtxCov const& collisions, aod::ReducedMCEvents const&, aod::ReducedMCTracks const&, MyMuonTracksWithCovWithAmbiguities const&, aod::ReducedMuonsAssoc const& assocs, MyMFTTrack const&)
  {
    //muon from pion
    //muon from D
    MCSignal* signalmufromD;
    MCProng mufromDprong(2, {13, 403}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromD = new MCSignal("signalmufromD", "Primary Muons", {mufromDprong}, {-1});
    //muon from B
    MCSignal* signalmufromB;
    MCProng mufromBprong(2, {13, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromB = new MCSignal("signalmufromB", "Primary Muons", {mufromBprong}, {-1});

    //Preslice<aod::ReducedMuonsAssoc> muonAssocsPerCollision = aod::reducedtrack_association::reducedevent;
    //Preslice<aod::ReducedMuonsAssoc> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;
    o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;

    bool propagateToPCA = true;
    float maxR = 200;
    float minParamChange = 1.0e-3;
    float minRelChi2Change = 0.9;
    bool useAbsDCA = true;
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

    if (!o2::base::GeometryManager::isGeometryLoaded())
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");


    for (auto& collision : collisions) {
      auto muonIdsThisCollision = assocs.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& [fwdtrackId1, fwdtrackId2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(muonIdsThisCollision, muonIdsThisCollision))) {
        auto fwdtrack1 = fwdtrackId1.template reducedmuon_as<MyMuonTracksWithCovWithAmbiguities>();
        auto fwdtrack2 = fwdtrackId2.template reducedmuon_as<MyMuonTracksWithCovWithAmbiguities>();
        //Check Ambiguity
        if (fwdtrack1.muonAmbiguityInBunch() > 1) continue;
        if (fwdtrack2.muonAmbiguityOutOfBunch() > 1) continue;

        //calculate DCA1
        float fwd0dcaX1 = fwdtrack1.fwdDcaX();
        float fwd0dcaY1 = fwdtrack1.fwdDcaY();
        float DCA01 = std::sqrt(fwd0dcaX1 * fwd0dcaX1 + fwd0dcaY1 * fwd0dcaY1);
        //set TrackParCovFwd for fwdtrack pairs
        double chi211 = fwdtrack1.chi2();
        SMatrix5 tpars11(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.phi(), fwdtrack1.tgl(), fwdtrack1.signed1Pt());
        std::vector<double> v111{fwdtrack1.cXX(), fwdtrack1.cXY(), fwdtrack1.cYY(), fwdtrack1.cPhiX(), fwdtrack1.cPhiY(), fwdtrack1.cPhiPhi(), fwdtrack1.cTglX(), fwdtrack1.cTglY(), fwdtrack1.cTglPhi(), fwdtrack1.cTglTgl(), fwdtrack1.c1PtX(), fwdtrack1.c1PtY(), fwdtrack1.c1PtPhi(), fwdtrack1.c1PtTgl(), fwdtrack1.c1Pt21Pt2()};
        SMatrix55 tcovs11(v111.begin(), v111.end());
        o2::track::TrackParCovFwd pars11{fwdtrack1.z(), tpars11, tcovs11, chi211};
        pars11.propagateToZlinear(collision.posZ());
        float fwd1dcaX1 = (pars11.getX() - collision.posX());
        float fwd1dcaY1 = (pars11.getY() - collision.posY());
        float DCA11 = std::sqrt(fwd1dcaX1 * fwd1dcaX1 + fwd1dcaY1 * fwd1dcaY1);

        //calculate DCA2
        float fwd0dcaX2 = fwdtrack2.fwdDcaX();
        float fwd0dcaY2 = fwdtrack2.fwdDcaY();
        float DCA02 = std::sqrt(fwd0dcaX2 * fwd0dcaX2 + fwd0dcaY2 * fwd0dcaY2);
        //set TrackParCovFwd for fwdtrack pairs
        double chi212 = fwdtrack2.chi2();
        SMatrix5 tpars12(fwdtrack2.x(), fwdtrack2.y(), fwdtrack2.phi(), fwdtrack2.tgl(), fwdtrack2.signed1Pt());
        std::vector<double> v112{fwdtrack2.cXX(), fwdtrack2.cXY(), fwdtrack2.cYY(), fwdtrack2.cPhiX(), fwdtrack2.cPhiY(), fwdtrack2.cPhiPhi(), fwdtrack2.cTglX(), fwdtrack2.cTglY(), fwdtrack2.cTglPhi(), fwdtrack2.cTglTgl(), fwdtrack2.c1PtX(), fwdtrack2.c1PtY(), fwdtrack2.c1PtPhi(), fwdtrack2.c1PtTgl(), fwdtrack2.c1Pt21Pt2()};
        SMatrix55 tcovs12(v112.begin(), v112.end());
        o2::track::TrackParCovFwd pars12{fwdtrack2.z(), tpars12, tcovs12, chi212};
        pars12.propagateToZlinear(collision.posZ());
        float fwd1dcaX2 = (pars12.getX() - collision.posX());
        float fwd1dcaY2 = (pars12.getY() - collision.posY());
        float DCA12 = std::sqrt(fwd1dcaX2 * fwd1dcaX2 + fwd1dcaY2 * fwd1dcaY2);

        //fwdtrack lenear propagation 1
        o2::track::TrackParCovFwd pars21{fwdtrack1.z(), tpars11, tcovs11, chi211};
        o2::dataformats::GlobalFwdTrack propmuon1;
        auto geoMan1 = o2::base::GeometryManager::meanMaterialBudget(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.z(), collision.posX(), collision.posY(), collision.posZ());
        auto x2x01 = static_cast<float>(geoMan1.meanX2X0);
        pars21.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, mMagField, x2x01);
        propmuon1.setParameters(pars21.getParameters());
        propmuon1.setZ(pars21.getZ());
        propmuon1.setCovariances(pars21.getCovariances());
        //calculate DCA1
        float fwd2dcaX1 = (pars21.getX() - collision.posX());
        float fwd2dcaY1 = (pars21.getY() - collision.posY());
        float DCA21 = std::sqrt(fwd2dcaX1 * fwd2dcaX1 + fwd2dcaY1 * fwd2dcaY1);

        //fwdtrack lenear propagation 2
        o2::track::TrackParCovFwd pars22{fwdtrack2.z(), tpars12, tcovs12, chi212};
        o2::dataformats::GlobalFwdTrack propmuon2;
        auto geoMan2 = o2::base::GeometryManager::meanMaterialBudget(fwdtrack2.x(), fwdtrack2.y(), fwdtrack2.z(), collision.posX(), collision.posY(), collision.posZ());
        auto x2x02 = static_cast<float>(geoMan2.meanX2X0);
        pars22.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, mMagField, x2x02);
        propmuon2.setParameters(pars22.getParameters());
        propmuon2.setZ(pars22.getZ());
        propmuon2.setCovariances(pars22.getCovariances());
        //calculate DCA2
        float fwd2dcaX2 = (pars22.getX() - collision.posX());
        float fwd2dcaY2 = (pars22.getY() - collision.posY());
        float DCA22 = std::sqrt(fwd2dcaX2 * fwd2dcaX2 + fwd2dcaY2 * fwd2dcaY2);

        //Get secondary vertex using DCAFitterN
        int procCode = fgFitterTwoProngFwd.process(pars21, pars22);
        double chi2PCA = -999;
        double VertexingSV = -999;
        double VertexingLxyz = -999;
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


        if (fwdtrack1.has_reducedMCTrack() && fwdtrack1.has_matchMFTTrack() && fwdtrack2.has_reducedMCTrack() && fwdtrack2.has_matchMFTTrack()) {
          auto mfttrack1 = fwdtrack1.matchMFTTrack_as<MyMFTTrack>();
          auto mfttrack2 = fwdtrack2.matchMFTTrack_as<MyMFTTrack>();
          if (mfttrack1.has_reducedMCTrack() && mfttrack2.has_reducedMCTrack()) {
            //auto fwdparticle = fwdtrack.mcParticle();
            auto fwdparticle1 = fwdtrack1.template reducedMCTrack_as<aod::ReducedMCTracks>();
            auto fwdparticle2 = fwdtrack2.template reducedMCTrack_as<aod::ReducedMCTracks>();
            auto mftparticle1 = mfttrack1.template reducedMCTrack_as<aod::ReducedMCTracks>();
            auto mftparticle2 = mfttrack2.template reducedMCTrack_as<aod::ReducedMCTracks>();
            bool isMatchCorrect = 0;
            //calculate mass
            TLorentzVector lv1, lv2, lv;
            lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
            lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
            lv = lv1 + lv2;

            if (lv.M() < 1.4 || lv.M() > 2.8) continue;

            if (fwdparticle1.globalIndex() == mftparticle1.globalIndex() && fwdparticle2.globalIndex() == mftparticle2.globalIndex()) isMatchCorrect = 1;
            //Fill histograms
            if (!signalmufromD->CheckSignal(true, fwdparticle1) && !signalmufromD->CheckSignal(true, fwdparticle2)) {
              if (!signalmufromB->CheckSignal(true, fwdparticle2) && !signalmufromB->CheckSignal(true, fwdparticle2)) {
                registry.fill(HIST("hDCANoRefitlinearNoHF"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearNoHF"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixNoHF"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixNoHF"), VertexingLxyz);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearNoHFCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearNoHFCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixNoHFCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixNoHFCorrectMatch"), VertexingLxyz);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearNoHFWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearNoHFWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixNoHFWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixNoHFWrongMatch"), VertexingLxyz);
                }
                if (lv.Pt() < 1.0) {
                  registry.fill(HIST("hDCANoRefitlinearNoHFLowPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearNoHFLowPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixNoHFLowPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixNoHFLowPt"), VertexingLxyz);
                  if (isMatchCorrect == 1) {
                    registry.fill(HIST("hDCANoRefitlinearNoHFLowPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                    registry.fill(HIST("hDCARefitlinearNoHFLowPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                    registry.fill(HIST("hDCARefithelixNoHFLowPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                    registry.fill(HIST("hLxyzRefithelixNoHFLowPtCorrectMatch"), VertexingLxyz);
                  } else {
                    registry.fill(HIST("hDCANoRefitlinearNoHFLowPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                    registry.fill(HIST("hDCARefitlinearNoHFLowPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                    registry.fill(HIST("hDCARefithelixNoHFLowPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                    registry.fill(HIST("hLxyzRefithelixNoHFLowPtWrongMatch"), VertexingLxyz);
                  }
                } else if (lv.Pt() >= 1.0 && lv.Pt() < 2.0) {
                  registry.fill(HIST("hDCANoRefitlinearNoHFMidPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearNoHFMidPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixNoHFMidPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixNoHFMidPt"), VertexingLxyz);
                  if (isMatchCorrect == 1) {
                    registry.fill(HIST("hDCANoRefitlinearNoHFMidPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                    registry.fill(HIST("hDCARefitlinearNoHFMidPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                    registry.fill(HIST("hDCARefithelixNoHFMidPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                    registry.fill(HIST("hLxyzRefithelixNoHFMidPtCorrectMatch"), VertexingLxyz);
                  } else {
                    registry.fill(HIST("hDCANoRefitlinearNoHFMidPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                    registry.fill(HIST("hDCARefitlinearNoHFMidPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                    registry.fill(HIST("hDCARefithelixNoHFMidPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                    registry.fill(HIST("hLxyzRefithelixNoHFMidPtWrongMatch"), VertexingLxyz);
                  }
                } else if (lv.Pt() >= 2.0) {
                  registry.fill(HIST("hDCANoRefitlinearNoHFHighPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearNoHFHighPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixNoHFHighPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixNoHFHighPt"), VertexingLxyz);
                  if (isMatchCorrect == 1) {
                    registry.fill(HIST("hDCANoRefitlinearNoHFHighPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                    registry.fill(HIST("hDCARefitlinearNoHFHighPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                    registry.fill(HIST("hDCARefithelixNoHFHighPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                    registry.fill(HIST("hLxyzRefithelixNoHFHighPtCorrectMatch"), VertexingLxyz);
                  } else {
                    registry.fill(HIST("hDCANoRefitlinearNoHFHighPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                    registry.fill(HIST("hDCARefitlinearNoHFHighPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                    registry.fill(HIST("hDCARefithelixNoHFHighPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                    registry.fill(HIST("hLxyzRefithelixNoHFHighPtWrongMatch"), VertexingLxyz);
                  }
                }
              }
            }
            if (signalmufromD->CheckSignal(true, fwdparticle1) && signalmufromD->CheckSignal(true, fwdparticle2)) {
              registry.fill(HIST("hDCANoRefitlinearFromD"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
              registry.fill(HIST("hDCARefitlinearFromD"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
              registry.fill(HIST("hDCARefithelixFromD"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
              registry.fill(HIST("hLxyzRefithelixFromD"), VertexingLxyz);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromDCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromDCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromDCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromDCorrectMatch"), VertexingLxyz);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromDWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromDWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromDWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromDWrongMatch"), VertexingLxyz);
              }
              if (lv.Pt() < 1.0) {
                registry.fill(HIST("hDCANoRefitlinearFromDLowPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromDLowPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromDLowPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromDLowPt"), VertexingLxyz);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearFromDLowPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromDLowPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromDLowPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromDLowPtCorrectMatch"), VertexingLxyz);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearFromDLowPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromDLowPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromDLowPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromDLowPtWrongMatch"), VertexingLxyz);
                }
              } else if (lv.Pt() >= 1.0 && lv.Pt() < 2.0) {
                registry.fill(HIST("hDCANoRefitlinearFromDMidPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromDMidPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromDMidPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromDMidPt"), VertexingLxyz);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearFromDMidPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromDMidPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromDMidPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromDMidPtCorrectMatch"), VertexingLxyz);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearFromDMidPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromDMidPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromDMidPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromDMidPtWrongMatch"), VertexingLxyz);
                }
              } else if (lv.Pt() >= 2.0) {
                registry.fill(HIST("hDCANoRefitlinearFromDHighPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromDHighPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromDHighPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromDHighPt"), VertexingLxyz);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearFromDHighPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromDHighPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromDHighPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromDHighPtCorrectMatch"), VertexingLxyz);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearFromDHighPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromDHighPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromDHighPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromDHighPtWrongMatch"), VertexingLxyz);
                }
              }
            }
            if (signalmufromB->CheckSignal(true, fwdparticle1) && signalmufromB->CheckSignal(true, fwdparticle2)) {
              registry.fill(HIST("hDCANoRefitlinearFromB"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
              registry.fill(HIST("hDCARefitlinearFromB"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
              registry.fill(HIST("hDCARefithelixFromB"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
              registry.fill(HIST("hLxyzRefithelixFromB"), VertexingLxyz);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromBCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromBCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromBCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromBCorrectMatch"), VertexingLxyz);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromBWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromBWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromBWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromBWrongMatch"), VertexingLxyz);
              }
              if (lv.Pt() < 1.0) {
                registry.fill(HIST("hDCANoRefitlinearFromBLowPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromBLowPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromBLowPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromBLowPt"), VertexingLxyz);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearFromBLowPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromBLowPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromBLowPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromBLowPtCorrectMatch"), VertexingLxyz);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearFromBLowPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromBLowPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromBLowPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromBLowPtWrongMatch"), VertexingLxyz);
                }
              } else if (lv.Pt() >= 1.0 && lv.Pt() < 2.0) {
                registry.fill(HIST("hDCANoRefitlinearFromBMidPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromBMidPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromBMidPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromBMidPt"), VertexingLxyz);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearFromBMidPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromBMidPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromBMidPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromBMidPtCorrectMatch"), VertexingLxyz);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearFromBMidPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromBMidPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromBMidPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromBMidPtWrongMatch"), VertexingLxyz);
                }
              } else if (lv.Pt() >= 2.0) {
                registry.fill(HIST("hDCANoRefitlinearFromBHighPt"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                registry.fill(HIST("hDCARefitlinearFromBHighPt"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                registry.fill(HIST("hDCARefithelixFromBHighPt"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                registry.fill(HIST("hLxyzRefithelixFromBHighPt"), VertexingLxyz);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearFromBHighPtCorrectMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromBHighPtCorrectMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromBHighPtCorrectMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromBHighPtCorrectMatch"), VertexingLxyz);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearFromBHighPtWrongMatch"), std::sqrt((DCA01 * DCA01 + DCA02 * DCA02)/2));
                  registry.fill(HIST("hDCARefitlinearFromBHighPtWrongMatch"), std::sqrt((DCA11 * DCA11 + DCA12 * DCA12)/2));
                  registry.fill(HIST("hDCARefithelixFromBHighPtWrongMatch"), std::sqrt((DCA21 * DCA21 + DCA22 * DCA22)/2));
                  registry.fill(HIST("hLxyzRefithelixFromBHighPtWrongMatch"), VertexingLxyz);
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
    adaptAnalysisTask<myQApair>(cfgc)
  };
}
