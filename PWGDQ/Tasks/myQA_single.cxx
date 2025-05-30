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

struct myQAsingle {

  float collisionZcut = 10.0f;
  Preslice<aod::McParticles> particlesIndicesPerCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  float mMagField = -5.0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{
    "registry", 
    {
      {"hFwdMatchCounter", "couneter;", {HistType::kTH1F, {{2, -0.5, 1.5}}}},

      {"hDCANoRefitlinearNoHF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPion", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaon", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromB", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsi", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPion", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaon", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromB", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsi", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPion", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaon", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromB", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsi", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLF", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},

      {"hDCANoRefitlinearNoHFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},

      {"hDCANoRefitlinearNoHFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},



      {"hDCANoRefitlinearNoHFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFLowPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},

      {"hDCANoRefitlinearNoHFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFLowPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},

      {"hDCANoRefitlinearNoHFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFLowPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},



      {"hDCANoRefitlinearNoHFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFMidPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},

      {"hDCANoRefitlinearNoHFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFMidPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},

      {"hDCANoRefitlinearNoHFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFMidPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},



      {"hDCANoRefitlinearNoHFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFHighPt", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},

      {"hDCANoRefitlinearNoHFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFHighPtCorrectMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},

      {"hDCANoRefitlinearNoHFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromPionHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromKaonHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromDHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromBHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromJpsiHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCANoRefitlinearFromLFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearNoHFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromPionHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromKaonHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromDHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromBHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromJpsiHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefitlinearFromLFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.5}}}},
      {"hDCARefithelixNoHFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromPionHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromKaonHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromDHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromBHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromJpsiHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}},
      {"hDCARefithelixFromLFHighPtWrongMatch", "DCA;DCA (cm)", {HistType::kTH1F, {{10000, 0, 0.005}}}}
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(MyEventsVtxCov const&, aod::ReducedMCEvents const&, aod::ReducedMCTracks const&, MyMuonTracksWithCovWithAmbiguities const&, aod::ReducedMuonsAssoc const& assocs, MyMFTTrack const&)
  {
    //muon from pion
    MCSignal* signalmufrompion;
    MCProng mufrompionprong(2, {13, 211}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufrompion = new MCSignal("signalmufrompion", "Primary Muons", {mufrompionprong}, {-1});
    //muon from kaon
    MCSignal* signalmufromkaon;
    MCProng mufromkaonprong(2, {13, 321}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromkaon = new MCSignal("signalmufromkaon", "Primary Muons", {mufromkaonprong}, {-1});
    //muon from D
    MCSignal* signalmufromD;
    MCProng mufromDprong(2, {13, 403}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromD = new MCSignal("signalmufromD", "Primary Muons", {mufromDprong}, {-1});
    //muon from B
    MCSignal* signalmufromB;
    MCProng mufromBprong(2, {13, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromB = new MCSignal("signalmufromB", "Primary Muons", {mufromBprong}, {-1});
    //muon from omega, eta, phi
    MCSignal* signalmufromphi;
    MCProng mufromphiprong(2, {13, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromphi = new MCSignal("signalmufromphi", "Primary Muons", {mufromphiprong}, {-1});
    MCSignal* signalmufromomega;
    MCProng mufromomegaprong(2, {13, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromomega = new MCSignal("signalmufromomega", "Primary Muons", {mufromomegaprong}, {-1});
    MCSignal* signalmufrometa;
    MCProng mufrometaprong(2, {13, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufrometa = new MCSignal("signalmufrometa", "Primary Muons", {mufrometaprong}, {-1});
    MCSignal* signalmufromrho;
    MCProng mufromrhoprong(2, {13, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromrho = new MCSignal("signalmufromrho", "Primary Muons", {mufromrhoprong}, {-1});
    //muon from J/psi
    MCSignal* signalmufromjpsi;
    MCProng mufromjpsiprong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); 
    signalmufromjpsi = new MCSignal("signalmufromjpsi", "Primary Muons", {mufromjpsiprong}, {-1});

    //Preslice<aod::ReducedMuonsAssoc> muonAssocsPerCollision = aod::reducedtrack_association::reducedevent;
    //Preslice<aod::ReducedMuonsAssoc> muonAssocsPerCollision = aod::reducedtrack_association::reducedeventId;

    if (!o2::base::GeometryManager::isGeometryLoaded())
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");

    for (auto& fwdtrackAssocId : assocs) {
      auto fwdtrack = fwdtrackAssocId.template reducedmuon_as<MyMuonTracksWithCovWithAmbiguities>();
      //Check Ambiguity
      if (fwdtrack.muonAmbiguityInBunch() > 1) continue;
      if (fwdtrack.muonAmbiguityOutOfBunch() > 1) continue;
      auto collision = fwdtrackAssocId.template reducedevent_as<MyEventsVtxCov>();

      //calculate DCA
      float fwd0dcaX = fwdtrack.fwdDcaX();
      float fwd0dcaY = fwdtrack.fwdDcaY();
      float DCA0 = std::sqrt(fwd0dcaX * fwd0dcaX + fwd0dcaY * fwd0dcaY);
      //set TrackParCovFwd for fwdtrack pairs
      double chi21 = fwdtrack.chi2();
      SMatrix5 tpars1(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
      std::vector<double> v11{fwdtrack.cXX(), fwdtrack.cXY(), fwdtrack.cYY(), fwdtrack.cPhiX(), fwdtrack.cPhiY(), fwdtrack.cPhiPhi(), fwdtrack.cTglX(), fwdtrack.cTglY(), fwdtrack.cTglPhi(), fwdtrack.cTglTgl(), fwdtrack.c1PtX(), fwdtrack.c1PtY(), fwdtrack.c1PtPhi(), fwdtrack.c1PtTgl(), fwdtrack.c1Pt21Pt2()};
      SMatrix55 tcovs1(v11.begin(), v11.end());
      o2::track::TrackParCovFwd pars1{fwdtrack.z(), tpars1, tcovs1, chi21};
      pars1.propagateToZlinear(collision.posZ());
      float fwd1dcaX = (pars1.getX() - collision.posX());
      float fwd1dcaY = (pars1.getY() - collision.posY());
      float DCA1 = std::sqrt(fwd1dcaX * fwd1dcaX + fwd1dcaY * fwd1dcaY);

      //fwdtrack lenear propagation
      o2::track::TrackParCovFwd pars2{fwdtrack.z(), tpars1, tcovs1, chi21};
      o2::dataformats::GlobalFwdTrack propmuon;
      auto geoMan = o2::base::GeometryManager::meanMaterialBudget(fwdtrack.x(), fwdtrack.y(), fwdtrack.z(), collision.posX(), collision.posY(), collision.posZ());
      auto x2x0 = static_cast<float>(geoMan.meanX2X0);
      pars2.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, mMagField, x2x0);
      propmuon.setParameters(pars2.getParameters());
      propmuon.setZ(pars2.getZ());
      propmuon.setCovariances(pars2.getCovariances());
      //calculate DCA
      float fwd2dcaX = (pars2.getX() - collision.posX());
      float fwd2dcaY = (pars2.getY() - collision.posY());
      float DCA2 = std::sqrt(fwd2dcaX * fwd2dcaX + fwd2dcaY * fwd2dcaY);

      //Get MFTTracks
      if (fwdtrack.has_reducedMCTrack() && fwdtrack.has_matchMFTTrack()) registry.fill(HIST("hFwdMatchCounter"), 0);
      if (fwdtrack.has_reducedMCTrack() && !fwdtrack.has_matchMFTTrack()) registry.fill(HIST("hFwdMatchCounter"), 1);

      if (fwdtrack.has_reducedMCTrack() && fwdtrack.has_matchMFTTrack()) {
        auto mfttrack = fwdtrack.matchMFTTrack_as<MyMFTTrack>();
        if (mfttrack.has_reducedMCTrack()) {
          auto mfttrack = fwdtrack.matchMFTTrack_as<MyMFTTrack>();
          //auto fwdparticle = fwdtrack.mcParticle();
          auto fwdparticle = fwdtrack.template reducedMCTrack_as<aod::ReducedMCTracks>();
          auto mftparticle = mfttrack.template reducedMCTrack_as<aod::ReducedMCTracks>();
          bool isMatchCorrect = 0;
          if (fwdparticle.globalIndex() == mftparticle.globalIndex()) isMatchCorrect = 1;
          //Fill histograms
          if (!signalmufromD->CheckSignal(true, fwdparticle)) {
            if (!signalmufromB->CheckSignal(true, fwdparticle)) {
              registry.fill(HIST("hDCANoRefitlinearNoHF"), DCA0);
              registry.fill(HIST("hDCARefitlinearNoHF"), DCA1);
              registry.fill(HIST("hDCARefithelixNoHF"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearNoHFCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearNoHFCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixNoHFCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearNoHFWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearNoHFWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixNoHFWrongMatch"), DCA2);
              }
              if (fwdtrack.pt() < 0.5) {
                registry.fill(HIST("hDCANoRefitlinearNoHFLowPt"), DCA0);
                registry.fill(HIST("hDCARefitlinearNoHFLowPt"), DCA1);
                registry.fill(HIST("hDCARefithelixNoHFLowPt"), DCA2);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearNoHFLowPtCorrectMatch"), DCA0);
                  registry.fill(HIST("hDCARefitlinearNoHFLowPtCorrectMatch"), DCA1);
                  registry.fill(HIST("hDCARefithelixNoHFLowPtCorrectMatch"), DCA2);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearNoHFLowPtWrongMatch"), DCA0);
                  registry.fill(HIST("hDCARefitlinearNoHFLowPtWrongMatch"), DCA1);
                  registry.fill(HIST("hDCARefithelixNoHFLowPtWrongMatch"), DCA2);
                }
              } else if (fwdtrack.pt() >= 0.5 && fwdtrack.pt() < 0.8) {
                registry.fill(HIST("hDCANoRefitlinearNoHFMidPt"), DCA0);
                registry.fill(HIST("hDCARefitlinearNoHFMidPt"), DCA1);
                registry.fill(HIST("hDCARefithelixNoHFMidPt"), DCA2);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearNoHFMidPtCorrectMatch"), DCA0);
                  registry.fill(HIST("hDCARefitlinearNoHFMidPtCorrectMatch"), DCA1);
                  registry.fill(HIST("hDCARefithelixNoHFMidPtCorrectMatch"), DCA2);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearNoHFMidPtWrongMatch"), DCA0);
                  registry.fill(HIST("hDCARefitlinearNoHFMidPtWrongMatch"), DCA1);
                  registry.fill(HIST("hDCARefithelixNoHFMidPtWrongMatch"), DCA2);
                }
              } else if (fwdtrack.pt() >= 0.8) {
                registry.fill(HIST("hDCANoRefitlinearNoHFHighPt"), DCA0);
                registry.fill(HIST("hDCARefitlinearNoHFHighPt"), DCA1);
                registry.fill(HIST("hDCARefithelixNoHFHighPt"), DCA2);
                if (isMatchCorrect == 1) {
                  registry.fill(HIST("hDCANoRefitlinearNoHFHighPtCorrectMatch"), DCA0);
                  registry.fill(HIST("hDCARefitlinearNoHFHighPtCorrectMatch"), DCA1);
                  registry.fill(HIST("hDCARefithelixNoHFHighPtCorrectMatch"), DCA2);
                } else {
                  registry.fill(HIST("hDCANoRefitlinearNoHFHighPtWrongMatch"), DCA0);
                  registry.fill(HIST("hDCARefitlinearNoHFHighPtWrongMatch"), DCA1);
                  registry.fill(HIST("hDCARefithelixNoHFHighPtWrongMatch"), DCA2);
                }
              }
            }
          }
          if (signalmufrompion->CheckSignal(true, fwdparticle)) {
            registry.fill(HIST("hDCANoRefitlinearFromPion"), DCA0);
            registry.fill(HIST("hDCARefitlinearFromPion"), DCA1);
            registry.fill(HIST("hDCARefithelixFromPion"), DCA2);
            if (isMatchCorrect == 1) {
              registry.fill(HIST("hDCANoRefitlinearFromPionCorrectMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromPionCorrectMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromPionCorrectMatch"), DCA2);
            } else {
              registry.fill(HIST("hDCANoRefitlinearFromPionWrongMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromPionWrongMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromPionWrongMatch"), DCA2);
            }
            if (fwdtrack.pt() < 0.5) {
              registry.fill(HIST("hDCANoRefitlinearFromPionLowPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromPionLowPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromPionLowPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromPionLowPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromPionLowPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromPionLowPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromPionLowPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromPionLowPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromPionLowPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.5 && fwdtrack.pt() < 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromPionMidPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromPionMidPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromPionMidPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromPionMidPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromPionMidPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromPionMidPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromPionMidPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromPionMidPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromPionMidPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromPionHighPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromPionHighPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromPionHighPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromPionHighPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromPionHighPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromPionHighPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromPionHighPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromPionHighPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromPionHighPtWrongMatch"), DCA2);
              }
            }
          }
          if (signalmufromkaon->CheckSignal(true, fwdparticle)) {
            registry.fill(HIST("hDCANoRefitlinearFromKaon"), DCA0);
            registry.fill(HIST("hDCARefitlinearFromKaon"), DCA1);
            registry.fill(HIST("hDCARefithelixFromKaon"), DCA2);
            if (isMatchCorrect == 1) {
              registry.fill(HIST("hDCANoRefitlinearFromKaonCorrectMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromKaonCorrectMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromKaonCorrectMatch"), DCA2);
            } else {
              registry.fill(HIST("hDCANoRefitlinearFromKaonWrongMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromKaonWrongMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromKaonWrongMatch"), DCA2);
            }
            if (fwdtrack.pt() < 0.5) {
              registry.fill(HIST("hDCANoRefitlinearFromKaonLowPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromKaonLowPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromKaonLowPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromKaonLowPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromKaonLowPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromKaonLowPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromKaonLowPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromKaonLowPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromKaonLowPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.5 && fwdtrack.pt() < 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromKaonMidPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromKaonMidPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromKaonMidPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromKaonMidPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromKaonMidPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromKaonMidPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromKaonMidPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromKaonMidPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromKaonMidPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromKaonHighPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromKaonHighPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromKaonHighPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromKaonHighPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromKaonHighPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromKaonHighPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromKaonHighPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromKaonHighPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromKaonHighPtWrongMatch"), DCA2);
              }
            }
          }
          if (signalmufromD->CheckSignal(true, fwdparticle)) {
            registry.fill(HIST("hDCANoRefitlinearFromD"), DCA0);
            registry.fill(HIST("hDCARefitlinearFromD"), DCA1);
            registry.fill(HIST("hDCARefithelixFromD"), DCA2);
            if (isMatchCorrect == 1) {
              registry.fill(HIST("hDCANoRefitlinearFromDCorrectMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromDCorrectMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromDCorrectMatch"), DCA2);
            } else {
              registry.fill(HIST("hDCANoRefitlinearFromDWrongMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromDWrongMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromDWrongMatch"), DCA2);
            }
            if (fwdtrack.pt() < 0.5) {
              registry.fill(HIST("hDCANoRefitlinearFromDLowPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromDLowPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromDLowPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromDLowPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromDLowPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromDLowPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromDLowPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromDLowPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromDLowPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.5 && fwdtrack.pt() < 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromDMidPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromDMidPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromDMidPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromDMidPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromDMidPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromDMidPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromDMidPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromDMidPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromDMidPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromDHighPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromDHighPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromDHighPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromDHighPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromDHighPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromDHighPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromDHighPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromDHighPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromDHighPtWrongMatch"), DCA2);
              }
            }
          }
          if (signalmufromB->CheckSignal(true, fwdparticle)) {
            registry.fill(HIST("hDCANoRefitlinearFromB"), DCA0);
            registry.fill(HIST("hDCARefitlinearFromB"), DCA1);
            registry.fill(HIST("hDCARefithelixFromB"), DCA2);
            if (isMatchCorrect == 1) {
              registry.fill(HIST("hDCANoRefitlinearFromBCorrectMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromBCorrectMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromBCorrectMatch"), DCA2);
            } else {
              registry.fill(HIST("hDCANoRefitlinearFromBWrongMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromBWrongMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromBWrongMatch"), DCA2);
            }
            if (fwdtrack.pt() < 0.5) {
              registry.fill(HIST("hDCANoRefitlinearFromBLowPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromBLowPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromBLowPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromBLowPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromBLowPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromBLowPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromBLowPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromBLowPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromBLowPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.5 && fwdtrack.pt() < 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromBMidPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromBMidPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromBMidPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromBMidPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromBMidPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromBMidPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromBMidPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromBMidPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromBMidPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromBHighPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromBHighPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromBHighPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromBHighPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromBHighPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromBHighPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromBHighPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromBHighPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromBHighPtWrongMatch"), DCA2);
              }
            }
          }
          if (signalmufromphi->CheckSignal(true, fwdparticle) || signalmufromomega->CheckSignal(true, fwdparticle) || signalmufrometa->CheckSignal(true, fwdparticle) || signalmufromrho->CheckSignal(true, fwdparticle)) {
            registry.fill(HIST("hDCANoRefitlinearFromLF"), DCA0);
            registry.fill(HIST("hDCARefitlinearFromLF"), DCA1);
            registry.fill(HIST("hDCARefithelixFromLF"), DCA2);
            if (isMatchCorrect == 1) {
              registry.fill(HIST("hDCANoRefitlinearFromLFCorrectMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromLFCorrectMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromLFCorrectMatch"), DCA2);
            } else {
              registry.fill(HIST("hDCANoRefitlinearFromLFWrongMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromLFWrongMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromLFWrongMatch"), DCA2);
            }
            if (fwdtrack.pt() < 0.5) {
              registry.fill(HIST("hDCANoRefitlinearFromLFLowPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromLFLowPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromLFLowPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromLFLowPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromLFLowPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromLFLowPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromLFLowPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromLFLowPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromLFLowPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.5 && fwdtrack.pt() < 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromLFMidPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromLFMidPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromLFMidPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromLFMidPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromLFMidPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromLFMidPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromLFMidPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromLFMidPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromLFMidPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromLFHighPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromLFHighPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromLFHighPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromLFHighPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromLFHighPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromLFHighPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromLFHighPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromLFHighPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromLFHighPtWrongMatch"), DCA2);
              }
            }
          }
          if (signalmufromjpsi->CheckSignal(true, fwdparticle)) {
            registry.fill(HIST("hDCANoRefitlinearFromJpsi"), DCA0);
            registry.fill(HIST("hDCARefitlinearFromJpsi"), DCA1);
            registry.fill(HIST("hDCARefithelixFromJpsi"), DCA2);
            if (isMatchCorrect == 1) {
              registry.fill(HIST("hDCANoRefitlinearFromJpsiCorrectMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromJpsiCorrectMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromJpsiCorrectMatch"), DCA2);
            } else {
              registry.fill(HIST("hDCANoRefitlinearFromJpsiWrongMatch"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromJpsiWrongMatch"), DCA1);
              registry.fill(HIST("hDCARefithelixFromJpsiWrongMatch"), DCA2);
            }
            if (fwdtrack.pt() < 0.5) {
              registry.fill(HIST("hDCANoRefitlinearFromJpsiLowPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromJpsiLowPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromJpsiLowPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromJpsiLowPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromJpsiLowPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromJpsiLowPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromJpsiLowPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromJpsiLowPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromJpsiLowPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.5 && fwdtrack.pt() < 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromJpsiMidPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromJpsiMidPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromJpsiMidPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromJpsiMidPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromJpsiMidPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromJpsiMidPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromJpsiMidPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromJpsiMidPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromJpsiMidPtWrongMatch"), DCA2);
              }
            } else if (fwdtrack.pt() >= 0.8) {
              registry.fill(HIST("hDCANoRefitlinearFromJpsiHighPt"), DCA0);
              registry.fill(HIST("hDCARefitlinearFromJpsiHighPt"), DCA1);
              registry.fill(HIST("hDCARefithelixFromJpsiHighPt"), DCA2);
              if (isMatchCorrect == 1) {
                registry.fill(HIST("hDCANoRefitlinearFromJpsiHighPtCorrectMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromJpsiHighPtCorrectMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromJpsiHighPtCorrectMatch"), DCA2);
              } else {
                registry.fill(HIST("hDCANoRefitlinearFromJpsiHighPtWrongMatch"), DCA0);
                registry.fill(HIST("hDCARefitlinearFromJpsiHighPtWrongMatch"), DCA1);
                registry.fill(HIST("hDCARefithelixFromJpsiHighPtWrongMatch"), DCA2);
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
    adaptAnalysisTask<myQAsingle>(cfgc)
  };
}
