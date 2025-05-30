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

//using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov, aod::ReducedMuonsLabels>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;

struct muon_ambi {
  Produces<aod::MuonAmbiguities> muonAmbiguities;
  std::map<int64_t, std::vector<int64_t>> fNAssocsInBunch;    // key: track global index, value: vector of global index for events associated in-bunch (events that have in-bunch pileup or splitting)
  std::map<int64_t, std::vector<int64_t>> fNAssocsOutOfBunch; // key: track global index, value: vector of global index for events associated out-of-bunch (events that have no in-bunch pileup)
                                                              //
  void process(aod::ReducedMuonsAssoc const& assocs, MyEventsVtxCov const&, MyMuonTracksWithCov const& muons) {
    fNAssocsInBunch.clear();
    fNAssocsOutOfBunch.clear();
    muonAmbiguities.reserve(muons.size());
    for (auto& assoc : assocs) {
      auto event = assoc.template reducedevent_as<MyEventsVtxCov>();
      auto track = assoc.template reducedmuon_as<MyMuonTracksWithCov>();
      if (fNAssocsInBunch.find(track.globalIndex()) == fNAssocsInBunch.end()) {
        std::vector<int64_t> evVector = {event.globalIndex()};
        fNAssocsInBunch[track.globalIndex()] = evVector;
      } else {
        auto& evVector = fNAssocsInBunch[track.globalIndex()];
        evVector.push_back(event.globalIndex());
      }
    }
    for (auto& track : muons) {
      int8_t nInBunch = 0;
      if (fNAssocsInBunch.find(track.globalIndex()) != fNAssocsInBunch.end()) {
        nInBunch = fNAssocsInBunch[track.globalIndex()].size();
      }
      int8_t nOutOfBunch = 0;
      if (fNAssocsOutOfBunch.find(track.globalIndex()) != fNAssocsOutOfBunch.end()) {
        nOutOfBunch = fNAssocsOutOfBunch[track.globalIndex()].size();
      }
      muonAmbiguities(nInBunch, nOutOfBunch);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<muon_ambi>(cfgc)
  };
}
