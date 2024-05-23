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
#include <CCDB/BasicCCDBManager.h>

#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "TDatabasePDG.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "Framework/ASoAHelpers.h"
#include <math.h>
#include <TLorentzVector.h>
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "DataFormatsParameters/GRPMagField.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTrkCompColls>;

struct checkeventmixing {

  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  Configurable<int> ndepth{"ndepth", 100, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType colBinning{{ConfVtxBins}, true};
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = -5.0;
  //o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  HistogramRegistry registry{
    "registry", 
    {
      {"hUSMassME", "mass (mumu);mass (GeV/c2)", {HistType::kTH1F, {{250, 0.0, 5}}}},
      {"hLSMassME", "mass (mumu);mass (GeV/c2)", {HistType::kTH1F, {{250, 0.0, 5}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }


  void process(aod::Collisions const& collisions, MyMuons const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.
      auto fwdtrackIndices1 = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision1.globalIndex());
      auto fwdtrackIndices2 = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision2.globalIndex());

      for (auto& [fwdtrackId1, fwdtrackId2] : combinations(o2::soa::CombinationsUpperIndexPolicy(fwdtrackIndices1, fwdtrackIndices2))) {
        for (auto& fwdtrack1 : fwdtracks) {
          if (fwdtrack1.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId1.fwdtrackId() == fwdtrack1.globalIndex() && fwdtrack1.chi2MatchMCHMFT() < 50 && fwdtrack1.compatibleCollIds().size() == 1) {
            for (auto& fwdtrack2 : fwdtracks) {
              if (fwdtrack2.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId2.fwdtrackId() == fwdtrack2.globalIndex() && fwdtrack2.chi2MatchMCHMFT() < 50 && fwdtrack2.compatibleCollIds().size() == 1) {
                if (fwdtrack1.eta() > -4 && fwdtrack1.eta() < -2.5 && (((17.6 < fwdtrack1.rAtAbsorberEnd()) && (fwdtrack1.rAtAbsorberEnd() < 26.5) && (fwdtrack1.pDca() < 594)) || ((26.5 < fwdtrack1.rAtAbsorberEnd()) && (fwdtrack1.rAtAbsorberEnd() < 89.5) && (fwdtrack1.pDca() < 324)))){
                  if (fwdtrack2.eta() > -4 && fwdtrack2.eta() < -2.5 && (((17.6 < fwdtrack2.rAtAbsorberEnd()) && (fwdtrack2.rAtAbsorberEnd() < 26.5) && (fwdtrack2.pDca() < 594)) || ((26.5 < fwdtrack2.rAtAbsorberEnd()) && (fwdtrack2.rAtAbsorberEnd() < 89.5) && (fwdtrack2.pDca() < 324)))){
                    //start from here
                    TLorentzVector lv1, lv2, lv;
                    lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
                    lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
                    lv = lv1 + lv2;
                    if (fwdtrack1.signed1Pt()*fwdtrack2.signed1Pt() < 0) {
                      registry.fill(HIST("hUSMassME"), lv.M());
                      LOGF(info, "unlike-sign");
                      LOGF(info, "mass = %f", lv.M());
                      LOGF(info, "fwdtrack1.gloalIndex = %d", fwdtrack1.globalIndex());
                      LOGF(info, "fwdtrack2.gloalIndex = %d", fwdtrack2.globalIndex());
                      LOGF(info, "fwdtrack1.matchMFTTrackId = %d", fwdtrack1.matchMFTTrackId());
                      LOGF(info, "fwdtrack2.matchMFTTrackId = %d", fwdtrack2.matchMFTTrackId());
                      LOGF(info, "fwdtrack1.matchMCHTrackId = %d", fwdtrack1.matchMCHTrackId());
                      LOGF(info, "fwdtrack2.matchMCHTrackId = %d", fwdtrack2.matchMCHTrackId());
                    } else {
                      registry.fill(HIST("hLSMassME"), lv.M());
                      LOGF(info, "like-sign");
                      LOGF(info, "mass = %f", lv.M());
                      LOGF(info, "fwdtrack1.gloalIndex = %d", fwdtrack1.globalIndex());
                      LOGF(info, "fwdtrack2.gloalIndex = %d", fwdtrack2.globalIndex());
                      LOGF(info, "fwdtrack1.matchMFTTrackId = %d", fwdtrack1.matchMFTTrackId());
                      LOGF(info, "fwdtrack2.matchMFTTrackId = %d", fwdtrack2.matchMFTTrackId());
                      LOGF(info, "fwdtrack1.matchMCHTrackId = %d", fwdtrack1.matchMCHTrackId());
                      LOGF(info, "fwdtrack2.matchMCHTrackId = %d", fwdtrack2.matchMCHTrackId());
                      if (fwdtrack1.matchMFTTrackId() == fwdtrack2.matchMFTTrackId()){
                        LOGF(info, "same MFT track!");
                      }
                      if (lv.M() <= 0.22){
                        LOGF(info, "low mass!");
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
    adaptAnalysisTask<checkeventmixing>(cfgc)
  };
}
