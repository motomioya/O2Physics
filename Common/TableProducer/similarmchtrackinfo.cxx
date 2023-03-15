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
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace similarmch
{
DECLARE_SOA_COLUMN(SIMILARid, similarid, int);
DECLARE_SOA_COLUMN(FIRSTMUONx, firstmuonx, double);
DECLARE_SOA_COLUMN(FIRSTMUONy, firstmuony, double);
DECLARE_SOA_COLUMN(FIRSTMUONeta, firstmuoneta, double);
DECLARE_SOA_COLUMN(FIRSTMUONphi, firstmuonphi, double);
DECLARE_SOA_COLUMN(FIRSTMUONpt, firstmuonpt, double);
DECLARE_SOA_COLUMN(SECONDMUONx, secondmuonx, double);
DECLARE_SOA_COLUMN(SECONDMUONy, secondmuony, double);
DECLARE_SOA_COLUMN(SECONDMUONeta, secondmuoneta, double);
DECLARE_SOA_COLUMN(SECONDMUONphi, secondmuonphi, double);
DECLARE_SOA_COLUMN(SECONDMUONpt, secondmuonpt, double);
}
DECLARE_SOA_TABLE(SimilarMch, "AOD", "SIMILARMCH",
                  similarmch::SIMILARid,
                  similarmch::FIRSTMUONx,
                  similarmch::FIRSTMUONy,
                  similarmch::FIRSTMUONeta,
                  similarmch::FIRSTMUONphi,
                  similarmch::FIRSTMUONpt,
                  similarmch::SECONDMUONx,
                  similarmch::SECONDMUONy,
                  similarmch::SECONDMUONeta,
                  similarmch::SECONDMUONphi,
                  similarmch::SECONDMUONpt);
}


struct simliarmchtrackinfo {

  Produces<aod::SimilarMch> similarmchTable;

  Configurable<float> slimilarThr{"similarThr", 0.1, "Threshold of similar event"};

  HistogramRegistry registry{
    "registry",
    {
      {"SimilarMCHTracks","Number of similar MCH tracks", {HistType::kTH1F, {{3,0.5,3.5}}}}
    }
  };
  void init(o2::framework::InitContext&)
  {
    auto simmch = registry.get<TH1>(HIST("SimilarMCHTracks"));
    auto* x = simmch->GetXaxis();
    x->SetBinLabel(1,"All");
    x->SetBinLabel(2,"Ambiguous");
    x->SetBinLabel(3,"In pileup event");

  }


  using CollisionSels = soa::Join<aod::Collisions, aod::EvSels>;
  void process(aod::FwdTracks const& fwdtracks, CollisionSels & collisions, aod::AmbiguousFwdTracks const & atracks)
  {
    int similarid = 0;

    for (auto const& fwdtrack : fwdtracks) {
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {

        for (auto const& secondfwdtrack : fwdtracks) {
          if (secondfwdtrack.has_collision() && secondfwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {

            if (secondfwdtrack.collisionId() != fwdtrack.collisionId()){
              if (fwdtrack.eta() - secondfwdtrack.eta() < slimilarThr && fwdtrack.eta() - secondfwdtrack.eta() > -slimilarThr){
                if (fwdtrack.phi() - secondfwdtrack.phi() < slimilarThr && fwdtrack.phi() - secondfwdtrack.phi() > -slimilarThr){
                  if (fwdtrack.pt() - secondfwdtrack.pt() < slimilarThr && fwdtrack.pt() - secondfwdtrack.pt() > -slimilarThr){

                    registry.fill(HIST("SimilarMCHTracks"),1.);
                    similarmchTable(similarid, fwdtrack.x(),fwdtrack.y(),fwdtrack.eta(),fwdtrack.phi(),fwdtrack.pt(),secondfwdtrack.x(),secondfwdtrack.y(),secondfwdtrack.eta(),secondfwdtrack.phi(),secondfwdtrack.pt());
                    similarid++;

                    for (auto& atrack : atracks)
                    {
                      if (fwdtrack.globalIndex() == atrack.fwdtrackId())
                      {
                        registry.fill(HIST("SimilarMCHTracks"),2.);
                      }
                    }

                    auto col = fwdtrack.collision_as<CollisionSels>();
                    if (col.selection()[evsel::kNoPileupFromSPD] == 0)
                    {
                      registry.fill(HIST("SimilarMCHTracks"),3.);
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
    adaptAnalysisTask<simliarmchtrackinfo>(cfgc)
  };
}
