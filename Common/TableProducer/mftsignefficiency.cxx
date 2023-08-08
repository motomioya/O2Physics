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
#include "GlobalTracking/MatchGlobalFwd.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "Framework/ASoAHelpers.h"
#include <math.h>
#include <TLorentzVector.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct mftsignefficiency {

  HistogramRegistry registry{
    "registry",
    {
      {"hMftPt", "pT distribution of MFT tracks", {HistType::kTH1F, {{10000, 0, 100}}}},
      {"hMftPtCorrectSign", "pT distribution of MFT tracks with correct sign", {HistType::kTH1F, {{10000, 0, 100}}}}
    }
  };


  void init(o2::framework::InitContext&)
  {
  }

  void process(soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&)
  {

    for (auto const& mfttrack : mfttracks){
      if (mfttrack.has_collision()){
        const auto mcParticle = mfttrack.mcParticle_as<aod::McParticles>();
        registry.fill(HIST("hMftPt"), mcParticle.pt());
        if (mfttrack.signed1Pt() > 0 && mcParticle.pdgCode() > 0){
          registry.fill(HIST("hMftPtCorrectSign"), mcParticle.pt());
        }
        else if (mfttrack.signed1Pt() < 0 && mcParticle.pdgCode() < 0){
          registry.fill(HIST("hMftPtCorrectSign"), mcParticle.pt());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftsignefficiency>(cfgc)
  };
}
