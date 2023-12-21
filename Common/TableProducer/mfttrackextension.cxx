// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//
// Task performing forward track DCA computation
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "CommonUtils/NameConf.h"

#include "Math/SMatrix.h"
#include "ReconstructionDataFormats/TrackFwd.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct MFTTrackExtension {

  HistogramRegistry registry{
    "registry",
    {
      {"hMFTDCAx", "hMFTDCAx", {HistType::kTH1F, {{2000, -10, 10}}}},
      {"hMFTDCAy", "hMFTDCAy", {HistType::kTH1F, {{2000, -10, 10}}}},
    }
  };

  void process(aod::MFTTracks const& mfttracks, aod::Collisions const&)
  {
    for (auto& track : mfttracks) {
      float dcaX = -999;
      float dcaY = -999;
      if (track.has_collision()) {

        auto const& collision = track.collision();
        double chi2 = track.chi2();
        SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
        std::vector<double> v1;
        SMatrix55 tcovs(v1.begin(), v1.end());
        o2::track::TrackParCovFwd pars1{track.z(), tpars, tcovs, chi2};
        pars1.propagateToZlinear(collision.posZ());

        dcaX = (pars1.getX() - collision.posX());
        dcaY = (pars1.getY() - collision.posY());
        registry.fill(HIST("hMFTDCAx"), dcaX);
        registry.fill(HIST("hMFTDCAy"), dcaY);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<MFTTrackExtension>(cfgc)};
  return workflow;
}
