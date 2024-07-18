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

///
/// \file   timestamp.cxx
/// \author Nicolò Jacazio
/// \since  2020-06-22
/// \brief  A task to fill the timestamp table from run number.
///         Uses headers from CCDB
///
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
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
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "Common/DataModel/colhasmuon.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace std;

struct ColhasMuon {
  Produces<aod::CHasMuon> colhasmuon;  /// Table with SOR timestamps produced by the task
  //float collisionZcut = 10.0f;

  //Filter collisionFilter = nabs(aod::collision::posZ) < collisionZcut;

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions const& collisions, aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    bool hasmuon = 0;
    for (auto& collision : collisions) {
      for (auto& fwdtrackIndice : fwdtrackIndices) {
        if (collision.globalIndex() == fwdtrackIndice.collisionId()) {
          hasmuon = 1;
        }
      }
      colhasmuon(hasmuon);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ColhasMuon>(cfgc)};
}
