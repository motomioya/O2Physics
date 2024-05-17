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
#include <math.h>
#include <TLorentzVector.h>
#include <string>
#include <regex>
#include <iostream>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
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

struct checkmftdca{

  HistogramRegistry registry{
    "registry", 
    {
      {"hDCA", "DCA;DCA (cm)", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"hDCAPrimary", "DCA;DCA (cm)", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"hDCAFromD", "DCA;DCA (cm)", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"hDCAPion", "DCA;DCA (cm)", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"hDCAKaon", "DCA;DCA (cm)", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"hRxyz", "Rxyz;DCA (cm)", {HistType::kTH1F, {{5000, 0, 10}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  void process(soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McCollisions const&, aod::McParticles const& particles)
  {
    //Everything From D
    MCSignal* signalEfromD;
    MCProng prongEfromD(2, {0, 403}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signalEfromD = new MCSignal("signalEfromD", "Everything from charm", {prongEfromD}, {-1});
    //Pion not from D
    MCSignal* signalpion;
    MCProng pionprong(2, {211, 403}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false}); 
    signalpion = new MCSignal("signalpion", "Primary Muons", {pionprong}, {-1});
    //Kaon not from D
    MCSignal* signalkaon;
    MCProng kaonprong(2, {321, 403}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false}); 
    signalkaon = new MCSignal("signalkaon", "Primary Muons", {kaonprong}, {-1});

    for (auto& mfttrack : mfttracks) {
      if (mfttrack.eta() > -3.6 && mfttrack.eta() < -2.5) {
        if (mfttrack.has_mcParticle()) {
          auto mftparticle = mfttrack.mcParticle();
          if (mftparticle.has_mcCollision()) {
            auto mccollision = mftparticle.mcCollision();

            float mftdcaX = -999;
            float mftdcaY = -999;

            //calculate DCA
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            std::vector<double> mftv1;
            SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
            o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
            //mftpars1.propagateToZlinear(MatchingPlaneZ);
            mftpars1.propagateToZlinear(mccollision.posZ());
            mftdcaX = (mftpars1.getX() - mccollision.posX());
            mftdcaY = (mftpars1.getY() - mccollision.posY());
            float MFTDCA = std::sqrt(mftdcaX * mftdcaX + mftdcaY * mftdcaY);
            registry.fill(HIST("hDCA"), MFTDCA);
            float rxyz = std::sqrt((mccollision.posX() - mftparticle.vx()) * (mccollision.posX() - mftparticle.vx()) + (mccollision.posX() - mftparticle.vx()) * (mccollision.posX() - mftparticle.vx()) + (mccollision.posX() - mftparticle.vx()) * (mccollision.posX() - mftparticle.vx()));
            registry.fill(HIST("hRxyz"), rxyz);
            if (rxyz < 0.01) {
              registry.fill(HIST("hDCAPrimary"), MFTDCA);
            }
            if (signalEfromD->CheckSignal(true, mftparticle)) {
              registry.fill(HIST("hDCAFromD"), MFTDCA);
            } else if (signalpion->CheckSignal(true, mftparticle)) {
              registry.fill(HIST("hDCAPion"), MFTDCA);
            } else if (signalkaon->CheckSignal(true, mftparticle)) {
              registry.fill(HIST("hDCAKaon"), MFTDCA);
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
    adaptAnalysisTask<checkmftdca>(cfgc)
  };
}
