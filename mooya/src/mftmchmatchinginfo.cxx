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
// Made by mooya

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"

namespace o2::aod
{
namespace matchedmuonmft
{
DECLARE_SOA_COLUMN(MFTx, mftx, double);
DECLARE_SOA_COLUMN(MFTy, mfty, double);
DECLARE_SOA_COLUMN(MUONx, muonx, double);
DECLARE_SOA_COLUMN(MUONy, muony, double);
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
}

DECLARE_SOA_TABLE(FwdMatchingInfo, "AOD", "FWDMATCHINGINFO", matchedmuonmft::MFTx,matchedmuonmft::MFTy, matchedmuonmft::MUONx, matchedmuonmft::MUONy, matchedmuonmft::MFTxMP, matchedmuonmft::MFTyMP, matchedmuonmft::MUONxMP, matchedmuonmft::MUONyMP);

}

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;

using FwdTracksLabeled = soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>;
using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;


struct mftmchmatchinginfo {

  Produces<aod::FwdMatchingInfo> matchedmuonmftTable;

  Configurable<bool> isMC{"isMC", false, "MC or not"};

  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true
  };

  void init(o2::framework::InitContext&)
  {
    AxisSpec trackXPos = {1000, -50, 50, "X in cm"};
    AxisSpec trackYPos = {1000, -50, 50, "Y in cm"};
    AxisSpec trackZPos = {1200, -100, 20, "Z in cm"};

    //HistogramConfigSpec HistVariable({HistType::kTHnSparseF, {ptRecoAxis, dcaxAxis, dcaxAxis, dcaAxis, zvtxAxis}});
    registry.add("MFTXYPos", "MFTTrack XYPosition", {HistType::kTH2F, {trackXPos, trackYPos}});
    registry.add("MFTZPos", "MFT ZPosition", {HistType::kTH1F, {trackZPos}});
    registry.add("MUONXYPos", "MUONTrack XYPosition", {HistType::kTH2F, {trackXPos, trackYPos}});
    registry.add("MUONZPos", "MUONTrack ZPosition", {HistType::kTH1F, {trackZPos}});
    registry.add("MFTXYPosProp", "MFTTrack XYPosition at Last MFT Disk", {HistType::kTH2F, {trackXPos, trackYPos}});
    registry.add("MFTZPosProp", "MFT ZPosition at Last MFT Disk", {HistType::kTH1F, {trackZPos}});
    registry.add("MUONXYPosProp", "MUONTrack XYPosition at Last MFT Disk", {HistType::kTH2F, {trackXPos, trackYPos}});
    registry.add("MUONZPosProp", "MUONTrack ZPosition at Last MFT Disk", {HistType::kTH1F, {trackZPos}});
    registry.add("MCHMFTDisp", "Displacement MUONTrack and MFTTrack in XY plane", {HistType::kTH2F, {trackXPos, trackYPos}});

    if (isMC){
      registry.add("MCHMFTDispX", "Displacement of MUONTrack and MFTTrack in X", {HistType::kTH1F, {trackXPos}});
      registry.add("MCHMFTDispY", "Displacement of MUONTrack and MFTTrack in Y", {HistType::kTH1F, {trackYPos}});
      registry.add("MCHMFTTrueDispX", "Displacement of true MUONTrack and MFTTrack in X", {HistType::kTH1F, {trackXPos}});
      registry.add("MCHMFTTrueDispY", "Displacement of true MUONTrack and MFTTrack in Y", {HistType::kTH1F, {trackYPos}});
      registry.add("MCHMFTBkgDispX", "Displacement of bkg MUONTrack and MFTTrack in X", {HistType::kTH1F, {trackXPos}});
      registry.add("MCHMFTBkgDispY", "Displacement of bkg MUONTrack and MFTTrack in Y", {HistType::kTH1F, {trackYPos}});
    }

  }

  void process(aod::Collisions::iterator const& collision, aod::FwdTracks const& tracks, aod::MFTTracks const& mfttracks)
  {
    static constexpr Double_t sLastMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[9];
    static const float mBz = -5.f;

    std::vector<std::array<int, 2>> mchmftpair;

    //create a pair of MFTTrack and MUONTrack in same globalmuontrack
    for (auto const& track : tracks) {
      if (track.has_collision() && track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        mchmftpair.push_back({track.matchMCHTrackId(), track.matchMFTTrackId()});
      }
    }


    for (auto const& track : tracks) {
      if (track.has_collision() && track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        for (auto const& pair : mchmftpair) {
	  if (track.globalIndex() == pair[0]) {
	    for (auto const& mfttrack : mfttracks) {
	      if (mfttrack.globalIndex() == pair[1]) {

		registry.fill(HIST("MUONXYPos"), track.x(), track.y());
		registry.fill(HIST("MUONZPos"), track.z());
		registry.fill(HIST("MFTXYPos"), mfttrack.x(), mfttrack.y());
		registry.fill(HIST("MFTZPos"), mfttrack.z());

                //propagate muontrack to matching position
                double muonchi2 = track.chi2();
                SMatrix5 muonpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
                std::vector<double> muonv1;
                SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
                o2::track::TrackParCovFwd muonpars1{track.z(), muonpars, muoncovs, muonchi2};
                muonpars1.propagateToZ(sLastMFTPlaneZ,mBz);
		registry.fill(HIST("MUONXYPosProp"), muonpars1.getX(), muonpars1.getY());
		registry.fill(HIST("MUONZPosProp"), muonpars1.getZ());

                //propagate mfttrack to matching position
                double mftchi2 = mfttrack.chi2();
                SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                std::vector<double> mftv1;
                SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                mftpars1.propagateToZ(sLastMFTPlaneZ,mBz);
		registry.fill(HIST("MFTXYPosProp"), mftpars1.getX(), mftpars1.getY());
		registry.fill(HIST("MFTZPosProp"), mftpars1.getZ());

		//Plot displacement
		double disX = muonpars1.getX() - mftpars1.getX();
		double disY = muonpars1.getY() - mftpars1.getY();
		registry.fill(HIST("MCHMFTDisp"), disX, disY);

		//update the talbe matchedmuonmft
		matchedmuonmftTable(track.x(), track.y(), mfttrack.x(), mfttrack.y(), muonpars1.getX(), muonpars1.getY(), mftpars1.getX(), mftpars1.getY());
	      }
	    }
	  }
	}
      }
    }
  }

  void processGen(aod::Collisions::iterator const& collision, FwdTracksLabeled const& fwdtracks, MFTTracksLabeled const& mfttracks)
  {
    static constexpr Double_t sLastMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[9];
    static const float mBz = -5.f;
    for (auto const& fwdtrack : fwdtracks) {
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        for (auto const& mfttrack : mfttracks) {
	  if (mfttrack.has_collision()){
            //propagate muontrack to matching position
            double muonchi2 = fwdtrack.chi2();
            SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
            std::vector<double> muonv1;
            SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
            o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
            muonpars1.propagateToZ(sLastMFTPlaneZ,mBz);

           //propagate mfttrack to matching position
           double mftchi2 = mfttrack.chi2();
           SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
           std::vector<double> mftv1;
           SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
           o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
           mftpars1.propagateToZ(sLastMFTPlaneZ,mBz);

	   //Fill displacement of all MCH and MFT tracks combinations
	   registry.fill(HIST("MCHMFTDispX"), muonpars1.getX() - mftpars1.getX());
	   registry.fill(HIST("MCHMFTDispY"), muonpars1.getY() - mftpars1.getY());
	   
	   //Fill displacement of true MCH and MFT tracks combinations
	   if (fwdtrack.mcParticleId() == mfttrack.mcParticleId())
	   {
	     registry.fill(HIST("MCHMFTTrueDispX"), muonpars1.getX() - mftpars1.getX());
	     registry.fill(HIST("MCHMFTTrueDispY"), muonpars1.getY() - mftpars1.getY());
	   }
	   //Fill displacement of true MCH and MFT tracks combinations
	  }
	}
      }
    }
  }
  PROCESS_SWITCH(mftmchmatchinginfo, processGen, "Show displacement of mft and mch tracks", isMC);

  void processGenMatchingBkg(FwdTracksLabeled const& fwdtracks, MFTTracksLabeled const& mfttracks, aod::McParticles const& particles)
  {
    static constexpr Double_t sLastMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[9];
    static const float mBz = -5.f;
    for (auto const& fwdtrack : fwdtracks) {
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        for (auto const& mfttrack : mfttracks) {
	  if (mfttrack.has_collision()){
	    auto mcfwdparticle = fwdtrack.mcParticle();
	    auto mcmftparticle = mfttrack.mcParticle();
	    {
               //propagate muontrack to matching position
               double muonchi2 = fwdtrack.chi2();
               SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
               std::vector<double> muonv1;
               SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
               o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
               muonpars1.propagateToZ(sLastMFTPlaneZ,mBz);             
               //propagate mfttrack to matching position
               double mftchi2 = mfttrack.chi2();
               SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
               std::vector<double> mftv1;
               SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
               o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
               mftpars1.propagateToZ(sLastMFTPlaneZ,mBz);
  
               //Fill displacement of all MCH and MFT tracks combinations
               registry.fill(HIST("MCHMFTBkgDispX"), muonpars1.getX() - mftpars1.getX());
               registry.fill(HIST("MCHMFTBkgDispY"), muonpars1.getY() - mftpars1.getY());
            }

	  }
	}
      }
    }
  }
  PROCESS_SWITCH(mftmchmatchinginfo, processGenMatchingBkg, "Show background of mch mft matching", isMC);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftmchmatchinginfo>(cfgc)
  };
}
