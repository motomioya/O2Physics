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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace evsel;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

namespace o2::aod
{
namespace mchmftpair
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
}
namespace mchmftpairtrue
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
}
namespace mchmftpairwrong
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
}
namespace mchmftpairbkg
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
}

DECLARE_SOA_TABLE(MchmftPair, "AOD", "MCHMFTPAIR",
                  mchmftpair::MFTxMP,
                  mchmftpair::MFTyMP,
                  mchmftpair::MFTetaMP,
                  mchmftpair::MFTphiMP,
                  mchmftpair::MFTptMP,
                  mchmftpair::MUONxMP,
                  mchmftpair::MUONyMP,
                  mchmftpair::MUONetaMP,
                  mchmftpair::MUONphiMP,
                  mchmftpair::MUONptMP,
                  mchmftpair::Deltax,
                  mchmftpair::Deltay,
                  mchmftpair::Deltaeta,
                  mchmftpair::Deltaphi,
                  mchmftpair::Deltapt);

DECLARE_SOA_TABLE(MchmftPairTrue, "AOD", "MCHMFTPAIRTRUE",
                  mchmftpairtrue::MFTxMP,
                  mchmftpairtrue::MFTyMP,
                  mchmftpairtrue::MFTetaMP,
                  mchmftpairtrue::MFTphiMP,
                  mchmftpairtrue::MFTptMP,
                  mchmftpairtrue::MUONxMP,
                  mchmftpairtrue::MUONyMP,
                  mchmftpairtrue::MUONetaMP,
                  mchmftpairtrue::MUONphiMP,
                  mchmftpairtrue::MUONptMP,
                  mchmftpairtrue::Deltax,
                  mchmftpairtrue::Deltay,
                  mchmftpairtrue::Deltaeta,
                  mchmftpairtrue::Deltaphi,
                  mchmftpairtrue::Deltapt);

DECLARE_SOA_TABLE(MchmftPairWrong, "AOD", "MCHMFTPAIRWRONG",
                  mchmftpairwrong::MFTxMP,
                  mchmftpairwrong::MFTyMP,
                  mchmftpairwrong::MFTetaMP,
                  mchmftpairwrong::MFTphiMP,
                  mchmftpairwrong::MFTptMP,
                  mchmftpairwrong::MUONxMP,
                  mchmftpairwrong::MUONyMP,
                  mchmftpairwrong::MUONetaMP,
                  mchmftpairwrong::MUONphiMP,
                  mchmftpairwrong::MUONptMP,
                  mchmftpairwrong::Deltax,
                  mchmftpairwrong::Deltay,
                  mchmftpairwrong::Deltaeta,
                  mchmftpairwrong::Deltaphi,
                  mchmftpairwrong::Deltapt);

DECLARE_SOA_TABLE(MchmftPairBkg, "AOD", "MCHMFTPAIRBKG",
                  mchmftpairbkg::MFTxMP,
                  mchmftpairbkg::MFTyMP,
                  mchmftpairbkg::MFTetaMP,
                  mchmftpairbkg::MFTphiMP,
                  mchmftpairbkg::MFTptMP,
                  mchmftpairbkg::MUONxMP,
                  mchmftpairbkg::MUONyMP,
                  mchmftpairbkg::MUONetaMP,
                  mchmftpairbkg::MUONphiMP,
                  mchmftpairbkg::MUONptMP,
                  mchmftpairbkg::Deltax,
                  mchmftpairbkg::Deltay,
                  mchmftpairbkg::Deltaeta,
                  mchmftpairbkg::Deltaphi,
                  mchmftpairbkg::Deltapt);
}


struct mftmchmatchinginfo {

  Produces<aod::MchmftPair> mchmftpairTable;
  Produces<aod::MchmftPairTrue> mchmftpairtrueTable;
  Produces<aod::MchmftPairWrong> mchmftpairwrongTable;
  Produces<aod::MchmftPairBkg> mchmftpairbkgTable;

  float etalow = -4;
  float etaup = -2.5;
  float pDCAcutrAtBsorberEndlow1 = 17.6;
  float pDCAcutrAtBsorberEndup1 = 26.5;
  float pDCAcutrAtBsorberEndlow2 = 26.5;
  float pDCAcutrAtBsorberEndup2 = 89.5;
  float pDCAcutdcaup1 = 594;
  float pDCAcutdcaup2 = 324;
  float chi2up = 1000000;
  float chi2MatchMCHMIDup = 1000000;
  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (etaup < aod::fwdtrack::eta ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) || (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) || (aod::fwdtrack::pDca < pDCAcutrAtBsorberEndup2)) && ((aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndlow2) || (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) || (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  Preslice<aod::FwdTracks> perCollision = aod::fwdtrack::collisionId;
  Preslice<aod::MFTTracks> perCollisionMFT = aod::fwdtrack::collisionId;

  using FwdTracksLabeled = soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>>;
  using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

  Configurable<bool> rejectSimilarTracks{"rejectSimilarTracks", false, "rejectSimilarTracks"};
  Configurable<float> slimilarThr{"similarThr", 0.1, "Threshold of similar event cut"};
  Configurable<bool> eventMixingCut{"eventMixingCut", false, "Event mixing cut using total pt of fwdtrack"};
  Configurable<double> eventMixingCutThr{"eventMixingCutThr", 0.5, "Range of total pt for binned event mixing"};

  void init(o2::framework::InitContext&)
  {
  }

  void process(FwdTracksLabeled const& fwdtracks, MFTTracksLabeled const& mfttracks, aod::Collisions const& collisions)
  {
    static constexpr Double_t MatchingPlaneZ = -505;
    for (auto& [fwdtrack, mfttrack] : combinations(CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {

        bool hasSimilarTrack = false;
        if (rejectSimilarTracks){
          for (auto const& secondfwdtrack : fwdtracks) {
            if (secondfwdtrack.has_collision() && secondfwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
              if (secondfwdtrack.collisionId() != fwdtrack.collisionId()){
                if (fwdtrack.eta() - secondfwdtrack.eta() < slimilarThr && fwdtrack.eta() - secondfwdtrack.eta() > -slimilarThr){
                  if (fwdtrack.phi() - secondfwdtrack.phi() < slimilarThr && fwdtrack.phi() - secondfwdtrack.phi() > -slimilarThr){
                    if (fwdtrack.pt() - secondfwdtrack.pt() < slimilarThr && fwdtrack.pt() - secondfwdtrack.pt() > -slimilarThr){
                      if (fwdtrack.globalIndex() > secondfwdtrack.globalIndex()) {
                        hasSimilarTrack = true;
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        if (hasSimilarTrack) continue;

        if (mfttrack.has_collision()){
          //propagate muontrack to matching position
          double muonchi2 = fwdtrack.chi2();
          SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
          std::vector<double> muonv1;
          SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
          o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
          muonpars1.propagateToZlinear(MatchingPlaneZ);

          //propagate mfttrack to matching position
          double mftchi2 = mfttrack.chi2();
          SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
          std::vector<double> mftv1;
          SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
          o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
          mftpars1.propagateToZlinear(MatchingPlaneZ);

          //update the talbe matchedmuonmft
          if (fwdtrack.collisionId() == mfttrack.collisionId()){
            mchmftpairTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
          }
          else
          {
            if (eventMixingCut){
              double totalPt1 = 0;
              double totalPt2 = 0;
              auto Col1 = fwdtrack.collision();
              auto Col2 = mfttrack.collision();
              auto groupedFwdTracks1 = fwdtracks.sliceBy(perCollision, Col1.globalIndex());
              auto groupedFwdTracks2 = fwdtracks.sliceBy(perCollision, Col2.globalIndex());
              auto groupedMFTTracks1 = mfttracks.sliceBy(perCollisionMFT, Col1.globalIndex());
              auto groupedMFTTracks2 = mfttracks.sliceBy(perCollisionMFT, Col2.globalIndex());

              if (groupedFwdTracks1.size() > 0 && groupedFwdTracks2.size() > 0){
                for (auto& track1 : groupedFwdTracks1) {
                  totalPt1 = totalPt1 + track1.pt();
                }
                for (auto& track2 : groupedFwdTracks2) {
                  totalPt2 = totalPt2 + track2.pt();
                }

                for (int i = 0; i < 20/eventMixingCutThr; i++){
                  if (eventMixingCutThr * i < totalPt1 && totalPt1 < eventMixingCutThr * (i + 1) && eventMixingCutThr * i < totalPt2 && totalPt2 < eventMixingCutThr * (i + 1) ){
                  //for (int i = 0; i < 200/20; i++){
                    //if (i * 5 < groupedMFTTracks1.size() && (i + 1) * 20 < groupedMFTTracks1.size() && i * 20 < groupedMFTTracks2.size() && (i + 1) * 20 < groupedMFTTracks2.size()){
                      mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    //}
                  //}
                  }
                }
              }
            }
            else
            {
              mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
            }
          }
        }
      }
    }
  }
  
  void processGen(FwdTracksLabeled const& fwdtracks, MFTTracksLabeled const& mfttracks)
  {
    static constexpr Double_t MatchingPlaneZ = -505;
    for (auto& [fwdtrack, mfttrack] : combinations(CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {

        bool hasSimilarTrack = false;
        if (rejectSimilarTracks){
          for (auto const& secondfwdtrack : fwdtracks) {
            if (secondfwdtrack.has_collision() && secondfwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
              if (secondfwdtrack.collisionId() != fwdtrack.collisionId()){
                if (fwdtrack.eta() - secondfwdtrack.eta() < slimilarThr && fwdtrack.eta() - secondfwdtrack.eta() > -slimilarThr){
                  if (fwdtrack.phi() - secondfwdtrack.phi() < slimilarThr && fwdtrack.phi() - secondfwdtrack.phi() > -slimilarThr){
                    if (fwdtrack.pt() - secondfwdtrack.pt() < slimilarThr && fwdtrack.pt() - secondfwdtrack.pt() > -slimilarThr){
                      if (fwdtrack.globalIndex() > secondfwdtrack.globalIndex()) {
                        hasSimilarTrack = true;
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        if (hasSimilarTrack) continue;

        if (mfttrack.has_collision()){
          if (fwdtrack.collisionId() == mfttrack.collisionId()){
            //propagate muontrack to matching position
            double muonchi2 = fwdtrack.chi2();
            SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
            std::vector<double> muonv1;
            SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
            o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
            muonpars1.propagateToZlinear(MatchingPlaneZ);

            //propagate mfttrack to matching position
            double mftchi2 = mfttrack.chi2();
            SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
            std::vector<double> mftv1;
            SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
            o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
            mftpars1.propagateToZlinear(MatchingPlaneZ);

            if (fwdtrack.mcParticleId() == mfttrack.mcParticleId())
            {
              mchmftpairtrueTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
            } else {
              mchmftpairwrongTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(mftmchmatchinginfo, processGen, "Show displacement of mft and mch tracks", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftmchmatchinginfo>(cfgc)
  };
}
