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
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}
namespace mchmftpairtrue
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}
namespace mchmftpairwrong
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}
namespace mchmftpairbkg
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}

DECLARE_SOA_TABLE(MchmftPair, "AOD", "MCHMFTPAIR",
                  mchmftpair::MFTxMP,
                  mchmftpair::MFTyMP,
                  mchmftpair::MFTetaMP,
                  mchmftpair::MFTphiMP,
                  mchmftpair::MFTptMP,
                  mchmftpair::MFTrMP,
                  mchmftpair::MUONxMP,
                  mchmftpair::MUONyMP,
                  mchmftpair::MUONetaMP,
                  mchmftpair::MUONphiMP,
                  mchmftpair::MUONptMP,
                  mchmftpair::MUONrMP,
                  mchmftpair::Deltax,
                  mchmftpair::Deltay,
                  mchmftpair::Deltaeta,
                  mchmftpair::Deltaphi,
                  mchmftpair::Deltapt,
                  mchmftpair::Deltar);

DECLARE_SOA_TABLE(MchmftPairTrue, "AOD", "MCHMFTPAIRTRUE",
                  mchmftpairtrue::MFTxMP,
                  mchmftpairtrue::MFTyMP,
                  mchmftpairtrue::MFTetaMP,
                  mchmftpairtrue::MFTphiMP,
                  mchmftpairtrue::MFTptMP,
                  mchmftpairtrue::MFTrMP,
                  mchmftpairtrue::MUONxMP,
                  mchmftpairtrue::MUONyMP,
                  mchmftpairtrue::MUONetaMP,
                  mchmftpairtrue::MUONphiMP,
                  mchmftpairtrue::MUONptMP,
                  mchmftpairtrue::MUONrMP,
                  mchmftpairtrue::Deltax,
                  mchmftpairtrue::Deltay,
                  mchmftpairtrue::Deltaeta,
                  mchmftpairtrue::Deltaphi,
                  mchmftpairtrue::Deltapt,
                  mchmftpairtrue::Deltar);

DECLARE_SOA_TABLE(MchmftPairWrong, "AOD", "MCHMFTPAIRWRONG",
                  mchmftpairwrong::MFTxMP,
                  mchmftpairwrong::MFTyMP,
                  mchmftpairwrong::MFTetaMP,
                  mchmftpairwrong::MFTphiMP,
                  mchmftpairwrong::MFTptMP,
                  mchmftpairwrong::MFTrMP,
                  mchmftpairwrong::MUONxMP,
                  mchmftpairwrong::MUONyMP,
                  mchmftpairwrong::MUONetaMP,
                  mchmftpairwrong::MUONphiMP,
                  mchmftpairwrong::MUONptMP,
                  mchmftpairwrong::MUONrMP,
                  mchmftpairwrong::Deltax,
                  mchmftpairwrong::Deltay,
                  mchmftpairwrong::Deltaeta,
                  mchmftpairwrong::Deltaphi,
                  mchmftpairwrong::Deltapt,
                  mchmftpairwrong::Deltar);

DECLARE_SOA_TABLE(MchmftPairBkg, "AOD", "MCHMFTPAIRBKG",
                  mchmftpairbkg::MFTxMP,
                  mchmftpairbkg::MFTyMP,
                  mchmftpairbkg::MFTetaMP,
                  mchmftpairbkg::MFTphiMP,
                  mchmftpairbkg::MFTptMP,
                  mchmftpairbkg::MFTrMP,
                  mchmftpairbkg::MUONxMP,
                  mchmftpairbkg::MUONyMP,
                  mchmftpairbkg::MUONetaMP,
                  mchmftpairbkg::MUONphiMP,
                  mchmftpairbkg::MUONptMP,
                  mchmftpairbkg::MUONrMP,
                  mchmftpairbkg::Deltax,
                  mchmftpairbkg::Deltay,
                  mchmftpairbkg::Deltaeta,
                  mchmftpairbkg::Deltaphi,
                  mchmftpairbkg::Deltapt,
                  mchmftpairbkg::Deltar);

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

  Configurable<bool> rejectSimilarTracks{"rejectSimilarTracks", false, "rejectSimilarTracks"};
  Configurable<float> slimilarThr{"similarThr", 0.1, "Threshold of similar event cut"};

  HistogramRegistry registry{
    "registry",
    {
      {"counter","Count manything", {HistType::kTH1F, {{4,0.5,4.5}}}},
      {"fwdinfo","FwdTracks Information", {HistType::kTH1F, {{3,0.5,3.5}}}},
      {"mftinfo","MFTTracks Information", {HistType::kTH1F, {{3,0.5,3.5}}}}
    }
  };

  void init(o2::framework::InitContext&)
  {
    auto count = registry.get<TH1>(HIST("counter"));
    auto* x = count->GetXaxis();
    x->SetBinLabel(1,"Event");
    x->SetBinLabel(2,"MCH tracks");
    x->SetBinLabel(3,"none");
    x->SetBinLabel(4,"none");
    auto fwd = registry.get<TH1>(HIST("fwdinfo"));
    auto* xfwd = fwd->GetXaxis();
    xfwd->SetBinLabel(1,"All");
    xfwd->SetBinLabel(2,"Is muon");
    xfwd->SetBinLabel(3,"Is generated muon");
    auto mft = registry.get<TH1>(HIST("mftinfo"));
    auto* xmft = mft->GetXaxis();
    xmft->SetBinLabel(1,"All");
    xmft->SetBinLabel(2,"Is muon");
    xmft->SetBinLabel(3,"Is generated muon");
  }

  void process(soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks, aod::Collisions const& collisions)
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
            mchmftpairTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            if (fwdtrack.sign() != mfttrack.sign()){
              mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            }
          }
          /*
          else
          {
            auto Col1 = fwdtrack.collision();
            auto Col2 = mfttrack.collision();
            auto groupedFwdTracks1 = fwdtracks.sliceBy(perCollision, Col1.globalIndex());
            auto groupedFwdTracks2 = fwdtracks.sliceBy(perCollision, Col2.globalIndex());
            auto groupedMFTTracks1 = mfttracks.sliceBy(perCollisionMFT, Col1.globalIndex());
            auto groupedMFTTracks2 = mfttracks.sliceBy(perCollisionMFT, Col2.globalIndex());

            if (groupedFwdTracks1.size() > 0 && groupedFwdTracks2.size() > 0){
              //for (int i = 0; i < 200/20; i++){
                //if (i * 5 < groupedMFTTracks1.size() && (i + 1) * 20 < groupedMFTTracks1.size() && i * 20 < groupedMFTTracks2.size() && (i + 1) * 20 < groupedMFTTracks2.size()){
                //}
              //}
              for (int i = 0; i< 8; i++){
                if (-20 + 5 * i < Col1.posZ() && Col1.posZ() < -20 + 5 * (i + 1) && -20 + 5 * i < Col2.posZ() && Col2.posZ() < -20 + 5 * (i + 1) ){
                  if (0 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 20 && 0 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 20){
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg0Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (20 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 40 && 20 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 40) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg20Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (40 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 60 && 40 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 60) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg40Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (60 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 80 && 60 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 80) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg60Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (80 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 90 && 80 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 90) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg80Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (90 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 100 && 90 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 100) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg90Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (100 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 110 && 100 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 110) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg100Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (110 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 120 && 110 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 120) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg110Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (120 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 130 && 120 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 130) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg120Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (130 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 140 && 130 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 140) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg130Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (140 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 150 && 140 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 150) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg140Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                  else if (150 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 160 && 150 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 160) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                    mchmftpairbkg150Table(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt());
                  }
                }
              }
            }
          }
          */
        }
      }
    }
  }
  
  void processGen(soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&)
  {
    static constexpr Double_t MatchingPlaneZ = -505;

    for(auto const& fwdtrack : fwdtracks) {
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        auto fwdparticle = fwdtrack.mcParticle();
        LOGF(info, "-------------------------fwdtrack-------------------------\n");
        LOGF(info, "mcParticleId = %d\n", fwdtrack.mcParticleId());
        LOGF(info, "pdgCode = %d\n", fwdparticle.pdgCode());
        LOGF(info, "vz = %g\n", fwdparticle.vz());
        LOGF(info, "producedByGenerator = %B\n", fwdparticle.producedByGenerator());
        registry.fill(HIST("counter"),2.);
        registry.fill(HIST("fwdinfo"),1.);
        if (fwdparticle.pdgCode() == 13 || fwdparticle.pdgCode() == -13) {
          registry.fill(HIST("fwdinfo"),2.);
          if (fwdparticle.producedByGenerator()) {
            registry.fill(HIST("fwdinfo"),3.);
          }
        }
      }
    }
    for(auto const& mfttrack : mfttracks) {
      if (mfttrack.has_collision()){
        auto mftparticle = mfttrack.mcParticle();
        LOGF(info, "-------------------------mfttrack-------------------------\n");
        LOGF(info, "mcParticleId = %d\n", mfttrack.mcParticleId());
        LOGF(info, "pdgCode = %d\n", mftparticle.pdgCode());
        LOGF(info, "vz = %g\n", mftparticle.vz());
        LOGF(info, "producedByGenerator = %d\n", mftparticle.producedByGenerator());
        registry.fill(HIST("mftinfo"),1.);
        if (mftparticle.pdgCode() == 13 || mftparticle.pdgCode() == -13) {
          registry.fill(HIST("mftinfo"),2.);
          if (mftparticle.producedByGenerator()) {
            registry.fill(HIST("mftinfo"),3.);
          }
        }
      }
    }

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
              mchmftpairtrueTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else {
              mchmftpairwrongTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
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
