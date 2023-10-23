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
#include <TLorentzVector.h>
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
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
DECLARE_SOA_COLUMN(MFTtanlMP, mfttanlmp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MFTchi2MP, mftchi2mp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONadjDeltaeta, muonadjdeltaeta, double);
DECLARE_SOA_COLUMN(MUONadjDeltaphi, muonadjdeltaphi, double);
DECLARE_SOA_COLUMN(MUONadjSimilarity, muonadjsimilarity, double);
DECLARE_SOA_COLUMN(MUONadjDeltasmallID, muonadjdeltasmallid, bool);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltatanl, deltatanl, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}
namespace mchmftpairtrue
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTtanlMP, mfttanlmp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MFTchi2MP, mftchi2mp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONadjDeltaeta, muonadjdeltaeta, double);
DECLARE_SOA_COLUMN(MUONadjDeltaphi, muonadjdeltaphi, double);
DECLARE_SOA_COLUMN(MUONadjSimilarity, muonadjsimilarity, double);
DECLARE_SOA_COLUMN(MUONadjDeltasmallID, muonadjdeltasmallid, bool);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltatanl, deltatanl, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}
namespace mchmftpairwrong
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTtanlMP, mfttanlmp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MFTchi2MP, mftchi2mp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONadjDeltaeta, muonadjdeltaeta, double);
DECLARE_SOA_COLUMN(MUONadjDeltaphi, muonadjdeltaphi, double);
DECLARE_SOA_COLUMN(MUONadjSimilarity, muonadjsimilarity, double);
DECLARE_SOA_COLUMN(MUONadjDeltasmallID, muonadjdeltasmallid, bool);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltatanl, deltatanl, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}
namespace mchmftpairbkg
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTtanlMP, mfttanlmp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MFTchi2MP, mftchi2mp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONadjDeltaeta, muonadjdeltaeta, double);
DECLARE_SOA_COLUMN(MUONadjDeltaphi, muonadjdeltaphi, double);
DECLARE_SOA_COLUMN(MUONadjSimilarity, muonadjsimilarity, double);
DECLARE_SOA_COLUMN(MUONadjDeltasmallID, muonadjdeltasmallid, bool);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltatanl, deltatanl, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}

namespace mchmftpairbkgem
{
DECLARE_SOA_COLUMN(MFTxMP, mftxmp, double);
DECLARE_SOA_COLUMN(MFTyMP, mftymp, double);
DECLARE_SOA_COLUMN(MFTetaMP, mftetamp, double);
DECLARE_SOA_COLUMN(MFTphiMP, mftphimp, double);
DECLARE_SOA_COLUMN(MFTtanlMP, mfttanlmp, double);
DECLARE_SOA_COLUMN(MFTptMP, mftptmp, double);
DECLARE_SOA_COLUMN(MFTrMP, mftrmp, double);
DECLARE_SOA_COLUMN(MFTchi2MP, mftchi2mp, double);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONadjDeltaeta, muonadjdeltaeta, double);
DECLARE_SOA_COLUMN(MUONadjDeltaphi, muonadjdeltaphi, double);
DECLARE_SOA_COLUMN(MUONadjSimilarity, muonadjsimilarity, double);
DECLARE_SOA_COLUMN(MUONadjDeltasmallID, muonadjdeltasmallid, bool);
DECLARE_SOA_COLUMN(Deltax, deltax, double);
DECLARE_SOA_COLUMN(Deltay, deltay, double);
DECLARE_SOA_COLUMN(Deltaeta, deltaeta, double);
DECLARE_SOA_COLUMN(Deltaphi, deltaphi, double);
DECLARE_SOA_COLUMN(Deltatanl, deltatanl, double);
DECLARE_SOA_COLUMN(Deltapt, deltapt, double);
DECLARE_SOA_COLUMN(Deltar, deltar, double);
}

DECLARE_SOA_TABLE(MchmftPair, "AOD", "MCHMFTPAIR",
                  mchmftpair::MFTxMP,
                  mchmftpair::MFTyMP,
                  mchmftpair::MFTetaMP,
                  mchmftpair::MFTphiMP,
                  mchmftpair::MFTtanlMP,
                  mchmftpair::MFTptMP,
                  mchmftpair::MFTrMP,
                  mchmftpair::MFTchi2MP,
                  mchmftpair::MUONxMP,
                  mchmftpair::MUONyMP,
                  mchmftpair::MUONetaMP,
                  mchmftpair::MUONphiMP,
                  mchmftpair::MUONtanlMP,
                  mchmftpair::MUONptMP,
                  mchmftpair::MUONrMP,
                  mchmftpair::MUONadjDeltaeta,
                  mchmftpair::MUONadjDeltaphi,
                  mchmftpair::MUONadjSimilarity,
                  mchmftpair::MUONadjDeltasmallID,
                  mchmftpair::Deltax,
                  mchmftpair::Deltay,
                  mchmftpair::Deltaeta,
                  mchmftpair::Deltaphi,
                  mchmftpair::Deltatanl,
                  mchmftpair::Deltapt,
                  mchmftpair::Deltar);

DECLARE_SOA_TABLE(MchmftPairTrue, "AOD", "MCHMFTPAIRTRUE",
                  mchmftpairtrue::MFTxMP,
                  mchmftpairtrue::MFTyMP,
                  mchmftpairtrue::MFTetaMP,
                  mchmftpairtrue::MFTphiMP,
                  mchmftpairtrue::MFTtanlMP,
                  mchmftpairtrue::MFTptMP,
                  mchmftpairtrue::MFTrMP,
                  mchmftpairtrue::MFTchi2MP,
                  mchmftpairtrue::MUONxMP,
                  mchmftpairtrue::MUONyMP,
                  mchmftpairtrue::MUONetaMP,
                  mchmftpairtrue::MUONphiMP,
                  mchmftpairtrue::MUONtanlMP,
                  mchmftpairtrue::MUONptMP,
                  mchmftpairtrue::MUONrMP,
                  mchmftpairtrue::MUONadjDeltaeta,
                  mchmftpairtrue::MUONadjDeltaphi,
                  mchmftpairtrue::MUONadjSimilarity,
                  mchmftpairtrue::MUONadjDeltasmallID,
                  mchmftpairtrue::Deltax,
                  mchmftpairtrue::Deltay,
                  mchmftpairtrue::Deltaeta,
                  mchmftpairtrue::Deltaphi,
                  mchmftpairtrue::Deltatanl,
                  mchmftpairtrue::Deltapt,
                  mchmftpairtrue::Deltar);

DECLARE_SOA_TABLE(MchmftPairWrong, "AOD", "MCHMFTPAIRWRONG",
                  mchmftpairwrong::MFTxMP,
                  mchmftpairwrong::MFTyMP,
                  mchmftpairwrong::MFTetaMP,
                  mchmftpairwrong::MFTphiMP,
                  mchmftpairwrong::MFTtanlMP,
                  mchmftpairwrong::MFTptMP,
                  mchmftpairwrong::MFTrMP,
                  mchmftpairwrong::MFTchi2MP,
                  mchmftpairwrong::MUONxMP,
                  mchmftpairwrong::MUONyMP,
                  mchmftpairwrong::MUONetaMP,
                  mchmftpairwrong::MUONphiMP,
                  mchmftpairwrong::MUONtanlMP,
                  mchmftpairwrong::MUONptMP,
                  mchmftpairwrong::MUONrMP,
                  mchmftpairwrong::MUONadjDeltaeta,
                  mchmftpairwrong::MUONadjDeltaphi,
                  mchmftpairwrong::MUONadjSimilarity,
                  mchmftpairwrong::MUONadjDeltasmallID,
                  mchmftpairwrong::Deltax,
                  mchmftpairwrong::Deltay,
                  mchmftpairwrong::Deltaeta,
                  mchmftpairwrong::Deltaphi,
                  mchmftpairwrong::Deltatanl,
                  mchmftpairwrong::Deltapt,
                  mchmftpairwrong::Deltar);

DECLARE_SOA_TABLE(MchmftPairBkg, "AOD", "MCHMFTPAIRBKG",
                  mchmftpairbkg::MFTxMP,
                  mchmftpairbkg::MFTyMP,
                  mchmftpairbkg::MFTetaMP,
                  mchmftpairbkg::MFTphiMP,
                  mchmftpairbkg::MFTtanlMP,
                  mchmftpairbkg::MFTptMP,
                  mchmftpairbkg::MFTrMP,
                  mchmftpairbkg::MFTchi2MP,
                  mchmftpairbkg::MUONxMP,
                  mchmftpairbkg::MUONyMP,
                  mchmftpairbkg::MUONetaMP,
                  mchmftpairbkg::MUONphiMP,
                  mchmftpairbkg::MUONtanlMP,
                  mchmftpairbkg::MUONptMP,
                  mchmftpairbkg::MUONrMP,
                  mchmftpairbkg::MUONadjDeltaeta,
                  mchmftpairbkg::MUONadjDeltaphi,
                  mchmftpairbkg::MUONadjSimilarity,
                  mchmftpairbkg::MUONadjDeltasmallID,
                  mchmftpairbkg::Deltax,
                  mchmftpairbkg::Deltay,
                  mchmftpairbkg::Deltaeta,
                  mchmftpairbkg::Deltaphi,
                  mchmftpairbkg::Deltatanl,
                  mchmftpairbkg::Deltapt,
                  mchmftpairbkg::Deltar);

DECLARE_SOA_TABLE(MchmftPairBkgem, "AOD", "MCHMFTPAIRBKGEM",
                  mchmftpairbkgem::MFTxMP,
                  mchmftpairbkgem::MFTyMP,
                  mchmftpairbkgem::MFTetaMP,
                  mchmftpairbkgem::MFTphiMP,
                  mchmftpairbkgem::MFTtanlMP,
                  mchmftpairbkgem::MFTptMP,
                  mchmftpairbkgem::MFTrMP,
                  mchmftpairbkgem::MFTchi2MP,
                  mchmftpairbkgem::MUONxMP,
                  mchmftpairbkgem::MUONyMP,
                  mchmftpairbkgem::MUONetaMP,
                  mchmftpairbkgem::MUONphiMP,
                  mchmftpairbkgem::MUONtanlMP,
                  mchmftpairbkgem::MUONptMP,
                  mchmftpairbkgem::MUONrMP,
                  mchmftpairbkgem::MUONadjDeltaeta,
                  mchmftpairbkgem::MUONadjDeltaphi,
                  mchmftpairbkgem::MUONadjSimilarity,
                  mchmftpairbkgem::MUONadjDeltasmallID,
                  mchmftpairbkgem::Deltax,
                  mchmftpairbkgem::Deltay,
                  mchmftpairbkgem::Deltaeta,
                  mchmftpairbkgem::Deltaphi,
                  mchmftpairbkgem::Deltatanl,
                  mchmftpairbkgem::Deltapt,
                  mchmftpairbkgem::Deltar);
}

struct mchmftmatchinginfoem {

  //Produce tables
  Produces<aod::MchmftPair> mchmftpairTable;
  Produces<aod::MchmftPairTrue> mchmftpairtrueTable;
  Produces<aod::MchmftPairWrong> mchmftpairwrongTable;
  Produces<aod::MchmftPairBkg> mchmftpairbkgTable;
  Produces<aod::MchmftPairBkgem> mchmftpairbkgemTable;

  //List of cut parameters
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
  
  //fwdtrack filtering
  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  Configurable<int> ndepth{"ndepth", 5, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType colBinning{{ConfVtxBins}, true};

  Preslice<aod::FwdTracks> perCollision = aod::fwdtrack::collisionId;
  Preslice<aod::MFTTracks> perCollisionMFT = aod::fwdtrack::collisionId;

  void init(o2::framework::InitContext&)
  {
  }

  void process(soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks, aod::Collisions const& collisions)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;

    for (auto const& fwdtrack : fwdtracks){
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        double adjDeltaeta = 10000000;
        double adjDeltaphi = 10000000;
        double adjSimilarity = 10000000;
        bool adjDeltasmallID = true;
        for (auto const& secondfwdtrack : fwdtracks) {
          if (secondfwdtrack.has_collision() && secondfwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
            if (secondfwdtrack.collisionId() == fwdtrack.collisionId()){
              if (secondfwdtrack.globalIndex() != fwdtrack.globalIndex()){
                if (adjSimilarity > std::sqrt(std::pow(fwdtrack.eta() - secondfwdtrack.eta(), 2) + std::pow(fwdtrack.phi() - secondfwdtrack.phi(), 2))){
                  adjDeltaeta = fwdtrack.eta() - secondfwdtrack.eta();
                  adjDeltaphi = fwdtrack.phi() - secondfwdtrack.phi();
                  adjSimilarity = std::sqrt(std::pow(adjDeltaeta, 2) + std::pow(adjDeltaphi, 2));
                  if (secondfwdtrack.globalIndex() < fwdtrack.globalIndex()){
                    adjDeltasmallID = false;
                  }
                }
              }
            }
          }
        }

        for (auto const& mfttrack : mfttracks){
          if (mfttrack.has_collision()){
            //update the talbe matchedmuonmft
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

              mchmftpairTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              if (fwdtrack.sign() != mfttrack.sign()){
                mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              }
            }
          }
        }
      }
    }
  }

  void processME(soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks, aod::Collisions const& collisions)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.
      auto groupedFwdTracks1 = fwdtracks.sliceBy(perCollision, collision1.globalIndex());
      auto groupedMFTTracks1 = mfttracks.sliceBy(perCollisionMFT, collision1.globalIndex());
      auto groupedFwdTracks2 = fwdtracks.sliceBy(perCollision, collision2.globalIndex());
      auto groupedMFTTracks2 = mfttracks.sliceBy(perCollisionMFT, collision2.globalIndex());

      for (auto const& fwdtrack : groupedFwdTracks1){
        if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
          double adjDeltaeta = 10000000;
          double adjDeltaphi = 10000000;
          double adjSimilarity = 10000000;
          bool adjDeltasmallID = true;
          for (auto const& secondfwdtrack : fwdtracks) {
            if (secondfwdtrack.has_collision() && secondfwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
              if (secondfwdtrack.collisionId() == fwdtrack.collisionId()){
                if (secondfwdtrack.globalIndex() != fwdtrack.globalIndex()){
                  if (adjSimilarity > std::sqrt(std::pow(fwdtrack.eta() - secondfwdtrack.eta(), 2) + std::pow(fwdtrack.phi() - secondfwdtrack.phi(), 2))){
                    adjDeltaeta = fwdtrack.eta() - secondfwdtrack.eta();
                    adjDeltaphi = fwdtrack.phi() - secondfwdtrack.phi();
                    adjSimilarity = std::sqrt(std::pow(adjDeltaeta, 2) + std::pow(adjDeltaphi, 2));
                    if (secondfwdtrack.globalIndex() < fwdtrack.globalIndex()){
                      adjDeltasmallID = false;
                    }
                  }
                }
              }
            }
          }

          for (auto const& mfttrack : groupedMFTTracks2){
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

            if (0 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 20 && 0 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 20){

              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else if (20 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 40 && 20 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 40){
              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else if (40 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 60 && 40 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 60){
              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else if (60 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 80 && 60 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 80){
              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else if (80 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 100 && 80 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 100){
              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else if (100 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 120 && 100 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 120){
              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else if (120 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 140 && 120 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 140){
              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else if (140 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 160 && 140 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 160){
              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else if (160 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 180 && 160 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 180){
              mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(mchmftmatchinginfoem, processME, "Show displacement of mft and mch tracks", false);
   
  void processGen(soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&, aod::Collisions const& colllisions)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;

    for (auto const& fwdtrack : fwdtracks){
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        double adjDeltaeta = 10000000;
        double adjDeltaphi = 10000000;
        double adjSimilarity = 10000000;
        bool adjDeltasmallID = true;
        for (auto const& secondfwdtrack : fwdtracks) {
          if (secondfwdtrack.has_collision() && secondfwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
            if (secondfwdtrack.collisionId() == fwdtrack.collisionId()){
              if (secondfwdtrack.globalIndex() != fwdtrack.globalIndex()){
                if (adjSimilarity > std::sqrt(std::pow(fwdtrack.eta() - secondfwdtrack.eta(), 2) + std::pow(fwdtrack.phi() - secondfwdtrack.phi(), 2))){
                  adjDeltaeta = fwdtrack.eta() - secondfwdtrack.eta();
                  adjDeltaphi = fwdtrack.phi() - secondfwdtrack.phi();
                  adjSimilarity = std::sqrt(std::pow(adjDeltaeta, 2) + std::pow(adjDeltaphi, 2));
                  if (secondfwdtrack.globalIndex() < fwdtrack.globalIndex()){
                    adjDeltasmallID = false;
                  }
                }
              }
            }
          }
        }

        for (auto const& mfttrack : mfttracks){
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
                //Write True MCH-MFT pair Table
              mchmftpairtrueTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else {
                //Write Wrong MCH-MFT pair Table
              mchmftpairwrongTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),adjDeltaeta,adjDeltaphi,adjSimilarity,adjDeltasmallID,muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(mchmftmatchinginfoem, processGen, "Show displacement of mft and mch tracks", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mchmftmatchinginfoem>(cfgc)
  };
}
