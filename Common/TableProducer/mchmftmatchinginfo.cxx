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
                  mchmftpairbkg::Deltax,
                  mchmftpairbkg::Deltay,
                  mchmftpairbkg::Deltaeta,
                  mchmftpairbkg::Deltaphi,
                  mchmftpairbkg::Deltatanl,
                  mchmftpairbkg::Deltapt,
                  mchmftpairbkg::Deltar);

}


struct mftmchmatchinginfo {

  //Produce tables
  Produces<aod::MchmftPair> mchmftpairTable;
  Produces<aod::MchmftPairTrue> mchmftpairtrueTable;
  Produces<aod::MchmftPairWrong> mchmftpairwrongTable;
  Produces<aod::MchmftPairBkg> mchmftpairbkgTable;

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
  
  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  //fwdtrack filtering
  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (etaup < aod::fwdtrack::eta ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  /*
  Preslice<aod::FwdTracks> perCollision = aod::fwdtrack::collisionId;
  Preslice<aod::MFTTracks> perCollisionMFT = aod::fwdtrack::collisionId;
  */

  HistogramRegistry registry{
    "registry",
    {
      //{"counter","Count manything", {HistType::kTH1F, {{3,0.5,3.5}}}},
      {"hInvariantMass", "Invariant Mass of MCH standalone track;Invariant Mass (GeV/c^{2});Counts", {HistType::kTH1F, {{5000, 0, 50}}}},
      //{"hInvariantMassMuon", "Invariant Mass of MCH standalone track (Muon);Invariant Mass (GeV/c^{2});Counts", {HistType::kTH1F, {{50, 0, 200}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
    /*
    auto count = registry.get<TH1>(HIST("counter"));
    auto* x = count->GetXaxis();
    x->SetBinLabel(1,"Event");
    x->SetBinLabel(2,"MCH tracks");
    x->SetBinLabel(3,"MFT tracks");
    */
  }

  //void process(aod::FwdTracks const& fwdtracks, aod::MFTTracks const& mfttracks, aod::Collisions const& collisions)
  void process(soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks, aod::Collisions const& collisions)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;
    //static constexpr Double_t MatchingPlaneZ = -300;
    //static constexpr Double_t MatchingPlaneZ = -505;

    /*
    for (auto& fwdtrack : fwdtracks){
      LOGF(info,"---------------------------------");
      LOGF(info,"fwdtrack.collisionId() = %d",fwdtrack.collisionId());
      LOGF(info,"fwdtrack.globalIndex() = %d",fwdtrack.globalIndex());
      LOGF(info,"fwdtrack.trackType() = %d",fwdtrack.trackType());
      LOGF(info,"fwdtrack.matchMCHTrackId() = %d",fwdtrack.matchMCHTrackId());
      LOGF(info,"fwdtrack.matchMFTTrackId() = %d",fwdtrack.matchMFTTrackId());
    }
    */
    for (auto& [fwdtrack, mfttrack] : combinations(CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
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

          //adding
          //TrackParam trackParamMuon(fwdtrack.z(), track.getParameters());
          //< 5 parameters: X (cm), SlopeX, Y (cm), SlopeY, q/pYZ ((GeV/c)^-1)
          //~adding

          //update the talbe matchedmuonmft
          if (fwdtrack.collisionId() == mfttrack.collisionId()){
            mchmftpairTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            if (fwdtrack.sign() != mfttrack.sign()){
              mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            }
          }
          /*
          else {
            auto Col1 = fwdtrack.collision();
            auto Col2 = mfttrack.collision();
            auto groupedFwdTracks1 = fwdtracks.sliceBy(perCollision, Col1.globalIndex());
            auto groupedFwdTracks2 = fwdtracks.sliceBy(perCollision, Col2.globalIndex());
            auto groupedMFTTracks1 = mfttracks.sliceBy(perCollisionMFT, Col1.globalIndex());
            auto groupedMFTTracks2 = mfttracks.sliceBy(perCollisionMFT, Col2.globalIndex());

            if (groupedFwdTracks1.size() > 0 && groupedFwdTracks2.size() > 0){
              for (int i = 0; i< 8; i++){
                if (-20 + 5 * i < Col1.posZ() && Col1.posZ() < -20 + 5 * (i + 1) && -20 + 5 * i < Col2.posZ() && Col2.posZ() < -20 + 5 * (i + 1) ){
                  if (0 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 20 && 0 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 20){
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (20 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 40 && 20 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 40) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (40 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 60 && 40 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 60) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (60 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 80 && 60 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 80) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (80 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 90 && 80 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 90) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (90 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 100 && 90 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 100) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (100 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 110 && 100 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 110) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (110 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 120 && 110 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 120) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (120 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 130 && 120 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 130) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (130 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 140 && 130 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 140) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (140 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 150 && 140 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 150) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
                  }
                  else if (150 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 160 && 150 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 160) {
                    mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
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
  
  void processGen(soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&, aod::Collisions const& colllisions)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;
    //static constexpr Double_t MatchingPlaneZ = -300;
    //static constexpr Double_t MatchingPlaneZ = -505;

    for (auto& [fwdtrack, mfttrack] : combinations(CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {

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
              mchmftpairtrueTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else {
              //Write Wrong MCH-MFT pair Table
              mchmftpairwrongTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            }


          }
        }
      }
    }
  }
  PROCESS_SWITCH(mftmchmatchinginfo, processGen, "Show displacement of mft and mch tracks", false);

  void processDimuonMass(aod::Collisions::iterator const& collision, soa::Filtered<o2::aod::FwdTracks> const& fwdtracks){
    for (auto& [fwdtrack1,fwdtrack2] : combinations(CombinationsStrictlyUpperIndexPolicy(fwdtracks, fwdtracks))) {
      if (fwdtrack1.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrack2.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        TLorentzVector lv1, lv2, lv;
        lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
        lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
        lv = lv1 + lv2;
        registry.fill(HIST("hInvariantMass"), lv.M());
      }
    }
  }
  PROCESS_SWITCH(mftmchmatchinginfo, processDimuonMass, "Show dimuon mass spectrum", true);

  /*
  void processTrueMCHMFTPair(aod::Collisions::iterator const& collision,soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&){

    static constexpr Double_t MatchingPlaneZ = -505;

    for(auto const& fwdtrack1 : fwdtracks) {
      if (fwdtrack1.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        for(auto const& fwdtrack2 : fwdtracks) {
          if (fwdtrack2.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
            if (fwdtrack1.globalIndex() < fwdtrack1.globalIndex()){
              TLorentzVector lv1, lv2, lv;
              lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
              lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
              lv = lv1 + lv2;
              registry.fill(HIST("hInvariantMass"), lv.M());
              auto fwdparticle1 = fwdtrack1.mcParticle();
              if (fwdparticle1.pdgCode() == 13 || fwdparticle1.pdgCode() == -13){
                registry.fill(HIST("hInvariantMassMuon"), lv.M());
              }
              if (2.697 < lv.M() && lv.M() < 3.697){
                //propagate muontrack to matching position
                double muonchi2 = fwdtrack1.chi2();
                SMatrix5 muonpars(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.phi(), fwdtrack1.tgl(), fwdtrack1.signed1Pt());
                std::vector<double> muonv1;
                SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
                o2::track::TrackParCovFwd muonpars1{fwdtrack1.z(), muonpars, muoncovs, muonchi2};
                muonpars1.propagateToZlinear(MatchingPlaneZ);

                double frgtableinput[18];
                int onlymftmcid;
                int mftnumber = 0;
                for(auto const& mfttrack : mfttracks) {
                  if (mfttrack.has_collision()){
                    //auto mftparticle = mfttrack.mcParticle();
                    //propagate mfttrack to matching position
                    double mftchi2 = mfttrack.chi2();
                    SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                    std::vector<double> mftv1;
                    SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                    o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                    mftpars1.propagateToZlinear(MatchingPlaneZ);
                    registry.fill(HIST("frgMCHandAllMFT"), mftpars1.getX() - muonpars1.getX());
                    if (fwdtrack1.mcParticleId() == mfttrack.mcParticleId()){
                      registry.fill(HIST("frgMCHandAllMFTTrue"), mftpars1.getX() - muonpars1.getX());
                    }
                    if (mftpars1.getX() - muonpars1.getX() < 5 && mftpars1.getX() - muonpars1.getX() > -5 && mftpars1.getX() - muonpars1.getX() < 5 && mftpars1.getX() - muonpars1.getX() > -5){
                      mftnumber++;
                      frgtableinput[0] = mftpars1.getX();
                      frgtableinput[1] = mftpars1.getY();
                      frgtableinput[2] = mftpars1.getEta();
                      frgtableinput[3] = mftpars1.getPhi();
                      frgtableinput[4] = mftpars1.getPt();
                      frgtableinput[5] = std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2));
                      frgtableinput[6] = muonpars1.getX();
                      frgtableinput[7] = muonpars1.getY();
                      frgtableinput[8] = muonpars1.getEta();
                      frgtableinput[9] = muonpars1.getPhi();
                      frgtableinput[10] = muonpars1.getPt();
                      frgtableinput[11] = std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2));
                      frgtableinput[12] = muonpars1.getX() - mftpars1.getX();
                      frgtableinput[13] = muonpars1.getY() - mftpars1.getY();
                      frgtableinput[14] = muonpars1.getEta() - mftpars1.getEta();
                      frgtableinput[15] = muonpars1.getPhi() - mftpars1.getPhi();
                      frgtableinput[16] = muonpars1.getPt() - mftpars1.getPt();
                      frgtableinput[17] = std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2));
                      onlymftmcid = mfttrack.mcParticleId();
                    }
                  }
                }

                if (mftnumber == 1){
                  mchmftpairfrgTable(frgtableinput[0], frgtableinput[1], frgtableinput[2], frgtableinput[3], frgtableinput[4],frgtableinput[5],frgtableinput[6],frgtableinput[7],frgtableinput[8],frgtableinput[9],frgtableinput[10],frgtableinput[11],frgtableinput[12],frgtableinput[13],frgtableinput[14],frgtableinput[15],frgtableinput[16],frgtableinput[17]);
                  registry.fill(HIST("frgpurity"), 1);
                  if (fwdtrack1.mcParticleId() == onlymftmcid){
                    registry.fill(HIST("frgpurity"), 2);
                  }
                }
              }
            }
          }
        }
      }
    }

  }
  PROCESS_SWITCH(mftmchmatchinginfo, processTrueMCHMFTPair, "Extract true MCH-MFT pair", false);
  */

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftmchmatchinginfo>(cfgc)
  };
}
