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
namespace mchmftpairfrg
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

DECLARE_SOA_TABLE(MchmftPairFrg, "AOD", "MCHMFTPAIRFRG",
                  mchmftpairfrg::MFTxMP,
                  mchmftpairfrg::MFTyMP,
                  mchmftpairfrg::MFTetaMP,
                  mchmftpairfrg::MFTphiMP,
                  mchmftpairfrg::MFTptMP,
                  mchmftpairfrg::MFTrMP,
                  mchmftpairfrg::MUONxMP,
                  mchmftpairfrg::MUONyMP,
                  mchmftpairfrg::MUONetaMP,
                  mchmftpairfrg::MUONphiMP,
                  mchmftpairfrg::MUONptMP,
                  mchmftpairfrg::MUONrMP,
                  mchmftpairfrg::Deltax,
                  mchmftpairfrg::Deltay,
                  mchmftpairfrg::Deltaeta,
                  mchmftpairfrg::Deltaphi,
                  mchmftpairfrg::Deltapt,
                  mchmftpairfrg::Deltar);
}


struct mftmchmatchinginfo {

  //Produce tables
  Produces<aod::MchmftPair> mchmftpairTable;
  Produces<aod::MchmftPairTrue> mchmftpairtrueTable;
  Produces<aod::MchmftPairWrong> mchmftpairwrongTable;
  Produces<aod::MchmftPairBkg> mchmftpairbkgTable;
  Produces<aod::MchmftPairFrg> mchmftpairfrgTable;

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
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) || (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) || (aod::fwdtrack::pDca < pDCAcutrAtBsorberEndup2)) && ((aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndlow2) || (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) || (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  HistogramRegistry registry{
    "registry",
    {
      {"counter","Count manything", {HistType::kTH1F, {{4,0.5,4.5}}}},
      {"fwdinfo","FwdTracks Information", {HistType::kTH1F, {{3,0.5,3.5}}}},
      {"mftinfo","MFTTracks Information", {HistType::kTH1F, {{3,0.5,3.5}}}},
      {"hInvariantMass", "Invariant Mass of MCH standalone track;Invariant Mass (GeV/c^{2});Counts", {HistType::kTH1F, {{50, 0, 200}}}},
      {"hInvariantMassMuon", "Invariant Mass of MCH standalone track (muon);Invariant Mass (GeV/c^{2});Counts", {HistType::kTH1F, {{50, 0, 200}}}},
      {"frgpurity","Purity of true pair extraction", {HistType::kTH1F, {{2,0.5,2.5}}}},
      {"frgMCHandAllMFT","frgMCHandAllMFT", {HistType::kTH1F, {{2000,-200,200}}}},
      {"frgMCHandAllMFTTrue","frgMCHandAllMFTTrue", {HistType::kTH1F, {{2000,-200,200}}}},
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
    auto frg = registry.get<TH1>(HIST("frgpurity"));
    auto* xfrg = frg->GetXaxis();
    xfrg->SetBinLabel(1,"Number of frg");
    xfrg->SetBinLabel(2,"True frg");
  }

  void process(soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks, aod::Collisions const& collisions)
  {
    static constexpr Double_t MatchingPlaneZ = -505;


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

          //update the talbe matchedmuonmft
          if (fwdtrack.collisionId() == mfttrack.collisionId()){
            mchmftpairTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            if (fwdtrack.sign() != mfttrack.sign()){
              mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            }
          }
        }
      }
    }
  }
  
  void processGen(aod::Collisions const& collisions,soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&)
  {
    static constexpr Double_t MatchingPlaneZ = -505;

    for(auto const& fwdtrack : fwdtracks) {
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        auto fwdparticle = fwdtrack.mcParticle();
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
              mchmftpairtrueTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            } else {
              //Write Wrong MCH-MFT pair Table
              mchmftpairwrongTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
            }


          }
        }
      }
    }
  }
  PROCESS_SWITCH(mftmchmatchinginfo, processGen, "Show displacement of mft and mch tracks", false);

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
                    auto mftparticle = mfttrack.mcParticle();
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

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftmchmatchinginfo>(cfgc)
  };
}
