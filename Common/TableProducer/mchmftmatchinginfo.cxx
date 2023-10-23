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
  
  //fwdtrack filtering
  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  void init(o2::framework::InitContext&)
  {
  }

  void process(aod::Collisions::iterator const& collision, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;
    for (auto& [fwdtrack, mfttrack] : combinations(CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
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
        mchmftpairTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
        if (fwdtrack.sign() != mfttrack.sign()){
          mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
        }
      }
    }
  }
   
  void processGen(aod::Collisions::iterator const& collision, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;

    for (auto& [fwdtrack, mfttrack] : combinations(CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {

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
  PROCESS_SWITCH(mftmchmatchinginfo, processGen, "Show displacement of mft and mch tracks", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftmchmatchinginfo>(cfgc)
  };
}
