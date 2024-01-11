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
DECLARE_SOA_COLUMN(MFTmother0ID, mftmother0id, int);
DECLARE_SOA_COLUMN(MFTmother0pdg, mftmother0pdg, int);
DECLARE_SOA_COLUMN(MFTmother1ID, mftmother1id, int);
DECLARE_SOA_COLUMN(MFTmother1pdg, mftmother1pdg, int);
DECLARE_SOA_COLUMN(MFTfromBE, mftfrombe, int);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONmother0ID, muonmother0id, int);
DECLARE_SOA_COLUMN(MUONmother1ID, muonmother1id, int);
DECLARE_SOA_COLUMN(MUONmother0pdg, muonmother0pdg, int);
DECLARE_SOA_COLUMN(MUONmother1pdg, muonmother1pdg, int);
DECLARE_SOA_COLUMN(MUONfromBE, muonfrombe, int);
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
DECLARE_SOA_COLUMN(MFTmother0ID, mftmother0id, int);
DECLARE_SOA_COLUMN(MFTmother0pdg, mftmother0pdg, int);
DECLARE_SOA_COLUMN(MFTmother1ID, mftmother1id, int);
DECLARE_SOA_COLUMN(MFTmother1pdg, mftmother1pdg, int);
DECLARE_SOA_COLUMN(MFTfromBE, mftfrombe, int);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONmother0ID, muonmother0id, int);
DECLARE_SOA_COLUMN(MUONmother1ID, muonmother1id, int);
DECLARE_SOA_COLUMN(MUONmother0pdg, muonmother0pdg, int);
DECLARE_SOA_COLUMN(MUONmother1pdg, muonmother1pdg, int);
DECLARE_SOA_COLUMN(MUONfromBE, muonfrombe, int);
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
DECLARE_SOA_COLUMN(MFTmother0ID, mftmother0id, int);
DECLARE_SOA_COLUMN(MFTmother0pdg, mftmother0pdg, int);
DECLARE_SOA_COLUMN(MFTmother1ID, mftmother1id, int);
DECLARE_SOA_COLUMN(MFTmother1pdg, mftmother1pdg, int);
DECLARE_SOA_COLUMN(MFTfromBE, mftfrombe, int);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONmother0ID, muonmother0id, int);
DECLARE_SOA_COLUMN(MUONmother1ID, muonmother1id, int);
DECLARE_SOA_COLUMN(MUONmother0pdg, muonmother0pdg, int);
DECLARE_SOA_COLUMN(MUONmother1pdg, muonmother1pdg, int);
DECLARE_SOA_COLUMN(MUONfromBE, muonfrombe, int);
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
DECLARE_SOA_COLUMN(MFTmother0ID, mftmother0id, int);
DECLARE_SOA_COLUMN(MFTmother0pdg, mftmother0pdg, int);
DECLARE_SOA_COLUMN(MFTmother1ID, mftmother1id, int);
DECLARE_SOA_COLUMN(MFTmother1pdg, mftmother1pdg, int);
DECLARE_SOA_COLUMN(MFTfromBE, mftfrombe, int);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONmother0ID, muonmother0id, int);
DECLARE_SOA_COLUMN(MUONmother1ID, muonmother1id, int);
DECLARE_SOA_COLUMN(MUONmother0pdg, muonmother0pdg, int);
DECLARE_SOA_COLUMN(MUONmother1pdg, muonmother1pdg, int);
DECLARE_SOA_COLUMN(MUONfromBE, muonfrombe, int);
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
DECLARE_SOA_COLUMN(MFTmother0ID, mftmother0id, int);
DECLARE_SOA_COLUMN(MFTmother0pdg, mftmother0pdg, int);
DECLARE_SOA_COLUMN(MFTmother1ID, mftmother1id, int);
DECLARE_SOA_COLUMN(MFTmother1pdg, mftmother1pdg, int);
DECLARE_SOA_COLUMN(MFTfromBE, mftfrombe, bool);
DECLARE_SOA_COLUMN(MUONxMP, muonxmp, double);
DECLARE_SOA_COLUMN(MUONyMP, muonymp, double);
DECLARE_SOA_COLUMN(MUONetaMP, muonetamp, double);
DECLARE_SOA_COLUMN(MUONphiMP, muonphimp, double);
DECLARE_SOA_COLUMN(MUONtanlMP, muontanlmp, double);
DECLARE_SOA_COLUMN(MUONptMP, muonptmp, double);
DECLARE_SOA_COLUMN(MUONrMP, muonrmp, double);
DECLARE_SOA_COLUMN(MUONmother0ID, muonmother0id, int);
DECLARE_SOA_COLUMN(MUONmother1ID, muonmother1id, int);
DECLARE_SOA_COLUMN(MUONmother0pdg, muonmother0pdg, int);
DECLARE_SOA_COLUMN(MUONmother1pdg, muonmother1pdg, int);
DECLARE_SOA_COLUMN(MUONfromBE, muonfrombe, int);
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
                  mchmftpair::MFTmother0ID,
                  mchmftpair::MFTmother0pdg,
                  mchmftpair::MFTmother1ID,
                  mchmftpair::MFTmother1pdg,
                  mchmftpair::MFTfromBE,
                  mchmftpair::MUONxMP,
                  mchmftpair::MUONyMP,
                  mchmftpair::MUONetaMP,
                  mchmftpair::MUONphiMP,
                  mchmftpair::MUONtanlMP,
                  mchmftpair::MUONptMP,
                  mchmftpair::MUONrMP,
                  mchmftpair::MUONmother0ID,
                  mchmftpair::MUONmother0pdg,
                  mchmftpair::MUONmother1ID,
                  mchmftpair::MUONmother1pdg,
                  mchmftpair::MUONfromBE,
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
                  mchmftpairtrue::MFTmother0ID,
                  mchmftpairtrue::MFTmother0pdg,
                  mchmftpairtrue::MFTmother1ID,
                  mchmftpairtrue::MFTmother1pdg,
                  mchmftpairtrue::MFTfromBE,
                  mchmftpairtrue::MUONxMP,
                  mchmftpairtrue::MUONyMP,
                  mchmftpairtrue::MUONetaMP,
                  mchmftpairtrue::MUONphiMP,
                  mchmftpairtrue::MUONtanlMP,
                  mchmftpairtrue::MUONptMP,
                  mchmftpairtrue::MUONrMP,
                  mchmftpairtrue::MUONmother0ID,
                  mchmftpairtrue::MUONmother0pdg,
                  mchmftpairtrue::MUONmother1ID,
                  mchmftpairtrue::MUONmother1pdg,
                  mchmftpairtrue::MUONfromBE,
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
                  mchmftpairwrong::MFTmother0ID,
                  mchmftpairwrong::MFTmother0pdg,
                  mchmftpairwrong::MFTmother1ID,
                  mchmftpairwrong::MFTmother1pdg,
                  mchmftpairwrong::MFTfromBE,
                  mchmftpairwrong::MUONxMP,
                  mchmftpairwrong::MUONyMP,
                  mchmftpairwrong::MUONetaMP,
                  mchmftpairwrong::MUONphiMP,
                  mchmftpairwrong::MUONtanlMP,
                  mchmftpairwrong::MUONptMP,
                  mchmftpairwrong::MUONrMP,
                  mchmftpairwrong::MUONmother0ID,
                  mchmftpairwrong::MUONmother0pdg,
                  mchmftpairwrong::MUONmother1ID,
                  mchmftpairwrong::MUONmother1pdg,
                  mchmftpairwrong::MUONfromBE,
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
                  mchmftpairbkg::MFTmother0ID,
                  mchmftpairbkg::MFTmother0pdg,
                  mchmftpairbkg::MFTmother1ID,
                  mchmftpairbkg::MFTmother1pdg,
                  mchmftpairbkg::MFTfromBE,
                  mchmftpairbkg::MUONxMP,
                  mchmftpairbkg::MUONyMP,
                  mchmftpairbkg::MUONetaMP,
                  mchmftpairbkg::MUONphiMP,
                  mchmftpairbkg::MUONtanlMP,
                  mchmftpairbkg::MUONptMP,
                  mchmftpairbkg::MUONrMP,
                  mchmftpairbkg::MUONmother0ID,
                  mchmftpairbkg::MUONmother0pdg,
                  mchmftpairbkg::MUONmother1ID,
                  mchmftpairbkg::MUONmother1pdg,
                  mchmftpairbkg::MUONfromBE,
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
                  mchmftpairbkgem::MFTmother0ID,
                  mchmftpairbkgem::MFTmother0pdg,
                  mchmftpairbkgem::MFTmother1ID,
                  mchmftpairbkgem::MFTmother1pdg,
                  mchmftpairbkgem::MFTfromBE,
                  mchmftpairbkgem::MUONxMP,
                  mchmftpairbkgem::MUONyMP,
                  mchmftpairbkgem::MUONetaMP,
                  mchmftpairbkgem::MUONphiMP,
                  mchmftpairbkgem::MUONtanlMP,
                  mchmftpairbkgem::MUONptMP,
                  mchmftpairbkgem::MUONrMP,
                  mchmftpairbkgem::MUONmother0ID,
                  mchmftpairbkgem::MUONmother0pdg,
                  mchmftpairbkgem::MUONmother1ID,
                  mchmftpairbkgem::MUONmother1pdg,
                  mchmftpairbkgem::MUONfromBE,
                  mchmftpairbkgem::Deltax,
                  mchmftpairbkgem::Deltay,
                  mchmftpairbkgem::Deltaeta,
                  mchmftpairbkgem::Deltaphi,
                  mchmftpairbkgem::Deltatanl,
                  mchmftpairbkgem::Deltapt,
                  mchmftpairbkgem::Deltar);
}

struct mchmftmatchinginfoemmc {

  //Produce tables
  Produces<aod::MchmftPair> mchmftpairTable;
  Produces<aod::MchmftPairTrue> mchmftpairtrueTable;
  Produces<aod::MchmftPairWrong> mchmftpairwrongTable;
  Produces<aod::MchmftPairBkg> mchmftpairbkgTable;
  Produces<aod::MchmftPairBkgem> mchmftpairbkgemTable;

  //List of cut parameters
  float etalow = -3.6;
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
  Configurable<int> cfgColWindow{"collision-window", 1, "Search window (collision ID) for MFT track"};
  Configurable<float> cfgXYWindow{"XY-window", 3, "Search window (delta XY) for MFT track"};

  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType colBinning{{ConfVtxBins}, true};

  Preslice<aod::FwdTracks> perCollision = aod::fwdtrack::collisionId;
  Preslice<aod::MFTTracks> perCollisionMFT = aod::fwdtrack::collisionId;


  void init(o2::framework::InitContext&)
  {
  }

  void process(soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const& mcparticles, aod::Collisions const& collisions)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;

    for (auto const& fwdtrack : fwdtracks){
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {

        for (auto const& mfttrack : mfttracks){
          if (mfttrack.has_collision()){
            //update the talbe matchedmuonmft
            if (0 <= fwdtrack.collisionId() - mfttrack.collisionId() && fwdtrack.collisionId() - mfttrack.collisionId() < cfgColWindow) {
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

              int mcfwdfirstmotherid = -1;
              int mcfwdlastmotherid = -1;
              int mcfwdfirstmotherpdg = -1;
              int mcfwdlastmotherpdg = -1;
              int mcfwdfrombkg = -1;
              int mcmftfirstmotherid = -1;
              int mcmftlastmotherid = -1;
              int mcmftfirstmotherpdg = -1;
              int mcmftlastmotherpdg = -1;
              int mcmftfrombkg = -1;

              if (fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()) {
                const auto mcfwdtrackmothers = fwdtrack.mcParticle().mothersIds();
                const auto mcmfttrackmothers = mfttrack.mcParticle().mothersIds();

                if (mcfwdtrackmothers.size() > 0) {
                  mcfwdfirstmotherid = mcfwdtrackmothers[0];
                  auto mcfwdfirstmother = mcparticles.iteratorAt(mcfwdfirstmotherid);
                  mcfwdfirstmotherpdg = mcfwdfirstmother.pdgCode();
                  mcfwdlastmotherid = mcfwdtrackmothers[mcfwdtrackmothers.size() - 1];
                  auto mcfwdlastmother = mcparticles.iteratorAt(mcfwdlastmotherid);
                  mcfwdlastmotherpdg = mcfwdlastmother.pdgCode();
                  if (mcfwdfirstmother.fromBackgroundEvent() == true) {
                    mcfwdfrombkg = 1;
                  } else {
                    mcfwdfrombkg = 0;
                  }
                }
                if (mcmfttrackmothers.size() > 0) {
                  mcmftfirstmotherid = mcmfttrackmothers[0];
                  auto mcmftfirstmother = mcparticles.iteratorAt(mcmftfirstmotherid);
                  mcmftfirstmotherpdg = mcmftfirstmother.pdgCode();
                  mcmftlastmotherid = mcmfttrackmothers[mcmfttrackmothers.size() - 1];
                  auto mcmftlastmother = mcparticles.iteratorAt(mcmftlastmotherid);
                  mcmftlastmotherpdg = mcmftlastmother.pdgCode();
                  mcmftfrombkg = mcmftfirstmother.fromBackgroundEvent();
                  if (mcmftfirstmother.fromBackgroundEvent() == true) {
                    mcmftfrombkg = 1;
                  } else {
                    mcmftfrombkg = 0;
                  }
                }
              }

              mchmftpairTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              if (fwdtrack.sign() != mfttrack.sign()){
                mchmftpairbkgTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              }
            }
          }
        }
      }
    }
  }

  void processME(soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const& mcparticles, aod::Collisions const& collisions)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.
      if (collision2.globalIndex() - collision1.globalIndex() > cfgColWindow || collision1.globalIndex() - collision2.globalIndex() > cfgColWindow) {
        auto groupedFwdTracks1 = fwdtracks.sliceBy(perCollision, collision1.globalIndex());
        auto groupedMFTTracks1 = mfttracks.sliceBy(perCollisionMFT, collision1.globalIndex());
        auto groupedFwdTracks2 = fwdtracks.sliceBy(perCollision, collision2.globalIndex());
        auto groupedMFTTracks2 = mfttracks.sliceBy(perCollisionMFT, collision2.globalIndex());

        for (auto const& fwdtrack : groupedFwdTracks1){
          if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {

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

                int mcfwdfirstmotherid = -1;
                int mcfwdlastmotherid = -1;
                int mcfwdfirstmotherpdg = -1;
                int mcfwdlastmotherpdg = -1;
                int mcfwdfrombkg = -1;
                int mcmftfirstmotherid = -1;
                int mcmftlastmotherid = -1;
                int mcmftfirstmotherpdg = -1;
                int mcmftlastmotherpdg = -1;
                int mcmftfrombkg = -1;

                if (fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()) {
                  const auto mcfwdtrackmothers = fwdtrack.mcParticle().mothersIds();
                  const auto mcmfttrackmothers = mfttrack.mcParticle().mothersIds();

                  if (mcfwdtrackmothers.size() > 0) {
                    mcfwdfirstmotherid = mcfwdtrackmothers[0];
                    auto mcfwdfirstmother = mcparticles.iteratorAt(mcfwdfirstmotherid);
                    mcfwdfirstmotherpdg = mcfwdfirstmother.pdgCode();
                    mcfwdlastmotherid = mcfwdtrackmothers[mcfwdtrackmothers.size() - 1];
                    auto mcfwdlastmother = mcparticles.iteratorAt(mcfwdlastmotherid);
                    mcfwdlastmotherpdg = mcfwdlastmother.pdgCode();
                    if (mcfwdfirstmother.fromBackgroundEvent() == true) {
                      mcfwdfrombkg = 1;
                    } else {
                      mcfwdfrombkg = 0;
                    }
                  }
                  if (mcmfttrackmothers.size() > 0) {
                    mcmftfirstmotherid = mcmfttrackmothers[0];
                    auto mcmftfirstmother = mcparticles.iteratorAt(mcmftfirstmotherid);
                    mcmftfirstmotherpdg = mcmftfirstmother.pdgCode();
                    mcmftlastmotherid = mcmfttrackmothers[mcmfttrackmothers.size() - 1];
                    auto mcmftlastmother = mcparticles.iteratorAt(mcmftlastmotherid);
                    mcmftlastmotherpdg = mcmftlastmother.pdgCode();
                    mcmftfrombkg = mcmftfirstmother.fromBackgroundEvent();
                    if (mcmftfirstmother.fromBackgroundEvent() == true) {
                      mcmftfrombkg = 1;
                    } else {
                      mcmftfrombkg = 0;
                    }
                  }
                }

              if (0 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 20 && 0 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 20){

                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else if (20 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 40 && 20 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 40){
                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else if (40 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 60 && 40 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 60){
                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else if (60 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 80 && 60 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 80){
                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else if (80 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 100 && 80 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 100){
                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else if (100 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 120 && 100 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 120){
                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else if (120 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 140 && 120 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 140){
                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else if (140 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 160 && 140 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 160){
                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else if (160 < groupedMFTTracks1.size() && groupedMFTTracks1.size() <= 180 && 160 < groupedMFTTracks2.size() && groupedMFTTracks2.size() <= 180){
                mchmftpairbkgemTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(mchmftmatchinginfoemmc, processME, "Show displacement of mft and mch tracks", false);
 
  void processGen(soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const& mcparticles, aod::Collisions const& collisions)
  {
    static constexpr Double_t MatchingPlaneZ = -77.5;

    for (auto const& fwdtrack : fwdtracks){
      if (fwdtrack.has_collision() && fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {

        for (auto const& mfttrack : mfttracks){
          if (mfttrack.has_collision()){
            if (0 <= fwdtrack.collisionId() - mfttrack.collisionId() && fwdtrack.collisionId() - mfttrack.collisionId() < cfgColWindow) {
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

              int mcfwdfirstmotherid = -1;
              int mcfwdlastmotherid = -1;
              int mcfwdfirstmotherpdg = -1;
              int mcfwdlastmotherpdg = -1;
              int mcfwdfrombkg = -1;
              int mcmftfirstmotherid = -1;
              int mcmftlastmotherid = -1;
              int mcmftfirstmotherpdg = -1;
              int mcmftlastmotherpdg = -1;
              int mcmftfrombkg = -1;

              if (fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()) {
                const auto mcfwdtrackmothers = fwdtrack.mcParticle().mothersIds();
                const auto mcmfttrackmothers = mfttrack.mcParticle().mothersIds();

                if (mcfwdtrackmothers.size() > 0) {
                  mcfwdfirstmotherid = mcfwdtrackmothers[0];
                  auto mcfwdfirstmother = mcparticles.iteratorAt(mcfwdfirstmotherid);
                  mcfwdfirstmotherpdg = mcfwdfirstmother.pdgCode();
                  mcfwdlastmotherid = mcfwdtrackmothers[mcfwdtrackmothers.size() - 1];
                  auto mcfwdlastmother = mcparticles.iteratorAt(mcfwdlastmotherid);
                  mcfwdlastmotherpdg = mcfwdlastmother.pdgCode();
                  if (mcfwdfirstmother.fromBackgroundEvent() == true) {
                    mcfwdfrombkg = 1;
                  } else {
                    mcfwdfrombkg = 0;
                  }
                }
                if (mcmfttrackmothers.size() > 0) {
                  mcmftfirstmotherid = mcmfttrackmothers[0];
                  auto mcmftfirstmother = mcparticles.iteratorAt(mcmftfirstmotherid);
                  mcmftfirstmotherpdg = mcmftfirstmother.pdgCode();
                  mcmftlastmotherid = mcmfttrackmothers[mcmfttrackmothers.size() - 1];
                  auto mcmftlastmother = mcparticles.iteratorAt(mcmftlastmotherid);
                  mcmftlastmotherpdg = mcmftlastmother.pdgCode();
                  mcmftfrombkg = mcmftfirstmother.fromBackgroundEvent();
                  if (mcmftfirstmother.fromBackgroundEvent() == true) {
                    mcmftfrombkg = 1;
                  } else {
                    mcmftfrombkg = 0;
                  }
                }
              }

              if (fwdtrack.mcParticleId() == mfttrack.mcParticleId())
              {
                //Write True MCH-MFT pair Table
              mchmftpairtrueTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              } else {
                //Write Wrong MCH-MFT pair Table
              mchmftpairwrongTable(mftpars1.getX(), mftpars1.getY(), mftpars1.getEta(), mftpars1.getPhi(),mftpars1.getTanl(), mftpars1.getPt(),std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)),mfttrack.chi2(), mcmftfirstmotherid,mcmftfirstmotherpdg,mcmftlastmotherid,mcmftlastmotherpdg,mcmftfrombkg, muonpars1.getX(), muonpars1.getY(), muonpars1.getEta(), muonpars1.getPhi(), muonpars1.getTanl(), muonpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)),mcfwdfirstmotherid,mcfwdfirstmotherpdg,mcfwdlastmotherid,mcfwdlastmotherpdg,mcfwdfrombkg, muonpars1.getX() - mftpars1.getX(),muonpars1.getY() - mftpars1.getY(),muonpars1.getEta() - mftpars1.getEta(),muonpars1.getPhi() - mftpars1.getPhi(),muonpars1.getTanl() - mftpars1.getTanl(), muonpars1.getPt() - mftpars1.getPt(),std::sqrt(std::pow(muonpars1.getX(), 2) + std::pow(muonpars1.getY(), 2)) - std::sqrt(std::pow(mftpars1.getX(), 2) + std::pow(mftpars1.getY(), 2)));
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(mchmftmatchinginfoemmc, processGen, "Show displacement of mft and mch tracks", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mchmftmatchinginfoemmc>(cfgc)
  };
}
