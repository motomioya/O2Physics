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
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include <math.h>
#include <TLorentzVector.h>
#include <string>
#include <regex>
#include <iostream>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"

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
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTrkCompColls, aod::McFwdTrackLabels>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;

struct myhfanalysis {

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
  float collisionZcut = 10.0f;

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);
  Filter collisionFilter = nabs(aod::collision::posZ) < collisionZcut;
  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();

  Preslice<aod::McParticles> particlesIndicesPerCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = aod::track_association::collisionId;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float mMagField = -5.0;
  //o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;


  HistogramRegistry registry{
    "registry", 
    {
      {"hRecMassvsPtPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtbbPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtbbMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtbbPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtccPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtccMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtccPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtjpsiPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtjpsiMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtjpsiPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtetaPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtetaMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtetaPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtetaprimePM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtetaprimeMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtetaprimePP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtrhoPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtrhoMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtrhoPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtomegaPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtomegaMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtomegaPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtphiPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtphiMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hRecMassvsPtphiPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},

      {"hGenMassvsPtPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtbbPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtccPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtccMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtccPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtjpsiPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtjpsiMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtjpsiPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtetaPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtetaMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtetaPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtetaprimePM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtetaprimeMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtetaprimePP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtrhoPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtrhoMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtrhoPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtomegaPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtomegaMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtomegaPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtphiPM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtphiMM", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
      {"hGenMassvsPtphiPP", "Invariant;Invariant Mass (GeV/c^{2});Pair Pt", {HistType::kTH2F, {{750, 0.0, 15.0}, {40, 0, 10}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, soa::Filtered<MyMuons> const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const& particles)
  {
    o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;

    bool propagateToPCA = true;
    float maxR = 200;
    float minParamChange = 1.0e-3;
    float minRelChi2Change = 0.9;
    bool useAbsDCA = false;

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    fgFitterTwoProngFwd.setTGeoMat(false);
    fgFitterTwoProngFwd.setMatLUT(lut);
    fgFitterTwoProngFwd.setBz(mMagField);
    fgFitterTwoProngFwd.setPropagateToPCA(propagateToPCA);
    fgFitterTwoProngFwd.setMaxR(maxR);
    fgFitterTwoProngFwd.setMinParamChange(minParamChange);
    fgFitterTwoProngFwd.setMinRelChi2Change(minRelChi2Change);
    fgFitterTwoProngFwd.setUseAbsDCA(useAbsDCA);

    //Configuring signals bb -> mumu
    MCSignal* signalbb;
    MCProng prongbmu(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongbmu.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb = new MCSignal("bbTomumu", "mumu pairs from b->mu and b->mu", {prongbmu, prongbmu}, {-1, -1});
    //Configuring signals b->cmu->mumu
    MCSignal* signalbb1;
    MCProng prongbmu1A(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongbmu1B(3, {13, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prongbmu1A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu1B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb1 = new MCSignal("bbcTomumusingle", "mumu pairs from b->e and b->c->e (single b)", {prongbmu1A, prongbmu1B}, {1, 2}); // signal at pair level
    //Configuring signals b->cmu->mumu, inversed
    MCSignal* signalbb2;
    MCProng prongbmu2A(3, {13, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongbmu2B(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongbmu2A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu2B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb2 = new MCSignal("bcbTomumusingle", "mumu pairs from b->e and b->c->e (single b), inversed", {prongbmu2A, prongbmu2B}, {2, 1}); // signal at pair level
    //Configuring signals b->cmu->mumu
    MCSignal* signalbb3;
    MCProng prongbmu3A(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongbmu3B(4, {13, 402, 402, 502}, {true, true, true, true}, {false, false, false, false}, {0, 0, 0, 0}, {0, 0, 0, 0}, {false, false, false, false});
    prongbmu3A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu3B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb3 = new MCSignal("bbcTomumusingle", "mumu pairs from b->e and b->c->e (single b)", {prongbmu3A, prongbmu3B}, {1, 3}); // signal at pair level
    //Configuring signals b->cmu->mumu, inversed
    MCSignal* signalbb4;
    MCProng prongbmu4A(4, {13, 402, 402, 502}, {true, true, true, true}, {false, false, false, false}, {0, 0, 0, 0}, {0, 0, 0, 0}, {false, false, false, false});
    MCProng prongbmu4B(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongbmu4A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu4B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb4 = new MCSignal("bcbTomumusingle", "mumu pairs from b->e and b->c->e (single b), inversed", {prongbmu4A, prongbmu4B}, {3, 1}); // signal at pair level
                                                                                                                                            //
    MCSignal* signalbb5;
    MCProng prongbmu5(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prongbmu5.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb5 = new MCSignal("bcbcTomumu", "mumu pairs from b->c->mu and b->c->mu", {prongbmu5, prongbmu5}, {-1, -1});

    MCSignal* signalbb6;
    MCProng prongbmu6A(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); // check if mother pdg code is in history
    MCProng prongbmu6B(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prongbmu6A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu6B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb6 = new MCSignal("bbcTomumu", "mumu pairs from b->mu and b->c->mu", {prongbmu6A, prongbmu6B}, {-1, -1}); // signal at pair level
                                                                                                                    
    MCSignal* signalbb7;
    MCProng prongbmu7A(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    MCProng prongbmu7B(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); // check if mother pdg code is in history
    prongbmu7A.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongbmu7B.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalbb7 = new MCSignal("bcbTomumu", "mumu pairs from b->c->mu and b->mu", {prongbmu7A, prongbmu7B}, {-1, -1}); // signal at pair level

    //Configuring signals cc -> mumu
    MCSignal* signalcc;
    MCProng prongcmu(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongcmu.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalcc = new MCSignal("ccTomumu", "mumu pairs from c->mu and c->mu", {prongcmu, prongcmu}, {-1, -1});

    //J/Psi
    MCSignal* signaljpsi;
    MCProng prongjpsi(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signaljpsi = new MCSignal("Jpsitomumu", "mumu pairs from j/psi decays", {prongjpsi, prongjpsi}, {1, 1}); // signal at pair level
    //LF

    MCSignal* signaleta;
    MCProng prongeta(2, {13, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongeta.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signaleta = new MCSignal("etatomumu", "mumu pairs from eta decays", {prongeta, prongeta}, {1, 1}); // signaleta at pair level

    MCSignal* signaletaprime;
    MCProng prongetaprime(2, {13, 331}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongetaprime.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signaletaprime = new MCSignal("etaprimetomumu", "mumu pairs from eta' decays", {prongetaprime, prongetaprime}, {1, 1}); // signal at pair level

    MCSignal* signalrho;
    MCProng prongrho(2, {13, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongrho.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalrho = new MCSignal("rhotomumu", "mumu pairs from rho decays", {prongrho, prongrho}, {1, 1}); // signal at pair level

    MCSignal* signalomega;
    MCProng prongomega(2, {13, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongomega.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalomega = new MCSignal("omegatomumu", "mumu pairs from omega decays", {prongomega, prongomega}, {1, 1}); // signal at pair level
                                                                                                                 //
    MCSignal* signalphi;
    MCProng prongphi(2, {13, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongphi.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalphi = new MCSignal("phitomumu", "mumu pairs from phi decays", {prongphi, prongphi}, {1, 1}); // signal at pair level


    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) { 
        continue;
      }
      //auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      //grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      //mMagField = grpmag->getNominalL3Field();
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& [fwdtrackId1, fwdtrackId2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(muonIdsThisCollision, muonIdsThisCollision))) {
        for (auto& fwdtrack1 : fwdtracks) {
          if (fwdtrack1.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId1.fwdtrackId() == fwdtrack1.globalIndex() && fwdtrack1.chi2MatchMCHMFT() < 50 && fwdtrack1.compatibleCollIds().size() == 1) {
            for (auto& fwdtrack2 : fwdtracks) {
              if (fwdtrack2.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId2.fwdtrackId() == fwdtrack2.globalIndex() && fwdtrack2.chi2MatchMCHMFT() < 50 && fwdtrack2.compatibleCollIds().size() == 1) {
                //calculate mass
                TLorentzVector lv1, lv2, lv;
                lv1.SetPtEtaPhiM(fwdtrack1.pt(), fwdtrack1.eta(), fwdtrack1.phi(), muonMass);
                lv2.SetPtEtaPhiM(fwdtrack2.pt(), fwdtrack2.eta(), fwdtrack2.phi(), muonMass);
                lv = lv1 + lv2;

                //set TrackParCovFwd for fwdtrack pairs
                double chi21 = fwdtrack1.chi2();
                SMatrix5 tpars1(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.phi(), fwdtrack1.tgl(), fwdtrack1.signed1Pt());
                std::vector<double> v11{fwdtrack1.cXX(), fwdtrack1.cXY(), fwdtrack1.cYY(), fwdtrack1.cPhiX(), fwdtrack1.cPhiY(), fwdtrack1.cPhiPhi(), fwdtrack1.cTglX(), fwdtrack1.cTglY(), fwdtrack1.cTglPhi(), fwdtrack1.cTglTgl(), fwdtrack1.c1PtX(), fwdtrack1.c1PtY(), fwdtrack1.c1PtPhi(), fwdtrack1.c1PtTgl(), fwdtrack1.c1Pt21Pt2()};
                SMatrix55 tcovs1(v11.begin(), v11.end());
                /*
                float detXY1 = fwdtrack1.cXX() * fwdtrack1.cYY() - fwdtrack1.cXY() * fwdtrack1.cXY();
                if (detXY1 <= 0) {
                  continue;
                }
                */
                o2::track::TrackParCovFwd pars1{fwdtrack1.z(), tpars1, tcovs1, chi21};

                double chi22 = fwdtrack2.chi2();
                SMatrix5 tpars2(fwdtrack2.x(), fwdtrack2.y(), fwdtrack2.phi(), fwdtrack2.tgl(), fwdtrack2.signed1Pt());
                std::vector<double> v12{fwdtrack2.cXX(), fwdtrack2.cXY(), fwdtrack2.cYY(), fwdtrack2.cPhiX(), fwdtrack2.cPhiY(), fwdtrack2.cPhiPhi(), fwdtrack2.cTglX(), fwdtrack2.cTglY(), fwdtrack2.cTglPhi(), fwdtrack2.cTglTgl(), fwdtrack2.c1PtX(), fwdtrack2.c1PtY(), fwdtrack2.c1PtPhi(), fwdtrack2.c1PtTgl(), fwdtrack2.c1Pt21Pt2()};
                SMatrix55 tcovs2(v12.begin(), v12.end());
                /*
                float detXY2 = fwdtrack2.cXX() * fwdtrack2.cYY() - fwdtrack2.cXY() * fwdtrack2.cXY();
                if (detXY2 <= 0) {
                  continue;
                }
                */
                o2::track::TrackParCovFwd pars2{fwdtrack2.z(), tpars2, tcovs2, chi22};

                //Get secondary vertex using DCAFitterN
                int procCode = fgFitterTwoProngFwd.process(pars1, pars2);
                //double chi2PCA = -999;
                //double VertexingSV = -999;
                //double VertexingLxyz = -999;
                if (procCode != 0 ) {
                  //Vec3D secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
                  fgFitterTwoProngFwd.calcPCACovMatrixFlat();
                  //chi2PCA = fgFitterTwoProngFwd.getChi2AtPCACandidate();
                  //auto VertexingLxy = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) + (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
                  //auto VertexingLz = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
                  //VertexingLxyz = VertexingLxy + VertexingLz;
                  //VertexingLxy = std::sqrt(VertexingLxy);
                  //VertexingLz = std::sqrt(VertexingLz);
                  //VertexingLxyz = std::sqrt(VertexingLxyz);
                  //VertexingSV = secondaryVertex[2];
                } 

                //Get secondary vertex using KFparticles
                /*
                KFPTrack kfpTrack0 = createKFPFwdTrackFromFwdTrack(fwdtrack1);
                trk0KF = KFParticle(kfpTrack0, -13 * fwdtrack1.sign());
                KFPTrack kfpTrack1 = createKFPFwdTrackFromFwdTrack(fwdtrack2);
                trk1KF = KFParticle(kfpTrack1, -13 * fwdtrack2.sign());
                KFGeoTwoProng.SetConstructMethod(2);
                KFGeoTwoProng.AddDaughter(trk0KF);
                KFGeoTwoProng.AddDaughter(trk1KF);
                
                KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
                KFParticle KFPV(kfpVertex);
                double dxPair2PV = KFGeoTwoProng.GetX() - KFPV.GetX();
                double dyPair2PV = KFGeoTwoProng.GetY() - KFPV.GetY();
                double dzPair2PV = KFGeoTwoProng.GetZ() - KFPV.GetZ();
                auto KFVertexingLxy = std::sqrt(dxPair2PV * dxPair2PV + dyPair2PV * dyPair2PV);
                auto KFVertexingLz = std::sqrt(dzPair2PV * dzPair2PV);
                auto KFVertexingLxyz = std::sqrt(dxPair2PV * dxPair2PV + dyPair2PV * dyPair2PV + dzPair2PV * dzPair2PV);
                auto KFVertexingSV = KFGeoTwoProng.GetZ();
                auto kKFDCAxyzBetweenProngs = trk0KF.GetDistanceFromParticle(trk1KF);
                */

                //fwdtrack propagation to collision
                pars1.propagateToZlinear(collision.posZ());
                pars2.propagateToZlinear(collision.posZ());

                //get MC particles
                if (fwdtrack1.has_mcParticle() == 1 && fwdtrack2.has_mcParticle() == 1) {
                  auto fwdparticle1 = fwdtrack1.mcParticle();
                  auto fwdparticle2 = fwdtrack2.mcParticle();
                  const auto mcmothersidlist1 = fwdparticle1.mothersIds();
                  const auto mcmothersidlist2 = fwdparticle2.mothersIds();
                  if (fwdparticle1.globalIndex() != fwdparticle2.globalIndex()) {
                    if (fwdtrack1.signed1Pt()*fwdtrack2.signed1Pt() < 0) {
                      //fill mass vs DCA
                      registry.fill(HIST("hRecMassvsPtPM"), lv.M(), lv.Pt());
                      //select muons from HF semileptonic decays
                      if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtbbPM"), lv.M(), lv.Pt());
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtccPM"), lv.M(), lv.Pt());
                      } else if (signaljpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtjpsiPM"), lv.M(), lv.Pt());
                      } else if (signaleta->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtetaPM"), lv.M(), lv.Pt());
                      } else if (signaletaprime->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtetaprimePM"), lv.M(), lv.Pt());
                      } else if (signalrho->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtrhoPM"), lv.M(), lv.Pt());
                      } else if (signalomega->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtomegaPM"), lv.M(), lv.Pt());
                      } else if (signalphi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtphiPM"), lv.M(), lv.Pt());
                      }
                    } else if (fwdtrack1.signed1Pt() >= 0 && fwdtrack2.signed1Pt() >= 0) {
                      //fill mass vs DCA
                      registry.fill(HIST("hRecMassvsPtPP"), lv.M(), lv.Pt());
                      //select muons from HF semileptonic decays
                      if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtbbPP"), lv.M(), lv.Pt());
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtccPP"), lv.M(), lv.Pt());
                      } else if (signaljpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtjpsiPP"), lv.M(), lv.Pt());
                      } else if (signaleta->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtetaPP"), lv.M(), lv.Pt());
                      } else if (signaletaprime->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtetaprimePP"), lv.M(), lv.Pt());
                      } else if (signalrho->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtrhoPP"), lv.M(), lv.Pt());
                      } else if (signalomega->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtomegaPP"), lv.M(), lv.Pt());
                      } else if (signalphi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtphiPP"), lv.M(), lv.Pt());
                      }
                    } else if (fwdtrack1.signed1Pt() < 0 && fwdtrack2.signed1Pt() < 0) {
                      //fill mass vs DCA
                      registry.fill(HIST("hRecMassvsPtMM"), lv.M(), lv.Pt());
                      //select muons from HF semileptonic decays
                      if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtbbMM"), lv.M(), lv.Pt());
                      } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtccMM"), lv.M(), lv.Pt());
                      } else if (signaljpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtjpsiMM"), lv.M(), lv.Pt());
                      } else if (signaleta->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtetaMM"), lv.M(), lv.Pt());
                      } else if (signaletaprime->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtetaprimeMM"), lv.M(), lv.Pt());
                      } else if (signalrho->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtrhoMM"), lv.M(), lv.Pt());
                      } else if (signalomega->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtomegaMM"), lv.M(), lv.Pt());
                      } else if (signalphi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                        registry.fill(HIST("hRecMassvsPtphiMM"), lv.M(), lv.Pt());
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      //auto mccolision = collision.mcCollision();
      auto particlessThisCollision = particles.sliceBy(particlesIndicesPerCollision, collision.globalIndex());
      for (auto& [fwdparticle1, fwdparticle2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(particlessThisCollision, particlessThisCollision))) {
        if (fwdparticle1.pdgCode() == 13 || fwdparticle1.pdgCode() == -13 ) {
          if (fwdparticle2.pdgCode() == 13 || fwdparticle2.pdgCode() == -13 ) {
            //calculate mass
            if (fwdparticle1.eta() > -4 && fwdparticle1.eta() < -2.5) {
              if (fwdparticle2.eta() > -4 && fwdparticle2.eta() < -2.5) {
                TLorentzVector lv1, lv2, lv;
                lv1.SetPtEtaPhiM(fwdparticle1.pt(), fwdparticle1.eta(), fwdparticle1.phi(), muonMass);
                lv2.SetPtEtaPhiM(fwdparticle2.pt(), fwdparticle2.eta(), fwdparticle2.phi(), muonMass);
                lv = lv1 + lv2;

                if (fwdparticle1.globalIndex() != fwdparticle2.globalIndex()) {
                  if (fwdparticle1.pdgCode()*fwdparticle2.pdgCode() < 0) {
                    //fill mass vs DCA
                    registry.fill(HIST("hGenMassvsPtPM"), lv.M(), lv.Pt());
                    //select muons from HF semileptonic decays
                    if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPM"), lv.M(), lv.Pt());
                    } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtccPM"), lv.M(), lv.Pt());
                    } else if (signaljpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtjpsiPM"), lv.M(), lv.Pt());
                    } else if (signaleta->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtetaPM"), lv.M(), lv.Pt());
                    } else if (signaletaprime->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtetaprimePM"), lv.M(), lv.Pt());
                    } else if (signalrho->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtrhoPM"), lv.M(), lv.Pt());
                    } else if (signalomega->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtomegaPM"), lv.M(), lv.Pt());
                    } else if (signalphi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtphiPM"), lv.M(), lv.Pt());
                    }
                  } else if (fwdparticle1.pdgCode() > 0 && fwdparticle2.pdgCode() > 0) {
                    //fill mass vs DCA
                    registry.fill(HIST("hGenMassvsPtPP"), lv.M(), lv.Pt());
                    //select muons from HF semileptonic decays
                    if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbPP"), lv.M(), lv.Pt());
                    } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtccPP"), lv.M(), lv.Pt());
                    } else if (signaljpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtjpsiPP"), lv.M(), lv.Pt());
                    } else if (signaleta->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtetaPP"), lv.M(), lv.Pt());
                    } else if (signaletaprime->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtetaprimePP"), lv.M(), lv.Pt());
                    } else if (signalrho->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtrhoPP"), lv.M(), lv.Pt());
                    } else if (signalomega->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtomegaPP"), lv.M(), lv.Pt());
                    } else if (signalphi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtphiPP"), lv.M(), lv.Pt());
                    }
                  } else if (fwdparticle1.pdgCode() < 0 && fwdparticle2.pdgCode() < 0) {
                    //fill mass vs DCA
                    registry.fill(HIST("hGenMassvsPtMM"), lv.M(), lv.Pt());
                    //select muons from HF semileptonic decays
                    if(signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb1->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb2->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb3->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb4->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb5->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb6->CheckSignal(true, fwdparticle1, fwdparticle2) || signalbb7->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtbbMM"), lv.M(), lv.Pt());
                    } else if (signalcc->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtccMM"), lv.M(), lv.Pt());
                    } else if (signaljpsi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtjpsiMM"), lv.M(), lv.Pt());
                    } else if (signaleta->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtetaMM"), lv.M(), lv.Pt());
                    } else if (signaletaprime->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtetaprimeMM"), lv.M(), lv.Pt());
                    } else if (signalrho->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtrhoMM"), lv.M(), lv.Pt());
                    } else if (signalomega->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtomegaMM"), lv.M(), lv.Pt());
                    } else if (signalphi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                      registry.fill(HIST("hGenMassvsPtphiMM"), lv.M(), lv.Pt());
                    }
                  }
                }
              }
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
    adaptAnalysisTask<myhfanalysis>(cfgc)
  };
}
