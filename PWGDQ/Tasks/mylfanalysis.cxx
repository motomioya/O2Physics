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
#include "PWGDQ/Core/VarManager.h"
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

struct myvmanalysis {

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
      {"hetaMassPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hetaMassMcPM", "Invariant;Invariant MassMc (GeV/c^{2});MassMc", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hetaDCAxyPM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"hetaLxyzPM", "Invariant;Invariant Mass (GeV/c^{2});Lxyz", {HistType::kTH1F, {{1000, 0, 8}}}},
      {"hetaPCAchiPM", "Invariant;Invariant Mass (GeV/c^{2});PCA(chi2)", {HistType::kTH1F, {{1000, 0., 6}}}},
      {"hetaprimeMassPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hetaprimeMassMcPM", "Invariant;Invariant MassMc (GeV/c^{2});MassMc", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hetaprimeDCAxyPM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"hetaprimeLxyzPM", "Invariant;Invariant Mass (GeV/c^{2});Lxyz", {HistType::kTH1F, {{1000, 0, 8}}}},
      {"hetaprimePCAchiPM", "Invariant;Invariant Mass (GeV/c^{2});PCA(chi2)", {HistType::kTH1F, {{1000, 0., 6}}}},
      {"hrhoMassPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hrhoMassMcPM", "Invariant;Invariant MassMc (GeV/c^{2});MassMc", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hrhoDCAxyPM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"hrhoLxyzPM", "Invariant;Invariant Mass (GeV/c^{2});Lxyz", {HistType::kTH1F, {{1000, 0, 8}}}},
      {"hrhoPCAchiPM", "Invariant;Invariant Mass (GeV/c^{2});PCA(chi2)", {HistType::kTH1F, {{1000, 0., 6}}}},
      {"homegaMassPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"homegaMassMcPM", "Invariant;Invariant MassMc (GeV/c^{2});MassMc", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"homegaDCAxyPM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"homegaLxyzPM", "Invariant;Invariant Mass (GeV/c^{2});Lxyz", {HistType::kTH1F, {{1000, 0, 8}}}},
      {"homegaPCAchiPM", "Invariant;Invariant Mass (GeV/c^{2});PCA(chi2)", {HistType::kTH1F, {{1000, 0., 6}}}},
      {"hphiMassPM", "Invariant;Invariant Mass (GeV/c^{2});Mass", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hphiMassMcPM", "Invariant;Invariant MassMc (GeV/c^{2});MassMc", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hphiDCAxyPM", "Invariant;Invariant Mass (GeV/c^{2});Pair DCA", {HistType::kTH1F, {{1000, 0.0, 2}}}},
      {"hphiLxyzPM", "Invariant;Invariant Mass (GeV/c^{2});Lxyz", {HistType::kTH1F, {{1000, 0, 8}}}},
      {"hphiPCAchiPM", "Invariant;Invariant Mass (GeV/c^{2});PCA(chi2)", {HistType::kTH1F, {{1000, 0., 6}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>("GLO/Config/GeometryAligned");
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, soa::Filtered<MyMuons> const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const& particles)
  {
    o2::vertexing::FwdDCAFitterN<2> fgFitterTwoProngFwd;

    bool propagateToPCA = true;
    float maxR = 200;
    float minParamChange = 1.0e-3;
    float minRelChi2Change = 0.9;
    bool useAbsDCA = false;
    //double singleptcut = 0.5;

    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    fgFitterTwoProngFwd.setTGeoMat(false);
    fgFitterTwoProngFwd.setMatLUT(lut);
    fgFitterTwoProngFwd.setBz(mMagField);
    fgFitterTwoProngFwd.setPropagateToPCA(propagateToPCA);
    fgFitterTwoProngFwd.setMaxR(maxR);
    fgFitterTwoProngFwd.setMinParamChange(minParamChange);
    fgFitterTwoProngFwd.setMinRelChi2Change(minRelChi2Change);
    fgFitterTwoProngFwd.setUseAbsDCA(useAbsDCA);

    MCSignal* signaleta;
    MCProng prongeta(2, {13, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongeta.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signaleta = new MCSignal("signaleta", "ee pairs from eta decays", {prongeta, prongeta}, {1, 1}); // signal at pair level
                                                                                              //
    MCSignal* signaletaprime;
    MCProng prongetaprime(2, {13, 331}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongetaprime.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signaletaprime = new MCSignal("signaletaprime", "ee pairs from etaprime decays", {prongetaprime, prongetaprime}, {1, 1}); // signal at pair level

    MCSignal* signalrho;
    MCProng prongrho(2, {13, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongrho.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalrho = new MCSignal("signalrho", "ee pairs from rho decays", {prongrho, prongrho}, {1, 1}); // signal at pair level
                                                                                              
    MCSignal* signalomega;
    MCProng prongomega(2, {13, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongomega.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalomega = new MCSignal("signalomega", "ee pairs from omega decays", {prongomega, prongomega}, {1, 1}); // signal at pair level
                                                                                                      
    MCSignal* signalphi;
    MCProng prongphi(2, {13, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongphi.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signalphi = new MCSignal("signalphi", "ee pairs from phi decays", {prongphi, prongphi}, {1, 1}); // signal at pair level

    for (auto& collision : collisions) {
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
                double chi2PCA = -999;
                double VertexingSV = -999;
                double VertexingLxyz = -999;
                /*
                if (procCode == 0 ) {
                  continue;
                }
                */
                if (procCode != 0 ) {
                  Vec3D secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
                  fgFitterTwoProngFwd.calcPCACovMatrixFlat();
                  chi2PCA = fgFitterTwoProngFwd.getChi2AtPCACandidate();
                  auto VertexingLxy = (collision.posX() - secondaryVertex[0]) * (collision.posX() - secondaryVertex[0]) + (collision.posY() - secondaryVertex[1]) * (collision.posY() - secondaryVertex[1]);
                  auto VertexingLz = (collision.posZ() - secondaryVertex[2]) * (collision.posZ() - secondaryVertex[2]);
                  VertexingLxyz = VertexingLxy + VertexingLz;
                  VertexingLxy = std::sqrt(VertexingLxy);
                  VertexingLz = std::sqrt(VertexingLz);
                  VertexingLxyz = std::sqrt(VertexingLxyz);
                  VertexingSV = secondaryVertex[2];
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
                //double centerMFT[3] = {0, 0, -61.4};
                //o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
                //auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
                //auto Bz = -5.;

                /*
                o2::dataformats::GlobalFwdTrack propmuon1;
                auto geoMan1 = o2::base::GeometryManager::meanMaterialBudget(fwdtrack1.x(), fwdtrack1.y(), fwdtrack1.z(), collision.posX(), collision.posY(), collision.posZ());
                auto x2x01 = static_cast<float>(geoMan1.meanX2X0);
                pars1.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x01);
                propmuon1.setParameters(pars1.getParameters());
                propmuon1.setZ(pars1.getZ());
                propmuon1.setCovariances(pars1.getCovariances());

                o2::dataformats::GlobalFwdTrack propmuon2;
                auto geoMan2 = o2::base::GeometryManager::meanMaterialBudget(fwdtrack2.x(), fwdtrack2.y(), fwdtrack2.z(), collision.posX(), collision.posY(), collision.posZ());
                auto x2x02 = static_cast<float>(geoMan2.meanX2X0);
                pars2.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x02);
                propmuon2.setParameters(pars2.getParameters());
                propmuon2.setZ(pars2.getZ());
                propmuon2.setCovariances(pars2.getCovariances());
                */

                //calculate pair DCA
                float fwd1dcaX = (pars1.getX() - collision.posX());
                float fwd1dcaY = (pars1.getY() - collision.posY());
                float fwd2dcaX = (pars2.getX() - collision.posX());
                float fwd2dcaY = (pars2.getY() - collision.posY());
                //float fwd1dcaX = (propmuon1.getX() - collision.posX());
                //float fwd1dcaY = (propmuon1.getY() - collision.posY());
                //float fwd2dcaX = (propmuon2.getX() - collision.posX());
                //float fwd2dcaY = (propmuon2.getY() - collision.posY());
                float DCA1 = std::sqrt(fwd1dcaX * fwd1dcaX + fwd1dcaY * fwd1dcaY);
                float DCA2 = std::sqrt(fwd2dcaX * fwd2dcaX + fwd2dcaY * fwd2dcaY);
                float DCAmumu = std::sqrt((DCA1 * DCA1 + DCA2 * DCA2)/2);

                //get MC particles
                if (fwdtrack1.has_mcParticle() == 1 && fwdtrack2.has_mcParticle() == 1) {
                  auto fwdparticle1 = fwdtrack1.mcParticle();
                  auto fwdparticle2 = fwdtrack2.mcParticle();
                  const auto mcmothersidlist1 = fwdparticle1.mothersIds();
                  const auto mcmothersidlist2 = fwdparticle2.mothersIds();
                  TLorentzVector lvmc1, lvmc2, lvmc;
                  lvmc1.SetPtEtaPhiM(fwdparticle1.pt(), fwdparticle1.eta(), fwdparticle1.phi(), muonMass);
                  lvmc2.SetPtEtaPhiM(fwdparticle2.pt(), fwdparticle2.eta(), fwdparticle2.phi(), muonMass);
                  lvmc = lvmc1 + lvmc2;
                  if (fwdparticle1.globalIndex() != fwdparticle2.globalIndex()) {
                    //svcut
                    if(VertexingSV > -20 || VertexingSV == -999) {
                      if (fwdtrack1.signed1Pt()*fwdtrack2.signed1Pt() < 0) {
                        if(signaleta->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hetaMassPM"), lv.M());
                          registry.fill(HIST("hetaMassMcPM"), lvmc.M());
                          registry.fill(HIST("hetaDCAxyPM"), DCAmumu);
                          registry.fill(HIST("hetaLxyzPM"), VertexingLxyz);
                          registry.fill(HIST("hetaPCAchiPM"), chi2PCA);
                        } 
                        if(signaletaprime->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hetaprimeMassPM"), lv.M());
                          registry.fill(HIST("hetaprimeMassMcPM"), lvmc.M());
                          registry.fill(HIST("hetaprimeDCAxyPM"), DCAmumu);
                          registry.fill(HIST("hetaprimeLxyzPM"), VertexingLxyz);
                          registry.fill(HIST("hetaprimePCAchiPM"), chi2PCA);
                        }
                        if(signalrho->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hrhoMassPM"), lv.M());
                          registry.fill(HIST("hrhoMassMcPM"), lvmc.M());
                          registry.fill(HIST("hrhoDCAxyPM"), DCAmumu);
                          registry.fill(HIST("hrhoLxyzPM"), VertexingLxyz);
                          registry.fill(HIST("hrhoPCAchiPM"), chi2PCA);
                        }
                        if(signalomega->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("homegaMassPM"), lv.M());
                          registry.fill(HIST("homegaMassMcPM"), lvmc.M());
                          registry.fill(HIST("hrhoDCAxyPM"), DCAmumu);
                          registry.fill(HIST("homegaDCAxyPM"), DCAmumu);
                          registry.fill(HIST("homegaLxyzPM"), VertexingLxyz);
                          registry.fill(HIST("homegaPCAchiPM"), chi2PCA);
                        }
                        if(signalphi->CheckSignal(true, fwdparticle1, fwdparticle2)) {
                          registry.fill(HIST("hphiMassPM"), lv.M());
                          registry.fill(HIST("hphiMassMcPM"), lvmc.M());
                          registry.fill(HIST("hphiDCAxyPM"), DCAmumu);
                          registry.fill(HIST("hphiLxyzPM"), VertexingLxyz);
                          registry.fill(HIST("hphiPCAchiPM"), chi2PCA);
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
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myvmanalysis>(cfgc)
  };
}
