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

struct checklikesign{

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
      {"hMassUS", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassLS", "Invariant Mass (++/--);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassUScorrectmatch", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassLScorrectmatch", "Invariant Mass (++/--);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassUSonesidematch", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassLSonesidematch", "Invariant Mass (++/--);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassUSwrongmatch", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassLSwrongmatch", "Invariant Mass (++/--);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassUSdifferentmft", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassLSdifferentmft", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassUSsamemft", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassLSsamemft", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassUSdifferentmftrec", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassLSdifferentmftrec", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassUSsamemftrec", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hMassLSsamemftrec", "Invariant Mass (+-);Invariant Mass (GeV/c^{2})", {HistType::kTH1F, {{750, 0.0, 15.0}}}},
      {"hSVUS", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVLS", "Secondary Vertex (++/--);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVUScorrectmatch", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVLScorrectmatch", "Secondary Vertex (++/--);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVUSonesidematch", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVLSonesidematch", "Secondary Vertex (++/--);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVUSwrongmatch", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVLSwrongmatch", "Secondary Vertex (++/--);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVUSdifferentmft", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVLSdifferentmft", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVUSsamemft", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVLSsamemft", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVUSdifferentmftrec", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVLSdifferentmftrec", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVUSsamemftrec", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
      {"hSVLSsamemftrec", "Secondary Vertex (+-);Secondary Vertex (cm)", {HistType::kTH1F, {{1100, -1000, 100}}}},
    },
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, MyMuons const& fwdtracks, aod::FwdTrackAssoc const& fwdtrackIndices, aod::McCollisions const&, aod::McParticles const& particles, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks)
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

    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) { 
        continue;
      }
      auto muonIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (auto& [fwdtrackId1, fwdtrackId2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(muonIdsThisCollision, muonIdsThisCollision))) {
        for (auto& fwdtrack1 : fwdtracks) {
          if (fwdtrack1.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId1.fwdtrackId() == fwdtrack1.globalIndex() && fwdtrack1.chi2MatchMCHMFT() < 50 && fwdtrack1.compatibleCollIds().size() == 1) {
            for (auto& fwdtrack2 : fwdtracks) {
              if (fwdtrack2.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack && fwdtrackId2.fwdtrackId() == fwdtrack2.globalIndex() && fwdtrack2.chi2MatchMCHMFT() < 50 && fwdtrack2.compatibleCollIds().size() == 1) {
                if (fwdtrack1.eta() > -4 && fwdtrack1.eta() < -2.5 && (((17.6 < fwdtrack1.rAtAbsorberEnd()) && (fwdtrack1.rAtAbsorberEnd() < 26.5) && (fwdtrack1.pDca() < 594)) || ((26.5 < fwdtrack1.rAtAbsorberEnd()) && (fwdtrack1.rAtAbsorberEnd() < 89.5) && (fwdtrack1.pDca() < 324)))){
                  if (fwdtrack2.eta() > -4 && fwdtrack2.eta() < -2.5 && (((17.6 < fwdtrack2.rAtAbsorberEnd()) && (fwdtrack2.rAtAbsorberEnd() < 26.5) && (fwdtrack2.pDca() < 594)) || ((26.5 < fwdtrack2.rAtAbsorberEnd()) && (fwdtrack2.rAtAbsorberEnd() < 89.5) && (fwdtrack2.pDca() < 324)))){
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
                    double VertexingSV = -999;
                    if (procCode != 0 ) {
                      Vec3D secondaryVertex = fgFitterTwoProngFwd.getPCACandidate();
                      fgFitterTwoProngFwd.calcPCACovMatrixFlat();
                      VertexingSV = secondaryVertex[2];
                    } 

                    //fwdtrack propagation to collision
                    pars1.propagateToZlinear(collision.posZ());
                    pars2.propagateToZlinear(collision.posZ());

                    //get MC particles
                    if (fwdtrack1.has_mcParticle() == 1 && fwdtrack2.has_mcParticle() == 1) {
                      auto fwdparticle1 = fwdtrack1.mcParticle();
                      auto fwdparticle2 = fwdtrack2.mcParticle();
                      //check match correction
                      bool correctmatch1 = 0;
                      bool correctmatch2 = 0;

                      int fwd1ID = 0;
                      int fwd2ID = 0;
                      int fwd1pdgcode = 0;
                      int fwd2pdgcode = 0;
                      double fwd1vz = 0;
                      double fwd2vz = 0;
                      int fwd1process = 0;
                      int fwd2process = 0;
                      int mft1recID = 0;
                      int mft1ID = 0;
                      int mft2recID = 0;
                      int mft2ID = 0;
                      int mft1pdgcode = 0;
                      int mft2pdgcode = 0;
                      double mft1vz = 0;
                      double mft2vz = 0;
                      int mft1process = 0;
                      int mft2process = 0;

                      const auto mcmothersidlist1 = fwdparticle1.mothersIds();
                      for (auto& matchMCH : fwdtracks) {
                        if (matchMCH.trackType() == 3) {
                          if (fwdtrack1.matchMCHTrackId() == matchMCH.globalIndex()) {
                            for (auto& mfttrack : mfttracks) {
                              if (fwdtrack1.matchMFTTrackId() == mfttrack.globalIndex()) {
                                if (mfttrack.has_mcParticle() == 1) {
                                  auto mftparticle = mfttrack.mcParticle();
                                  fwd1ID = fwdparticle1.globalIndex();
                                  fwd1pdgcode = fwdparticle1.pdgCode();
                                  fwd1vz = fwdparticle1.vz();
                                  fwd1process = fwdparticle1.getProcess();
                                  mft1ID = mftparticle.globalIndex();
                                  mft1recID = mfttrack.globalIndex();
                                  mft1pdgcode = mftparticle.pdgCode();
                                  mft1vz = mftparticle.vz();
                                  mft1process = mftparticle.getProcess();
                                  if (mftparticle.globalIndex() == fwdparticle1.globalIndex()) correctmatch1 = 1;
                                }
                              }
                            }
                          }
                          if (fwdtrack2.matchMCHTrackId() == matchMCH.globalIndex()) {
                            for (auto& mfttrack : mfttracks) {
                              if (fwdtrack2.matchMFTTrackId() == mfttrack.globalIndex()) {
                                if (mfttrack.has_mcParticle() == 1) {
                                  auto mftparticle = mfttrack.mcParticle();
                                  fwd2ID = fwdparticle2.globalIndex();
                                  fwd2pdgcode = fwdparticle2.pdgCode();
                                  fwd2vz = fwdparticle2.vz();
                                  fwd2process = fwdparticle2.getProcess();
                                  mft2ID = mftparticle.globalIndex();
                                  mft2recID = mfttrack.globalIndex();
                                  mft2pdgcode = mftparticle.pdgCode();
                                  mft2vz = mftparticle.vz();
                                  mft2process = mftparticle.getProcess();
                                  if (mftparticle.globalIndex() == fwdparticle2.globalIndex()) correctmatch2 = 1;
                                }
                              }
                            }
                          }
                        }
                      }
                      //end check match correction
                      if (fwdtrack1.signed1Pt()*fwdtrack2.signed1Pt() < 0) {
                        registry.fill(HIST("hMassUS"), lv.M());
                        registry.fill(HIST("hSVUS"), VertexingSV);
                        if (mft1ID == mft2ID) {
                          registry.fill(HIST("hMassUSsamemft"), lv.M());
                          registry.fill(HIST("hSVUSsamemft"), VertexingSV);
                        } else {
                          registry.fill(HIST("hMassUSdifferentmft"), lv.M());
                          registry.fill(HIST("hSVUSdifferentmft"), VertexingSV);
                        }
                        if (mft1recID == mft2recID) {
                          registry.fill(HIST("hMassUSsamemftrec"), lv.M());
                          registry.fill(HIST("hSVUSsamemftrec"), VertexingSV);
                        } else {
                          registry.fill(HIST("hMassUSdifferentmftrec"), lv.M());
                          registry.fill(HIST("hSVUSdifferentmftrec"), VertexingSV);
                        }
                        if (correctmatch1 == 1 && correctmatch2 ==1) {
                          registry.fill(HIST("hMassUScorrectmatch"), lv.M());
                          registry.fill(HIST("hSVUScorrectmatch"), VertexingSV);
                        } else if (correctmatch1 == 0 && correctmatch2 ==0) {
                          registry.fill(HIST("hMassUSwrongmatch"), lv.M());
                          registry.fill(HIST("hSVUSwrongmatch"), VertexingSV);
                          /*
                          if (VertexingSV < -40 && VertexingSV != -999) {
                            LOGF(info, "--- Unlike-sign, Wrongmatch, VertexingSV < 40 ---");
                            LOGF(info, "VertexingSV = %f", VertexingSV);
                            LOGF(info, "fwd1ID = %d", fwd1ID);
                            LOGF(info, "fwd1pdgcode = %d", fwd1pdgcode);
                            LOGF(info, "fwd1vz = %f", fwd1vz);
                            LOGF(info, "fwd1process = %d", fwd1process);
                            LOGF(info, "mft1ID = %d", mft1ID);
                            LOGF(info, "mft1recID = %d", mft1recID);
                            LOGF(info, "mft1pdgcode = %d", mft1pdgcode);
                            LOGF(info, "mft1vz = %f", mft1vz);
                            LOGF(info, "mft1process = %d", mft1process);
                            LOGF(info, "fwd2ID = %d", fwd2ID);
                            LOGF(info, "fwd2pdgcode = %d", fwd2pdgcode);
                            LOGF(info, "fwd2vz = %f", fwd2vz);
                            LOGF(info, "fwd2process = %d", fwd2process);
                            LOGF(info, "mft2ID = %d", mft2ID);
                            LOGF(info, "mft2recID = %d", mft2recID);
                            LOGF(info, "mft2pdgcode = %d", mft2pdgcode);
                            LOGF(info, "mft2vz = %f", mft2vz);
                            LOGF(info, "mft2process = %d", mft2process);
                          }
                          */
                        } else {
                          registry.fill(HIST("hMassUSonesidematch"), lv.M());
                          registry.fill(HIST("hSVUSonesidematch"), VertexingSV);
                          /*
                          if (VertexingSV < -40 && VertexingSV != -999) {
                            LOGF(info, "--- Unlike-sign, Onematch, VertexingSV < 40 ---");
                            LOGF(info, "VertexingSV = %f", VertexingSV);
                            LOGF(info, "fwd1ID = %d", fwd1ID);
                            LOGF(info, "fwd1pdgcode = %d", fwd1pdgcode);
                            LOGF(info, "fwd1vz = %f", fwd1vz);
                            LOGF(info, "fwd1process = %d", fwd1process);
                            LOGF(info, "mft1ID = %d", mft1ID);
                            LOGF(info, "mft1pdgcode = %d", mft1pdgcode);
                            LOGF(info, "mft1vz = %f", mft1vz);
                            LOGF(info, "mft1process = %d", mft1process);
                            LOGF(info, "fwd2ID = %d", fwd2ID);
                            LOGF(info, "fwd2pdgcode = %d", fwd2pdgcode);
                            LOGF(info, "fwd2vz = %f", fwd2vz);
                            LOGF(info, "fwd2process = %d", fwd2process);
                            LOGF(info, "mft2ID = %d", mft2ID);
                            LOGF(info, "mft2pdgcode = %d", mft2pdgcode);
                            LOGF(info, "mft2vz = %f", mft2vz);
                            LOGF(info, "mft2process = %d", mft2process);
                          }
                          */
                        }
                      } else if ((fwdtrack1.signed1Pt() >= 0 && fwdtrack2.signed1Pt() >= 0) || (fwdtrack1.signed1Pt() < 0 && fwdtrack2.signed1Pt() < 0)) {
                        registry.fill(HIST("hMassLS"), lv.M());
                        registry.fill(HIST("hSVLS"), VertexingSV);
                        if (mft1ID == mft2ID) {
                          registry.fill(HIST("hMassLSsamemft"), lv.M());
                          registry.fill(HIST("hSVLSsamemft"), VertexingSV);
                        } else {
                          registry.fill(HIST("hMassLSdifferentmft"), lv.M());
                          registry.fill(HIST("hSVLSdifferentmft"), VertexingSV);
                        }
                        if (mft1recID == mft2recID) {
                          registry.fill(HIST("hMassLSsamemftrec"), lv.M());
                          registry.fill(HIST("hSVLSsamemftrec"), VertexingSV);
                        } else {
                          registry.fill(HIST("hMassLSdifferentmftrec"), lv.M());
                          registry.fill(HIST("hSVLSdifferentmftrec"), VertexingSV);
                        }
                        if (correctmatch1 == 1 && correctmatch2 ==1) {
                          registry.fill(HIST("hMassLScorrectmatch"), lv.M());
                          registry.fill(HIST("hSVLScorrectmatch"), VertexingSV);
                        } else if (correctmatch1 == 0 && correctmatch2 ==0) {
                          registry.fill(HIST("hMassLSwrongmatch"), lv.M());
                          registry.fill(HIST("hSVLSwrongmatch"), VertexingSV);
                          /*
                          if (VertexingSV < -40 && VertexingSV != -999) {
                            LOGF(info, "--- Like-sign, Wrongmatch, VertexingSV < 40 ---");
                            LOGF(info, "VertexingSV = %f", VertexingSV);
                            LOGF(info, "fwd1ID = %d", fwd1ID);
                            LOGF(info, "fwd1pdgcode = %d", fwd1pdgcode);
                            LOGF(info, "fwd1vz = %f", fwd1vz);
                            LOGF(info, "fwd1process = %d", fwd1process);
                            LOGF(info, "mft1ID = %d", mft1ID);
                            LOGF(info, "mft1pdgcode = %d", mft1pdgcode);
                            LOGF(info, "mft1vz = %f", mft1vz);
                            LOGF(info, "mft1process = %d", mft1process);
                            LOGF(info, "fwd2ID = %d", fwd2ID);
                            LOGF(info, "fwd2pdgcode = %d", fwd2pdgcode);
                            LOGF(info, "fwd2vz = %f", fwd2vz);
                            LOGF(info, "fwd2process = %d", fwd2process);
                            LOGF(info, "mft2ID = %d", mft2ID);
                            LOGF(info, "mft2pdgcode = %d", mft2pdgcode);
                            LOGF(info, "mft2vz = %f", mft2vz);
                            LOGF(info, "mft2process = %d", mft2process);
                          }
                          */
                        } else {
                          registry.fill(HIST("hMassLSonesidematch"), lv.M());
                          registry.fill(HIST("hSVLSonesidematch"), VertexingSV);
                          /*
                          if (VertexingSV < -40 && VertexingSV != -999) {
                            LOGF(info, "--- Like-sign, Onematch, VertexingSV < 40 ---");
                            LOGF(info, "VertexingSV = %f", VertexingSV);
                            LOGF(info, "fwd1ID = %d", fwd1ID);
                            LOGF(info, "fwd1pdgcode = %d", fwd1pdgcode);
                            LOGF(info, "fwd1vz = %f", fwd1vz);
                            LOGF(info, "fwd1process = %d", fwd1process);
                            LOGF(info, "mft1ID = %d", mft1ID);
                            LOGF(info, "mft1pdgcode = %d", mft1pdgcode);
                            LOGF(info, "mft1vz = %f", mft1vz);
                            LOGF(info, "mft1process = %d", mft1process);
                            LOGF(info, "fwd2ID = %d", fwd2ID);
                            LOGF(info, "fwd2pdgcode = %d", fwd2pdgcode);
                            LOGF(info, "fwd2vz = %f", fwd2vz);
                            LOGF(info, "fwd2process = %d", fwd2process);
                            LOGF(info, "mft2ID = %d", mft2ID);
                            LOGF(info, "mft2pdgcode = %d", mft2pdgcode);
                            LOGF(info, "mft2vz = %f", mft2vz);
                            LOGF(info, "mft2process = %d", mft2process);
                          }
                          */
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
    adaptAnalysisTask<checklikesign>(cfgc)
  };
}
