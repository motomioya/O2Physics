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
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/FwdDCAFitterN.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/mftmchMatchingML.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include <TLorentzVector.h>
#include <CCDB/BasicCCDBManager.h>
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "TDatabasePDG.h"
#include "Math/Vector3D.h"
#include <vector>
#include <cmath>
#include <utility>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using MyMuonsColl = soa::Join<aod::FwdTracks, aod::FwdTrkCompColls>;
using MyMFTsColl = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;
using Vec3D = ROOT::Math::SVector<double, 3>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps>;

struct hfcheck {
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

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  int muonPDGCode = 13;
  TParticlePDG* muonParticle = TDatabasePDG::Instance()->GetParticle(muonPDGCode);
  double muonMass = muonParticle->Mass();
  int kaonPDGCode = 311;
  TParticlePDG* kaonParticle = TDatabasePDG::Instance()->GetParticle(kaonPDGCode);
  double kaonMass = kaonParticle->Mass();

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
                                                        //
  o2::parameters::GRPMagField* grpmag = nullptr;
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  HistogramRegistry registry{
    "registry",
    {
      {"hMuInclusivePt", "hMuInclusivePt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuInclusiveP", "hMuInclusiveP", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuFromHfPt", "hMuFromHfPt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuFromHfP", "hMuFromHfP", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuFromHadronPt", "hMuFromHadronPt", {HistType::kTH1F, {{2000, 0, 20}}}},
      {"hMuFromHadronP", "hMuFromHadronP", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuFromHadronZpos", "hMuFromHadronZpos", {HistType::kTH1F, {{1000,-1000,0}}}},
      {"hConter", "hConter", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
      {"hInvariantMassHF", "hInvariantMassHF", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hInvariantMassLF", "hInvariantMassLF", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hMuMufromHfLxyz", "hMuMfromHfLxyz", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuMufromHfPCA", "hMuMfromHfPCA", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuMufromLfLxyz", "hMuMfromLfLxyz", {HistType::kTH1F, {{3000, 0, 30}}}},
      {"hMuMufromLfPCA", "hMuMfromLfPCA", {HistType::kTH1F, {{3000, 0, 30}}}},
    }
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    o2::base::Propagator::initFieldFromGRP(grpmag);

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    Bz = field->getBz(centerMFT);
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz;
  }

std::vector<double> PCA_Cal(std::vector<double> &vec_mftdca)
{
  std::vector<double> vec_pca;
  double mumftX = vec_mftdca[0];
  double mumftY = vec_mftdca[1];
  double mumftZ = vec_mftdca[2];
  double mudcaX = vec_mftdca[3];
  double mudcaY = vec_mftdca[4];
  double canmftX = vec_mftdca[5];
  double canmftY = vec_mftdca[6];
  double canmftZ = vec_mftdca[7];
  double candcaX = vec_mftdca[8];
  double candcaY = vec_mftdca[9];
  double Col_x = vec_mftdca[10];
  double Col_y = vec_mftdca[11];
  double Col_z = vec_mftdca[12];

  auto unit_Na = sqrt(pow(mumftX-mudcaX,2)+pow(mumftY-mudcaY,2)+pow(mumftZ-Col_z,2));
  auto unit_Nc = sqrt(pow(canmftX-candcaX,2)+pow(canmftY-candcaY,2)+pow(canmftZ-Col_z,2));
  auto Nax = (mumftX-mudcaX)/unit_Na;
  auto Nay = (mumftY-mudcaY)/unit_Na;
  auto Naz = (mumftZ-Col_z)/unit_Na;
  auto Ncx = (canmftX-candcaX)/unit_Nc;
  auto Ncy = (canmftY-candcaY)/unit_Nc;
  auto Ncz = (canmftZ-Col_z)/unit_Nc;
  auto A1 = Nax*Nax + Nay*Nay + Naz*Naz;
  auto A2 = -(Nax*Ncx + Nay*Ncy + Naz*Ncz);
  auto A3 = (mudcaX-candcaX)*Nax + (mudcaY-candcaY)*Nay + (Col_z-Col_z)*Naz;
  auto B1 = A2;
  auto B2 = Ncx*Ncx + Ncy*Ncy + Ncz*Ncz;
  auto B3 = (candcaX-mudcaX)*Ncx + (candcaY-mudcaY)*Ncy + (Col_z-Col_z)*Ncz;
  auto t = (A1*B3-A3*B1)/(A2*B1-A1*B2);
  auto s = -((A2*t+A3)/A1);

  double predict_mux = mudcaX + s*Nax;
  double predict_muy = mudcaY + s*Nay;
  double predict_muz = Col_z + s*Naz;
  double predict_canx = candcaX + t*Ncx;
  double predict_cany = candcaY + t*Ncy;
  double predict_canz = Col_z + t*Ncz;
  double r_xyz = sqrt(pow(predict_canx-predict_mux,2) + pow(predict_cany-predict_muy,2) + pow(predict_canz-predict_muz,2));
  auto vecx_mu = mumftX - mudcaX;
  auto vecy_mu = mumftY - mudcaY;
  auto vecz_mu = mumftZ - Col_z;
  auto vecx_can = canmftX - candcaX;
  auto vecy_can = canmftY - candcaY;
  auto vecz_can = canmftZ - Col_z;
  auto cosxy = (vecx_mu*vecx_can + vecy_mu*vecy_can + vecz_mu*vecz_can)/((sqrt(pow(vecx_mu,2)+pow(vecy_mu,2)+pow(vecz_mu,2))) * (sqrt(pow(vecx_can,2)+pow(vecy_can,2)+pow(vecz_can,2))));

  auto pcaX = ((predict_mux+predict_canx)/2)-Col_x;
  auto pcaY = ((predict_muy+predict_cany)/2)-Col_y;
  auto pcaZ = ((predict_muz+predict_canz)/2)-Col_z;
  auto pcaD = sqrt(pow(pcaX,2)+pow(pcaY,2)+pow(pcaZ,2));
  vec_pca = {r_xyz, cosxy, pcaX, pcaY, pcaZ, pcaD};
  return vec_pca;
}

  void process(aod::Collisions::iterator const& collision, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, aod::MFTTracks const& mfttracks, aod::McParticles const& mcparticles)
  {
    for (auto& fwdtrack : fwdtracks) {

      if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        if (fwdtrack.has_mcParticle()){
          registry.fill(HIST("hMuInclusivePt"), fwdtrack.pt());
          registry.fill(HIST("hMuInclusiveP"), fwdtrack.p());
          auto fwdparticle = fwdtrack.mcParticle();
          const auto mcfwdtrackmothers = fwdtrack.mcParticle().mothersIds();
          registry.fill(HIST("hConter"), 0);
          if (mcfwdtrackmothers.size() > 0) {
            registry.fill(HIST("hConter"), 1);
            int mcfwdfirstmotherid = mcfwdtrackmothers[0];
            auto mcfwdfirstmother = mcparticles.iteratorAt(mcfwdfirstmotherid);
            //int mcfwdfirstmotherpdg = mcfwdfirstmother.pdgCode();
            int mcfwdlastmotherid = mcfwdtrackmothers[mcfwdtrackmothers.size() - 1];
            auto mcfwdlastmother = mcparticles.iteratorAt(mcfwdlastmotherid);
            int mcfwdlastmotherpdg = mcfwdlastmother.pdgCode();
            if (mcfwdfirstmother.fromBackgroundEvent() == false) {
              registry.fill(HIST("hConter"), 2);
              if ( (mcfwdlastmotherpdg >= 400 && mcfwdlastmotherpdg <= 439) || (mcfwdlastmotherpdg <= -400 && mcfwdlastmotherpdg >= -439) ) {
                registry.fill(HIST("hMuFromHfPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHfP"), fwdtrack.p());
                registry.fill(HIST("hConter"), 3);
                for (auto& mfttrack : mfttracks) {

                  TLorentzVector lv1, lv2, lv;
                  lv1.SetPtEtaPhiM(fwdtrack.pt(), fwdtrack.eta(), fwdtrack.phi(), muonMass);
                  lv2.SetPtEtaPhiM(mfttrack.pt(), mfttrack.eta(), mfttrack.phi(), kaonMass);
                  lv = lv1 + lv2;
                  registry.fill(HIST("hInvariantMassHF"), lv.M());

                  //propagate to 
                  double muonchi2 = fwdtrack.chi2();
                  SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
                  std::vector<double> muonv1;
                  SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
                  o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
                  o2::track::TrackParCovFwd muonpars2{fwdtrack.z(), muonpars, muoncovs, muonchi2};
                  muonpars1.propagateToZ(-77.5, Bz);
                  muonpars2.propagateToZ(collision.posZ(), Bz);

                  double mftchi2 = mfttrack.chi2();
                  SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                  std::vector<double> mftv1;
                  SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                  o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                  o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars, mftcovs, mftchi2};
                  mftpars1.propagateToZ(-77.5, Bz);
                  mftpars2.propagateToZ(collision.posZ(), Bz);

                  std::vector<double> pca_cal_input = {muonpars1.getX(),muonpars1.getY(),muonpars1.getZ(),muonpars2.getX(), muonpars2.getY(),mftpars1.getX(), mftpars1.getY(),mftpars1.getZ(),mftpars2.getX(), mftpars2.getY(),collision.posX(), collision.posY(), collision.posZ()};
                  std::vector<double> pca_cal_result = PCA_Cal(pca_cal_input);
                  registry.fill(HIST("hMuMufromHfLxyz"), pca_cal_result[5]);
                  registry.fill(HIST("hMuMufromHfPCA"), pca_cal_result[0]);
                }
              }
            } else {
              registry.fill(HIST("hConter"), 4);
              if ( (mcfwdlastmotherpdg >= 100 && mcfwdlastmotherpdg <= 119) || (mcfwdlastmotherpdg <= -100 && mcfwdlastmotherpdg >= -119) || (mcfwdlastmotherpdg >= 1000 && mcfwdlastmotherpdg <= 1999) || (mcfwdlastmotherpdg <= -1000 && mcfwdlastmotherpdg >= -1999) || (mcfwdlastmotherpdg >= 200 && mcfwdlastmotherpdg <= 299) || (mcfwdlastmotherpdg <= -200 && mcfwdlastmotherpdg >= -299) || (mcfwdlastmotherpdg >= 2000 && mcfwdlastmotherpdg <= 2999) || (mcfwdlastmotherpdg <= -2000 && mcfwdlastmotherpdg >= -2999) || (mcfwdlastmotherpdg >= 300 && mcfwdlastmotherpdg <= 399) || (mcfwdlastmotherpdg <= -300 && mcfwdlastmotherpdg >= -399) || (mcfwdlastmotherpdg >= 3000 && mcfwdlastmotherpdg <= 3999) || (mcfwdlastmotherpdg <= -3000 && mcfwdlastmotherpdg >= -3999) ) {
                registry.fill(HIST("hMuFromHadronPt"), fwdtrack.pt());
                registry.fill(HIST("hMuFromHadronP"), fwdtrack.p());
                registry.fill(HIST("hMuFromHadronZpos"), fwdparticle.vz());
                registry.fill(HIST("hConter"), 5);
                for (auto& mfttrack : mfttracks) {
                  TLorentzVector lv1, lv2, lv;
                  lv1.SetPtEtaPhiM(fwdtrack.pt(), fwdtrack.eta(), fwdtrack.phi(), muonMass);
                  lv2.SetPtEtaPhiM(mfttrack.pt(), mfttrack.eta(), mfttrack.phi(), kaonMass);
                  lv = lv1 + lv2;
                  registry.fill(HIST("hInvariantMassLF"), lv.M());

                  //propagate to 
                  double muonchi2 = fwdtrack.chi2();
                  SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
                  std::vector<double> muonv1;
                  SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
                  o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
                  o2::track::TrackParCovFwd muonpars2{fwdtrack.z(), muonpars, muoncovs, muonchi2};
                  muonpars1.propagateToZ(-77.5, Bz);
                  muonpars2.propagateToZ(collision.posZ(), Bz);

                  double mftchi2 = mfttrack.chi2();
                  SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                  std::vector<double> mftv1;
                  SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                  o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                  o2::track::TrackParCovFwd mftpars2{mfttrack.z(), mftpars, mftcovs, mftchi2};
                  mftpars1.propagateToZ(-77.5, Bz);
                  mftpars2.propagateToZ(collision.posZ(), Bz);

                  std::vector<double> pca_cal_input = {muonpars1.getX(),muonpars1.getY(),muonpars1.getZ(),muonpars2.getX(), muonpars2.getY(),mftpars1.getX(), mftpars1.getY(),mftpars1.getZ(),mftpars2.getX(), mftpars2.getY(),collision.posX(), collision.posY(), collision.posZ()};
                  std::vector<double> pca_cal_result = PCA_Cal(pca_cal_input);
                  registry.fill(HIST("hMuMufromLfLxyz"), pca_cal_result[5]);
                  registry.fill(HIST("hMuMufromLfPCA"), pca_cal_result[0]);
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
    adaptAnalysisTask<hfcheck>(cfgc)
  };
}
