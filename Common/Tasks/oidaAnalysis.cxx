#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <utility>

#include <TObject.h>
#include <TString.h>
#include "TRandom.h"
#include "TH3F.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/Boost.h"

#include "Framework/DataTypes.h"
//#include "MCHTracking/TrackExtrap.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "DCAFitter/DCAFitterN.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Math/SMatrix.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "DCAFitter/FwdDCAFitterN.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "CommonConstants/PhysicsConstants.h"

/*
#include "/Users/masahirolaptop/alice/sw/osx_arm64/KFParticle/v1.1-5-local1/include/KFParticle.h"
#include "/Users/masahirolaptop/alice/sw/osx_arm64/KFParticle/v1.1-5-local1/include/KFPTrack.h"
#include "/Users/masahirolaptop/alice/sw/osx_arm64/KFParticle/v1.1-5-local1/include/KFPVertex.h"
#include "/Users/masahirolaptop/alice/sw/osx_arm64/KFParticle/v1.1-5-local1/include/KFParticleBase.h"
#include "/Users/masahirolaptop/alice/sw/osx_arm64/KFParticle/v1.1-5-local1/include/KFVertex.h"
*/


// my Analysis include
#include "Framework/runDataProcessing.h" 
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using Vec3D = ROOT::Math::SVector<double, 3>;

struct myAnalysis{
// Histogram registry: an object to hold your histograms 
HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<double> totalbinCollid{"totalbinCollid",100000,"bin number of Collisionid histo"};
  Configurable<int> binwidthCollid{"binwidthCollid",100000,"bin width of Collisionid histo"};
  Configurable<double> totalbinNClusters{"totalbinNClusters",20,"total bin Number of clusters"};
  Configurable<int> binwidthNClusters{"binwidthNClusters",20,"bin width of number of Clasters"};
  Configurable<double> totalbinx{"totalbinx",0.1,"total bin x"};
  Configurable<int> binwidthx{"binwidthx",1000,"bin width x"};
  Configurable<double> totalbiny{"totalbiny",0.1,"total bin y"};
  Configurable<int> binwidthy{"binwidthy",1000,"bin width y"};
  Configurable<double> totalbinz{"totalbinz",1000,"total bin z"};
  Configurable<int> binwidthz{"binwidthz",100000,"bin width z"};

  void init(InitContext const&){
    // define axes you want to use
    const AxisSpec axisEta{300, -5, +5, "#eta"};
    const AxisSpec axisCollisionID{binwidthCollid,0,totalbinCollid,"CollisionID"};
    const AxisSpec axisNClusters{binwidthNClusters,0,totalbinNClusters,"NClusters"};
    const AxisSpec axisx{binwidthx,-totalbinx,totalbinx,"x"};
    const AxisSpec axisy{binwidthy,-totalbiny,totalbiny,"y"};
    const AxisSpec axisz{binwidthz,-totalbinz,totalbinz,"z"};
    const AxisSpec axisTrackType{11,-1,10,"trackType"};
    const AxisSpec axisPt{600,0,30};
    const AxisSpec axisPosX{2000,-0.1,0.1};
    const AxisSpec axisPosY{2000,-0.1,0.1};
    const AxisSpec axisPosZ{400,-20,20};
    const AxisSpec axisMFTTracksNumber{251,-1,250};
    const AxisSpec axisFwdTracksNumber{6,-1,5};
    const AxisSpec axisIndex{1000000,0,1000000};
    const AxisSpec axisBCId{binwidthCollid,0,totalbinCollid};
    const AxisSpec axisglobalBC{binwidthCollid,0,totalbinCollid};
    const AxisSpec axisglobalIndex{binwidthCollid,0,totalbinCollid};

    // create histograms for FwdTracks
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta}); 
    histos.add("CollisionID","CollisionID",kTH1F,{axisCollisionID});
    histos.add("NClusters","NClusters",kTH1F,{axisNClusters});
    histos.add("x","x",kTH1F,{axisx});
    histos.add("y","y",kTH1F,{axisy});
    histos.add("z","z",kTH1F,{axisz});
    histos.add("xy","xy",kTH2F,{axisx,axisy});
    histos.add("trackType","trackType",kTH1F,{axisTrackType});
    histos.add("pt","pt",kTH1F,{axisPt});
    histos.add("pt_2","pt_2",kTH1F,{axisPt});
    histos.add("fwdglobalIndex","fwdglobalIndex",kTH1F,{axisIndex});

    // create histograms for MFTTracks
    histos.add("CollisionIDforMFT","CollisionIDforMFT",kTH1F,{axisCollisionID});
    histos.add("NClustersforMFT","NClustersforMFT",kTH1F,{axisNClusters});
    histos.add("xforMFT","xforMFT",kTH1F,{axisx});
    histos.add("yforMFT","yforMFT",kTH1F,{axisy});
    histos.add("zforMFT","zforMFT",kTH1F,{axisz});
    histos.add("xyforMFT","xyforMFT",kTH2F,{axisx,axisy});
    histos.add("ptforMFT","ptforMFT",kTH1F,{axisPt});
    histos.add("etaHistogramforMFT", "etaHistogramforMFT", kTH1F, {axisEta}); 
    histos.add("mftglobalIndex","mftglobalIndex",kTH1F,{axisIndex});
    histos.add("ptforMFT_2","ptforMFT_2",kTH1F,{axisPt});

    // create histogram for collisions
    histos.add("PosX","PosX",kTH1F,{axisPosX});
    histos.add("PosY","PosY",kTH1F,{axisPosY});
    histos.add("PosZ","PosZ",kTH1F,{axisPosZ});
    histos.add("PosXY","PosXY",kTH2F,{axisPosX,axisPosY});
    histos.add("PosZX","PosZX",kTH2F,{axisPosZ,axisPosX});
    histos.add("PosZY","PosZY",kTH2F,{axisPosZ,axisPosY});
    histos.add("BCId","BCId",kTH1F,{axisBCId});
    histos.add("globalIndex","globalIndex",kTH1F,{axisglobalIndex});
    // create histogram for BCs
    histos.add("globalBC","globalBC",kTH1F,{axisglobalBC});

    // create histogram for Track Number
    histos.add("MFTTracksNumber","MFTTracksNumber",kTH1F,{axisMFTTracksNumber});
    histos.add("FwdTracksNumber","FwdTracksNumber",kTH1F,{axisFwdTracksNumber});
  }
     
  //void process(aod::FwdTracks const& fwdtracks,aod::MFTTracks const& mfttracks,aod::Collisions const& collisions,aod::BCs const& bcs,soa::Filtered<MyMuonsWithCov> const& tracksMuon)
  void process(aod::FwdTracks const& fwdtracks,aod::MFTTracks const& mfttracks,aod::Collisions const& collisions,aod::BCs const& bcs,MyMuonsWithCov const& tracksMuon)
  {
   //static o2::globaltracking::MatchGlobalFwd mMatching;
    for (auto& fwdtrack : fwdtracks) {
      if (fwdtrack.trackType() == 3 ){
      histos.fill(HIST("etaHistogram"),fwdtrack.eta());
      histos.fill(HIST("CollisionID"),fwdtrack.collisionId());
      histos.fill(HIST("NClusters"),fwdtrack.nClusters());
      histos.fill(HIST("x"),fwdtrack.x());
      histos.fill(HIST("y"),fwdtrack.y());
      histos.fill(HIST("z"),fwdtrack.z());
      histos.fill(HIST("xy"),fwdtrack.x(),fwdtrack.y());
      histos.fill(HIST("trackType"),fwdtrack.trackType());
      histos.fill(HIST("fwdglobalIndex"),fwdtrack.globalIndex()); 
      histos.fill(HIST("pt"),fwdtrack.pt());
      }
    }

    for (auto& mfttrack : mfttracks){
      histos.fill(HIST("CollisionIDforMFT"),mfttrack.collisionId());
      histos.fill(HIST("NClustersforMFT"),mfttrack.nClusters());
      histos.fill(HIST("xforMFT"),mfttrack.x());
      histos.fill(HIST("yforMFT"),mfttrack.y());
      histos.fill(HIST("zforMFT"),mfttrack.z());
      histos.fill(HIST("xyforMFT"),mfttrack.x(),mfttrack.y());
      histos.fill(HIST("ptforMFT"),mfttrack.pt());
      histos.fill(HIST("etaHistogramforMFT"),mfttrack.eta());
      histos.fill(HIST("mftglobalIndex"),mfttrack.globalIndex());
    }

    for (auto& bc : bcs){
      histos.fill(HIST("globalBC"),bc.globalBC());
    }

    int allFwdTracknumber;
    int allMFTTracknumber;
    for (auto& collision : collisions){
      histos.fill(HIST("BCId"),collision.bcId());
      histos.fill(HIST("PosX"),collision.posX());
      histos.fill(HIST("PosY"),collision.posY());
      histos.fill(HIST("PosZ"),collision.posZ());
      histos.fill(HIST("PosXY"),collision.posX(),collision.posY());
      histos.fill(HIST("PosZX"),collision.posZ(),collision.posX());
      histos.fill(HIST("PosZY"),collision.posZ(),collision.posY());
      histos.fill(HIST("globalIndex"),collision.globalIndex());

      int NumberofFwdTracks = 0;

      for (auto& trackMuon : tracksMuon){
        if (trackMuon.trackType() == 3 ){
          if(trackMuon.collisionId() == collision.globalIndex()){
            NumberofFwdTracks = NumberofFwdTracks + 1;
            allFwdTracknumber = allFwdTracknumber + 1;
            histos.fill(HIST("pt_2"),trackMuon.pt()); // test entry fwd tracks 

              double chi2 = trackMuon.chi2();
              SMatrix5 tpars(trackMuon.x(), trackMuon.y(), trackMuon.phi(), trackMuon.tgl(), trackMuon.signed1Pt());
              std::vector<double> v1{trackMuon.cXX(), trackMuon.cXY(), trackMuon.cYY(), trackMuon.cPhiX(), trackMuon.cPhiY(),
                           trackMuon.cPhiPhi(), trackMuon.cTglX(), trackMuon.cTglY(), trackMuon.cTglPhi(), trackMuon.cTglTgl(),
                         trackMuon.c1PtX(), trackMuon.c1PtY(), trackMuon.c1PtPhi(), trackMuon.c1PtTgl(), trackMuon.c1Pt21Pt2()};
              SMatrix55 tcovs(v1.begin(), v1.end());
              o2::track::TrackParCovFwd fwdtrack{trackMuon.z(), tpars, tcovs, chi2};
              o2::mch::TrackExtrap::setField();
              o2::dataformats::GlobalFwdTrack track;
              track.setParameters(tpars);
              track.setZ(fwdtrack.getZ());
              track.setCovariances(tcovs);
              //auto mchtrack = mMatching.FwdtoMCH(track);
              //o2::mch::TrackExtrap::extrapToVertex(mchtrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

          }
        }
      }
      histos.fill(HIST("FwdTracksNumber"),NumberofFwdTracks);

      int NumberofMFTTracks = 0;
      for (auto& mfttrack : mfttracks){
        if(mfttrack.collisionId() == collision.globalIndex()){
          NumberofMFTTracks = NumberofMFTTracks + 1;
          allMFTTracknumber = allMFTTracknumber + 1;
          histos.fill(HIST("ptforMFT_2"),mfttrack.pt()); // test entry MFT tracks
        }
      }
      histos.fill(HIST("MFTTracksNumber"),NumberofMFTTracks);
    }
    //std::cout << "MFTTracknumber" << allMFTTracknumber << std::endl;
    //std::cout << "FwdTracknumber" << allFwdTracknumber << std::endl;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myAnalysis>(cfgc)};
}


