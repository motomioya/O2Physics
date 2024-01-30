// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <math.h>
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include <string>
#include <regex>
#include <TLorentzVector.h>
#include "Common/DataModel/MftmchMatchingML.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "CCDB/CcdbApi.h"
#include "Tools/ML/model.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "DataFormatsParameters/GRPMagField.h"
#include <CCDB/BasicCCDBManager.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::ml;
using o2::globaltracking::MatchingFunc_t;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps>;

struct checkMlEff {
  Produces<aod::FwdTracksML> fwdtrackml;

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

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  HistogramRegistry registry{
    "registry",
    {
      {"MuonPt", "ParticlePt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MuonHasColPt", "ParticleHasColPt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MuonHasGlobalPt", "ParticleHasGlobalPt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"Chi2Pt", "Chi2Pt", {HistType::kTH2F, {{200, 0, 20},{200,0,200}}}},
      {"Chi2TruePt", "Chi2TruePt", {HistType::kTH2F, {{200, 0, 20},{200,0,200}}}},
      {"MlPt", "MlPt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MlOnePt", "MlOnePt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MlTruePt", "MlOnePt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MlOneTruePt", "MlOnePt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MuonX", "ParticleX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MuonHasColX", "ParticleX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MuonHasGlobalX", "ParticleHasGlobalX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"Chi2X", "Chi2X", {HistType::kTH2F, {{500, 0, 50},{200,0,200}}}},
      {"Chi2TrueX", "Chi2TrueX", {HistType::kTH2F, {{500, 0, 50},{200,0,200}}}},
      {"MlX", "MlX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlOneX", "MlOneX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlTrueX", "MlOneX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlOneTrueX", "MlOneX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MuonY", "ParticleY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MuonHasColY", "ParticleHasColY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MuonHasGlobalY", "ParticleHasGlobalY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"Chi2Y", "Chi2Y", {HistType::kTH2F, {{500, 0, 50},{200,0,200}}}},
      {"Chi2TrueY", "Chi2TrueY", {HistType::kTH2F, {{500, 0, 50},{200,0,200}}}},
      {"MlY", "MlY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlOneY", "MlOneY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlTrueY", "MlOneY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlOneTrueY", "MlOneY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MuonEta", "ParticleEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MuonHasColEta", "ParticleHasColEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MuonHasGlobalEta", "ParticleHasGlobalEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"Chi2Eta", "Chi2Eta", {HistType::kTH2F, {{500, -5, 0},{200,0,200}}}},
      {"Chi2TrueEta", "Chi2TrueEta", {HistType::kTH2F, {{500, -5, 0},{200,0,200}}}},
      {"MlEta", "MlEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MlOneEta", "MlOneEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MlTrueEta", "MlOneEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MlOneTrueEta", "MlOneEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MuonPhi", "ParticlePhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MuonHasColPhi", "ParticleHasColPhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MuonHasGlobalPhi", "ParticleHasGlobalPhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"Chi2Phi", "Chi2Phi", {HistType::kTH2F, {{1000, -5, 5},{200,0,200}}}},
      {"Chi2TruePhi", "Chi2TruePhi", {HistType::kTH2F, {{1000, -5, 5},{200,0,200}}}},
      {"MlPhi", "MlPhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MlOnePhi", "MlOnePhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MlTruePhi", "MlOnePhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MlOneTruePhi", "MlOnePhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MldeltaColId", "MldeltaColId", {HistType::kTH1F, {{40, -20, 20}}}},
      {"MlOnedeltaColId", "MlOnedeltaColId", {HistType::kTH1F, {{40, -20, 20}}}},
      {"MlTruedeltaColId", "MlOnedeltaColId", {HistType::kTH1F, {{40, -20, 20}}}},
      {"MlOneTruedeltaColId", "MlOnedeltaColId", {HistType::kTH1F, {{40, -20, 20}}}},
      {"MlScore", "MlScore", {HistType::kTH1F, {{1000, 0, 1}}}},
      {"MlScoreTrue", "MlScoreTrue", {HistType::kTH1F, {{1000, 0, 1}}}},
      {"MlScoreFalse", "MlScoreFalse", {HistType::kTH1F, {{1000, 0, 1}}}},
      {"Input_MFT_X", "Input_MFT_X", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MFT_Y", "Input_MFT_Y", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MFT_Phi", "Input_MFT_Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MFT_Tanl", "Input_MFT_Tanl", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_X", "Input_MCH_X", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MCH_Y", "Input_MCH_Y", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MCH_Phi", "Input_MCH_Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_Tanl", "Input_MCH_Tanl", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_Pt", "Input_MCH_Pt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"Input_Delta_XY", "Input_Delta_XY", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_X", "Input_Delta_X", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_Y", "Input_Delta_Y", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_Phi", "Input_Delta_Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_Delta_Tanl", "Input_Delta_Tanl", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MFT_X_TRUE", "Input_MFT_X_TRUE", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MFT_Y_TRUE", "Input_MFT_Y_TRUE", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MFT_Phi_TRUE", "Input_MFT_Phi_TRUE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MFT_Tanl_TRUE", "Input_MFT_Tanl_TRUE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_X_TRUE", "Input_MCH_X_TRUE", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MCH_Y_TRUE", "Input_MCH_Y_TRUE", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MCH_Phi_TRUE", "Input_MCH_Phi_TRUE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_Tanl_TRUE", "Input_MCH_Tanl_TRUE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_Pt_TRUE", "Input_MCH_Pt_TRUE", {HistType::kTH1F, {{200, 0, 20}}}},
      {"Input_Delta_XY_TRUE", "Input_Delta_XY_TRUE", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_X_TRUE", "Input_Delta_X_TRUE", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_Y_TRUE", "Input_Delta_Y_TRUE", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_Phi_TRUE", "Input_Delta_Phi_TRUE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_Delta_Tanl_TRUE", "Input_Delta_Tanl_TRUE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MFT_X_FALSE", "Input_MFT_X_FALSE", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MFT_Y_FALSE", "Input_MFT_Y_FALSE", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MFT_Phi_FALSE", "Input_MFT_Phi_FALSE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MFT_Tanl_FALSE", "Input_MFT_Tanl_FALSE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_X_FALSE", "Input_MCH_X_FALSE", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MCH_Y_FALSE", "Input_MCH_Y_FALSE", {HistType::kTH1F, {{100, -50, 50}}}},
      {"Input_MCH_Phi_FALSE", "Input_MCH_Phi_FALSE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_Tanl_FALSE", "Input_MCH_Tanl_FALSE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_MCH_Pt_FALSE", "Input_MCH_Pt_FALSE", {HistType::kTH1F, {{200, 0, 20}}}},
      {"Input_Delta_XY_FALSE", "Input_Delta_XY_FALSE", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_X_FALSE", "Input_Delta_X_FALSE", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_Y_FALSE", "Input_Delta_Y_FALSE", {HistType::kTH1F, {{200, -50, 50}}}},
      {"Input_Delta_Phi_FALSE", "Input_Delta_Phi_FALSE", {HistType::kTH1F, {{200, -10, 10}}}},
      {"Input_Delta_Tanl_FALSE", "Input_Delta_Tanl_FALSE", {HistType::kTH1F, {{200, -10, 10}}}}
    }
  };

  Configurable<std::string> cfgModelFile{"ccdb-file", "model.onnx", "name of ONNX model file"};
  Configurable<float> cfgThrScore{"threshold-score", 0.5, "Threshold value for matching score"};
  Configurable<int> cfgColWindow{"collision-window", 5, "Search window (collision ID) for MFT track"};
  Configurable<float> cfgXYWindow{"XY-window", 2.21, "Search window (delta XY) for MFT track"};
  Configurable<float> cfgPhiTanlWindow{"PhiTanl-window", 2.24, "Search window (delta PhiTanl) for MFT track"};
  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "model-explorer"};
  Ort::SessionOptions session_options;
  std::shared_ptr<Ort::Experimental::Session> onnx_session = nullptr;
  OnnxModel model;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT

  template <typename F, typename M>
  std::vector<float> getVariables(F const& muonpars1, M const& mftpars1, bool isTrue)
  {

    Float_t MFT_X = mftpars1.getX();
    Float_t MFT_Y = mftpars1.getY();
    Float_t MFT_Phi = mftpars1.getPhi();
    Float_t MFT_Tanl = mftpars1.getTanl();

    Float_t MCH_X = muonpars1.getX();
    Float_t MCH_Y = muonpars1.getY();
    Float_t MCH_Phi = muonpars1.getPhi();
    Float_t MCH_Tanl = muonpars1.getTanl();
    Float_t MCH_Pt = muonpars1.getPt();

    //Float_t Ratio_X = MFT_X / MCH_X;
    //Float_t Ratio_Y = MFT_Y / MCH_Y;
    //Float_t Ratio_Phi = MFT_Phi / MCH_Phi;
    //Float_t Ratio_Tanl = MFT_Tanl / MCH_Tanl;

    Float_t Delta_X = MFT_X - MCH_X;
    Float_t Delta_Y = MFT_Y - MCH_Y;
    Float_t Delta_Phi = MFT_Phi - MCH_Phi;
    Float_t Delta_Tanl = MFT_Tanl - MCH_Tanl;

    Float_t Delta_XY = sqrt(Delta_X * Delta_X + Delta_Y * Delta_Y);

    registry.fill(HIST("Input_MFT_X"), MFT_X);
    registry.fill(HIST("Input_MFT_Y"), MFT_Y);
    registry.fill(HIST("Input_MFT_Phi"), MFT_Phi);
    registry.fill(HIST("Input_MFT_Tanl"), MFT_Tanl);
    registry.fill(HIST("Input_MCH_X"), MCH_X);
    registry.fill(HIST("Input_MCH_Y"), MCH_Y);
    registry.fill(HIST("Input_MCH_Phi"), MCH_Phi);
    registry.fill(HIST("Input_MCH_Tanl"), MCH_Tanl);
    registry.fill(HIST("Input_MCH_Pt"), MCH_Pt);
    registry.fill(HIST("Input_Delta_XY"), Delta_XY);
    registry.fill(HIST("Input_Delta_X"), Delta_X);
    registry.fill(HIST("Input_Delta_Y"), Delta_Y);
    registry.fill(HIST("Input_Delta_Phi"), Delta_Phi);
    registry.fill(HIST("Input_Delta_Tanl"), Delta_Tanl);

    if (isTrue == 1) {
      registry.fill(HIST("Input_MFT_X_TRUE"), MFT_X);
      registry.fill(HIST("Input_MFT_Y_TRUE"), MFT_Y);
      registry.fill(HIST("Input_MFT_Phi_TRUE"), MFT_Phi);
      registry.fill(HIST("Input_MFT_Tanl_TRUE"), MFT_Tanl);
      registry.fill(HIST("Input_MCH_X_TRUE"), MCH_X);
      registry.fill(HIST("Input_MCH_Y_TRUE"), MCH_Y);
      registry.fill(HIST("Input_MCH_Phi_TRUE"), MCH_Phi);
      registry.fill(HIST("Input_MCH_Tanl_TRUE"), MCH_Tanl);
      registry.fill(HIST("Input_MCH_Pt_TRUE"), MCH_Pt);
      registry.fill(HIST("Input_Delta_XY_TRUE"), Delta_XY);
      registry.fill(HIST("Input_Delta_X_TRUE"), Delta_X);
      registry.fill(HIST("Input_Delta_Y_TRUE"), Delta_Y);
      registry.fill(HIST("Input_Delta_Phi_TRUE"), Delta_Phi);
      registry.fill(HIST("Input_Delta_Tanl_TRUE"), Delta_Tanl);
    } else {
      registry.fill(HIST("Input_MFT_X_FALSE"), MFT_X);
      registry.fill(HIST("Input_MFT_Y_FALSE"), MFT_Y);
      registry.fill(HIST("Input_MFT_Phi_FALSE"), MFT_Phi);
      registry.fill(HIST("Input_MFT_Tanl_FALSE"), MFT_Tanl);
      registry.fill(HIST("Input_MCH_X_FALSE"), MCH_X);
      registry.fill(HIST("Input_MCH_Y_FALSE"), MCH_Y);
      registry.fill(HIST("Input_MCH_Phi_FALSE"), MCH_Phi);
      registry.fill(HIST("Input_MCH_Tanl_FALSE"), MCH_Tanl);
      registry.fill(HIST("Input_MCH_Pt_FALSE"), MCH_Pt);
      registry.fill(HIST("Input_Delta_XY_FALSE"), Delta_XY);
      registry.fill(HIST("Input_Delta_X_FALSE"), Delta_X);
      registry.fill(HIST("Input_Delta_Y_FALSE"), Delta_Y);
      registry.fill(HIST("Input_Delta_Phi_FALSE"), Delta_Phi);
      registry.fill(HIST("Input_Delta_Tanl_FALSE"), Delta_Tanl);
    }
    
    std::vector<float> input_tensor_values{
      MFT_X,
      MFT_Y,
      MFT_Phi,
      MFT_Tanl,
      MCH_X,
      MCH_Y,
      MCH_Phi,
      MCH_Tanl,
      MCH_Pt,
      Delta_XY,
      Delta_X,
      Delta_Y,
      Delta_Phi,
      Delta_Tanl,
      //Ratio_X,
      //Ratio_Y,
      //Ratio_Phi,
      //Ratio_Tanl,
    };
    return input_tensor_values;
  }

  template <typename F, typename M>
  double matchONNX(F const& muonpars1, M const& mftpars1, bool isTrue)
  {
    std::vector<std::string> input_names;
    std::vector<std::vector<int64_t>> input_shapes;
    std::vector<std::string> output_names;
    std::vector<std::vector<int64_t>> output_shapes;

    input_names = onnx_session->GetInputNames();
    input_shapes = onnx_session->GetInputShapes();
    output_names = onnx_session->GetOutputNames();
    output_shapes = onnx_session->GetOutputShapes();

    auto input_shape = input_shapes[0];
    input_shape[0] = 1;

    std::vector<float> input_tensor_values;
    input_tensor_values = getVariables(muonpars1, mftpars1, isTrue);

    std::vector<Ort::Value> input_tensors;
    input_tensors.push_back(Ort::Experimental::Value::CreateTensor<float>(input_tensor_values.data(), input_tensor_values.size(), input_shape));

    std::vector<Ort::Value> output_tensors = onnx_session->Run(input_names, input_tensors, output_names);

    const float* output_value = output_tensors[0].GetTensorData<float>();

    auto score = output_value[0];
    return score;
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    o2::base::Propagator::initFieldFromGRP(grpmag);

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    Bz = field->getBz(centerMFT);
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz;
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    model.initModel(cfgModelFile, false, 1);
    onnx_session = model.getSession();
  }

  void process(aod::Collisions const& collisions, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks,soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels> const& fwdtracksnonfiltered, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&, ExtBCs const&)
  {
    auto bc = collisions.begin().bc_as<ExtBCs>();
    initCCDB(bc);
    static constexpr Double_t MatchingPlaneZ = -77.5;

    for (auto& fwdtrack : fwdtracks) {
      double bestscore = 0;
      int bestmfttrackid = -1;
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        registry.fill(HIST("MuonX"), fwdtrack.x());
        registry.fill(HIST("MuonY"), fwdtrack.y());
        registry.fill(HIST("MuonPhi"), fwdtrack.phi());
        registry.fill(HIST("MuonEta"), fwdtrack.eta());
        registry.fill(HIST("MuonPt"), fwdtrack.pt());
        if (fwdtrack.has_collision() && fwdtrack.has_mcParticle()) {
          auto fwdparticle = fwdtrack.mcParticle();
          registry.fill(HIST("MuonHasColX"), fwdtrack.x());
          registry.fill(HIST("MuonHasColY"), fwdtrack.y());
          registry.fill(HIST("MuonHasColPhi"), fwdtrack.phi());
          registry.fill(HIST("MuonHasColEta"), fwdtrack.eta());
          registry.fill(HIST("MuonHasColPt"), fwdtrack.pt());
          for (auto& fwdtrackglobal : fwdtracksnonfiltered) {
            if (fwdtrackglobal.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
              if (fwdtrackglobal.matchMCHTrackId() == fwdtrack.globalIndex() ){
                registry.fill(HIST("MuonHasGlobalX"), fwdtrack.x());
                registry.fill(HIST("MuonHasGlobalY"), fwdtrack.y());
                registry.fill(HIST("MuonHasGlobalPhi"), fwdtrack.phi());
                registry.fill(HIST("MuonHasGlobalEta"), fwdtrack.eta());
                registry.fill(HIST("MuonHasGlobalPt"), fwdtrack.pt());
                for (auto& mfttrack : mfttracks) {
                  if (mfttrack.has_mcParticle() && mfttrack.has_collision()){
                    auto mftparticle = mfttrack.mcParticle();
                    if (fwdtrackglobal.matchMFTTrackId() == mfttrack.globalIndex() ){
                      registry.fill(HIST("Chi2X"), fwdtrack.x(), fwdtrackglobal.chi2MatchMCHMFT());
                      registry.fill(HIST("Chi2Y"), fwdtrack.y(), fwdtrackglobal.chi2MatchMCHMFT());
                      registry.fill(HIST("Chi2Phi"), fwdtrack.phi(), fwdtrackglobal.chi2MatchMCHMFT());
                      registry.fill(HIST("Chi2Eta"), fwdtrack.eta(), fwdtrackglobal.chi2MatchMCHMFT());
                      registry.fill(HIST("Chi2Pt"), fwdtrack.pt(), fwdtrackglobal.chi2MatchMCHMFT());
                      if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                        registry.fill(HIST("Chi2TrueX"), fwdtrack.x(), fwdtrackglobal.chi2MatchMCHMFT());
                        registry.fill(HIST("Chi2TrueY"), fwdtrack.y(), fwdtrackglobal.chi2MatchMCHMFT());
                        registry.fill(HIST("Chi2TruePhi"), fwdtrack.phi(), fwdtrackglobal.chi2MatchMCHMFT());
                        registry.fill(HIST("Chi2TrueEta"), fwdtrack.eta(), fwdtrackglobal.chi2MatchMCHMFT());
                        registry.fill(HIST("Chi2TruePt"), fwdtrack.pt(), fwdtrackglobal.chi2MatchMCHMFT());
                      }
                    }
                  }
                }
              }
            }
          }
          for (auto& mfttrack : mfttracks) {
            if (mfttrack.has_collision() && mfttrack.has_mcParticle()) {
              auto mftparticle = mfttrack.mcParticle();
              if (0 <= fwdtrack.collisionId() - mfttrack.collisionId() && fwdtrack.collisionId() - mfttrack.collisionId() <= cfgColWindow) {
                // propagate muontrack to matching position
                double muonchi2 = fwdtrack.chi2();
                SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
                std::vector<double> muonv1;
                SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
                o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
                //muonpars1.propagateToZlinear(MatchingPlaneZ);
                muonpars1.propagateToZ(MatchingPlaneZ, Bz);

                // propagate mfttrack to matching position
                double mftchi2 = mfttrack.chi2();
                SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
                std::vector<double> mftv1;
                SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
                o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
                //mftpars1.propagateToZlinear(MatchingPlaneZ);
                mftpars1.propagateToZ(MatchingPlaneZ, Bz);

                Float_t Delta_XY = sqrt(((muonpars1.getX() - mftpars1.getX()) * (muonpars1.getX() - mftpars1.getX())) + ((muonpars1.getY() - mftpars1.getY()) * (muonpars1.getY() - mftpars1.getY())));
                Float_t Delta_PhiTanl = sqrt(((muonpars1.getPhi() - mftpars1.getPhi()) * (muonpars1.getPhi() - mftpars1.getPhi())) + ((muonpars1.getTanl() - mftpars1.getTanl()) * (muonpars1.getTanl() - mftpars1.getTanl())));

                if (Delta_XY <= cfgXYWindow && Delta_PhiTanl <= cfgPhiTanlWindow) {
                  double result = 0;
                  if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                    result = matchONNX(muonpars1, mftpars1, 1);
                  } else {
                    result = matchONNX(muonpars1, mftpars1, 0);
                  }

                  registry.fill(HIST("MlScore"), result);
                  if (result > cfgThrScore) {
                    registry.fill(HIST("MlX"), fwdtrack.x());
                    registry.fill(HIST("MlY"), fwdtrack.y());
                    registry.fill(HIST("MlPhi"), fwdtrack.phi());
                    registry.fill(HIST("MlEta"), fwdtrack.eta());
                    registry.fill(HIST("MlPt"), fwdtrack.pt());
                    registry.fill(HIST("MldeltaColId"), fwdtrack.collisionId() - mfttrack.collisionId());
                    if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                      registry.fill(HIST("MlTrueX"), fwdtrack.x());
                      registry.fill(HIST("MlTrueY"), fwdtrack.y());
                      registry.fill(HIST("MlTruePhi"), fwdtrack.phi());
                      registry.fill(HIST("MlTrueEta"), fwdtrack.eta());
                      registry.fill(HIST("MlTruePt"), fwdtrack.pt());
                      registry.fill(HIST("MlTruedeltaColId"), fwdtrack.collisionId() - mfttrack.collisionId());
                    }
                    if (result > bestscore) {
                      bestscore = result;
                      bestmfttrackid = mfttrack.globalIndex();
                    }
                  }
                }
              }
            }
          }
          if (bestmfttrackid != -1) {
            for (auto& mfttrack : mfttracks) {
              if (mfttrack.globalIndex() == bestmfttrackid) {
                registry.fill(HIST("MlOneX"), fwdtrack.x());
                registry.fill(HIST("MlOneY"), fwdtrack.y());
                registry.fill(HIST("MlOnePhi"), fwdtrack.phi());
                registry.fill(HIST("MlOneEta"), fwdtrack.eta());
                registry.fill(HIST("MlOnePt"), fwdtrack.pt());
                registry.fill(HIST("MlOnedeltaColId"), fwdtrack.collisionId() - mfttrack.collisionId());
                if (fwdtrack.has_mcParticle() && mfttrack.has_mcParticle()){
                  auto fwdparticle = fwdtrack.mcParticle();
                  auto mftparticle = mfttrack.mcParticle();
                  if (fwdparticle.globalIndex() == mftparticle.globalIndex()){
                    registry.fill(HIST("MlOneTrueX"), fwdtrack.x());
                    registry.fill(HIST("MlOneTrueY"), fwdtrack.y());
                    registry.fill(HIST("MlOneTruePhi"), fwdtrack.phi());
                    registry.fill(HIST("MlOneTrueEta"), fwdtrack.eta());
                    registry.fill(HIST("MlOneTruePt"), fwdtrack.pt());
                    registry.fill(HIST("MlOneTruedeltaColId"), fwdtrack.collisionId() - mfttrack.collisionId());
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
    adaptAnalysisTask<checkMlEff>(cfgc)};
}
