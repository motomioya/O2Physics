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
#include "Common/DataModel/EventSelection.h"
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
#include "CCDB/CcdbApi.h"
#include "Tools/ML/model.h"
#include <string>
#include <regex>
#include <math.h>
#include <TLorentzVector.h>
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::ml;
using namespace o2::aod::evsel;
using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using o2::globaltracking::MatchingFunc_t;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct compareMlMatching {
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

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup ));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://ccdb-test.cern.ch:8080", "URL of the CCDB repository"};
  Configurable<std::string> cfgModelDir{"ccdb-path", "Users/m/mooya/models", "base path to the ONNX models"};
  Configurable<std::string> cfgModelName{"ccdb-file", "model_LHC22o.onnx", "name of ONNX model file"};
  Configurable<float> cfgThrScore{"threshold-score", 0.5, "Threshold value for matching score"};
  Configurable<int> cfgColWindow{"collision-window", 6, "Search window (collision ID) for MFT track"};
  Configurable<float> cfgXYWindow{"XY-window", 3, "Search window (delta XY) for MFT track"};

  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "model-explorer"};
  Ort::SessionOptions session_options;
  std::shared_ptr<Ort::Experimental::Session> onnx_session = nullptr;
  OnnxModel model;

  HistogramRegistry registry{
    "registry",
    {
      {"hmlscore", "mlscore", {HistType::kTH1F, {{10000, 0, 1}}}},
      {"hXmuon", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYmuon", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZmuon", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTmuon", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhimuon", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlmuon", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtamuon", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXmuonwithmft", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYmuonwithmft", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZmuonwithmft", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTmuonwithmft", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhimuonwithmft", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlmuonwithmft", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtamuonwithmft", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXchi2", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYchi2", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZchi2", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTchi2", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhichi2", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlchi2", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtachi2", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXchi2cut10", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYchi2cut10", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZchi2cut10", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTchi2cut10", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhichi2cut10", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlchi2cut10", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtachi2cut10", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXchi2cut20", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYchi2cut20", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZchi2cut20", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTchi2cut20", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhichi2cut20", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlchi2cut20", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtachi2cut20", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXchi2cut30", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYchi2cut30", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZchi2cut30", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTchi2cut30", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhichi2cut30", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlchi2cut30", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtachi2cut30", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXchi2cut40", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYchi2cut40", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZchi2cut40", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTchi2cut40", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhichi2cut40", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlchi2cut40", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtachi2cut40", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXchi2cut50", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYchi2cut50", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZchi2cut50", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTchi2cut50", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhichi2cut50", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlchi2cut50", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtachi2cut50", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXml", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYml", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZml", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTml", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhiml", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlml", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtaml", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXmlfwdtrue", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYmlfwdtrue", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZmlfwdtrue", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTmlfwdtrue", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhimlfwdtrue", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlmlfwdtrue", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtamlfwdtrue", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXmlmfttrue", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYmlmfttrue", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZmlmfttrue", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTmlmfttrue", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhimlmfttrue", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlmlmfttrue", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtamlmfttrue", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXmldeltatrue", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYmldeltatrue", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZmldeltatrue", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTmldeltatrue", "PT", {HistType::kTH1F, {{10000, -50, 50}}}},
      {"hPhimldeltatrue", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlmldeltatrue", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtamldeltatrue", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXYmldeltatrue", "ZY", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hXmlfwdfalse", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYmlfwdfalse", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZmlfwdfalse", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTmlfwdfalse", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhimlfwdfalse", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlmlfwdfalse", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtamlfwdfalse", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXmlmftfalse", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYmlmftfalse", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZmlmftfalse", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTmlmftfalse", "PT", {HistType::kTH1F, {{5000, 0, 50}}}},
      {"hPhimlmftfalse", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlmlmftfalse", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtamlmftfalse", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXmldeltafalse", "X", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hYmldeltafalse", "Y", {HistType::kTH1F, {{1000, -50, 50}}}},
      {"hZmldeltafalse", "Z", {HistType::kTH1F, {{2000, -100, 100}}}},
      {"hPTmldeltafalse", "PT", {HistType::kTH1F, {{5000, -50, 50}}}},
      {"hPhimldeltafalse", "Phi", {HistType::kTH1F, {{200, -10, 10}}}},
      {"hTanlmldeltafalse", "Tanl", {HistType::kTH1F, {{600, -30, 30}}}},
      {"hEtamldeltafalse", "Eta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"hXYmldeltafalse", "XY", {HistType::kTH1F, {{200, -10, 10}}}},
    }
  };

  template <typename F, typename M>
  std::vector<float> getVariables(F const& fwdtrack, M const& mfttrack){

    static constexpr Double_t MatchingPlaneZ = -77.5;

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

    Float_t MFT_X      = mftpars1.getX();
    Float_t MFT_Y      = mftpars1.getY();
    Float_t MFT_Phi    = mftpars1.getPhi();
    Float_t MFT_Tanl   = mftpars1.getTanl();

    Float_t MCH_X      = muonpars1.getX();
    Float_t MCH_Y      = muonpars1.getY();
    Float_t MCH_Phi    = muonpars1.getPhi();
    Float_t MCH_Tanl   = muonpars1.getTanl();

    Float_t Ratio_X      = MFT_X      / MCH_X;
    Float_t Ratio_Y      = MFT_Y      / MCH_Y;
    Float_t Ratio_Phi    = MFT_Phi    / MCH_Phi;
    Float_t Ratio_Tanl   = MFT_Tanl   / MCH_Tanl;
          
    Float_t Delta_X      = MFT_X      - MCH_X;
    Float_t Delta_Y      = MFT_Y      - MCH_Y;
    Float_t Delta_Phi    = MFT_Phi    - MCH_Phi;
    Float_t Delta_Tanl   = MFT_Tanl   - MCH_Tanl;

    Float_t Delta_XY     = sqrt(Delta_X*Delta_X + Delta_Y*Delta_Y);
    
    std::vector<float> input_tensor_values{
      MFT_X,
      MFT_Y,
      MFT_Phi,
      MFT_Tanl,
      MCH_X,
      MCH_Y,
      MCH_Phi,
      MCH_Tanl,
      Delta_XY,
      Delta_X,
      Delta_Y,
      Delta_Phi,
      Delta_Tanl,
      Ratio_X,
      Ratio_Y,
      Ratio_Phi,
      Ratio_Tanl,
    };
    return input_tensor_values;  
  }

  template <typename F, typename M>
  double matchONNX(F const& fwdtrack, M const& mfttrack)
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
    input_tensor_values = getVariables(fwdtrack,mfttrack);

    if (input_tensor_values[8] < 3) {
      std::vector<Ort::Value> input_tensors;
      input_tensors.push_back(Ort::Experimental::Value::CreateTensor<float> (input_tensor_values.data(), input_tensor_values.size(), input_shape));  

      std::vector<Ort::Value> output_tensors = onnx_session->Run(input_names, input_tensors, output_names);
       
      const float* output_value = output_tensors[0].GetTensorData<float>();

      auto score = output_value[0];
      
      return score;
    } else {
      auto score = 0;
      return score;
    }
  };

  void init(o2::framework::InitContext&)
  {
    o2::ccdb::CcdbApi ccdbApi;
    std::map<std::string, std::string> metadata;

    ccdbApi.init(cfgCCDBURL);
    //retrieving onnx file from ccdb
    std::string modelFile = cfgModelDir.value;
    bool retrieveSuccess = ccdbApi.retrieveBlob(modelFile, ".", metadata, 1642502592629, false, cfgModelName.value);

    //start session
    if (retrieveSuccess){
      std::map<std::string, std::string> headers = ccdbApi.retrieveHeaders(modelFile, metadata, -1);
      LOG(info) << "Network file downloaded from: " << modelFile << " to: " << "." << "/" << cfgModelName.value;
      model.initModel(cfgModelName, false, 1, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0)); //temporary
      onnx_session = model.getSession();
    }else{
      LOG(info) << "Failed to retrieve Network file";
    }
  }

  void process(aod::Collisions const& collisions, soa::Filtered<aod::FwdTracks> const& fwdtracks, aod::MFTTracks const& mfttracks)
  {
    for (auto const& fwdtrack : fwdtracks){
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack){
        registry.fill(HIST("hXchi2"), fwdtrack.x());
        registry.fill(HIST("hYchi2"), fwdtrack.y());
        registry.fill(HIST("hZchi2"), fwdtrack.z());
        registry.fill(HIST("hPTchi2"), fwdtrack.pt());
        registry.fill(HIST("hPhichi2"), fwdtrack.phi());
        registry.fill(HIST("hTanlchi2"), fwdtrack.tgl());
        registry.fill(HIST("hEtachi2"), fwdtrack.eta());
        if (fwdtrack.chi2MatchMCHMFT() < 50){
          registry.fill(HIST("hXchi2cut50"), fwdtrack.x());
          registry.fill(HIST("hYchi2cut50"), fwdtrack.y());
          registry.fill(HIST("hZchi2cut50"), fwdtrack.z());
          registry.fill(HIST("hPTchi2cut50"), fwdtrack.pt());
          registry.fill(HIST("hPhichi2cut50"), fwdtrack.phi());
          registry.fill(HIST("hTanlchi2cut50"), fwdtrack.tgl());
          registry.fill(HIST("hEtachi2cut50"), fwdtrack.eta());
          if (fwdtrack.chi2MatchMCHMFT() < 40){
            registry.fill(HIST("hXchi2cut40"), fwdtrack.x());
            registry.fill(HIST("hYchi2cut40"), fwdtrack.y());
            registry.fill(HIST("hZchi2cut40"), fwdtrack.z());
            registry.fill(HIST("hPTchi2cut40"), fwdtrack.pt());
            registry.fill(HIST("hPhichi2cut40"), fwdtrack.phi());
            registry.fill(HIST("hTanlchi2cut40"), fwdtrack.tgl());
            registry.fill(HIST("hEtachi2cut40"), fwdtrack.eta());
            if (fwdtrack.chi2MatchMCHMFT() < 30){
              registry.fill(HIST("hXchi2cut30"), fwdtrack.x());
              registry.fill(HIST("hYchi2cut30"), fwdtrack.y());
              registry.fill(HIST("hZchi2cut30"), fwdtrack.z());
              registry.fill(HIST("hPTchi2cut30"), fwdtrack.pt());
              registry.fill(HIST("hPhichi2cut30"), fwdtrack.phi());
              registry.fill(HIST("hTanlchi2cut30"), fwdtrack.tgl());
              registry.fill(HIST("hEtachi2cut30"), fwdtrack.eta());
              if (fwdtrack.chi2MatchMCHMFT() < 20){
                registry.fill(HIST("hXchi2cut20"), fwdtrack.x());
                registry.fill(HIST("hYchi2cut20"), fwdtrack.y());
                registry.fill(HIST("hZchi2cut20"), fwdtrack.z());
                registry.fill(HIST("hPTchi2cut20"), fwdtrack.pt());
                registry.fill(HIST("hPhichi2cut20"), fwdtrack.phi());
                registry.fill(HIST("hTanlchi2cut20"), fwdtrack.tgl());
                registry.fill(HIST("hEtachi2cut20"), fwdtrack.eta());
                if (fwdtrack.chi2MatchMCHMFT() < 10){
                  registry.fill(HIST("hXchi2cut10"), fwdtrack.x());
                  registry.fill(HIST("hYchi2cut10"), fwdtrack.y());
                  registry.fill(HIST("hZchi2cut10"), fwdtrack.z());
                  registry.fill(HIST("hPTchi2cut10"), fwdtrack.pt());
                  registry.fill(HIST("hPhichi2cut10"), fwdtrack.phi());
                  registry.fill(HIST("hTanlchi2cut10"), fwdtrack.tgl());
                  registry.fill(HIST("hEtachi2cut10"), fwdtrack.eta());
                }
              }
            }
          }
        }
      }
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack){
        registry.fill(HIST("hXmuon"), fwdtrack.x());
        registry.fill(HIST("hYmuon"), fwdtrack.y());
        registry.fill(HIST("hZmuon"), fwdtrack.z());
        registry.fill(HIST("hPTmuon"), fwdtrack.pt());
        registry.fill(HIST("hPhimuon"), fwdtrack.phi());
        registry.fill(HIST("hTanlmuon"), fwdtrack.tgl());
        registry.fill(HIST("hEtamuon"), fwdtrack.eta());
        for (auto& fwdtrackglobal : fwdtracks) {
          if (fwdtrackglobal.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
            if (fwdtrackglobal.matchMCHTrackId() == fwdtrack.globalIndex()) {
              registry.fill(HIST("hXmuonwithmft"), fwdtrack.x());
              registry.fill(HIST("hYmuonwithmft"), fwdtrack.y());
              registry.fill(HIST("hZmuonwithmft"), fwdtrack.z());
              registry.fill(HIST("hPTmuonwithmft"), fwdtrack.pt());
              registry.fill(HIST("hPhimuonwithmft"), fwdtrack.phi());
              registry.fill(HIST("hTanlmuonwithmft"), fwdtrack.tgl());
              registry.fill(HIST("hEtamuonwithmft"), fwdtrack.eta());
            }
          }
        }
      }
    }

    for (auto& [fwdtrack, mfttrack] : combinations(CombinationsFullIndexPolicy(fwdtracks, mfttracks))) {

      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack){
        if (0 <= fwdtrack.collisionId() - mfttrack.collisionId() && fwdtrack.collisionId() - mfttrack.collisionId() < cfgColWindow) {
          double result = matchONNX(fwdtrack, mfttrack);
          registry.fill(HIST("hmlscore"), result);

          static constexpr Double_t MatchingPlaneZ = -77.5;

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

          double px = fwdtrack.p() * sin(M_PI/2 - atan(mfttrack.tgl())) * cos(mfttrack.phi());
          double py = fwdtrack.p() * sin(M_PI/2 - atan(mfttrack.tgl())) * sin(mfttrack.phi());

          Float_t Delta_XY     = sqrt((muonpars1.getX() - mftpars1.getX())*(muonpars1.getX() - mftpars1.getX()) + (muonpars1.getY() - mftpars1.getY())*(muonpars1.getY() - mftpars1.getY()));
          if (result > cfgThrScore){
            registry.fill(HIST("hXml"), mfttrack.x());
            registry.fill(HIST("hYml"), mfttrack.y());
            registry.fill(HIST("hZml"), mfttrack.z());
            registry.fill(HIST("hPTml"), std::sqrt(std::pow(px,2) + std::pow(py,2)));
            registry.fill(HIST("hPhiml"), mfttrack.phi());
            registry.fill(HIST("hTanlml"), mfttrack.tgl());
            registry.fill(HIST("hEtaml"), mfttrack.eta());

            registry.fill(HIST("hXmlfwdtrue"), muonpars1.getX());
            registry.fill(HIST("hYmlfwdtrue"), muonpars1.getY());
            registry.fill(HIST("hZmlfwdtrue"), muonpars1.getY());
            registry.fill(HIST("hPTmlfwdtrue"), muonpars1.getPt());
            registry.fill(HIST("hPhimlfwdtrue"), muonpars1.getPhi());
            registry.fill(HIST("hTanlmlfwdtrue"), muonpars1.getTanl());
            registry.fill(HIST("hEtamlfwdtrue"), muonpars1.getEta());
            registry.fill(HIST("hXmlmfttrue"), mftpars1.getX());
            registry.fill(HIST("hYmlmfttrue"), mftpars1.getY());
            registry.fill(HIST("hZmlmfttrue"), mftpars1.getZ());
            registry.fill(HIST("hPTmlmfttrue"), mftpars1.getPt());
            registry.fill(HIST("hPhimlmfttrue"), mftpars1.getPhi());
            registry.fill(HIST("hTanlmlmfttrue"), mftpars1.getTanl());
            registry.fill(HIST("hEtamlmfttrue"), mftpars1.getEta());
            registry.fill(HIST("hXmldeltatrue"), muonpars1.getX() - mftpars1.getX());
            registry.fill(HIST("hYmldeltatrue"), muonpars1.getY() - mftpars1.getY());
            registry.fill(HIST("hZmldeltatrue"), muonpars1.getZ() - mftpars1.getZ());
            registry.fill(HIST("hPTmldeltatrue"), muonpars1.getPt() - mftpars1.getPt());
            registry.fill(HIST("hPhimldeltatrue"), muonpars1.getPhi() - mftpars1.getPhi());
            registry.fill(HIST("hTanlmldeltatrue"), muonpars1.getTanl() - mftpars1.getTanl());
            registry.fill(HIST("hEtamldeltatrue"), muonpars1.getEta() - mftpars1.getEta());
            registry.fill(HIST("hXYmldeltatrue"), Delta_XY);
          } else {
            if (Delta_XY < 3) {
              registry.fill(HIST("hXmlfwdfalse"), muonpars1.getX());
              registry.fill(HIST("hYmlfwdfalse"), muonpars1.getY());
              registry.fill(HIST("hZmlfwdfalse"), muonpars1.getY());
              registry.fill(HIST("hPTmlfwdfalse"), muonpars1.getPt());
              registry.fill(HIST("hPhimlfwdfalse"), muonpars1.getPhi());
              registry.fill(HIST("hTanlmlfwdfalse"), muonpars1.getTanl());
              registry.fill(HIST("hEtamlfwdfalse"), muonpars1.getEta());
              registry.fill(HIST("hXmlmftfalse"), mftpars1.getX());
              registry.fill(HIST("hYmlmftfalse"), mftpars1.getY());
              registry.fill(HIST("hZmlmftfalse"), mftpars1.getZ());
              registry.fill(HIST("hPTmlmftfalse"), mftpars1.getPt());
              registry.fill(HIST("hPhimlmftfalse"), mftpars1.getPhi());
              registry.fill(HIST("hTanlmlmftfalse"), mftpars1.getTanl());
              registry.fill(HIST("hEtamlmftfalse"), mftpars1.getEta());
              registry.fill(HIST("hXmldeltafalse"), muonpars1.getX() - mftpars1.getX());
              registry.fill(HIST("hYmldeltafalse"), muonpars1.getY() - mftpars1.getY());
              registry.fill(HIST("hZmldeltafalse"), muonpars1.getZ() - mftpars1.getZ());
              registry.fill(HIST("hPTmldeltafalse"), muonpars1.getPt() - mftpars1.getPt());
              registry.fill(HIST("hPhimldeltafalse"), muonpars1.getPhi() - mftpars1.getPhi());
              registry.fill(HIST("hTanlmldeltafalse"), muonpars1.getTanl() - mftpars1.getTanl());
              registry.fill(HIST("hEtamldeltafalse"), muonpars1.getEta() - mftpars1.getEta());
              registry.fill(HIST("hXYmldeltafalse"), Delta_XY);
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
    adaptAnalysisTask<compareMlMatching>(cfgc)
  };
}
