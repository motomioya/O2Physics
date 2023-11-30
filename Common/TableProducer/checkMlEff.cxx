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

struct checkMlEff {
  Produces<aod::FwdTracksML> fwdtrackml;

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

  Filter etaFilter = ((etalow < aod::fwdtrack::eta) && (aod::fwdtrack::eta < etaup));
  Filter pDcaFilter = (((pDCAcutrAtBsorberEndlow1 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup1) && (aod::fwdtrack::pDca < pDCAcutdcaup1)) || ((pDCAcutrAtBsorberEndlow2 < aod::fwdtrack::rAtAbsorberEnd) && (aod::fwdtrack::rAtAbsorberEnd < pDCAcutrAtBsorberEndup2) && (aod::fwdtrack::pDca < pDCAcutdcaup2)));
  Filter chi2Filter = (aod::fwdtrack::chi2 < chi2up);
  Filter chi2MatchFilter = (aod::fwdtrack::chi2MatchMCHMID < chi2MatchMCHMIDup);

  HistogramRegistry registry{
    "registry",
    {
      {"MuonPt", "ParticlePt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"Chi2Pt", "Chi2Pt", {HistType::kTH2F, {{200, 0, 20},{200,0,200}}}},
      {"Chi2TruePt", "Chi2TruePt", {HistType::kTH2F, {{200, 0, 20},{200,0,200}}}},
      {"MlPt", "MlPt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MlOnePt", "MlOnePt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MlTruePt", "MlOnePt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MlOneTruePt", "MlOnePt", {HistType::kTH1F, {{200, 0, 20}}}},
      {"MuonX", "ParticleX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"Chi2X", "Chi2X", {HistType::kTH2F, {{500, 0, 50},{200,0,200}}}},
      {"Chi2TrueX", "Chi2TrueX", {HistType::kTH2F, {{500, 0, 50},{200,0,200}}}},
      {"MlX", "MlX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlOneX", "MlOneX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlTrueX", "MlOneX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlOneTrueX", "MlOneX", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MuonY", "ParticleY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"Chi2Y", "Chi2Y", {HistType::kTH2F, {{500, 0, 50},{200,0,200}}}},
      {"Chi2TrueY", "Chi2TrueY", {HistType::kTH2F, {{500, 0, 50},{200,0,200}}}},
      {"MlY", "MlY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlOneY", "MlOneY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlTrueY", "MlOneY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MlOneTrueY", "MlOneY", {HistType::kTH1F, {{500, 0, 50}}}},
      {"MuonEta", "ParticleEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"Chi2Eta", "Chi2Eta", {HistType::kTH2F, {{500, -5, 0},{200,0,200}}}},
      {"Chi2TrueEta", "Chi2TrueEta", {HistType::kTH2F, {{500, -5, 0},{200,0,200}}}},
      {"MlEta", "MlEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MlOneEta", "MlOneEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MlTrueEta", "MlOneEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MlOneTrueEta", "MlOneEta", {HistType::kTH1F, {{500, -5, 0}}}},
      {"MuonPhi", "ParticlePhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"Chi2Phi", "Chi2Phi", {HistType::kTH2F, {{1000, -5, 5},{200,0,200}}}},
      {"Chi2TruePhi", "Chi2TruePhi", {HistType::kTH2F, {{1000, -5, 5},{200,0,200}}}},
      {"MlPhi", "MlPhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MlOnePhi", "MlOnePhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MlTruePhi", "MlOnePhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MlOneTruePhi", "MlOnePhi", {HistType::kTH1F, {{1000, -5, 5}}}},
      {"MldeltaColId", "MldeltaColId", {HistType::kTH1F, {{40, -20, 20}}}},
      {"MlOnedeltaColId", "MlOnedeltaColId", {HistType::kTH1F, {{40, -20, 20}}}},
      {"MlTruedeltaColId", "MlOnedeltaColId", {HistType::kTH1F, {{40, -20, 20}}}},
      {"MlOneTruedeltaColId", "MlOnedeltaColId", {HistType::kTH1F, {{40, -20, 20}}}}
    }
  };

  Configurable<std::string> cfgModelFile{"ccdb-file", "model.onnx", "name of ONNX model file"};
  Configurable<float> cfgThrScore{"threshold-score", 0.5, "Threshold value for matching score"};
  Configurable<int> cfgColWindow{"collision-window", 1, "Search window (collision ID) for MFT track"};
  Configurable<float> cfgXYWindow{"XY-window", 3, "Search window (delta XY) for MFT track"};

  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "model-explorer"};
  Ort::SessionOptions session_options;
  std::shared_ptr<Ort::Experimental::Session> onnx_session = nullptr;
  OnnxModel model;

  template <typename F, typename M>
  std::vector<float> getVariables(F const& fwdtrack, M const& mfttrack)
  {

    static constexpr Double_t MatchingPlaneZ = -77.5;

    // propagate muontrack to matching position
    double muonchi2 = fwdtrack.chi2();
    SMatrix5 muonpars(fwdtrack.x(), fwdtrack.y(), fwdtrack.phi(), fwdtrack.tgl(), fwdtrack.signed1Pt());
    std::vector<double> muonv1;
    SMatrix55 muoncovs(muonv1.begin(), muonv1.end());
    o2::track::TrackParCovFwd muonpars1{fwdtrack.z(), muonpars, muoncovs, muonchi2};
    muonpars1.propagateToZlinear(MatchingPlaneZ);

    // propagate mfttrack to matching position
    double mftchi2 = mfttrack.chi2();
    SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
    std::vector<double> mftv1;
    SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
    o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
    mftpars1.propagateToZlinear(MatchingPlaneZ);

    Float_t MFT_X = mftpars1.getX();
    Float_t MFT_Y = mftpars1.getY();
    Float_t MFT_Phi = mftpars1.getPhi();
    Float_t MFT_Tanl = mftpars1.getTanl();

    Float_t MCH_X = muonpars1.getX();
    Float_t MCH_Y = muonpars1.getY();
    Float_t MCH_Phi = muonpars1.getPhi();
    Float_t MCH_Tanl = muonpars1.getTanl();

    Float_t Ratio_X = MFT_X / MCH_X;
    Float_t Ratio_Y = MFT_Y / MCH_Y;
    Float_t Ratio_Phi = MFT_Phi / MCH_Phi;
    Float_t Ratio_Tanl = MFT_Tanl / MCH_Tanl;

    Float_t Delta_X = MFT_X - MCH_X;
    Float_t Delta_Y = MFT_Y - MCH_Y;
    Float_t Delta_Phi = MFT_Phi - MCH_Phi;
    Float_t Delta_Tanl = MFT_Tanl - MCH_Tanl;

    Float_t Delta_XY = sqrt(Delta_X * Delta_X + Delta_Y * Delta_Y);

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
    input_tensor_values = getVariables(fwdtrack, mfttrack);

    if (input_tensor_values[8] < cfgXYWindow) {
      std::vector<Ort::Value> input_tensors;
      input_tensors.push_back(Ort::Experimental::Value::CreateTensor<float>(input_tensor_values.data(), input_tensor_values.size(), input_shape));

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
    model.initModel(cfgModelFile, false, 1);
    onnx_session = model.getSession();
  }

  void process(aod::Collisions const& collisions, soa::Filtered<soa::Join<o2::aod::FwdTracks, aod::McFwdTrackLabels>> const& fwdtracks, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks, aod::McParticles const&)
  {
    for (auto& fwdtrack : fwdtracks) {
      double bestscore = 0;
      int bestmfttrackid = -1;
      auto fwdparticle = fwdtrack.mcParticle();
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        registry.fill(HIST("MuonX"), fwdtrack.x());
        registry.fill(HIST("MuonY"), fwdtrack.y());
        registry.fill(HIST("MuonPhi"), fwdtrack.phi());
        registry.fill(HIST("MuonEta"), fwdtrack.eta());
        registry.fill(HIST("MuonPt"), fwdtrack.pt());
        for (auto& fwdtrackglobal : fwdtracks) {
          if (fwdtrackglobal.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
            if (fwdtrackglobal.matchMCHTrackId() == fwdtrack.globalIndex() ){
              for (auto& mfttrack : mfttracks) {
                if (fwdtrackglobal.matchMFTTrackId() == mfttrack.globalIndex() ){
                  registry.fill(HIST("Chi2X"), fwdtrack.x(), fwdtrackglobal.chi2MatchMCHMFT());
                  registry.fill(HIST("Chi2Y"), fwdtrack.y(), fwdtrackglobal.chi2MatchMCHMFT());
                  registry.fill(HIST("Chi2Phi"), fwdtrack.phi(), fwdtrackglobal.chi2MatchMCHMFT());
                  registry.fill(HIST("Chi2Eta"), fwdtrack.eta(), fwdtrackglobal.chi2MatchMCHMFT());
                  registry.fill(HIST("Chi2Pt"), fwdtrack.pt(), fwdtrackglobal.chi2MatchMCHMFT());
                  auto mftparticle = mfttrack.mcParticle();
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
      if (fwdtrack.trackType() == aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
        for (auto& mfttrack : mfttracks) {
          auto mftparticle = mfttrack.mcParticle();
          if (fwdtrack.has_collision() && mfttrack.has_collision()) {
            if (0 <= fwdtrack.collisionId() - mfttrack.collisionId() && fwdtrack.collisionId() - mfttrack.collisionId() < cfgColWindow) {
              double result = matchONNX(fwdtrack, mfttrack);
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
        if (bestmfttrackid != -1) {
          for (auto& mfttrack : mfttracks) {
            if (mfttrack.globalIndex() == bestmfttrackid) {
              registry.fill(HIST("MlOneX"), fwdtrack.x());
              registry.fill(HIST("MlOneY"), fwdtrack.y());
              registry.fill(HIST("MlOnePhi"), fwdtrack.phi());
              registry.fill(HIST("MlOneEta"), fwdtrack.eta());
              registry.fill(HIST("MlOnePt"), fwdtrack.pt());
              registry.fill(HIST("MlOnedeltaColId"), fwdtrack.collisionId() - mfttrack.collisionId());
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<checkMlEff>(cfgc)};
}
