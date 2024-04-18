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
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "DetectorsBase/Propagator.h"
#include "MFTTracking/Tracker.h"
#include <math.h>
#include <TLorentzVector.h>
#include <string>
#include <regex>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace std;

struct mccheck {

  //MCSignal* mySignal;
  HistogramRegistry registry{
    "registry", //
    {
    },
  };

  void init(o2::framework::InitContext&)
  {
  }

  //void process(soa::Filtered<aod::McParticles> const& mcTracks)
  void process(aod::McCollisions const& mcCollisions, aod::McParticles const& mcTracks)
  {
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

                                                                                                                //
    for (auto& collision : mcCollisions) {
      LOGF(info, "---------collision---------");
      for (auto& mctrackHF : mcTracks) {
        if (mctrackHF.mcCollisionId() == collision.globalIndex()){
          if ( (mctrackHF.pdgCode() >= 500 && mctrackHF.pdgCode() <= 549) || (mctrackHF.pdgCode() <= -500 && mctrackHF.pdgCode() >= -549) || (mctrackHF.pdgCode() >= 5000 && mctrackHF.pdgCode() <= 5490) || (mctrackHF.pdgCode() <= -5000 && mctrackHF.pdgCode() >= -5490)) {
            LOGF(info, "------B------");
            LOGF(info, "mctrackHF.globalIndex() = %d", mctrackHF.globalIndex());
            LOGF(info, "mctrackHF.pdgCode() = %d", mctrackHF.pdgCode());
            LOGF(info, "mctrackHF.fromBackgroundEvent() = %d", mctrackHF.fromBackgroundEvent());
          }
        }
      }
      for (auto& mctrack1 : mcTracks) {
        if (mctrack1.mcCollisionId() == collision.globalIndex()){
          if (mctrack1.pdgCode() == 13 || mctrack1.pdgCode() == -13){
            for (auto& mctrack2 : mcTracks) {
              if (mctrack2.mcCollisionId() == collision.globalIndex()){
                if (mctrack1.globalIndex() != mctrack2.globalIndex()){
                  if (mctrack2.pdgCode() == 13 || mctrack2.pdgCode() == -13){
                    if (signalbb->CheckSignal(true, mctrack1, mctrack2) || signalbb1->CheckSignal(true, mctrack1, mctrack2) || signalbb2->CheckSignal(true, mctrack1, mctrack2) || signalbb3->CheckSignal(true, mctrack1, mctrack2) || signalbb4->CheckSignal(true, mctrack1, mctrack2) || signalbb5->CheckSignal(true, mctrack1, mctrack2) || signalbb6->CheckSignal(true, mctrack1, mctrack2) || signalbb7->CheckSignal(true, mctrack1, mctrack2)){
                      if (signalbb->CheckSignal(true, mctrack1, mctrack2)){
                        LOGF(info, "---------signalbb---------");
                      } else if (signalbb1->CheckSignal(true, mctrack1, mctrack2)){
                        LOGF(info, "---------signalbb1---------");
                      } else if (signalbb2->CheckSignal(true, mctrack1, mctrack2)){
                        LOGF(info, "---------signalbb2---------");
                      } else if (signalbb3->CheckSignal(true, mctrack1, mctrack2)){
                        LOGF(info, "---------signalbb3---------");
                      } else if (signalbb4->CheckSignal(true, mctrack1, mctrack2)){
                        LOGF(info, "---------signalbb4---------");
                      } else if (signalbb5->CheckSignal(true, mctrack1, mctrack2)){
                        LOGF(info, "---------signalbb5---------");
                      } else if (signalbb6->CheckSignal(true, mctrack1, mctrack2)){
                        LOGF(info, "---------signalbb6---------");
                      } else if (signalbb7->CheckSignal(true, mctrack1, mctrack2)){
                        LOGF(info, "---------signalbb7---------");
                      }
                      const auto mcmothersidlist1 = mctrack1.mothersIds();
                      const auto mcmothersidlist2 = mctrack2.mothersIds();
                      LOGF(info, "------mu1------");
                      LOGF(info, "mctrack.globalIndex() = %d", mctrack1.globalIndex());
                      LOGF(info, "mctrack.pdgCode() = %d", mctrack1.pdgCode());
                      if (mcmothersidlist1.size() > 0) {
                        for (int i = 0; i <= mcmothersidlist1.size() - 1; i++) {
                          auto mcmother = mcTracks.iteratorAt(mcmothersidlist1[i]);
                          if ( (mcmother.pdgCode() >= 500 && mcmother.pdgCode() <= 549) || (mcmother.pdgCode() <= -500 && mcmother.pdgCode() >= -549) || (mcmother.pdgCode() >= 5000 && mcmother.pdgCode() <= 5490) || (mcmother.pdgCode() <= -5000 && mcmother.pdgCode() >= -5490) || (mcmother.pdgCode() >= 400 && mcmother.pdgCode() <= 449) || (mcmother.pdgCode() <= -400 && mcmother.pdgCode() >= -449) || (mcmother.pdgCode() >= 4000 && mcmother.pdgCode() <= 4490) || (mcmother.pdgCode() <= -4000 && mcmother.pdgCode() >= -4490) ) {
                            LOGF(info, "muon1, mcmother.globalIndex() = %d", mcmother.globalIndex());
                            LOGF(info, "muon1, mcmother.pdgCode() = %d", mcmother.pdgCode());
                            const auto mcgmothersidlist = mcmother.mothersIds();
                            if (mcgmothersidlist.size() > 0) {
                              for (int i = 0; i <= mcgmothersidlist.size() - 1; i++) {
                                auto mcgmother = mcTracks.iteratorAt(mcgmothersidlist[i]);
                                LOGF(info, "muon1, mcgmother.globalIndex() = %d", mcgmother.globalIndex());
                                LOGF(info, "muon1, mcgmother.pdgCode() = %d", mcgmother.pdgCode());
                                const auto mcggmothersidlist = mcgmother.mothersIds();
                                if (mcggmothersidlist.size() > 0) {
                                  for (int i = 0; i <= mcggmothersidlist.size() - 1; i++) {
                                    auto mcggmother = mcTracks.iteratorAt(mcggmothersidlist[i]);
                                    LOGF(info, "muon1, mcggmother.globalIndex() = %d", mcggmother.globalIndex());
                                    LOGF(info, "muon1, mcggmother.pdgCode() = %d", mcggmother.pdgCode());
                                    const auto mcgggmothersidlist = mcggmother.mothersIds();
                                    if (mcgggmothersidlist.size() > 0) {
                                      for (int i = 0; i <= mcgggmothersidlist.size() - 1; i++) {
                                        auto mcgggmother = mcTracks.iteratorAt(mcgggmothersidlist[i]);
                                        LOGF(info, "muon1, mcgggmother.globalIndex() = %d", mcgggmother.globalIndex());
                                        LOGF(info, "muon1, mcgggmother.pdgCode() = %d", mcgggmother.pdgCode());
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      LOGF(info, "------mu2------");
                      LOGF(info, "mctrack.globalIndex() = %d", mctrack2.globalIndex());
                      LOGF(info, "mctrack.pdgCode() = %d", mctrack2.pdgCode());
                      if (mcmothersidlist2.size() > 0) {
                        for (int i = 0; i <= mcmothersidlist2.size() - 1; i++) {
                          auto mcmother = mcTracks.iteratorAt(mcmothersidlist2[i]);
                          if ( (mcmother.pdgCode() >= 500 && mcmother.pdgCode() <= 549) || (mcmother.pdgCode() <= -500 && mcmother.pdgCode() >= -549) || (mcmother.pdgCode() >= 5000 && mcmother.pdgCode() <= 5490) || (mcmother.pdgCode() <= -5000 && mcmother.pdgCode() >= -5490) || (mcmother.pdgCode() >= 400 && mcmother.pdgCode() <= 449) || (mcmother.pdgCode() <= -400 && mcmother.pdgCode() >= -449) || (mcmother.pdgCode() >= 4000 && mcmother.pdgCode() <= 4490) || (mcmother.pdgCode() <= -4000 && mcmother.pdgCode() >= -4490) ) {
                            LOGF(info, "muon2, mcmother.globalIndex() = %d", mcmother.globalIndex());
                            LOGF(info, "muon2, mcmother.pdgCode() = %d", mcmother.pdgCode());
                            const auto mcgmothersidlist = mcmother.mothersIds();
                            if (mcgmothersidlist.size() > 0) {
                              for (int i = 0; i <= mcgmothersidlist.size() - 1; i++) {
                                auto mcgmother = mcTracks.iteratorAt(mcgmothersidlist[i]);
                                LOGF(info, "muon2, mcgmother.globalIndex() = %d", mcgmother.globalIndex());
                                LOGF(info, "muon2, mcgmother.pdgCode() = %d", mcgmother.pdgCode());
                                const auto mcggmothersidlist = mcgmother.mothersIds();
                                if (mcggmothersidlist.size() > 0) {
                                  for (int i = 0; i <= mcggmothersidlist.size() - 1; i++) {
                                    auto mcggmother = mcTracks.iteratorAt(mcggmothersidlist[i]);
                                    LOGF(info, "muon2, mcggmother.globalIndex() = %d", mcggmother.globalIndex());
                                    LOGF(info, "muon2, mcggmother.pdgCode() = %d", mcggmother.pdgCode());
                                    const auto mcgggmothersidlist = mcggmother.mothersIds();
                                    if (mcgggmothersidlist.size() > 0) {
                                      for (int i = 0; i <= mcgggmothersidlist.size() - 1; i++) {
                                        auto mcgggmother = mcTracks.iteratorAt(mcgggmothersidlist[i]);
                                        LOGF(info, "muon2, mcgggmother.globalIndex() = %d", mcgggmother.globalIndex());
                                        LOGF(info, "muon2, mcgggmother.pdgCode() = %d", mcgggmother.pdgCode());
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
              }
            }
          }
        }
      }
      /*
      for (auto& mctrack : mcTracks) {
        if (mctrack.mcCollisionId() == collision.globalIndex()){
          const auto mcmothersidlist = mctrack.mothersIds();
          if ( (mctrack.pdgCode() >= 500 && mctrack.pdgCode() <= 549) || (mctrack.pdgCode() <= -500 && mctrack.pdgCode() >= -549) || (mctrack.pdgCode() >= 5000 && mctrack.pdgCode() <= 5490) || (mctrack.pdgCode() <= -5000 && mctrack.pdgCode() >= -5490) ) {
            LOGF(info, "------B------");
            LOGF(info, "mctrack.globalIndex() = %d", mctrack.globalIndex());
            LOGF(info, "mctrack.pdgCode() = %d", mctrack.pdgCode());
            LOGF(info, "mctrack.FromBackgroundEvent() = %d", mctrack.fromBackgroundEvent());
            if (mcmothersidlist.size() > 0) {
              for (int i = 0; i <= mcmothersidlist.size() - 1; i++) {
                auto mcmother = mcTracks.iteratorAt(mcmothersidlist[i]);
                LOGF(info, "---this B has mother---");
                LOGF(info, "mcmother.globalIndex() = %d", mcmother.globalIndex());
                LOGF(info, "mcmother.pdgCode() = %d", mcmother.pdgCode());
              }
            }
          } else if (mctrack.pdgCode() == 13 || mctrack.pdgCode() == -13){
            if (mcmothersidlist.size() == 0) {
            } else if (mcmothersidlist.size() > 0) {
              for (int i = 0; i <= mcmothersidlist.size() - 1; i++) {
                auto mcmother = mcTracks.iteratorAt(mcmothersidlist[i]);
                if ( (mcmother.pdgCode() >= 500 && mcmother.pdgCode() <= 549) || (mcmother.pdgCode() <= -500 && mcmother.pdgCode() >= -549) || (mcmother.pdgCode() >= 5000 && mcmother.pdgCode() <= 5490) || (mcmother.pdgCode() <= -5000 && mcmother.pdgCode() >= -5490) ) {
                  LOGF(info, "------mu, mother is B------");
                  LOGF(info, "mctrack.globalIndex() = %d", mctrack.globalIndex());
                  LOGF(info, "mctrack.pdgCode() = %d", mctrack.pdgCode());
                  LOGF(info, "mcmother.globalIndex() = %d", mcmother.globalIndex());
                  LOGF(info, "mcmother.pdgCode() = %d", mcmother.pdgCode());
                }
              }
            }
          }
        }
      }
      */
    }
  }
};
  


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mccheck>(cfgc)
  };
}
