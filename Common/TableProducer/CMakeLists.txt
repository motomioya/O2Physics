# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

add_subdirectory(Converters)
add_subdirectory(PID)

o2physics_add_dpl_workflow(trackextension
                    SOURCES trackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2::DetectorsBase
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(trackselection
                    SOURCES trackselection.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(event-selection
                    SOURCES eventSelection.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCCDB
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multiplicity-table
                    SOURCES multiplicityTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(multiplicity-extra-table
                    SOURCES multiplicityExtraTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(centrality-table
                    SOURCES centralityTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(timestamp
                    SOURCES timestamp.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(weak-decay-indices
                    SOURCES weakDecayIndices.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(ft0-corrected-table
                    SOURCES ft0CorrectedTable.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(track-propagation
                    SOURCES trackPropagation.cxx
                    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(track-propagation-tester
                    SOURCES trackPropagationTester.cxx
                    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore O2Physics::trackSelectionRequest
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(calo-clusters
                    SOURCES caloClusterProducer.cxx
                    PUBLIC_LINK_LIBRARIES O2::DataFormatsPHOS O2::PHOSBase O2::PHOSReconstruction O2Physics::DataModel
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fwdtrackextension
                    SOURCES fwdtrackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2::DetectorsBase
                                          O2::DetectorsCommonDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mfttrackextension
                    SOURCES mfttrackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2::DetectorsBase
                                          O2::DetectorsCommonDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mchmftmatching-info
                    SOURCES mchmftmatchinginfo.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mchmftmatching-infoem
                    SOURCES mchmftmatchinginfoem.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mchmftmatching-infoemmc
                    SOURCES mchmftmatchinginfoemmc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(similarmchtrack-info
                    SOURCES similarmchtrackinfo.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(globalfwdtrackinfo
                    SOURCES globalFwdtrackInfo.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mccheck
                    SOURCES mccheck.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mccheckpi
                    SOURCES mccheckpi.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2Physics::PWGDQCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(quickcheck
                    SOURCES quickcheck.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::GlobalTracking
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(jpsicheck
                    SOURCES jpsicheck.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

                  #o2physics_add_dpl_workflow(mlcheck
                  #SOURCES mlcheck.cxx
                  #PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  #O2::ReconstructionDataFormats
                  #COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(masscheck
                    SOURCES masscheck.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(phicheck
                    SOURCES phicheck.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(hfcheck
                    SOURCES hfcheck.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(checkambiguous
                    SOURCES checkambiguous.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mftsignefficiency
                    SOURCES mftsignefficiency.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(track-to-collision-associator
                    SOURCES trackToCollisionAssociator.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(fwdtrack-to-collision-associator
                    SOURCES fwdtrackToCollisionAssociator.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(mccollisionextra
                    SOURCES mcCollsExtra.cxx
                    PUBLIC_LINK_LIBRARIES O2::DetectorsBase O2Physics::AnalysisCore
                    COMPONENT_NAME Analysis)

                  #o2physics_add_dpl_workflow(comparemlmatching
                  #  SOURCES compareMlMatching.cxx
                  #  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ONNXRuntime::ONNXRuntime O2::CCDB O2Physics::MLCore O2::ReconstructionDataFormats
                  #  COMPONENT_NAME Analysis)

                  #o2physics_add_dpl_workflow(checkmleff
                  #SOURCES checkMlEff.cxx
                  #PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ONNXRuntime::ONNXRuntime O2::CCDB O2Physics::MLCore O2::ReconstructionDataFormats
                  #COMPONENT_NAME Analysis)

                  #o2physics_add_dpl_workflow(checkmlparameters
                  #SOURCES checkMlParameters.cxx
                  #PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ONNXRuntime::ONNXRuntime O2::CCDB O2Physics::MLCore O2::ReconstructionDataFormats
                  #COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qvector-table
                    SOURCES qVectorsTable.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                                          O2Physics::AnalysisCore
                                          O2::DetectorsCommonDataFormats
                                          O2Physics::DataModel
                                          O2::CCDB
                    COMPONENT_NAME Analysis)

                  #o2physics_add_dpl_workflow(mftmchmatchingml
                  #SOURCES mftmchMatchingML.cxx
                  #PUBLIC_LINK_LIBRARIES O2::Framework
                  #O2Physics::AnalysisCore
                  #ONNXRuntime::ONNXRuntime
                  #O2::CCDB O2Physics::MLCore
                  #O2::ReconstructionDataFormats
                  #COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(runwithmchwithmft
                    SOURCES runwithmchwithmft.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework
                                          O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(checkeventmixing
                    SOURCES checkeventmixing.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(colhasmuon
                    SOURCES colhasmuon.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(etacheck
                    SOURCES etacheck.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                    COMPONENT_NAME Analysis)
