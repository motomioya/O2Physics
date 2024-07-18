// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef COMMON_DATAMODEL_COLHASMUON_H_
#define COMMON_DATAMODEL_COLHASMUON_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace collisionhasmuon
{
  DECLARE_SOA_COLUMN(ColHasMuon, colhasmuon, bool);     //! stores N times this PV was recoed
}
DECLARE_SOA_TABLE(CHasMuon, "AOD", "CHASMUON",
                  o2::soa::Index<>,
                  collisionhasmuon::ColHasMuon)
}

#endif // COMMON_DATAMODEL_COLLISIONASSOCIATIONTABLES_H_
