/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

///This file contains the Interpolated Pressure Boundary
///This is a new version of the Boundary, which only contains free floating functions

#ifndef SET_INTERPOLATED_PRESSURE_BOUNDARY_2D_H
#define SET_INTERPOLATED_PRESSURE_BOUNDARY_2D_H


#include <vector>
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"
#include "geometry/superGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "geometry/blockGeometry2D.h"
#include "core/superLattice2D.h"
#include "core/blockLatticeStructure2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "dynamics/dynamics.h"
#include "momentaOnBoundaries2D.h"
#include "boundaryPostProcessors2D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setInterpolatedPressureBoundary function on the superLattice domain
///Interpolated Boundaries use the BGKdynamics collision-operator
template<typename T, typename DESCRIPTOR, class MixinDynamics=BGKdynamics<T,DESCRIPTOR>>
void setInterpolatedPressureBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice,T omega, SuperGeometry2D<T>& superGeometry, int material);

///Initialising the setInterpolatedPressureBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedPressureBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF2D<T>>&& indicator);

////////// BlockLattice Domain  /////////////////////////////////////////

/// Add interpolated pressure boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedPressureBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block,T omega, BlockIndicatorF2D<T>& indicator, bool includeOuterCells=false);

}//namespace olb

#endif
