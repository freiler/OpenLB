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

///This file contains the Regularized Heat Flux Boundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_REGULARIZED_HEAT_FLUX_BOUNDARY_H
#define SET_REGULARIZED_HEAT_FLUX_BOUNDARY_H

#include <vector>
#include <list>
#include "geometry/blockGeometryStatistics2D.h"
#include "core/superLattice2D.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry2D.h"
#include "utilities/functorPtr.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "core/blockLatticeStructure2D.h"
#include "advectionDiffusionMomentaOnBoundaries.h"
#include "boundaryPostProcessors2D.h"
#include "dynamics/dynamics.h"
#include "dynamics/advectionDiffusionDynamics.h"
#include "advectionDiffusionBoundaries.h"
#include "geometry/blockGeometry2D.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"

namespace olb {
///Initialising the RegularizedHeatFluxBoundary on the superLattice domain
///This is an advection diffusion boundary -->MixinDynamics = AdvectionDiffusionRLBdynamics
template<typename T, typename DESCRIPTOR, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,DESCRIPTOR>>
void setRegularizedHeatFluxBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice,T omega, SuperGeometry2D<T>& superGeometry, int material, T *heatFlux=nullptr);

///Initialising the RegularizedHeatFluxBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setRegularizedHeatFluxBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice,T omega, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,  T *heatFlux=nullptr);

/// Set RegularizedHeatFluxBoundary for indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setRegularizedHeatFluxBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, T omega, BlockIndicatorF2D<T>& indicator, T *heatFlux);


}//namespace olb


#endif
