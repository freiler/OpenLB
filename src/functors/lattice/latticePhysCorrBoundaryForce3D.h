/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef LATTICE_PHYS_CORR_BOUNDARY_FORCE_3D_H
#define LATTICE_PHYS_CORR_BOUNDARY_FORCE_3D_H

#include<vector>

#include "superBaseF3D.h"
#include "superCalcF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "core/superLattice3D.h"
#include "blockBaseF3D.h"
#include "geometry/blockGeometry3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/blockIndicatorBaseF3D.h"
#include "dynamics/smagorinskyBGKdynamics.h"
#include "dynamics/porousBGKdynamics.h"


/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
/// see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysCorrBoundaryForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperLatticePhysCorrBoundaryForce3D(SuperLattice3D<T,DESCRIPTOR>&      sLattice,
                                      FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysCorrBoundaryForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                      SuperGeometry3D<T>& superGeometry, const int material,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
};

/**
 *  functor returns pointwise phys force acting on a indicated boundary on local lattice
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysCorrBoundaryForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockIndicatorF3D<T>&        _indicatorF;
  const BlockGeometryStructure3D<T>& _blockGeometry;
public:
  BlockLatticePhysCorrBoundaryForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                      BlockIndicatorF3D<T>& indicatorF,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
