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

#ifndef LATTICE_PHYS_EXTERNAL_3D_H
#define LATTICE_PHYS_EXTERNAL_3D_H

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

template <typename T, typename DESCRIPTOR, typename FIELD>
class SuperLatticePhysExternal3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternal3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             T convFactorToPhysUnits);
};

template <typename T, typename DESCRIPTOR, typename FIELD>
class BlockLatticePhysExternal3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternal3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             T convFactorToPhysUnits);
  bool operator() (T output[], const int input[]) override;
private:
  const T _convFactorToPhysUnits;
};

}
#endif
