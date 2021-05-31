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

#ifndef LATTICE_PHYS_CROPPED_PERMEABILITY_3D_HH
#define LATTICE_PHYS_CROPPED_PERMEABILITY_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysCroppedPermeability3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry3D.h"
#include "blockBaseF3D.h"
#include "core/blockLatticeStructure3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::SuperLatticePhysCroppedPermeability3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "cropped_permeability";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); iC++ ) {
    this->_blockF.emplace_back( new BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>(
                                  this->_sLattice.getBlockLattice(iC), this->getConverter() ) );
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::BlockLatticePhysCroppedPermeability3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1)
{
  this->getName() = "cropped_permeability";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T value;
  // get porosity
  this->_blockLattice.get(input[0], input[1], input[2]).template computeField<descriptors::POROSITY>(
    &value);
  // convert to physPermeability
//  output[0] = _converter.physPermeability(value); // TODO converter MG
  // \todo Why do we need this???
  if (output[0] >= 42 && output[0] <= 42 && output[0] != 42) {
    output[0] = 1;
  }
  if (std::isinf(output[0])) {
    output[0] = 1;
  }
  if (output[0] > 1.) {
    output[0] = 1.;
  }
  return true;
}

}
#endif
