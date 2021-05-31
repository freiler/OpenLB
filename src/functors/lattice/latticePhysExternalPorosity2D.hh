/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_PHYS_EXTERNAL_POROSITY_2D_HH
#define LATTICE_PHYS_EXTERNAL_POROSITY_2D_HH

#include <vector>
#include <cmath>
#include <limits>

#include "latticePhysExternalPorosity2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry2D.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "core/blockLattice2D.h"
#include "communication/mpiManager.h"
#include "core/blockLatticeStructure2D.h"


namespace olb {

template <typename T,typename DESCRIPTOR>
SuperLatticePhysExternalPorosity2D<T,DESCRIPTOR>::SuperLatticePhysExternalPorosity2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "ExtPorosityField";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalPorosity2D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysExternalPorosity2D<T,DESCRIPTOR>::BlockLatticePhysExternalPorosity2D(
  BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2),
    _overlap(overlap)
{
  this->getName() = "ExtPorosityField";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternalPorosity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap
  ).template computeField<descriptors::POROSITY>(output);
  return true;
}

}
#endif
