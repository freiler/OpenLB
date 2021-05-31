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

#ifndef LATTICE_DENSITY_2D_HH
#define LATTICE_DENSITY_2D_HH

#include <vector>
#include <cmath>
#include <limits>

#include "latticeDensity2D.h"
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

template<typename T,typename DESCRIPTOR>
SuperLatticeDensity2D<T,DESCRIPTOR>::SuperLatticeDensity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "density";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeDensity2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)) );
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticeDensity2D<T,DESCRIPTOR>::BlockLatticeDensity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{
  this->getName() = "density";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticeDensity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = this->_blockLattice.get( input[0], input[1] ).computeRho();
  return true;
}

}
#endif
