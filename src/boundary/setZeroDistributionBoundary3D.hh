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

//This file contains the Zero Distribution Boundary.
//This is a new version of the Boundary, which only contains free floating functions
//This boundary is an Advection Diffusion Boundary
//All functions, which are contained in this file are transferred and modified from:
//superBoundaryCondition3D.h/hh
//advectionDiffusionBoundaries.h/hh
//advectionDiffusionBoundaryInstantiator3D.h/hh
#ifndef SET_ZERO_DISTRIBUTION_BOUNDARY_3D_HH
#define SET_ZERO_DISTRIBUTION_BOUNDARY_3D_HH

#include "setZeroDistributionBoundary3D.h"

namespace olb {


//setZeroDistributionBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setZeroDistributionBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice,SuperGeometry3D<T>& superGeometry, int material)
{
  setZeroDistributionBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, superGeometry.getMaterialIndicator(material));


}

//setZeroDistributionBoundary function on the superLattice domain
//depending on the application, the first function can be skipped and this function can be called directly in the app
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setZeroDistributionBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  bool includeOuterCells = false;
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 1;
  OstreamManager clout(std::cout, "setZeroDistributionBoundary");
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    //sets ZeroDistributionBoundary on the block level
    setZeroDistributionBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc),
        indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  //the addPoints2CommBC is currently located inside setLocalVelocityBoundary3D.h
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);


}

//sets the ZeroDistributionBoundary on the block level
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setZeroDistributionBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setZeroDistributionBoundary");
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4, 0);
  int x0 = margin;
  int y0 = margin;
  int z0 = margin;
  int x1 = blockGeometryStructure.getNx()-1 -margin;
  int y1 = blockGeometryStructure.getNy()-1 -margin;
  int z1 = blockGeometryStructure.getNz()-1 -margin;
  //sets the boundary on cells
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
            PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new ZeroDistributionBoundaryProcessorGenerator3D<T, DESCRIPTOR>(iX, iX, iY, iY, iZ, iZ, -discreteNormal[1], -discreteNormal[2], -discreteNormal[3]);
            if (postProcessor) {
              _block.addPostProcessor(*postProcessor);
            }
          }
          else {
            clout << "Warning: Could not setZeroDistributionBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] << "," << discreteNormal[3] << ")" << std::endl;
          }
        }
      }
    }
  }

}


} //namespace olb


#endif

