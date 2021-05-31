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

#ifndef SET_INTERPOLATED_PRESSURE_BOUNDARY_2D_HH
#define SET_INTERPOLATED_PRESSURE_BOUNDARY_2D_HH

#include "setInterpolatedPressureBoundary2D.h"

namespace olb {
///Initialising the setInterpolatedPressureBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedPressureBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry2D<T>& superGeometry, int material)
{
  setInterpolatedPressureBoundary<T,DESCRIPTOR, MixinDynamics>(sLattice,omega, superGeometry.getMaterialIndicator(material));
}

///Initialising the setInterpolatedPressureBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedPressureBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setLocalPressureBoundary2D");
  bool includeOuterCells = false;
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 1;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setInterpolatedPressureBoundary<T,DESCRIPTOR, MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc),omega,
        indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

////////// BlockLattice Domain  /////////////////////////////////////////

/// Set interpolated pressure boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedPressureBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block,T omega, BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  /*
   *x0,x1,y0,y1 Range of cells to be traversed
   **/
  int x0 = margin;
  int y0 = margin;
  int x1 = blockGeometryStructure.getNx()-1 -margin;
  int y1 = blockGeometryStructure.getNy()-1 -margin;
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      //momenta vector provisionally inside src/core/blockLatticeStructure2D.h
      Momenta<T, DESCRIPTOR>* momenta = nullptr;
      //dynamics vector provisionally inside src/core/blockLatticeStructure2D.h
      Dynamics<T, DESCRIPTOR>* dynamics = nullptr;
      PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = nullptr;
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          //set momenta,dynamics and post processors on indicated pressure boundary cells
          if (discreteNormal[1] == -1) {
            momenta = new BasicDirichletBM<T,DESCRIPTOR,PressureBM,0,-1>;
            dynamics = new MixinDynamics(omega, *momenta);
            postProcessor = new StraightFdBoundaryProcessorGenerator2D < T,DESCRIPTOR,0,-1 >  (iX, iX, iY, iY);
          }
          else if (discreteNormal[1] == 1) {
            momenta = new BasicDirichletBM<T,DESCRIPTOR,PressureBM,0,1>;
            dynamics = new MixinDynamics(omega, *momenta);
            postProcessor = new StraightFdBoundaryProcessorGenerator2D < T,DESCRIPTOR,0,1 >  (iX, iX, iY, iY);
          }
          else if (discreteNormal[2] == -1) {
            momenta = new BasicDirichletBM<T,DESCRIPTOR,PressureBM,1,-1>;
            dynamics = new MixinDynamics(omega, *momenta);
            postProcessor = new StraightFdBoundaryProcessorGenerator2D < T,DESCRIPTOR,1,-1 >  (iX, iX, iY, iY);
          }
          else if (discreteNormal[2] == 1) {
            momenta = new BasicDirichletBM<T,DESCRIPTOR,PressureBM,1,1>;
            dynamics = new MixinDynamics(omega, *momenta);
            postProcessor = new StraightFdBoundaryProcessorGenerator2D < T,DESCRIPTOR,1,1 >  (iX, iX, iY, iY);
          }
          //initialised in boundary/setLocalVelocityBoundary2D.h/hh
          //sets boundary on indicated cellss
          setBoundary<T, DESCRIPTOR, MixinDynamics>(block, omega,iX,iY, momenta, dynamics, postProcessor);

        }
      }

    }
  }
}


}//namespace olb

#endif

