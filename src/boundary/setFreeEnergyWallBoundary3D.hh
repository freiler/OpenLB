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
#ifndef SET_FREE_ENERGY_WALL_BOUNDARY_3D_HH
#define SET_FREE_ENERGY_WALL_BOUNDARY_3D_HH

#include "setFreeEnergyWallBoundary3D.h"

namespace olb {
/// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
/// dynamics and an FreeEnergyWall PostProcessor.
/// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
/// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
/// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly

////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry, int material,
                               T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber)
{

  setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material),
                                          alpha, kappa1, kappa2, h1, h2, latticeNumber);

}

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                               T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber)
{
  OstreamManager clout(std::cout, "setFreeEnergyWallBoundary");
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
  T addend = 0;
  if (latticeNumber==1) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) );
  }
  else if (latticeNumber==2) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (-h2/kappa2) );
  }
  else if (latticeNumber==3) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) );
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice.getExtendedBlockLattice(iCloc),
                                            indicator->getExtendedBlockIndicatorF(iCloc), addend, latticeNumber, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}




/// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
/// dynamics and an FreeEnergyWall PostProcessor.
/// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
/// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa3_ - Parameter related to surface tension. [lattice units]
/// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h3_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly

////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry, int material,
                               T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2,T h3, int latticeNumber)
{
  setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material),
                                          alpha, kappa1, kappa2, kappa3, h1, h2, h3, latticeNumber);
}


///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                               T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber)
{
  OstreamManager clout(std::cout, "setFreeEnergyWallBoundary");
  bool includeOuterCells = false;
  int _overlap = indicator->getSuperGeometry().getOverlap();
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  T addend = 0;
  if (latticeNumber==1) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) + (h3/kappa3) );
  }
  else if (latticeNumber==2) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (-h2/kappa2) );
  }
  else if (latticeNumber==3) {
    addend = 1./(alpha*alpha) * ( (h3/kappa3) );
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice.getExtendedBlockLattice(iCloc),
                                            indicator->getExtendedBlockIndicatorF(iCloc), addend, latticeNumber, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

////////// BlockLattice Domain  /////////////////////////////////////////


//set FreeEnergyWallBoundary on block domain.
//This function works for the setFreeEnergyWallBoundary with h1,h2,h3 Parameters and h1,h2 Parameters
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator,
                               T addend, int latticeNumber, bool includeOuterCells)
{
  bool _output = false;
  OstreamManager clout(std::cout, "setFreeEnergyWallBoundary");
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  /*
   *x0,x1,y0,y1, z0, z1 Range of cells to be traversed
   **/
  int x0 = margin;
  int y0 = margin;
  int z0 = margin;
  int x1 = blockGeometryStructure.getNx()-1 -margin;
  int y1 = blockGeometryStructure.getNy()-1 -margin;
  int z1 = blockGeometryStructure.getNz()-1 -margin;
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        ///set dynamics and post Processors on indicated cells
        if (indicator(iX,iY,iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ, true);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {

            Dynamics<T, DESCRIPTOR>* dynamics = NULL;
            if (latticeNumber == 1) {
              dynamics = &instances::getBounceBack<T, DESCRIPTOR>();
            }
            else {
              dynamics = new FreeEnergyWallDynamics<T, DESCRIPTOR>;
              _block.dynamicsVector.push_back(dynamics);
            }
            _block.get(iX,iY,iZ).defineDynamics(dynamics);

            PostProcessorGenerator3D<T, DESCRIPTOR>* wettingPostProcessor =
              new FreeEnergyWallProcessorGenerator3D<T, DESCRIPTOR> ( iX, iX, iY, iY, iZ, iZ,
                  discreteNormal[1], discreteNormal[2], discreteNormal[3], addend );
            PostProcessorGenerator3D<T, DESCRIPTOR>* chemPotPostProcessor =
              new FreeEnergyChemPotBoundaryProcessorGenerator3D<T, DESCRIPTOR> ( iX, iX, iY, iY, iZ, iZ,
                  discreteNormal[1], discreteNormal[2], discreteNormal[3], latticeNumber );
            if (wettingPostProcessor) {
              _block.addPostProcessor(*wettingPostProcessor);
            }
            if (chemPotPostProcessor) {
              _block.addPostProcessor(*chemPotPostProcessor);
            }
          }
          if (_output) {
            clout << "setFreeEnergyWallBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << ")" << std::endl;
          }
        }
      }
    }
  }


}
}//namespace olb

#endif
