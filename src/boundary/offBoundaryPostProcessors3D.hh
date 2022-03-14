/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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

#ifndef OFF_BOUNDARY_POST_PROCESSORS_3D_HH
#define OFF_BOUNDARY_POST_PROCESSORS_3D_HH

#include "offBoundaryPostProcessors3D.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "core/cell.h"

namespace olb {

/////////// LinearBouzidiPostProcessor3D /////////////////////////////////////

/* Bouzidi Interpolation scheme of first order
 *
 * fluid nodes               wall  solid node
 * --o-------<-o->-----<-o->--|----x----
 *            xB         x        xN
 * directions: --> iPop
 *             <-- opp
 *
*/

template<typename T, typename DESCRIPTOR>
ZeroVelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
ZeroVelocityBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
  this->getName() = "ZeroVelocityBouzidiLinearPostProcessor3D";
  this->_priority = -1;
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];

  if (dist >= 0.5) {
    xB = x - c[0];
    yB = y - c[1];
    zB = z - c[2];
    q = 1/(2*dist);
    iPop2 = opp;
  }
  else {
    xB = x;
    yB = y;
    zB = z;
    q = 2*dist;
    iPop2 = iPop;
  }
  /*
  if ( x == 100 && z == 28)
    std::cout << "ZeroVelocityLinear (" << x << "," << y << "," << z <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," << zB <<
      "), dist: " << dist << ", q: " << q << std::endl;
      */

}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  blockLattice.get(x, y, z)[opp] = q*blockLattice.get(xN, yN, zN)[iPop] +
                                   (1-q)*blockLattice.get(xB, yB, zB)[iPop2];
}

template<typename T, typename DESCRIPTOR>
VelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
VelocityBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
  this->getName() = "VelocityBouzidiLinearPostProcessor3D";
  this->_priority = -1;

#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];

  if (dist >= 0.5) {
    xB = x - c[0];
    yB = y - c[1];
    zB = z - c[2];
    q = 1/(2*dist);
    ufrac = q;
    iPop2 = opp;
  }
  else {
    xB = x;
    yB = y;
    zB = z;
    q = 2*dist;
    iPop2 = iPop;
    ufrac = 1;
  }
  /*
    std::cout << "VelocityLinear (" << x << "," << y << "," << z <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," << zB <<
      "), dist: " << dist << ", q: " << q << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void VelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void VelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN, zN);
  T u = ufrac*dynamics->getVelocityCoefficient(iPop);
  auto cellN = blockLattice.get(xN, yN, zN);
  auto cell = blockLattice.get(x, y, z);
  dynamics->defineRho(cellN, cell.computeRho());
  cell[opp] = q*cellN[iPop] + (1-q)*blockLattice.get(xB, yB, zB)[iPop2] + u;
}


//////// CornerBouzidiPostProcessor3D ///////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
ZeroVelocityBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
  this->getName() = "ZeroVelocityBounceBackPostProcessor3D";
  this->_priority = -1;

#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<L>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];
  /*
    std::cout << "Corner (" << x << "," << y << "," << z <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  blockLattice.get(x, y, z)[opp] = blockLattice.get(xN, yN, zN)[iPop];
}

//////// CornerBouzidiPostProcessor3D ///////////////////

template<typename T, typename DESCRIPTOR>
VelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
VelocityBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
  this->getName() = "VelocityBounceBackPostProcessor3D";
  this->_priority = -1;

#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<L>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];

  /*
    std::cout << "Corner (" << x << "," << y << "," << z <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void VelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void VelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN, zN);
  T u = dynamics->getVelocityCoefficient(iPop);
  auto cellN = blockLattice.get(xN, yN, zN);
  auto cell  = blockLattice.get(x, y, z);
  dynamics->defineRho(cellN, cell.computeRho());
  cell[opp] = cellN[iPop] + u;
}

//////// EnthalpyBouzidi ///////////////////

template<typename T, typename DESCRIPTOR>
EnthalpyBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
EnthalpyBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
  this->getName() = "EnthalpyBouzidiLinearPostProcessor3D";
  this->_priority = -1;

#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];

  if (dist >= 0.5) {
    xB = x - c[0];
    yB = y - c[1];
    zB = z - c[2];
    q = 1/(2*dist);
    ufrac = q;
    sign = 1.;
    iPop2 = opp;
  }
  else {
    xB = x;
    yB = y;
    zB = z;
    q = 2*dist;
    ufrac = 1.;
    sign = -1.;
    iPop2 = iPop;
  }
}

template<typename T, typename DESCRIPTOR>
void EnthalpyBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void EnthalpyBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  //const Vector<int,3> c_1 = descriptors::c<DESCRIPTOR>(opp);
  //const Vector<int,3> c_2 = descriptors::c<DESCRIPTOR>(iPop);

  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN, zN);

  T concentration = dynamics->getConcentrationCoefficient();
  auto cellN = blockLattice.get(xN, yN, zN);
  auto cell = blockLattice.get(x, y, z);

  //dynamics->defineRho(cellN, cell.computeRho());
  //if (x==100 && y == 3 && z==27)
    //std::cout << "f5+f6_before: " <<  4.*(cell[opp] + cell[iPop])  << std::endl; //4 is 1/cs^2 --> 4.*(cell[opp] + cell[iPop]) should be rho

  cell[opp] = -q*cellN[iPop] + sign*(1-q)*blockLattice.get(xB, yB, zB)[iPop2] + ufrac*concentration; //Boundary Conditions for thermal lattice Boltzann equation method, Li et al.


  //if (x==200 && y == 3 && z==26){
    //std::cout << "iPop: " << iPop << "iOpp: "<< opp << std::endl;
    //std::cout << "c inward: (" << c_1[0] << "," << c_1[1] << "," << c_1[2] << ")" << "c outward: (" << c_2[0] << "," << c_2[1] << "," << c_2[2] << ")" <<std::endl;
    //std::cout << "Temperature: " << cell.computeRho() << "\n" <<std::endl;
  //}

  //cell[opp] = -1.*cellN[iPop] + ufrac*concentration;

  //std::cout << "Manual Enthalpy: " << manual_enthalpy << " Computed Enthalpy: " << cell.computeRho() << std::endl;
  //std::cout << "f5+f6: " <<  cell[opp] + cell[iPop]  << std::endl;
}

//////// CornerBouzidiPostProcessor3D ///////////////////

template<typename T, typename DESCRIPTOR>
EnthalpyBounceBackPostProcessor3D<T,DESCRIPTOR>::
EnthalpyBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
  this->getName() = "EnthalpyBounceBackPostProcessor3D";
  this->_priority = -1;

#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<L>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];

  /*
    std::cout << "Corner (" << x << "," << y << "," << z <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void EnthalpyBounceBackPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void EnthalpyBounceBackPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN, zN);
  T concentration = dynamics->getConcentrationCoefficient();
  auto cellN = blockLattice.get(xN, yN, zN);
  auto cell = blockLattice.get(x, y, z);
  dynamics->defineRho(cellN, cell.computeRho());
  cell[opp] = -cellN[iPop] + concentration;
}






////////  LinearBouzidiBoundaryPostProcessorGenerator ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::
ZeroVelocityBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroVelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroVelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
VelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::
VelocityBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
VelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new VelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
VelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new VelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
EnthalpyBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::
EnthalpyBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
EnthalpyBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new EnthalpyBouzidiLinearPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
EnthalpyBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new EnthalpyBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

/////////// CornerBouzidiBoundaryPostProcessorGenerator /////////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::
ZeroVelocityBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroVelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityBounceBackPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroVelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
VelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::
VelocityBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
VelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new VelocityBounceBackPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
VelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new VelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
EnthalpyBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::
EnthalpyBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
EnthalpyBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new EnthalpyBounceBackPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
EnthalpyBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new EnthalpyBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

}  // namespace olb

#endif
