/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcekt, Davide Dapelo
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

/** \file
 * BGK Dynamics with adjusted omega -- generic implementation.
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 */
#ifndef POWER_LAW_BGK_DYNAMICS_HH
#define POWER_LAW_BGK_DYNAMICS_HH

#include "descriptorAlias.h"
#include "powerLawBGKdynamics.h"
#include "core/cell.h"
#include "core/util.h"
#include "lbHelpers.h"
#include "latticeDescriptors.h"
#include "math.h"

namespace olb {

////////////////////// Class PowerLawDynamics //////////////////////////

template<typename T, typename DESCRIPTOR>
PowerLawDynamics<T,DESCRIPTOR>::PowerLawDynamics (T m, T n , T nuMin, T nuMax)
  : _m(m),
    _n(n)
{
  setViscoLimits(nuMin, nuMax);
}

template<typename T, typename DESCRIPTOR>
void PowerLawDynamics<T,DESCRIPTOR>::setViscoLimits (T nuMin, T nuMax)
{
  _omegaMin = 1./(nuMax*descriptors::invCs2<T,DESCRIPTOR>() + 0.5);
  _omegaMax = 1./(nuMin*descriptors::invCs2<T,DESCRIPTOR>() + 0.5);
}

template<typename T, typename DESCRIPTOR>
void PowerLawDynamics<T,DESCRIPTOR>::setM (T m)
{
  _m = m;
}

template<typename T, typename DESCRIPTOR>
void PowerLawDynamics<T,DESCRIPTOR>::setN (T n)
{
  _n = n;
}

template<typename T, typename DESCRIPTOR>
T PowerLawDynamics<T,DESCRIPTOR>::getM ()
{
  return _m;
}

template<typename T, typename DESCRIPTOR>
T PowerLawDynamics<T,DESCRIPTOR>::getN ()
{
  return _n;
}

template<typename T, typename DESCRIPTOR>
T PowerLawDynamics<T,DESCRIPTOR>::getNuMin ()
{
  return (1./_omegaMax -0.5) / descriptors::invCs2<T,DESCRIPTOR>();
}

template<typename T, typename DESCRIPTOR>
T PowerLawDynamics<T,DESCRIPTOR>::getNuMax ()
{
  return (1./_omegaMin -0.5) / descriptors::invCs2<T,DESCRIPTOR>();
}

template<typename T, typename DESCRIPTOR>
T PowerLawDynamics<T,DESCRIPTOR>::computeOmegaPL( Cell<T,DESCRIPTOR>& cell, T omega0,
           T rho, T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  T pre2 = pow(descriptors::invCs2<T,DESCRIPTOR>()/2.* omega0/rho,2.); // strain rate tensor prefactor
  //T gamma = sqrt(2.*pre2*PiNeqNormSqr(cell)); // shear rate
  T gamma = sqrt(2.*pre2*lbHelpers<T,DESCRIPTOR>::computePiNeqNormSqr(cell)); // shear rate

  T nuNew = _m*pow(gamma,_n-1.); //nu for non-Newtonian fluid
  T newOmega = 1./(nuNew*descriptors::invCs2<T,DESCRIPTOR>() + 0.5);

  if (newOmega>_omegaMax) {
    newOmega = _omegaMax;
  }
  if (newOmega<_omegaMin) {
    newOmega = _omegaMin;
  }

  //std::cout << nuNew << " " << (1./newOmega -0.5) / descriptors::invCs2<T,DESCRIPTOR>() << "        "
            //<< omega0 << " " << _omegaMin << " " << _omegaMax << " "  << newOmega << std::endl; // <---
  return newOmega;
  //return omega0;
}


////////////////////// Class PowerLawBGKdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
PowerLawBGKdynamics<T,DESCRIPTOR>::PowerLawBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T m, T n , T nuMin, T nuMax)
  : BGKdynamics<T,DESCRIPTOR>(omega, momenta),
    PowerLawDynamics<T,DESCRIPTOR>(m, n, nuMin, nuMax)
{ }

template<typename T, typename DESCRIPTOR>
void PowerLawBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);

  // Computation of the power-law omega.
  // An external is used in place of BGKdynamics::_omega to keep generality and flexibility.
  const auto oldOmega = cell.template getField<descriptors::OMEGA>();
  const auto newOmega = this->computeOmegaPL(cell, oldOmega, rho, pi);

  const T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  cell.template setField<descriptors::OMEGA>(newOmega);
  statistics.incrementStats(rho, uSqr);
}


////////////////////// Class ForcedPowerLawBGKdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
PowerLawForcedBGKdynamics<T,DESCRIPTOR>::PowerLawForcedBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T m, T n , T nuMin, T nuMax)
  : ForcedBGKdynamics<T,DESCRIPTOR>(omega, momenta),
    PowerLawDynamics<T,DESCRIPTOR>(m, n, nuMin, nuMax)
{ }

template<typename T, typename DESCRIPTOR>
void PowerLawForcedBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, pi[util::TensorVal<DESCRIPTOR >::n];
  FieldD<T,DESCRIPTOR,descriptors::VELOCITY> u;
  this->_momenta.computeAllMomenta(cell, rho, u.data(), pi);

  // Computation of the power-law omega.
  // An external is used in place of BGKdynamics::_omega to keep generality and flexibility.
  T oldOmega = cell.template getField<descriptors::OMEGA>();
  T newOmega = this->computeOmegaPL(cell, oldOmega, rho, pi);

  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  u += 0.5 * force;

  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u.data(), this->getOmega());
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u.data(), newOmega, rho);
  cell.template setField<descriptors::OMEGA>(newOmega);
  statistics.incrementStats(rho, uSqr);
}

}

#endif
