/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Local boundary cell 2D dynamics -- header file.
 */
#ifndef MOMENTA_ON_BOUNDARIES_2D_H
#define MOMENTA_ON_BOUNDARIES_2D_H

#include "momentaOnBoundaries.h"

namespace olb {

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
class InnerCornerVelBM2D : public DirichletBoundaryMomenta<T,DESCRIPTOR> {
public:
  /// Default Constructor: initialization to zero
  InnerCornerVelBM2D();
  /// Constructor with boundary initialization
  InnerCornerVelBM2D(const T u_[DESCRIPTOR::d]);

  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override;
  void computeU (
    ConstCell<T,DESCRIPTOR>& cell,
    T u[DESCRIPTOR::d] ) const override;
  void computeJ (
    ConstCell<T,DESCRIPTOR>& cell,
    T j[DESCRIPTOR::d] ) const override;
  void computeU(T u[DESCRIPTOR::d]) const;
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override ;
  void defineU(Cell<T,DESCRIPTOR>& cell,
               const T u[DESCRIPTOR::d]) override ;
  void defineU(const T u[DESCRIPTOR::d]);
  void defineAllMomenta (
    Cell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d],
    const T pi[util::TensorVal<DESCRIPTOR >::n] ) override;
  /// Stress tensor
  void computeStress (
    ConstCell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const override;
private:
  T _u[DESCRIPTOR::d];   ///< value of the velocity on the boundary
};


}

#endif
