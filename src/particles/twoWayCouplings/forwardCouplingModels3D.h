/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Models for Lagrangian forward-coupling methods -- header file.
 */

#ifndef LB_FORWARD_COUPLING_MODELS_H
#define LB_FORWARD_COUPLING_MODELS_H

#include "functors/lattice/reductionF3D.h"
#include "twoWayHelperFunctionals.h"
#include "smoothingFunctionals3D.h"

namespace olb {

/** Abstact base class for all the forward-coupling models
  * Its raison d'etre consists of not being templetized in Lattice.
  */
template<typename T, template<typename V> class Particle>
class ForwardCouplingModel {
public:
  /// Class operator to apply the coupling, for overload.
  virtual bool operator() (Particle<T>* p, int globic)=0;
protected:
  /// Constructor
  ForwardCouplingModel( SuperGeometry3D<T>& sGeometry, DragModel<T, Particle>& dragModel );
  SuperGeometry3D<T>& _sGeometry;
  DragModel<T, Particle>& _dragModel; // reference to a drag model
};

/** Abstact base class for all the local/non-local forward-coupling models.
  * Adding the templatization in Lattice.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class BaseForwardCouplingModel : public ForwardCouplingModel<T,Particle> {
protected:
  /// Constructor
  BaseForwardCouplingModel ( UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry,
                         DragModel<T, Particle>& dragModel );
  UnitConverter<T, Lattice>& _converter; // reference to a UnitConverter
  SuperLattice3D<T, Lattice>& _sLattice; // reference to a lattice
  // Functional to interpolate lattice density at particle's location
  std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > _interpLatticeDensity;
  // Functional to interpolate lattice velocity at particle's location
  std::shared_ptr<SuperLatticeInterpPhysVelocity3D<T, Lattice> > _interpLatticeVelocity;
  // Momentum-exchange helper functional
  std::shared_ptr<TwoWayHelperFunctional<T, Lattice> > _momentumExchange;
};

/** Abstact class for all the local forward-coupling models,
  * viz., momentum coupling from fluid to particle.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class LocalBaseForwardCouplingModel : public BaseForwardCouplingModel<T,Lattice,Particle> {
public:
  /// Class operator to apply the coupling, for overload.
  virtual bool operator() (Particle<T>* p, int globic) override;
  //virtual bool operator() (SuperParticleSystem3D<T,Particle>& spSys) override;
protected:
  /// Constructor
  LocalBaseForwardCouplingModel ( UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry,
                         DragModel<T, Particle>& dragModel );
};

/** Class for a naive forward-coupling model.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class NaiveForwardCouplingModel : public LocalBaseForwardCouplingModel<T,Lattice,Particle> {
public:
  /// Constructor
  NaiveForwardCouplingModel ( UnitConverter<T, Lattice>& converter,
                              SuperLattice3D<T, Lattice>& sLattice,
                              SuperGeometry3D<T>& sGeometry,
                              DragModel<T, Particle>& dragModel );
};

/** Class for a forward-coupling model following the Ladd's mechanism.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class LaddForwardCouplingModel : public LocalBaseForwardCouplingModel<T,Lattice,Particle> {
public:
  /// Constructor
  LaddForwardCouplingModel ( UnitConverter<T, Lattice>& converter,
                             SuperLattice3D<T, Lattice>& sLattice,
                             SuperGeometry3D<T>& sGeometry,
                             DragModel<T, Particle>& dragModel );
};

/** Abstact class for all the non-local forward-coupling models,
  * viz., momentum coupling from fluid to particle.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class NonLocalBaseForwardCouplingModel : public BaseForwardCouplingModel<T,Lattice,Particle> {
public:
  /// Class operator to apply the coupling, for overload.
  virtual bool operator() (Particle<T>* p, int globic) override;
//  virtual bool operator() (SuperParticleSystem3D<T,Particle>& spSys) override;
protected:
  /// Constructor
  NonLocalBaseForwardCouplingModel ( UnitConverter<T, Lattice>& converter,
                              SuperLattice3D<T, Lattice>& sLattice,
                              SuperGeometry3D<T>& sGeometry,
                              DragModel<T, Particle>& dragModel,
                              SmoothingFunctional<T, Lattice>& smoothingFunctional );
  SmoothingFunctional<T, Lattice>& _smoothingFunctional; // Functional to treat non-local smoothing
};

/** Class for a naive, non-local forward-coupling model as in Sungkorn et al. (2011), but with
  * an extra-normalization of the smoothing function.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class NaiveNonLocalForwardCouplingModel : public NonLocalBaseForwardCouplingModel<T,Lattice,Particle> {
public:
  /// Constructor
  NaiveNonLocalForwardCouplingModel ( UnitConverter<T, Lattice>& converter,
                                 SuperLattice3D<T, Lattice>& sLattice,
                                 SuperGeometry3D<T>& sGeometry,
                                 DragModel<T, Particle>& dragModel,
                                 SmoothingFunctional<T, Lattice>& smoothingFunctional );
};

/** Class for a forward-coupling model as in Evrard, Denner and van Wachem (2019), but with
  * an extra-normalization of the smoothing function.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class vanWachemForwardCouplingModel : public BaseForwardCouplingModel<T,Lattice,Particle> {
public:
  /// Constructor
  vanWachemForwardCouplingModel ( UnitConverter<T, Lattice>& converter,
                                 SuperLattice3D<T, Lattice>& sLattice,
                                 SuperGeometry3D<T>& sGeometry,
                                 DragModel<T, Particle>& dragModel,
                                 SmoothingFunctional<T, Lattice>& smoothingFunctional, int nVoxelInterpPoints );
  /// Class operator to apply the coupling, for overload.
  virtual bool operator() (Particle<T>* p, int globic) override;
protected:
  SmoothingFunctional<T, Lattice>& _smoothingFunctional; // Functional to treat non-local smoothing
};


}

#endif
