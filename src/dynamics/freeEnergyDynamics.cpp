/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Sam Avis
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#include "freeEnergyDynamics.h"
#include "freeEnergyDynamics.hh"
#include "latticeDescriptors.h"
#include "core/cell.hh"

namespace olb {

template class FreeEnergyWallDynamics<double, descriptors::D2Q9<>>;
template class FreeEnergyWallDynamics<double, descriptors::D3Q19<>>;

template class FreeEnergyInletOutletDynamics<double, descriptors::D2Q9<>, 0, -1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D2Q9<>, 0,  1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D2Q9<>, 1, -1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D2Q9<>, 1,  1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D3Q19<>, 0, -1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D3Q19<>, 0,  1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D3Q19<>, 1, -1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D3Q19<>, 1,  1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D3Q19<>, 2, -1>;
template class FreeEnergyInletOutletDynamics<double, descriptors::D3Q19<>, 2,  1>;

}
