/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Jakob Mangold, Mathias J. Krause
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

#ifndef SUPER_STATISTIC_F3D_HH
#define SUPER_STATISTIC_F3D_HH

#include "superStatisticF3D.h"
#include "indicator/superIndicatorF3D.h"

namespace olb {


template <typename T, typename W>
SuperVarianceF3D<T,W>::SuperVarianceF3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                                    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                              	  T expectedValue)
  : SuperF3D<T,W>(f->getSuperStructure(), f->getTargetDim()+1),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF)),
    _expectedValue(expectedValue)
{
  this->getName() = "Variance("+_f->getName()+")";

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();
  CuboidGeometry3D<T>& cuboid = _f->getSuperStructure().getCuboidGeometry();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockVarianceF3D<T,W>(_f->getBlockF(iC),
                                _indicatorF->getBlockIndicatorF(iC),
                                cuboid.get(load.glob(iC)),
                                _expectedValue)
      );
    }
  }
}

template <typename T, typename W>
SuperVarianceF3D<T,W>::SuperVarianceF3D(FunctorPtr<SuperF3D<T,W>>&& f,
                                    SuperGeometry3D<T>& superGeometry,
                                    const int material,
                              	  T expectedValue)
  : SuperVarianceF3D(
      std::forward<decltype(f)>(f),
      superGeometry.getMaterialIndicator(material),
      expectedValue)
{ }

template <typename T, typename W>
bool SuperVarianceF3D<T,W>::operator() (W output[], const int input[])
{
  _f->getSuperStructure().communicate();
  CuboidGeometry3D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&     load     = _f->getSuperStructure().getLoadBalancer();

  std::size_t voxels(0);

  if (this->_blockF.empty()) {
    W outputTmp[_f->getTargetDim()];
    int inputTmp[_f->getSourceDim()];

    for (int i = 0; i < _f->getTargetDim(); ++i) {
      output[i] = W(0);
    }

    for (int iC = 0; iC < load.size(); ++iC) {
      const Cuboid3D<T> cuboid = geometry.get(load.glob(iC));
      inputTmp[0] = load.glob(iC);
      for (inputTmp[1] = 0; inputTmp[1] < cuboid.getNx(); ++inputTmp[1]) {
        for (inputTmp[2] = 0; inputTmp[2] < cuboid.getNy(); ++inputTmp[2]) {
          for (inputTmp[3] = 0; inputTmp[3] < cuboid.getNz(); ++inputTmp[3]) {
            if (_indicatorF(inputTmp)) {
              _f(outputTmp,inputTmp);
              for (int i = 0; i < _f->getTargetDim(); ++i) {
                output[i] += std::pow(outputTmp[i], 2);
              }
              voxels += 1;
            }
          }
        }
      }
      output[_f->getTargetDim()] = voxels;
    }
  }
  else{
    for (int iC = 0; iC < load.size(); ++iC) {
      this->getBlockF(iC)(output, input);
    }
  }
    
#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim() + 1; ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
#endif


  for (int i = 0; i < _f->getTargetDim(); ++i) {
  	output[i] = output[i] - std::pow(_expectedValue * output[_f->getTargetDim()], 2) / output[_f->getTargetDim()];
  }



  return true;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template <typename T, typename W>
SuperStdDeviationF3D<T,W>::SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                                    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                              	  T expectedValue)
  : SuperF3D<T,W>(f->getSuperStructure(), f->getTargetDim()+1),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF)),
    _expectedValue(expectedValue)
{
  this->getName() = "Variance("+_f->getName()+")";

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();
  CuboidGeometry3D<T>& cuboid = _f->getSuperStructure().getCuboidGeometry();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockStdDeviationF3D<T,W>(_f->getBlockF(iC),
                                _indicatorF->getBlockIndicatorF(iC),
                                cuboid.get(load.glob(iC)),
                                _expectedValue)
      );
    }
  }
}

template <typename T, typename W>
SuperStdDeviationF3D<T,W>::SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&& f,
                                    SuperGeometry3D<T>& superGeometry,
                                    const int material,
                              	  T expectedValue)
  : SuperStdDeviationF3D(
      std::forward<decltype(f)>(f),
      superGeometry.getMaterialIndicator(material),
      expectedValue)
{ }

template <typename T, typename W>
bool SuperStdDeviationF3D<T,W>::operator() (W output[], const int input[])
{
  _f->getSuperStructure().communicate();
  CuboidGeometry3D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&     load     = _f->getSuperStructure().getLoadBalancer();

  std::size_t voxels(0);

  W outputTmp[_f->getTargetDim()];
  int inputTmp[_f->getSourceDim()];

  for (int i = 0; i < _f->getTargetDim(); ++i) {
    output[i] = W(0);
  }
  /*
  for (int iC = 0; iC < load.size(); ++iC) {
    const Cuboid3D<T> cuboid = geometry.get(load.glob(iC));
    inputTmp[0] = load.glob(iC);
    for (inputTmp[1] = 0; inputTmp[1] < cuboid.getNx(); ++inputTmp[1]) {
      for (inputTmp[2] = 0; inputTmp[2] < cuboid.getNy(); ++inputTmp[2]) {
        for (inputTmp[3] = 0; inputTmp[3] < cuboid.getNz(); ++inputTmp[3]) {
          if (_indicatorF(inputTmp)) {
            _f(outputTmp,inputTmp);
            for (int i = 0; i < _f->getTargetDim(); ++i) {
              output[i] += std::pow(outputTmp[i] - _expectedValue, 2);
            }
            voxels += 1;
          }
        }
      }
    }
    output[_f->getTargetDim()] = voxels;
  }

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim() + 1; ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
#endif
    */
  SuperVarianceF3D<T,W>::operator()(output, input, _expectedValue);
  
  
  for (int i = 0; i < _f->getTargetDim(); ++i) {
  	output[i] = std::sqrt(output[i]);
  }



  return true;
}



}

#endif

