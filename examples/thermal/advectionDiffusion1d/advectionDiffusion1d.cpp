
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Louis Kronberg, Stephan Simonis
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

/*  advectionDiffusion1d:
 *  The solution to a linear, scalar, one-dimensional advection-diffusion
 *  equation is approximated. The computation takes place in 2D and is
 *  projected to 1D by slicing the domain along a centerline.
 *  The numerical setup and the analytical solution are taken from
 *  [Simonis, S., Frank, M., and Krause, M. J. 2020. Phil. Trans. R. Soc. A378:
 *  20190400. DOI: 10.1098/rsta.2019.0400]
 *  Error norms are calculated for three subsequent resolutions and stored
 *  in the respective /tmp folders. A python script is provided to calculate
 *  the experimental order of convergence towards the analytical solution.
 */


#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

typedef D2Q5<VELOCITY> TDESCRIPTOR;

const int runs = 3;                // # simulations with increasing resolution

const int N0 = 50;                 // initial # discrete points per dimension
const int statIter0 = N0;          // initial # lattice output timesteps
const T mue0 = 1.5;                // physical diffusivity
const T peclet0 = 40./3.;          // Peclet number (Pe = u*L/mue)
const T physLength0 = 2.;          // physical domain x length



template <typename T>
class AdePhysTemp1D : public AnalyticalF2D<T,T> {

protected:
  T t;
  T mue;
  T uMag;
  T res;
public:
  AdePhysTemp1D(T time, AdeUnitConverter<T, TDESCRIPTOR> converter) : AnalyticalF2D<T,T>(1),
    t(time),
    mue(converter.getPhysDiffusivity()),
    uMag(converter.getCharPhysVelocity()),
    res(converter.getResolution())

  {};

  bool operator()(T output[], const T input[]) override
  {
    T x = input[0];
    T gf = res/(res+1.);

    // initial condition (pseudo 1D)
    output[0] = sin(M_PI*(x-uMag*t)*gf) * exp(-mue*pow(M_PI,2)*t*pow(gf,2));

    return true;
  };
};


void prepareGeometry(SuperGeometry2D<T>& superGeometry,
                     IndicatorF2D<T>& indicator )
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 1); // , indicator);
  superGeometry.communicate();

  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice(  SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                      Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                      SuperGeometry2D<T>& superGeometry,
                      AdeUnitConverter<T, TDESCRIPTOR> converter )
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // buffer layer
  ADlattice.defineDynamics(superGeometry.getMaterialIndicator({0}), &instances::getNoDynamics<T, TDESCRIPTOR>());
  // bulk
  ADlattice.defineDynamics(superGeometry.getMaterialIndicator({1}), &advectionDiffusionBulkDynamics);

  AnalyticalConst2D<T,T> u0( converter.getCharLatticeVelocity(), 0.);
  AdePhysTemp1D<T> Tinit( 0.0, converter );

  auto bulkIndicator = superGeometry.getMaterialIndicator({0,1});

  ADlattice.defineU( bulkIndicator, u0 );
  ADlattice.defineRho( bulkIndicator, Tinit );
  ADlattice.iniEquilibrium( bulkIndicator, Tinit, u0 );

  /// Make the lattice ready for simulation
  ADlattice.initialize();
  clout << "Prepare Lattice ... OK" << std::endl;
}

T errorOverLine( SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                 SuperGeometry2D<T>& superGeometry,
                 int iT,
                 AdeUnitConverter<T, TDESCRIPTOR> converter )
{
  OstreamManager clout(std::cout,"error");
  Gnuplot<T> plt("UNUSED_FILE");

  T lx = converter.getCharPhysLength();
  T dist = converter.getPhysDeltaX();

  SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
  AnalyticalFfromSuperF2D<T> aTemp( temperature, true, 1);
  AdePhysTemp1D<T> temperatureSol( converter.getPhysTime(iT), converter);

  T result[1] = {0.};    // will contain result of analytical solution
  T resultSim[1] = {0.}; // will contain result of lattice evaluation
  T input[2] = {T(), 0}; // lattice node at which to evaluate the temparature

  int ndatapoints = converter.getResolution(); // number of data points on line
  T tempMSE = .0;                              // mean squared error
  T tempSquaredSum = .0;                       // sum of analytical temperature
  for (int i = 0; i <= ndatapoints; i++) {
    input[0] = -.5*lx + i*dist;
    aTemp(resultSim, input);
    temperatureSol(result, input);

    tempMSE += (result[0]  - resultSim[0]) * (result[0]  - resultSim[0]);
    // remove division by zero
    tempSquaredSum += (result[0]+1) * (result[0]+1);
    plt.datFileOut( input[0], resultSim[0], "simulation" + to_string(converter.getPhysTime(iT)));
    plt.datFileOut( input[0], result[0], "analytical" + to_string(converter.getPhysTime(iT)));
  }
  plt.datFileOut(converter.getPhysTime(iT), sqrt(tempMSE/tempSquaredSum), "averageL2RelError");
  return sqrt(tempMSE/tempSquaredSum);
}


T getResults( SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
              int iT,
              int statIter,
              AdeUnitConverter<T, TDESCRIPTOR> converter,
              SuperGeometry2D<T>& superGeometry,
              Timer<T>& timer)
{
  OstreamManager clout(std::cout,"getResults");
  SuperVTMwriter2D<T> vtkWriter("advectionDiffusion1d");

  AdePhysTemp1D<T> temperatureSol( converter.getPhysTime(iT), converter );

  T avg = 0;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, TDESCRIPTOR> geometry(ADlattice, superGeometry);
    SuperLatticeCuboid2D<T, TDESCRIPTOR> cuboid(ADlattice);
    SuperLatticeRank2D<T, TDESCRIPTOR> rank(ADlattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();

    SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> solution(temperatureSol, ADlattice);

    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( solution );
    vtkWriter.write( iT );

    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print( iT, converter.getPhysTime(iT) );

    /// Output data for plot of inital temperatre value of lattice
    Gnuplot<T> plt("UNUSED_FILE");
    T result[1] = {0.};
    T input[2] = {T(), 0};
    for (int i = 0; i <= converter.getResolution(); i++) {
      input[0] = -.5*converter.getCharPhysLength() + i*converter.getPhysDeltaX();
      temperatureSol(result, input);
      plt.datFileOut(input[0],result[0], "initalValue");
    }
  }
  else if (iT % statIter == 0 ) {
    /// Writes the VTK files
    SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> solution(temperatureSol, ADlattice);

    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( solution );
    vtkWriter.write( iT );

    /// ADLattice statistics console output
    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print(iT, converter.getPhysTime( iT ));

    // compute relative error
    avg = errorOverLine(ADlattice, superGeometry, iT, converter);
    clout << "Relative L2-error norm over centerline: "  << avg << std::endl;
  }

  return avg;
}


void simulate(int N, int statIter, T mue, T peclet, T physLength)
{

  OstreamManager clout(std::cout,"simulate");
  clout << "Executing the simulation with N=" << to_string(N) << std::endl;

  AdeUnitConverter<T,TDESCRIPTOR> converter(
    physLength/N,            // physDeltaX
    pow(physLength/N, 2),    // physDeltaT (diffusive scaling)
    physLength,              // charPhysLength
    peclet*mue/physLength,   // charPhysVelocity from Peclet
    mue,                     // physDiffusivity
    1                        // physDensity,
  );

  converter.print();

  // switch outdirectory only if there are multiple simulation runs
  if (runs > 1) {
    singleton::directories().setOutputDir("./tmp/N" + to_string(N) + "/");
  }

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  std::vector<T> origin(2,T());

  extend[0] = converter.getCharPhysLength();
  extend[1] = converter.getCharPhysLength();

  origin[0] = - converter.getCharPhysLength()/2;
  origin[1] = - converter.getCharPhysLength()/2;

  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidGeometry
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), noOfCuboids);
  cuboidGeometry.setPeriodicity(true, true);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(superGeometry, cuboid);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, TDESCRIPTOR> ADlattice(superGeometry);

  AdvectionDiffusionBGKdynamics<T, TDESCRIPTOR> TbulkDynamics(
    converter.getLatticeAdeRelaxationFrequency(),
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>());

  prepareLattice( ADlattice,
                  TbulkDynamics,
                  superGeometry,
                  converter );

  /// === 4th Step: Main Loop with Timer ===
  Timer<T> timer( superGeometry.getStatistics().getNvoxel() );
  timer.start();

  T pulseDiffBound = 1e-3;
  int timeCount = ceil( -1./converter.getPhysDeltaT() * log(pulseDiffBound)/(converter.getPhysDiffusivity()*pow(M_PI,2)) /statIter );
  int iTmax = timeCount * statIter;

  int iT;
  T simulationAverage = .0;
  for (iT = 0; iT < iTmax; ++iT) {
    simulationAverage += getResults(ADlattice, iT, statIter, converter, superGeometry, timer);

    /// === 6th Step: Collide and Stream Execution ===
    ADlattice.collideAndStream();
  }

  simulationAverage /= timeCount;

  // this outputs into ./tmp/gnuplotData/data/averageSimL2RelErr
  singleton::directories().setOutputDir("./tmp/");
  Gnuplot<T> plt("UNUSED");
  plt.datFileOut(N, simulationAverage, "averageSimL2RelErr");

  clout << "Simulation Average Relative L2 Error: " << simulationAverage << std::endl;

  timer.stop();
  timer.printSummary();
}


int main(int argc, char *argv[])
{
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);

  singleton::directories().setOutputDir("./tmp/");

  for (int i = 0; i < runs; ++i) {
    simulate( pow(2,i) * N0,
              pow(4,i) * statIter0,
              mue0,
              peclet0,
              physLength0 );
  }
}
