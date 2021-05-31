/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2006-2019 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod, Adrian Kummerländer
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

/* cylinder2d.cpp:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
typedef D2Q9<> DESCRIPTOR;


// Parameters for the simulation setup
const int N = 10;       // resolution of the model
const T Re = 20.;       // Reynolds number
const T maxPhysT = 16;  // max. simulation time in s, SI unit
const T L = 0.1/N;      // latticeL
const T lengthX = 2.2;
const T lengthY = .41+L;
const T centerCylinderX = 0.2;
const T centerCylinderY = 0.2+L/2.;
const T radiusCylinder = 0.05;


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter,
                      SuperGeometry2D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend( lengthX,lengthY );
  Vector<T,2> center( centerCylinderX,centerCylinderY );
  Vector<T,2> origin;
  IndicatorCircle2D<T> circle( center, radiusCylinder );

  superGeometry.rename( 0,2 );

  superGeometry.rename( 2,1,1,1 );

  // Set material number for inflow
  extend[0] = 2.*L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );
  // Set material number for outflow
  origin[0] = lengthX-L;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );
  // Set material number for cylinder
  superGeometry.rename( 1,5,circle );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice2D<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T, DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     SuperGeometry2D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});
  sLattice.defineDynamics( bulkIndicator, &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  // Setting of the boundary conditions

  //if boundary conditions are chosen to be local
	setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
	setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);

	//if boundary conditions are chosen to be interpolated
	//setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, converter.getLatticeRelaxationFrequency(), superGeometry, 3);
	//setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, converter.getLatticeRelaxationFrequency(), superGeometry, 4);

	// Material=5 -->bounce back
  //sLattice.defineDynamics(superGeometry, 5, &instances::getBounceBack<T, DESCRIPTOR>());

  // Material=5 -->bouzidi

  Vector<T,2> center( centerCylinderX,centerCylinderY );
  IndicatorCircle2D<T> circle ( center, radiusCylinder );

  sLattice.defineDynamics( superGeometry, 5, &instances::getNoDynamics<T,DESCRIPTOR>() );
  setBouzidiZeroVelocityBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 5, circle);

  // Initial conditions
  AnalyticalConst2D<T,T> rhoF( 1 );
  std::vector<T> velocity( 2,T( 0 ) );
  AnalyticalConst2D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                        SuperGeometry2D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
  int iTupdate = 5;

  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = L/2.;
    Poiseuille2D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall );

    sLattice.defineU( superGeometry, 3, poiseuilleU );
  }

}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry2D<T>& superGeometry, Timer<T>& timer,
                 CircularBuffer<T>& buffer )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "cylinder2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR> quality( sLattice, converter);
  SuperRoundingF2D<T> roundedQuality ( quality, RoundingMode::NearestInteger);
  SuperDiscretizationF2D<T> discretization ( roundedQuality, 0., 2. );

  vtmWriter.addFunctor( quality );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int vtkIter  = converter.getLatticeTime( .3 );
  const int statIter = converter.getLatticeTime( .1 );

  T point[2] = {};
  point[0] = centerCylinderX + 3*radiusCylinder;
  point[1] = centerCylinderY;
  AnalyticalFfromSuperF2D<T> intpolateP( pressure, true );
  T p;
  intpolateP( &p,point );
  buffer.insert(p);

  if ( iT == 0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtk files
  if ( iT%vtkIter == 0 && iT > 0 ) {
    vtmWriter.write( iT );

    {
      SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
      BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }

    {
      BlockReduction2D2D<T> planeReduction( discretization, 600, BlockDataSyncMode::ReduceOnly );
      heatmap::plotParam<T> jpeg_scale;
      jpeg_scale.name = "quality";
      jpeg_scale.colour = "blackbody";
      heatmap::write( planeReduction, iT, jpeg_scale );
    }

  }

  // Gnuplot constructor (must be static!)
  // for real-time plotting: gplot("name", true) // experimental!
  static Gnuplot<T> gplot( "drag" );

  // write pdf at last time step
  if ( iT == converter.getLatticeTime( maxPhysT )-1 ) {
    // writes pdf
    gplot.writePDF();
  }

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    clout << "Circular buffer test: moving average pointwise value=" << buffer.average() << std::endl;

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T> intpolatePressure( pressure, true );
    SuperLatticePhysDrag2D<T,DESCRIPTOR> drag( sLattice, superGeometry, 5, converter );


    T point1[2] = {};
    T point2[2] = {};

    point1[0] = centerCylinderX - radiusCylinder;
    point1[1] = centerCylinderY;

    point2[0] = centerCylinderX + radiusCylinder;
    point2[1] = centerCylinderY;

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    int input[3] = {};
    T _drag[drag.getTargetDim()];
    drag( _drag,input );
    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << endl;

    // set data for gnuplot: input={xValue, yValue(s), names (optional), position of key (optional)}
    gplot.setData( converter.getPhysTime( iT ), {_drag[0], 5.58}, {"drag(openLB)", "drag(schaeferTurek)"}, "bottom right", {'l','l'} );
    // writes a png in one file for every timestep, if the file is open it can be used as a "liveplot"
    gplot.writePNG();

    // every (iT%vtkIter) write an png of the plot
    if ( iT%( vtkIter ) == 0 ) {
      // writes pngs: input={name of the files (optional), x range for the plot (optional)}
      gplot.writePNG( iT, maxPhysT );
    }
  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},                        // resolution: number of voxels per charPhysL
    (T)   0.56,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   2.0*radiusCylinder,       // charPhysLength: reference length of simulation geometry
    (T)   0.2,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.2*2.*radiusCylinder/Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("cylinder2d");

  // === 2rd Step: Prepare Geometry ===
  Vector<T,2> extend( lengthX,lengthY );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry2D<T> cuboidGeometry( cuboid, L, noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice( superGeometry );

  BGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() );

  //prepareLattice and set boundaryConditions
  prepareLattice( sLattice, converter, bulkDynamics, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  CircularBuffer<T> buffer(converter.getLatticeTime(.2));
  clout << "starting simulation..." << endl;
  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, buffer );
  }

  timer.stop();
  timer.printSummary();
}
