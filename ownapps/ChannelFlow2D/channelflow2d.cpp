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
const int N = 50.;       // resolution of the model
const T Re = 100.;       // Reynolds number
const T maxPhysT = 100.;  // max. simulation time in s, SI unit
const T L = 1./N;      // latticeL
const T lengthX = 10.;
const T lengthY = 1.;
const T hydraulicD = 2*lengthY;

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter,
                      SuperGeometry2D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend( lengthX,lengthY );
  Vector<T,2> origin;

  superGeometry.rename( 0,2 );

  superGeometry.rename( 2,1,1,1 );

  // Set material number for
  extend[0] = lengthX+2.*L;
  extend[1] = L;
  origin[0] = 0.;
  origin[1] = lengthY;
  IndicatorCuboid2D<T> wall( extend, origin );
  //superGeometry.rename( 2,3,1,wall );

  extend[0] = L;
  extend[1] = lengthY;
  origin[0] = 0.;
  origin[1] = -L;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,4,1,inflow );
  // Set material number for
  extend[0] = L;
  extend[1] = lengthY;
  origin[0] = lengthX;
  origin[1] = 0.;

  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,5,1,outflow );

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
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 4, 5});
  //auto noslipIndicator = superGeometry.getMaterialIndicator({3});
  sLattice.defineDynamics( bulkIndicator, &bulkDynamics );

  //Material=3 --> Symmetry
  //sLattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>());
  //setSlipBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 2);

  // Setting of the boundary conditions
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  //if boundary conditions are chosen to be local
	setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);
	setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 5);

  // Initial conditions
  AnalyticalConst2D<T,T> rhoF( 1 );
  std::vector<T> velocity( 2,T( 0 ) );
  AnalyticalConst2D<T,T> uzero( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  //sLattice.defineRhoU( superGeometry, 3, rhoF, uzero );
  //sLattice.iniEquilibrium( superGeometry, 3, rhoF, uzero );

  velocity[0] = converter.getCharLatticeVelocity();
  AnalyticalConst2D<T,T> u( velocity );

  sLattice.iniEquilibrium( bulkIndicator, rhoF, u );
  sLattice.defineRhoU( bulkIndicator, rhoF, u );
  sLattice.defineRhoU( superGeometry, 2, rhoF, uzero );
  sLattice.iniEquilibrium( superGeometry, 2, rhoF, uzero );
	//if boundary conditions are chosen to be interpolated
	//setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, converter.getLatticeRelaxationFrequency(), superGeometry, 3);
	//setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, converter.getLatticeRelaxationFrequency(), superGeometry, 4);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                        SuperGeometry2D<T>& superGeometry )
{


}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry2D<T>& superGeometry, Timer<T>& timer,
                 CircularBuffer<T>& buffer )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "channel2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );

  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int vtkIter  = converter.getLatticeTime( .5 );
  const int statIter = converter.getLatticeTime( .3 );
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

  }
  // Writes output on the console
  if ( iT%statIter == 0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},                        // resolution: number of voxels per charPhysL
    (T)   0.56,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1.,       // charPhysLength: reference length of simulation geometry
    (T)   1.0,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.0*hydraulicD/Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("channel2d");

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
