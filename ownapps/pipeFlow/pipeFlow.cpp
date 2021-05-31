#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

typedef D3Q19<> NSDESCRIPTOR;

//Hydrodynamic Sector
const T length = 0.5;
const T diameter = 0.2;
const T radius = diameter/2.;
const T entrance_length = 0.05;
const T exit_length = 0.05;
int N = 50.;       // resolution of the model
const T physDeltaX = diameter/N;      // latticeL
T       tau = 0.51;         // relaxation time
const T hydraulicD = diameter;
const T Re = 500.;
const T rho = 11000.;
const T mu = 2.22e-3;
const T nu = mu/rho;
const T u_inlet = Re*nu/hydraulicD;
//const T nu = 1.0*hydraulicD/Re;

//const T Ste = 0.039;        // Stephan number
const T maxPhysT = 500.;   // simulation time
const T output_frequency = 10.2;

//Thermal Sector
const T Tcold_real = 600.59; //in K
const T Tmelt_real = 600.6; //in K
const T Thot_real = 640.0; //in K
const T Tcold = 0.5; //in K
const T Tmelt = (Tmelt_real - Tcold_real)/(Thot_real - Tcold_real) + 0.5; //in K
const T Thot = 1.5; //in K

//const T lambda_s = 30.; // W / m K
const T lambda_s = 30.; // W / m K
const T lambda_l = 16.6; // W / m K
const T R_lambda =  lambda_s/lambda_l;

const T cp_real = 138.8; // J / kg K
//const T cp_l = 138.8; // J / kg K
const T cp_s = 1.; // J / kg K
const T cp_l = 1.; // J / kg K
//onst T R_cp = cp_s/cp_l;

T lattice_Hcold = cp_s * Tcold;
T lattice_Hhot = cp_l * Thot;

//for this case, the harmonic mean (cp_ref) is applicable
const T cp_ref = 2.0 * cp_s * cp_l / (cp_s + cp_l); // J / kg K

const T L_real = 23648.6; // J/kg
const T Ste = cp_real*(Thot_real-Tmelt_real)/L_real;
const T L = cp_s * (Thot - Tcold) / Ste; // J / kg
//const T L = 1.5;

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry3D<T>& superGeometry,
                     UnitConverter<T, NSDESCRIPTOR> const& converter)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

/*
  //Vector<T,2> extend( lengthX,lengthY );
  //Vector<T,2> origin;
  std::vector<T> center1(3,T());
  std::vector<T> center2(3,T());
  center2[0] = length;
  IndicatorCylinder3D<T> pipe(center1,center2,radius);

  superGeometry.rename( 0,2,pipe );
  superGeometry.rename( 2,1,1,1,1 );

  center1[0] = 0.;
  center2[0] = 2*physDeltaX;
  IndicatorCylinder3D<T> inflow (center1,center2,radius);
  superGeometry.rename( 2,3,1, inflow );

  center1[0] = -physDeltaX+length;
  center2[0] = physDeltaX+length;
  IndicatorCylinder3D<T> outflow (center1,center2,radius);
  superGeometry.rename( 2,4,1, outflow );
*/


  //Vector<T,2> extend( lengthX,lengthY );
  //Vector<T,2> origin;
  std::vector<T> center1(3,T());
  std::vector<T> center2(3,T());
  center2[0] = length;
  IndicatorCylinder3D<T> pipe (center1,center2,radius+converter.getPhysDeltaX());

  superGeometry.rename( 0,2,pipe );
  superGeometry.rename( 2,1,1,1,1 );

  center1[0] = 0.;
  center2[0] = 2*physDeltaX;
  IndicatorCylinder3D<T> inflow (center1,center2,radius);
  superGeometry.rename( 2,3,1, inflow );

  center1[0] = -physDeltaX+length;
  center2[0] = physDeltaX+length;
  IndicatorCylinder3D<T> outflow (center1,center2,radius);
  superGeometry.rename( 2,4,1, outflow );

  std::vector<T> extend(3,T());
  extend[0] = length-exit_length-entrance_length;
  extend[1] = diameter+physDeltaX;
  extend[2] = diameter+physDeltaX;
  std::vector<T> origin(3,T());
  origin[0] = entrance_length;
  origin[1] = -radius - converter.getPhysLength(1);
  origin[2] = -radius - converter.getPhysLength(1);
  IndicatorCuboid3D<T> cooledsection(extend, origin);
  superGeometry.rename(2,5,cooledsection); //changes from 4 to 1 only in region of cuboid2

  //center1[0] = entrance_length;
  //center2[0] = length-exit_length;
  //const T radius1 = radius+physDeltaX;
  //IndicatorCylinder3D<T> cooledsection (center1,center2,radius1);
  //superGeometry.rename( 2,5,1, cooledsection );


  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.innerClean();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;

}


void prepareLattice( UnitConverter<T, NSDESCRIPTOR> const& converter,
                     SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                     Dynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     SuperGeometry3D<T>& superGeometry )
{

    OstreamManager clout(std::cout,"prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    T omega  =  converter.getLatticeRelaxationFrequency();

    // Hydrodynamic begins here
    NSlattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>() );

    //auto bulkIndicator = superGeometry.getMaterialIndicator({1, 4, 5});
    auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});
    NSlattice.defineDynamics( bulkIndicator, &bulkDynamics );


    // Regular bounce back
    //NSlattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, NSDESCRIPTOR>() );
    //NSlattice.defineDynamics( superGeometry, 5, &instances::getBounceBack<T, NSDESCRIPTOR>() );

    setInterpolatedVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry,3);
    //setInterpolatedPressureBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry, 4);


    // Bouzidi Block
    NSlattice.defineDynamics( superGeometry,2,&instances::getNoDynamics<T,NSDESCRIPTOR>() );
    NSlattice.defineDynamics( superGeometry,5,&instances::getNoDynamics<T,NSDESCRIPTOR>() );

    std::vector<T> center1(3,T());
    std::vector<T> center2(3,T());
    center1[0] = 0.;
    center2[0] = length;
    IndicatorCylinder3D<T> pipe (center1,center2,radius);
    setBouzidiZeroVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, superGeometry, 2, pipe);
    setBouzidiZeroVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, superGeometry, 5, pipe);

    //NSlattice.defineDynamics( superGeometry,3,&instances::getNoDynamics<T,NSDESCRIPTOR>() );
    //setBouzidiVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, superGeometry, 3, pipe);


    setInterpolatedPressureBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry, 4);

    // Initial conditions
    AnalyticalConst3D<T,T> rhoF( 1. );
    std::vector<T> velocity( 3,T( 0 ) );
    AnalyticalConst3D<T,T> uzero( velocity );

    //std::vector<T> origin = { length, 0., 0.};
    //std::vector<T> axis = { 1, 0, 0 };
    //CirclePoiseuille3D<T> u(origin, axis, converter.getCharLatticeVelocity(), radius);

    velocity[0] = converter.getCharLatticeVelocity();
    AnalyticalConst3D<T,T> u( velocity );

    NSlattice.iniEquilibrium( bulkIndicator, rhoF, u );
    NSlattice.defineRhoU( bulkIndicator, rhoF, u );
    //NSlattice.defineRhoU( superGeometry,2 , rhoF, uzero );
    //NSlattice.iniEquilibrium( superGeometry,2 , rhoF, uzero );
    //NSlattice.defineRhoU( superGeometry,5 , rhoF, uzero );
    //NSlattice.iniEquilibrium( superGeometry,5 , rhoF, uzero );

    /// Make the lattice ready for simulation
    NSlattice.initialize();

    clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( UnitConverter<T, NSDESCRIPTOR> const& converter,
                        SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                        int iT, SuperGeometry3D<T>& superGeometry)
{

}

void getResults( UnitConverter<T, NSDESCRIPTOR> const& converter,
                 SuperLattice3D<T, NSDESCRIPTOR>& NSlattice, int iT,
                 SuperGeometry3D<T>& superGeometry,
                 Timer<T>& timer,
                 bool converged)
{

    OstreamManager clout(std::cout,"getResults");

    SuperVTMwriter3D<T> vtkWriter("pipeFlow");
    SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( NSlattice, converter );
    SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure(NSlattice, converter);

    vtkWriter.addFunctor( geometry );
    vtkWriter.addFunctor( pressure );
    vtkWriter.addFunctor( velocity );

    //const int vtkIter = converter.getLatticeTime(output_frequency);
    const int vtkIter = int (output_frequency/converter.getPhysDeltaT() );
    if (iT == 0) {
        /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
        SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
        SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSlattice);
        SuperLatticeRank3D<T, NSDESCRIPTOR> rank(NSlattice);
        vtkWriter.write(geometry);
        vtkWriter.write(cuboid);
        vtkWriter.write(rank);

        vtkWriter.createMasterFile();
    }

    /// Writes the VTK files
    if (iT % vtkIter == 0 || converged) {

        timer.update(iT);
        timer.printStep();

        /// NSLattice statistics console output
        NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));

        if ( NSlattice.getStatistics().getAverageRho() != NSlattice.getStatistics().getAverageRho()  ) {
            clout << "simulation diverged! stopping now." << endl;
            exit(1);
        }
        vtkWriter.write(iT);
    }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp3/");

  clout << "L is equal to :" << L << std::endl;

  const T physDeltaT = 1 / nu / descriptors::invCs2<T,NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX ;

  UnitConverterFromResolutionAndRelaxationTime<T, NSDESCRIPTOR> const converter(
    int {N},                        // resolution: number of voxels per charPhysL
    (T)   tau,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   diameter,       // charPhysLength: reference length of simulation geometry
    (T)   u_inlet,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   nu, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   rho                       // physDensity: physical density in __kg / m^3__
  );

  converter.print();

  std::vector<T> center1(3,T());
  std::vector<T> center2(3,T());
  center2[0] = length;
  IndicatorCylinder3D<T> pipe(center1,center2,radius+converter.getPhysDeltaX());

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry3D<T> cuboidGeometry( pipe, converter.getPhysDeltaX(), singleton::mpi().getSize() );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( superGeometry, converter );

    /// === 3rd Step: Prepare Lattice ===

    SuperLattice3D<T, NSDESCRIPTOR> NSlattice(superGeometry);

    //ForcedPSMBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    BGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
        converter.getLatticeRelaxationFrequency(),
        instances::getBulkMomenta<T,NSDESCRIPTOR>());

    //prepareLattice and setBoundaryConditions
		 prepareLattice(converter,
										NSlattice,
										NSbulkDynamics,
										superGeometry );

    /// === 4th Step: Main Loop with Timer ===
    Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
    timer.start();

    for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT)+1; ++iT) {

        /// === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(converter, NSlattice, iT, superGeometry);

        /// === 6th Step: Collide and Stream Execution ===
        NSlattice.collideAndStream();

        /// === 7th Step: Computation and Output of the Results ===
        getResults(converter, NSlattice, iT, superGeometry, timer, false);
    }


  timer.stop();
  timer.printSummary();

}
