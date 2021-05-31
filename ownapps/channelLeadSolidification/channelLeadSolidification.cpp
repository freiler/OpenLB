#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

#define NSDESCRIPTOR D2Q9<POROSITY,VELOCITY_SOLID,FORCE>
#define TDESCRIPTOR D2Q5<VELOCITY,TEMPERATURE>

//Hydrodynamic Sector
const T lengthX = 4.;
const T lengthY = 1.;
const T entrance_length = 0.25;
const T exit_length = 1.0;
int N = 100.;       // resolution of the model
const T physDeltaX = 1./N;      // latticeL
T       tau = 0.51;         // relaxation time
const T hydraulicD = 2.*lengthY;
const T Re = 100.;
const T rho = 11000.;
const T mu = 2.22e-3;
const T nu = mu/rho;
const T u_inlet = Re*nu/hydraulicD;
//const T nu = 1.0*hydraulicD/Re;

//const T Ste = 0.039;        // Stephan number
const T maxPhysT = 20000.;   // simulation time

//Thermal Sector
const T Tcold_real = 580; //in K
const T Tmelt_real = 600.6; //in K
const T Thot_real = 650.0; //in K
const T Tcold = 0.5; //in K
const T Tmelt = (Tmelt_real - Tcold_real)/(Thot_real - Tcold_real) + 0.5; //in K
const T Thot = 1.5; //in K

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
//const T L = cp_s * (Thot - Tmelt) / Ste; // J / kg
const T L = 1.5;

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry2D<T>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  //Vector<T,2> extend( lengthX,lengthY );
  //Vector<T,2> origin;

  superGeometry.rename( 0,2 );

  std::vector<T> extend(2,T());
  extend[0] = lengthX;
  extend[1] = lengthY;
  std::vector<T> origin(2,T());
  origin[0] = converter.getPhysLength(1);
  origin[1] = 0.5*converter.getPhysLength(1);
  IndicatorCuboid2D<T> bulk(extend, origin);
  superGeometry.rename(2,1,bulk); //changes from 4 to 1 only in region of cuboid2

  extend[0] = lengthX - entrance_length -  exit_length;
  extend[1] = physDeltaX;
  origin[0] = physDeltaX+entrance_length;
  origin[1] = 0.5*physDeltaX+lengthY;
  IndicatorCuboid2D<T> upper_wall( extend, origin );
  superGeometry.rename(2,3,1,upper_wall);

  extend[0] = physDeltaX;
  extend[1] = lengthY;
  origin[0] = 0.;
  origin[1] = 0.;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,4,1,inflow );

  extend[0] = physDeltaX;
  extend[1] = lengthY;
  origin[0] = lengthX+2.*physDeltaX;
  origin[1] = physDeltaX*0.5;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,5,1,outflow );

  extend[0] = lengthX - entrance_length  -  exit_length;
  extend[1] = physDeltaX;
  origin[0] = physDeltaX+entrance_length;
  origin[1] = -physDeltaX*0.5;
  IndicatorCuboid2D<T> lower_wall( extend, origin );
  superGeometry.rename( 2,6,1,lower_wall);

/*
  superGeometry.rename( 2,1,1,1 );

  extend[0] = lengthX-entrance_length+physDeltaX  -  exit_length;
  extend[1] = physDeltaX;
  origin[0] = entrance_length;
  origin[1] = lengthY;
  IndicatorCuboid2D<T> upper_wall( extend, origin );
  superGeometry.rename( 2,3,1,upper_wall);

  extend[0] = physDeltaX;
  extend[1] = lengthY;
  origin[0] = 0.;
  origin[1] = -physDeltaX;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,4,1,inflow );

  // Set material number for
  extend[0] = physDeltaX*1.1;
  extend[1] = lengthY;
  origin[0] = lengthX-physDeltaX*0.1;
  origin[1] = 0.;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,5,1,outflow );

  extend[0] = lengthX-entrance_length+physDeltaX  -  exit_length;
  extend[1] = physDeltaX;
  origin[0] = entrance_length;
  origin[1] = 0.;
  IndicatorCuboid2D<T> lower_wall( extend, origin );
  superGeometry.rename( 2,6,1,lower_wall);
*/

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.innerClean();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;

}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                     ForcedBGKdynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     SuperGeometry2D<T>& superGeometry )
{

    OstreamManager clout(std::cout,"prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    T omega  =  converter.getLatticeRelaxationFrequency();
    T Tomega  =  converter.getLatticeThermalRelaxationFrequency();

    // Hydrodynamic begins here
    NSlattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>() );

    //auto bulkIndicator = superGeometry.getMaterialIndicator({1, 4, 5});
    auto bulkIndicator = superGeometry.getMaterialIndicator({1, 4, 5});
    NSlattice.defineDynamics( bulkIndicator, &bulkDynamics );

    // Setting of the boundary conditions
    NSlattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, NSDESCRIPTOR>() );
    NSlattice.defineDynamics( superGeometry, 3, &instances::getBounceBack<T, NSDESCRIPTOR>() );
    NSlattice.defineDynamics( superGeometry, 6, &instances::getBounceBack<T, NSDESCRIPTOR>() );
    //NSlattice.defineDynamics( superGeometry, 7, &instances::getBounceBack<T, NSDESCRIPTOR>() );
    //NSlattice.defineDynamics( superGeometry, 8, &instances::getBounceBack<T, NSDESCRIPTOR>() );

    //if boundary conditions are chosen to be local
    setInterpolatedVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry,4);
    setLocalPressureBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry, 5);

    // Initial conditions
    AnalyticalConst2D<T,T> rhoF( 1. );
    std::vector<T> velocity( 2,T( 0 ) );
    AnalyticalConst2D<T,T> uzero( velocity );

    velocity[0] = converter.getCharLatticeVelocity();
    AnalyticalConst2D<T,T> u( velocity );

    NSlattice.iniEquilibrium( bulkIndicator, rhoF, u );
    NSlattice.defineRhoU( bulkIndicator, rhoF, u );
    NSlattice.defineRhoU( superGeometry.getMaterialIndicator({2, 3,6 }), rhoF, uzero );
    NSlattice.iniEquilibrium( superGeometry.getMaterialIndicator({2, 3,6}), rhoF, uzero );

    ////Thermal Begins here
    ADlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, TDESCRIPTOR>());

    ADlattice.defineDynamics(superGeometry.getMaterialIndicator({1,3,4,6}), &advectionDiffusionBulkDynamics);
    //ADlattice.defineDynamics(superGeometry,1, &advectionDiffusionBulkDynamics);
    //Adiabatic Outlet
    //ADlattice.defineDynamics(superGeometry.getMaterialIndicator({2,3, 4,5,6}), &instances::getBounceBack<T, TDESCRIPTOR>());
    ADlattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, TDESCRIPTOR>());
    ADlattice.defineDynamics(superGeometry, 5, &instances::getBounceBack<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 2, &advectionDiffusionBulkDynamics);

    //setRegularizedHeatFluxBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 2);
    //setRegularizedHeatFluxBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 5);
    setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 4);
    setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 3);
    setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 6);

    //ADlattice.defineDynamics(superGeometry, 6, &instances::getBounceBack<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 7, &instances::getBounceBack<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 8, &instances::getBounceBack<T, TDESCRIPTOR>());
    /// sets boundary
    //setRegularizedTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry.getMaterialIndicator({3, 4, 6}));

    /// define initial conditions
    AnalyticalConst2D<T,T> rho(1.);
    AnalyticalConst2D<T,T> u0(0.0, 0.0);
    AnalyticalConst2D<T,T> T_cold(Tcold);
    //AnalyticalConst2D<T,T> T_hot(lattice_Hhot+L);
    AnalyticalConst2D<T,T> T_hot_BC(Thot);
    AnalyticalConst2D<T,T> T_dummy_BC(1.3);
    AnalyticalConst2D<T,T> H_dummy_bulk(1.+L);
    AnalyticalConst2D<T,T> T_hot(lattice_Hhot+L);
    //AnalyticalConst2D<T,T> T_hot2(lattice_Hhot+10.);

    //ADlattice.defineRho(superGeometry.getMaterialIndicator({1,4}), T_hot);
    ADlattice.defineRho(superGeometry, 4 , T_hot_BC);
    ADlattice.defineRho(superGeometry.getMaterialIndicator({3, 6}), T_cold);
    //ADlattice.iniEquilibrium(superGeometry.getMaterialIndicator({4,5}), T_hot_BC, u);//????
    ADlattice.iniEquilibrium(superGeometry, 4 , H_dummy_bulk, u); //inlet
    ADlattice.iniEquilibrium(superGeometry, 5 , H_dummy_bulk, u); // outlet
    ADlattice.iniEquilibrium(superGeometry, 2 , H_dummy_bulk, u0); // adiabatic walls
    ADlattice.iniEquilibrium(superGeometry.getMaterialIndicator({3, 6}), H_dummy_bulk, u0); //cold walls
    ADlattice.iniEquilibrium(superGeometry, 1 , H_dummy_bulk, u); //bulk

    /// Make the lattice ready for simulation
    NSlattice.initialize();
    ADlattice.initialize();

    clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                        SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                        SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                        int iT, SuperGeometry2D<T>& superGeometry)
{

// nothing to do here

}

void getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                 SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice2D<T, TDESCRIPTOR>& ADlattice, int iT,
                 SuperGeometry2D<T>& superGeometry,
                 Timer<T>& timer,
                 bool converged)
{

    OstreamManager clout(std::cout,"getResults");

    SuperVTMwriter2D<T> vtkWriter("leadsolidification");
    SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperLatticeField2D<T, TDESCRIPTOR, VELOCITY> velocity(ADlattice);
    SuperLatticePhysPressure2D<T, NSDESCRIPTOR> pressure(NSlattice, converter);

    SuperLatticeDensity2D<T, TDESCRIPTOR> enthalpy(ADlattice);
    enthalpy.getName() = "enthalpy";
    SuperLatticeField2D<T, NSDESCRIPTOR, POROSITY> liquid_frac(NSlattice);
    liquid_frac.getName() = "liquid fraction";
    SuperLatticeField2D<T, TDESCRIPTOR, TEMPERATURE> temperature(ADlattice);
    temperature.getName() = "temperature";
    SuperLatticeField2D<T, NSDESCRIPTOR, FORCE> force(NSlattice);
    force.getName() = "force";
    vtkWriter.addFunctor( geometry );
    vtkWriter.addFunctor( pressure );
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( enthalpy );
    vtkWriter.addFunctor( liquid_frac );
    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( force );

    const int vtkIter = converter.getLatticeTime(500.);

    if (iT == 0) {
        /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
        SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
        SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
        SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
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

        /// ADLattice statistics console output
        ADlattice.getStatistics().print(iT,converter.getPhysTime(iT));

        if ( NSlattice.getStatistics().getAverageRho() != NSlattice.getStatistics().getAverageRho() or ADlattice.getStatistics().getAverageRho() != ADlattice.getStatistics().getAverageRho() ) {
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
  singleton::directories().setOutputDir("./tmp/");

  clout << "L is equal to :" << L << std::endl;

  const T physDeltaT = 1 / nu / descriptors::invCs2<T,NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX ;

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const converter(
      (T) physDeltaX, // physDeltaX
      (T) physDeltaT, // physDeltaT
      (T) hydraulicD, // charPhysLength
      (T) u_inlet, // charPhysVelocity
      (T) nu, // physViscosity
      (T) rho, // physDensity
      (T) lambda_l, // physThermalConductivity
      (T) 138.8, // physSpecificHeatCapacity
      (T) 1.0e-4, // physThermalExpansionCoefficient
      (T) Tcold, // charPhysLowTemperature
      (T) Thot // charPhysHighTemperature
  );
  converter.print();
  clout << "lattice cp " << converter.getLatticeSpecificHeatCapacity(cp_l) << endl;
  // Writes the converter log in a file
  //converter.write("channelLeadSolidification");

  // === 2rd Step: Prepare Geometry ===
  /*
  Vector<T,2> extend( lengthX,lengthY );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );
  */
  std::vector<T> extend(2,T());
  extend[0] = lengthX + 2*converter.getPhysLength(1);
  extend[1] = lengthY + converter.getPhysLength(1);
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize() );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( superGeometry, converter );


    /// === 3rd Step: Prepare Lattice ===

    SuperLattice2D<T, TDESCRIPTOR> ADlattice(superGeometry);
    SuperLattice2D<T, NSDESCRIPTOR> NSlattice(superGeometry);

    ForcedPSMBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
        converter.getLatticeRelaxationFrequency(),
        instances::getBulkMomenta<T,NSDESCRIPTOR>());

    TotalEnthalpyAdvectionDiffusionTRTdynamics<T, TDESCRIPTOR> TbulkDynamics (
        converter.getLatticeThermalRelaxationFrequency(),
        instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>(),
        1./4.,
        Tmelt,
        Tmelt,
        cp_s,
        cp_l,
        cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5) * R_lambda,
        cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5),
        L
    );

    std::vector<T> dir{0.0, 1.0};
    //T boussinesqForcePrefactor = Ra / pow(T(N),3) * Pr * pow(cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5), 2);
    //clout << "boussinesq " << Ra / pow(T(N), 3) * Pr * lambda_l * lambda_l << endl;

    TotalEnthalpyPhaseChangeCouplingGenerator2D<T,NSDESCRIPTOR>
    coupling(0, converter.getLatticeLength(lengthX), 0, converter.getLatticeLength(lengthY),
             0.0, Tcold, 1., dir);

    NSlattice.addLatticeCoupling(superGeometry, 1, coupling, ADlattice);

    //prepareLattice and setBoundaryConditions
		 prepareLattice(converter,
										NSlattice, ADlattice,
										NSbulkDynamics, TbulkDynamics,
										superGeometry );

    /// === 4th Step: Main Loop with Timer ===
    Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
    timer.start();

    for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT)+1; ++iT) {

        /// === 5th Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

        /// === 6th Step: Collide and Stream Execution ===
        NSlattice.executeCoupling();
        NSlattice.collideAndStream();
        ADlattice.collideAndStream();

        /// === 7th Step: Computation and Output of the Results ===
        getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, false);
    }


  timer.stop();
  timer.printSummary();
}
