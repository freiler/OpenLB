#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

//#define NSDESCRIPTOR D3Q19<POROSITY,VELOCITY_SOLID,FORCE>
//#define TDESCRIPTOR D3Q7<VELOCITY,TEMPERATURE>
typedef D3Q7<VELOCITY> TDESCRIPTOR;
typedef D3Q19<FORCE> NSDESCRIPTOR;

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
const T output_frequency = 40.2;

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
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
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
  IndicatorCylinder3D<T> pipe(center1,center2,radius+converter.getPhysDeltaX());

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


void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice3D<T, TDESCRIPTOR>& ADlattice,
                     ForcedBGKdynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     SuperGeometry3D<T>& superGeometry )
{

    OstreamManager clout(std::cout,"prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    T omega  =  converter.getLatticeRelaxationFrequency();
    T Tomega  =  converter.getLatticeThermalRelaxationFrequency();

    // Hydrodynamic begins here
    NSlattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>() );

    //auto bulkIndicator = superGeometry.getMaterialIndicator({1, 4, 5});
    auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});
    NSlattice.defineDynamics( bulkIndicator, &bulkDynamics );

    // Setting of the boundary conditions
    NSlattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, NSDESCRIPTOR>() );
    NSlattice.defineDynamics( superGeometry, 5, &instances::getBounceBack<T, NSDESCRIPTOR>() );

    //if boundary conditions are chosen to be local
    setInterpolatedVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry,3);
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
    NSlattice.defineRhoU( superGeometry,2 , rhoF, uzero );
    NSlattice.iniEquilibrium( superGeometry,2 , rhoF, uzero );
    NSlattice.defineRhoU( superGeometry,5 , rhoF, uzero );
    NSlattice.iniEquilibrium( superGeometry,5 , rhoF, uzero );

    ////Thermal Begins here
    ADlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>());
    ADlattice.defineDynamics(superGeometry, 4, &instances::getNoDynamics<T, TDESCRIPTOR>());

    ADlattice.defineDynamics(superGeometry.getMaterialIndicator({1,3,5}), &advectionDiffusionBulkDynamics);
    //ADlattice.defineDynamics(superGeometry,1, &advectionDiffusionBulkDynamics);
    //Adiabatic Outlet
    //ADlattice.defineDynamics(superGeometry.getMaterialIndicator({2,3, 4,5,6}), &instances::getBounceBack<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 4, &instances::getBounceBack<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 5, &instances::getBounceBack<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 2, &advectionDiffusionBulkDynamics);

    //setRegularizedHeatFluxBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 4);
    //setRegularizedHeatFluxBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 5);
    //setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 2);
    setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 5);
    setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry, 3);
    setAdvectionDiffusionConvectionBoundary<T,TDESCRIPTOR>(ADlattice, superGeometry, 4);
    ADlattice.defineDynamics(superGeometry,2, &instances::getBounceBack<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 7, &instances::getBounceBack<T, TDESCRIPTOR>());
    //ADlattice.defineDynamics(superGeometry, 8, &instances::getBounceBack<T, TDESCRIPTOR>());
    /// sets boundary
    //setRegularizedTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry.getMaterialIndicator({3, 4, 6}));

    /// define initial conditions
    AnalyticalConst3D<T,T> T_cold(Tcold);
    //AnalyticalConst2D<T,T> T_hot(lattice_Hhot+L);
    AnalyticalConst3D<T,T> T_hot(Thot);
    AnalyticalConst3D<T,T> H_hot(lattice_Hhot+L);
    AnalyticalConst2D<T,T> H_dummy(2.0);

    //ADlattice.defineRho(superGeometry.getMaterialIndicator({1,4}), T_hot);
    //ADlattice.defineRho(superGeometry, 2 , T_hot);
    ADlattice.defineRho(superGeometry, 3 , T_hot);
    ADlattice.defineRho(superGeometry, 5 , T_cold);
    //ADlattice.iniEquilibrium(superGeometry.getMaterialIndicator({4,5}), T_hot_BC, u);//????
    ADlattice.iniEquilibrium(superGeometry, 1 , T_hot, u);
    ADlattice.iniEquilibrium(superGeometry, 5 , T_cold, uzero);
    ADlattice.iniEquilibrium(superGeometry, 3 , T_hot, u);
    ADlattice.iniEquilibrium(superGeometry, 4 , T_hot, u);
    //ADlattice.iniEquilibrium(superGeometry, 2 , T_hot, uzero);

    /// Make the lattice ready for simulation
    NSlattice.initialize();
    ADlattice.initialize();

    clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                        SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                        SuperLattice3D<T, TDESCRIPTOR>& ADlattice,
                        int iT, SuperGeometry3D<T>& superGeometry)
{
/*
  int iTupdate = 100;

  if ( iT%iTupdate == 0 ) {
    AnalyticalConst3D<T,T> T_hot(Thot);
    ADlattice.defineRho(superGeometry, 3 , T_hot);
  }
*/
}

void getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                 SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice3D<T, TDESCRIPTOR>& ADlattice, int iT,
                 SuperGeometry3D<T>& superGeometry,
                 Timer<T>& timer,
                 bool converged)
{

    OstreamManager clout(std::cout,"getResults");

    SuperVTMwriter3D<T> vtkWriter("solidification3D");
    SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( NSlattice, converter );
    SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
    SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);

    SuperLatticeField3D<T, TDESCRIPTOR, VELOCITY> velocity_coupled(ADlattice);
    velocity_coupled.getName() = "velocity_coupled";
    //SuperLatticeDensity3D<T, TDESCRIPTOR> enthalpy(ADlattice);
    //enthalpy.getName() = "enthalpy";
    //SuperLatticeField3D<T, NSDESCRIPTOR, POROSITY> liquid_frac(NSlattice);
    //liquid_frac.getName() = "liquid fraction";
    //SuperLatticeField3D<T, TDESCRIPTOR, TEMPERATURE> temperature(ADlattice);
    //temperature.getName() = "temperature";
    //SuperLatticeField3D<T, NSDESCRIPTOR, FORCE> force(NSlattice);
    //force.getName() = "force";
    vtkWriter.addFunctor( geometry );
    vtkWriter.addFunctor( pressure );
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( velocity_coupled );
    //vtkWriter.addFunctor( enthalpy );
    //vtkWriter.addFunctor( liquid_frac );
    vtkWriter.addFunctor( temperature );
    //vtkWriter.addFunctor( force );

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
  singleton::directories().setOutputDir("./tmp2/");

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
      (T) cp_real, // physSpecificHeatCapacity
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

    SuperLattice3D<T, TDESCRIPTOR> ADlattice(superGeometry);
    SuperLattice3D<T, NSDESCRIPTOR> NSlattice(superGeometry);

    //ForcedPSMBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    ForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
        converter.getLatticeRelaxationFrequency(),
        instances::getBulkMomenta<T,NSDESCRIPTOR>()
    );

    AdvectionDiffusionBGKdynamics<T, TDESCRIPTOR> TbulkDynamics (
      converter.getLatticeThermalRelaxationFrequency(),
      instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>()
    );

/*
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
*/
    std::vector<T> dir{0.0, 1.0, 0.0};
    //T boussinesqForcePrefactor = Ra / pow(T(N),3) * Pr * pow(cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5), 2);
    //clout << "boussinesq " << Ra / pow(T(N), 3) * Pr * lambda_l * lambda_l << endl;
/*
    TotalEnthalpyPhaseChangeCouplingGenerator3D<T,NSDESCRIPTOR>
    coupling(0, converter.getLatticeLength(length),
             -1 *converter.getLatticeLength(radius),converter.getLatticeLength(radius),
             -1 *converter.getLatticeLength(radius),converter.getLatticeLength(radius), 0.0, Tcold, 1., dir);
    //TotalEnthalpyPhaseChangeCouplingGenerator3D<T,NSDESCRIPTOR> coupling(0, converter.getLatticeLength(length),-1 * converter.getLatticeLength(radius), converter.getLatticeLength(radius),-1 * converter.getLatticeLength(radius), converter.getLatticeLength(radius),0.0, Tcold, 1., dir);
*/

    NavierStokesAdvectionDiffusionCouplingGenerator3D<T,NSDESCRIPTOR>
    coupling(0, converter.getLatticeLength(length),
             -1 *converter.getLatticeLength(radius),converter.getLatticeLength(radius),
             -1 *converter.getLatticeLength(radius),converter.getLatticeLength(radius), 0.0, Tcold, 1., dir);

    //NSlattice.addLatticeCoupling(superGeometry, 1, coupling, ADlattice);
    NSlattice.addLatticeCoupling(coupling, ADlattice);
    //SuperVTMwriter3D<T> vtkWriter("solidification3d");
    //SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    //vtkWriter.addFunctor( geometry );
    //vtkWriter.write(geometry);

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
