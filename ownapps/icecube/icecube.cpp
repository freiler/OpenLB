#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

#define NSDESCRIPTOR D2Q9<POROSITY,VELOCITY_SOLID,FORCE>
#define TDESCRIPTOR D2Q5<VELOCITY,TEMPERATURE>

// Parameters for the simulation setup
const T lx  = 0.06;      // length of the channel
const T ly  = 0.12;      // height of the channel
const T l_icecube  = 0.03;      // ice cube side
int     N = 80;            // resolution of the model
T       tau = 0.51;         // relaxation time
const T maxPhysT = 250.;   // simulation time

const T Tcold = 0.5;
const T Tmelt = 0.65;
const T Thot = 1.5;

const T lambda_s = 5.56; // W / m K Water at 0.01 C 1 bar
const T lambda_l = 5.79; // W / m K Water at 10 C 1 bar
const T R_lambda =  lambda_s/lambda_l;

const T mu = 1.52e-3;
const T rho = 1000.;
const T beta = 4.0e-5;
const T cp = 4184.;
const T nu = mu/rho;
const T alpha = lambda_l/(rho*cp);

const T cp_s = 1.0; // J / kg K
const T cp_l = 1.0; // J / kg K
const T R_cp = cp_s/cp_l;

//for this case, the harmonic mean (cp_ref) is applicable
const T cp_ref = 2.0 * cp_s * cp_l / (cp_s + cp_l); // J / kg K

//const T Ra = 9.81*beta*4.*ly*ly*ly/nu/alpha;           // Rayleigh number
const T Ra =10000.;
const T Pr = nu/alpha;        // Prandtl number
const T Ste = 0.12;        // Stefan number calculated as cpl*(Thot-Tmelt)/L = 4.184*(10-0)/334.4

const T R_alpha = lambda_s / lambda_l * cp_l / cp_s;
const T L = cp_l * (Thot - Tmelt) / Ste; // J / kg

T lattice_Hcold, lattice_Hhot;
T physDeltaX, physDeltaT;

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry2D<T>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
{

    OstreamManager clout(std::cout,"prepareGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    superGeometry.rename(0,2);

    std::vector<T> extend(2,T());
    extend[0] = lx;
    extend[1] = ly;
    std::vector<T> origin(2,T());
    IndicatorCuboid2D<T> glassIndicator(extend, origin);

    superGeometry.rename(2,1,glassIndicator);

    /// Removes all not needed boundary voxels outside the surface
    superGeometry.clean();
    /// Removes all not needed boundary voxels inside the surface
    superGeometry.innerClean();
    superGeometry.checkForErrors();

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

    ADlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>());
    NSlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>());

    ADlattice.defineDynamics(superGeometry, 1, &advectionDiffusionBulkDynamics);
    ADlattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, TDESCRIPTOR>());

    NSlattice.defineDynamics(superGeometry.getMaterialIndicator({1, 2}), &bulkDynamics);

    /// sets boundary
    //setRegularizedTemperatureBoundary<T,TDESCRIPTOR>(ADlattice, Tomega, superGeometry.getMaterialIndicator({2, 3}));
    setInterpolatedVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry,2);

    /// define initial conditions
    AnalyticalConst2D<T,T> rho(1.);
    AnalyticalConst2D<T,T> u0(0.0, 0.0);

    /// for each material set Rho, U and the Equilibrium
    NSlattice.defineRhoU(superGeometry.getMaterialIndicator({1, 2}), rho, u0);
    NSlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2}), rho, u0);

    //IceCube
    Vector<T,2> cubeCenter = {lx/2.,ly-0.5*l_icecube};

    SmoothIndicatorCuboid2D<T,T> icecube(cubeCenter, l_icecube, l_icecube, converter.getPhysLength(1));

    AnalyticalConst2D<T,T> H_cold(lattice_Hcold);
    //AnalyticalConst2D<T,T> H_hot(lattice_Hhot+L);
    AnalyticalConst2D<T,T> H_hot(lattice_Hhot+1.2*L);
    //AnalyticalConst2D<T,T> H_L(L);
    AnalyticalIdentity2D<T,T> H_0(H_hot - icecube*(H_hot-H_cold));

    ADlattice.defineRhoU(superGeometry,1 , H_0, u0);
    ADlattice.iniEquilibrium(superGeometry,1 , H_0, u0);

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

    SuperVTMwriter2D<T> vtkWriter("thermalNaturalConvection2D");
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

    const int vtkIter = converter.getLatticeTime(0.5);

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

int main(int argc, char *argv[])
{

    /// === 1st Step: Initialization ===
    OstreamManager clout(std::cout,"main");
    olbInit(&argc, &argv);
    singleton::directories().setOutputDir("./tmp/");

    T char_lattice_u = 0.2;

    const T char_u = sqrt( 9.81 * beta * (4. - 0.) * rho );
    const T conversion_u = char_u / char_lattice_u;

    physDeltaX = lx / N;
    physDeltaT = physDeltaX / conversion_u;
    physDeltaT = rho / mu / descriptors::invCs2<T,NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX;

    lattice_Hcold = cp_s * Tcold;
    lattice_Hhot = cp_l * Thot;

    clout << "H_cold " << lattice_Hcold << " H_hot " << lattice_Hhot << endl;

    ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const converter(
        (T) physDeltaX, // physDeltaX
        (T) physDeltaT, // physDeltaT
        (T) lx, // charPhysLength
        (T) char_u, // charPhysVelocity
        (T) nu, // physViscosity
        (T) rho, // physDensity
        (T) lambda_l, // physThermalConductivity
        (T) cp, // physSpecificHeatCapacity
        (T) beta, // physThermalExpansionCoefficient
        (T) Tcold, // charPhysLowTemperature
        (T) Thot // charPhysHighTemperature
    );
    converter.print();
    clout << "lattice cp " << converter.getLatticeSpecificHeatCapacity(cp_l) << endl;

    /// === 2nd Step: Prepare Geometry ===
    std::vector<T> extend(2,T());
    extend[0] = lx + converter.getPhysLength(1);
    extend[1] = ly + converter.getPhysLength(1);
    std::vector<T> origin(2,T());
    origin[0] = -0.5*converter.getPhysLength(1);
    origin[1] = -0.5*converter.getPhysLength(1);
    IndicatorCuboid2D<T> cuboid(extend, origin);

    /// Instantiation of an empty cuboidGeometry
    CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());

    /// Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    /// Instantiation of a superGeometry
    SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

    prepareGeometry(superGeometry, converter);


    /// === 3rd Step: Prepare Lattice ===

    SuperLattice2D<T, TDESCRIPTOR> ADlattice(superGeometry);
    SuperLattice2D<T, NSDESCRIPTOR> NSlattice(superGeometry);

    //Check geometry
    SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperVTMwriter2D<T> vtkWriter("iceCube2D");
    vtkWriter.addFunctor( geometry );
    vtkWriter.write(geometry);

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
    T boussinesqForcePrefactor = Ra / pow(T(N),3) * Pr * pow(cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5), 2);
    clout << "boussinesq " << Ra / pow(T(N), 3) * Pr * lambda_l * lambda_l << endl;

    TotalEnthalpyPhaseChangeCouplingGenerator2D<T,NSDESCRIPTOR>
    coupling(0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(ly),
             boussinesqForcePrefactor, Thot, 1., dir);

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
    //for (std::size_t iT = 0; iT < 1; ++iT) {
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
