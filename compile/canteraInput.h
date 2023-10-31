using namespace std;

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>			// Geometry 
#include <stdio.h>      
#include <math.h>     
#include <cmath>
#include <memory>

// Include Boost libraries
#include <boost/config.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/lu.hpp> 
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric::ublas;

// Include Cantera files
#include <cantera/thermo.h>
#include <cantera/kinetics.h>
#include <cantera/transport.h>
#include <cantera/transport/DustyGasTransport.h>
#include <cantera/kinetics/GasKinetics.h>
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include <cantera/kinetics/ImplicitSurfChem.h>
#include "cantera/kinetics/solveSP.h"
#include <cantera/numerics/DenseMatrix.h>
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/Solution.h"

using namespace Cantera;
typedef std::vector<double> vector_double;

class inputClass
{
public:

    // Functions
    void readInputFile(std::string filename);
    void updateEosType(int);
    void getVarPointers(int);
    void resizeArrays(int);

    // Variables 

    // Mechanisms for the Lumen and the support structures:
    string m_filenameGas[2];				// Gas-phase Mechanism file name
    string m_filenameSurf[2];				// Surface Mechanism file name
    string m_gasPhase[2];					// Name of the gas phase
    string m_surfPhase[2];					// Name of the surface phase

    // Lumen and Support parameters
    int m_nPoints = -1;						// Total number of grid points	
    vector_double m_phi_gas; 				// Gas-phase porosity 
    vector_double m_tau;					// Gas-phase tortuosity	 
    vector_double m_As;						// Reactive specific surface area [1/m]
    vector_double m_Dp;                     // Particle diameter [m]
    vector_double m_Bg;						// Permeability
    vector_double m_rPore;                  // Pore radius
    vector_double m_phi_micro;              // Micro-porosities []
    vector_double m_rMicroPore;             // Micro-pore radius [m]

    // Support:
    int is_support = 0;                     // Flag to solve equations inside the support
    double m_thickness_2 = 0.0;             // Thickness of the support
    
    // Species related properties
    std::vector<int> m_nsp, m_nspSurf;
    std::vector<int> offset_Ygas, offset_theta, offset_T, offset_Ts, offset_U, offset_V, offset_P, m_Neq;
    int m_NeqTotal = 0;
    int nParts = 1;

    // Cantera objects
    ThermoPhase* m_gas[2] = {NULL, NULL};
    SurfPhase* m_surf[2] = {NULL, NULL};
    GasKinetics* m_kin[2] = {NULL, NULL};
    InterfaceKinetics* m_kin_surf[2] = {NULL, NULL};
    Transport* m_tran[2] = {NULL, NULL};
    DustyGasTransport* m_tranDGM[2] = {NULL, NULL};

    // Add flag for input mass or mole fraction
    int is_moleFrac = -1;                   // 0 for mole fractions, 1 for mass fractions

    // Flag to specify inlet mass flow rate or velocity
	int is_sccm = -1;

    // Annular channel flows
    int is_channel = 0;
    double m_rChannel;
    double u0_ch;

    // Boundary conditions
    vector_double m_yInlet = {};			    // Inlet mass fractions
    vector_double m_xInlet = {};			    // Inlet mole fractions
    vector_double m_thetaInlet = {};		    // Inlet surface coverages
    double m_outletPres = 0.0;					// Outlet pressure 
    double m_inletVel = 0.0;				    // Inlet velocity
    double m_inletTemp = 0.0;				    // Inlet temperature
    double m_Twall = 0.0;                       // Wall temperature
    double m_hcoeff;                            // Convective heat transfer coefficient (Wall)
    double flowRate_sccm = 0;                   // Inlet flowrate [SCCM]
	double m_Vel = 0;                           // Inlet velocity [m/s] 
    double m_inletFlux = 0;                     // Inlet mass flux [kg/m2-s]
    double m_Temp;
    int m_do_energy = 0;

    // Flag to solve radial transport
    int m_meshPoints = 0;						// Number of mesh points to solve Dusty gas transport model
    int m_lumenPoints = 1;                      // Radial mesh points inside the Lumen
    int m_supportPoints = 1;                    // Radial mesh points inside the support
    double m_axialLength = 0.0;					// Length of the bed
    int m_tranModel = 0;                        // 0 for Fickian diffusion, 1 for DGM
    int randPoreModel = 0;                      // 1 for random-pore model
    int is_radVel = 0;                          // 1 to include radial velocity

    // Flag to solve solid transport
    int m_solveSolid = 0;
    double m_condSolid, m_rhoSolid, m_CpSolid, m_emissivity, T_env, m_Aenv;
        
    //Membrane related parameters
    int m_solveMem = 0;                         // Flag to include membrane
    double m_memThickness;						// Membrane thickness
    double m_presMem;							// Permeate (sweep channel) pressure
    double m_perm0, m_Ea_R, m_perm;				// Membrane permeability parameters k0 and (Ea/R)
    double m_memArea2Vol;						// membrane area to volume ratio
    double m_Rin;								// Inner radius of the membrane tube
    double m_Rout;								// Outer radius of the membrane tube
    int m_memSpeciesIndex;						// Index of the permeable species through membrane 

    //Solver tolerances
    double m_Atol_ss, m_Rtol_ss, m_Atol_ts, m_Rtol_ts;
    int m_loglevel = 1;
    int m_refine = 0;
    int m_maxGridPoints = 100;
    int m_maxTimeSteps = 500;                   // Cantera default
    int m_minTimeSteps = 10;                    // Cantera default
    double m_minTimeStepSize = 1e-5;              // Cantera default
    int m_maxTimeStepSize = 1e+08;              // Cantera default

    // Grif refinement
    double max_grid_ratio, max_delta, prune, max_delta_slope;

    // Restart
    int m_restart = 0;

    // Initial guess from solution of diffusion-free code
    int m_initGuess = 0;
    string m_initGuessFile;					// Name of the initial guess file
};