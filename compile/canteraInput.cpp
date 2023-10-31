#include "canteraInput.h"

/*
    Read the input file
*/
void inputClass::readInputFile(std::string filename)
{
	std::ifstream input(filename);
	std::map<std::string, double> vars;
	std::string line;
	std::string varName;
	double varValue;
    int iSpecies;
	std::string memSpeciesName;
	
	cout<<"\nReading from the input file \n";

    resizeArrays(2);

	// Open file
	while(input.good())
	{
		getline(input, line);
		istringstream ss(line);
        string file;
        ss >> varName;
        if (varName == "GASCTIFILE_LUMEN")
        {
            ss >> m_filenameGas[0];
        } 
        else if (varName == "GASPHASE_LUMEN")
        {
            ss >> m_gasPhase[0];
        } 
        else if (varName == "SURFCTIFILE_LUMEN")
		{
			ss >> m_filenameSurf[0];
		}
        else if (varName == "SURFPHASE_LUMEN")
        {
            ss >> m_surfPhase[0];
        }
        if (varName == "GASCTIFILE_SUPPORT")
        {
            ss >> m_filenameGas[1];
        } 
        else if (varName == "GASPHASE_SUPPORT")
        {
            ss >> m_gasPhase[1];
        } 
        else if (varName == "SURFCTIFILE_SUPPORT")
		{
			ss >> m_filenameSurf[1];
		}
        else if (varName == "SURFPHASE_SUPPORT")
        {
            ss >> m_surfPhase[1];
        }
		else if (varName == "MEMSPNAME")
		{
			ss >> memSpeciesName;
		}        
        else if (varName == "SPECIESFRAC")
        {
            std::string spFrac;
            ss >> spFrac;
            if(spFrac == "MOLE")
            {
                is_moleFrac = 0;
            } 
            else if (spFrac == "MASS")
            {
                is_moleFrac = 1;
            }
        }
        else if (varName == "TRAN_MODEL")
        {
            std::string tranModel;
            ss >> tranModel;
            if(tranModel == "FICK")
            {
                m_tranModel = 0;
            }
            else if (tranModel == "DGM")
            {
                m_tranModel = 1;
            }
        }
        else if (varName == "VEL")
        {
            is_sccm = 0;
        }
        else if (varName == "SCCM")
        {
            is_sccm = 1;      
        }
        else if (varName == "INITGUESS_FILENAME")
        {
            ss >> m_initGuessFile;
        }
        ss >> varValue;
        // Read the parameter value
		vars[varName] = varValue;
	}

    // Lumen Geometry
    m_phi_gas[0] = vars["POROSITY_LUMEN"];
    m_As[0] = vars["SPECIFICAREA_LUMEN"];
    m_tau[0] = vars["TORTUOSITY_LUMEN"];
	m_Rin = vars["CHANNELRAD_LUMEN"];
    m_Dp[0] = vars["PARTICLEDIA_LUMEN"];
    m_rPore[0] = vars["PORERAD_LUMEN"];
    is_radVel = vars["RADIAL_VEL"];
    m_phi_micro[0] = vars["MICROPOROSITY_LUMEN"];
    m_rMicroPore[0] = vars["MICROPORERAD_LUMEN"];
    randPoreModel = vars["RAND_PORE_MODEL"];

    // Support geometry
    is_support = vars["SOLVE_SUPPORT"];
    m_phi_gas[1] = vars["POROSITY_SUPPORT"];
    m_As[1] = vars["SPECIFICAREA_SUPPORT"];
    m_tau[1] = vars["TORTUOSITY_SUPPORT"];
	m_Rout = vars["CHANNELRAD_SUPPORT"];
    m_Dp[1] = vars["PARTICLEDIA_SUPPORT"];
    m_rPore[1] = vars["PORERAD_SUPPORT"];
    m_thickness_2 = vars["THICKNESS_SUPPORT"];
    m_phi_micro[1] = vars["MICROPOROSITY_SUPPORT"];
    m_rMicroPore[1] = vars["MICROPORERAD_SUPPORT"];
    
    //EoS type
	updateEosType(0);
    if(is_support)
    {
        updateEosType(1);
    }

	// Boundary condition
	m_yInlet.resize(m_nsp[0],0.0);
    m_xInlet.resize(m_nsp[0], 0.0);
	m_thetaInlet.resize(m_nspSurf[0],0.0);
	m_inletTemp = vars["INLETTEMP"];
	m_outletPres = vars["OUTLETPRES"];

    // Energy equation
    m_Twall = vars["TWALL"];
    m_hcoeff = vars["HCOEFF"];
    
    // Membrane
    m_solveMem = vars["MEMSOLVE"];

    // Energy equation for Solid particle
    m_solveSolid = vars["SOLVE_SOLID"];
    m_condSolid = vars["COND_SOLID"];
    m_rhoSolid = vars["DENSITY_SOLID"];
    m_CpSolid = vars["CP_SOLID"];
    T_env = vars["T_ENV"];
    m_emissivity = vars["EMISSIVITY"];
    m_Aenv = vars["A_ENV"];
    
    // Apply initial conditions
    for (iSpecies = 0; iSpecies < m_nsp[0]; iSpecies++)
    {
        if (is_moleFrac == 0)  // Specify inlet mole fractions
        {
            m_xInlet[iSpecies] = 0.0; 
            m_xInlet[iSpecies] = vars[m_gas[0]->speciesName(iSpecies)];
        }
        else if (is_moleFrac == 1) //Specify inlet mass fractions
        {
            m_yInlet[iSpecies] = 0.0; 
            m_yInlet[iSpecies] = vars[m_gas[0]->speciesName(iSpecies)];
        }
    }
    if (is_moleFrac == 0)
    {
        m_gas[0]->setMoleFractions(m_xInlet.data());
        m_gas[0]->getMassFractions(m_yInlet.data());
        m_gas[0]->setState_TP(m_inletTemp, m_outletPres);
    }
    else
    {
        m_gas[0]->setMassFractions(m_yInlet.data());
        m_gas[0]->setState_TP(m_inletTemp, m_outletPres);
    }
    m_gas[0]->getMassFractions(m_yInlet.data());
    m_gas[0]->getMoleFractions(m_xInlet.data());
    
    // Surface coverages
    for (iSpecies = 0; iSpecies < m_nspSurf[0]; iSpecies++)
	{
		m_thetaInlet[iSpecies] = 0.0;
        m_thetaInlet[iSpecies] = vars[m_surf[0]->speciesName(iSpecies)];
	}

    if (is_sccm == 0)
	{
		m_inletVel = vars["VEL"];
	}
	else if (is_sccm == 1)
	{
		flowRate_sccm = vars["SCCM"];
		double area_out = M_PI*m_Rin*m_Rin;
        
		// Convert to standard velocity
		m_inletVel = flowRate_sccm * 1e-6 / 60 / area_out;	
        
        // Add pressure and temperature correction
		m_inletVel = m_inletVel * m_inletTemp / 273 * 1e5/m_outletPres;
	}

    // Print initial conditions
    cout << "\n Inlet velocity = " << m_inletVel << " m/s"; 
    cout << "\n Outlet pressure = " << m_outletPres << " Pa"; 
    cout << "\n Inlet temperature = " << m_inletTemp << " K"; 

    // Calculate mass flux
    m_inletFlux = m_inletVel * (m_gas[0]->density());

    // Membrane related input parameters
    if(m_solveMem)
    {
        if(is_support)
        {
            m_memSpeciesIndex = m_gas[1]->speciesIndex(memSpeciesName);
        }
        else{
            m_memSpeciesIndex = m_gas[0]->speciesIndex(memSpeciesName);
        }
        m_memThickness = vars["MEMTHICKNESS"];
        m_presMem = vars["SWEEPPRES"];
        m_perm0 = vars["k0"];
        m_Ea_R = vars["Ea_R"];
    }
    
    // Axial transport
    m_meshPoints = vars["MESHPOINTS"];
    m_axialLength = vars["LENGTH_LUMEN"];

    // Radial transport
    m_lumenPoints = vars["LUMEN_POINTS"];
    m_supportPoints = vars["SUPPORT_POINTS"];

    // Solve energy
    m_do_energy = vars["SOLVEENERGY"];

	//If the temperature is constant, this is called only once
	m_perm = m_perm0 * exp(-m_Ea_R / m_inletTemp);

    // Solve annular channel flow
    is_channel = vars["SOLVE_CHANNEL"];
    m_rChannel = vars["CHANNELRAD_ANNULUS"];
    u0_ch = vars["CHANNEL_VEL"];
    if(!m_solveMem)
    {
        is_channel = 0;
        cout<<"\n Since the membrane is not present, flow through the channel is not solved.";
    }

    //Solver tolerances
    m_Atol_ss = vars["ATOL_SS"];
    m_Rtol_ss = vars["RTOL_SS"];
    m_Atol_ts = vars["ATOL_TS"];
    m_Rtol_ts = vars["RTOL_TS"];
    m_loglevel = vars["LOGLEVEL"];
    m_refine = vars["REFINE"];
    m_maxGridPoints = vars["MAXGRIDPOINTS"];
    m_maxTimeSteps = vars["MAXTIMESTEPS"];
    m_minTimeSteps = vars["MINTIMESTEPS"];
    m_minTimeStepSize = vars["MINTIMESTEPSIZE"];
    m_maxTimeStepSize = vars["MAXTIMESTEPSIZE"];
  
    cout<<"\n refine = "<<m_refine;

    max_grid_ratio = vars["MAX_GRID_RATIO"];
    max_delta = vars["MAX_SLOPE"];
    max_delta_slope = vars["MAX_CURVE"];
    prune = vars["PRUNE"];
    m_restart = vars["RESTART"];
    m_initGuess = vars["INITGUESS_FROMFILE"];
    
    // Variable pointers
    m_NeqTotal = 0;
    if(is_support)
    {
        nParts = 2;
    }
    for (int i = 0; i<nParts; i++)
    {
        getVarPointers(i);
    }  
    if(!is_support || !m_supportPoints) 
    {
        m_supportPoints = 0;
        is_support = 0;
        m_NeqTotal =  m_Neq[0] * m_lumenPoints;
    }
    else{
        m_NeqTotal = m_Neq[0] * m_lumenPoints + m_Neq[1] * m_supportPoints;
    }

    if(is_channel)
    {
        m_NeqTotal += 3;  // Additional variables for rho, u and T inside the external flow channel
    }

    cout << "\nFinished reading the input file \n";
	// close file
	input.close();
} 

/*
    This subroutine reads the CTI file and creates gas and surf objects based on the given equation of state
*/
void inputClass::updateEosType(int i)
{
	// Read Gas-Phase
    if(m_filenameGas[i].empty())
    {
        cout<<"\n Error: No gas phase mechanism file is specified. Program will exit now...\n";
        exit(-1);
    } 
    //XML_Node* xgas = get_XML_File(m_filenameGas[i]);
	
    // Read Surface-Phase
    if(m_filenameSurf[i].empty())
    {
        cout<<"\n Error: No surface phase mechanism file is specified. Program will exit now... \n";
        exit(-1);
    } 
		
    // Read the gas and surface phase

	//Ideal gas EoS			
	cout << "Using Ideal gas equation of state.... \n";      
    cout << "Gas-phase name = "<< m_gasPhase[i] <<endl; 
    if(m_gasPhase[i].empty())
    {
        cout<<"\n Error: No gas phase is specified. Program will exit now... \n";
        exit(-1);
    }

    cout<<"\nReading gas-phase mechanism file \t" << m_filenameGas[i] << endl;	
    if (boost::algorithm::ends_with(m_filenameGas[i], "yaml"))
    {
        m_gas[i] = newPhase(m_filenameGas[i], m_gasPhase[i]);
        auto kin = newKinetics({m_gas[i]}, m_filenameGas[i], m_gasPhase[i]);
        m_kin[i] = dynamic_cast<GasKinetics*>(kin.release());
        m_tran[i] = newDefaultTransportMgr(m_gas[i]);
    } else {
        cout<<"\n Error: Unsupported file format. Specify YAML file. Program will exit now... \n";
        exit(-1);
    }
    m_nsp[i] = m_gas[i]->nSpecies();
	cout << "Number of gas species = " << m_nsp[i] << endl;
	cout << "Number of gas reactions = " << m_kin[i]->nReactions() << endl;

    if(m_tranModel ==1)
    {
        m_tranDGM[i] = dynamic_cast<DustyGasTransport*>(m_tran[i]);
    }
	//Surface phase (uses only ideal gas EoS)
    cout<<"\nReading surface-phase mechanism file \t" << m_filenameSurf[i] << endl;	
    if(m_surfPhase[i].empty())
    {
        cout<<"\n Error: No surface phase is specified. Program will exit now... \n";
        exit(-1);
    }
    if (boost::algorithm::ends_with(m_filenameSurf[i], "yaml"))
    {
        cout << "surface-phase name = "<< m_surfPhase[i] <<endl;
        m_surf[i] = dynamic_cast<SurfPhase*> (newPhase(m_filenameSurf[i], m_surfPhase[i]));
        auto kinSurf = newKinetics({m_surf[i], m_gas[i]},m_filenameSurf[i], m_surfPhase[i]);
        m_kin_surf[i] = dynamic_cast<InterfaceKinetics*>(kinSurf.release());
    } else {
        cout<<"\n Error: Unsupported file format. Specify YAML file. Program will exit now... \n";
        exit(-1);
    }
    m_nspSurf[i] = m_surf[i]->nSpecies();
    cout << "Number of species in surface phase = " << m_nspSurf[i] << endl;
	cout << "Number of surface reactions = " << m_kin_surf[i]->nReactions() << endl;
}

/* Define variable pointers
    U = 0, T = 1, Ygas = 2 to nsp+1
	and theta = (nsp + 2) to (nsp + m_nspSurf + 2) 
*/
void inputClass::getVarPointers(int i)
{
    offset_U[i] = 0;
    if(is_radVel)
    {
        offset_V[i] = 1;
        offset_T[i] = offset_V[i] + 1;
    } else{
        offset_T[i] = offset_U[i] + 1;
    }	
    offset_Ts[i] = offset_T[i] + 1;
    offset_Ygas[i] = offset_Ts[i] + 1;				        // offset_Ygas to (offset_Ygas + m_nsp-1)   
    offset_theta[i] = offset_Ygas[i] + m_nsp[i];	    	// offset_theta to (offset_theta + m_nspSurf)    
    m_Neq[i] = m_nspSurf[i] + offset_theta[i];
}

/*  
    Resize arrays
*/
void inputClass::resizeArrays(int n)
{
    m_phi_gas.resize(n, 0.0);
    m_tau.resize(n, 0.0);
    m_As.resize(n, 0.0);
    m_Dp.resize(n, 0.0);
    m_Bg.resize(n, 0.0);
    m_rPore.resize(n, 0.0);
    m_phi_micro.resize(n, 0.0);
    m_rMicroPore.resize(n, 0.0);

    m_nsp.resize(n, 0.0);
    m_nspSurf.resize(n, 0.0);
    
    offset_Ygas.resize(n, 0.0);
    offset_theta.resize(n, 0.0);
    offset_T.resize(n, 0.0);
    offset_Ts.resize(n, 0.0);
    offset_U.resize(n, 0.0);
    offset_P.resize(n, 0.0);
    m_Neq.resize(n, 0.0);
}