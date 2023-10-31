#include "BVP_solver.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"

using namespace std;
using namespace Cantera;
using namespace BVP;

/* 
    Constructor for the bvpQ2D class
*/
bvpQ2D::bvpQ2D(std::string filename, inputClass *ip_ptr, int nv, int points, double L, int energyFlag):
    BVP::BoundaryValueProblem(nv, points, 0.0, L),
    m_nsp(ip_ptr->m_nsp[0]),
	m_nspSurf(ip_ptr->m_nspSurf[0])
{
    // Crate user_func class instance
    // user_func write_instance();
    //m_func = &write_instance;

    // Input data
    m_axialPoints = points;
    m_nv = nv;
    m_nsp = ip_ptr->m_nsp[0];
    m_nspSurf = ip_ptr->m_nspSurf[0];
    m_do_energy = energyFlag;
    m_lumenPoints = ip_ptr->m_lumenPoints;
    m_supportPoints = ip_ptr->m_supportPoints;
    
    //Annular channel flow
    m_do_channel = ip_ptr->is_channel;
    m_rChannel = ip_ptr->m_rChannel;
    u0_ch = ip_ptr->u0_ch;
    
    // Total number of mesh points
    m_radialPoints = m_lumenPoints + m_supportPoints;
    m_totRadialPoints = m_radialPoints;
    if(m_do_channel)
    {
        m_totRadialPoints += 1;
    }
    m_points = m_totRadialPoints * m_axialPoints;

    // Support geometry
    m_solveSupport = ip_ptr->is_support;

    // Radial velocity
    m_solveRadialVel = ip_ptr->is_radVel;
    
    // Solid phase
    m_solveSolid = ip_ptr->m_solveSolid;
    m_condSolid = ip_ptr->m_condSolid;
    m_rhoSolid = ip_ptr->m_rhoSolid;
    m_CpSolid = ip_ptr->m_CpSolid;
    m_emissivity = ip_ptr->m_emissivity;
    T_env = ip_ptr->T_env;
    m_Aenv = ip_ptr->m_Aenv;

    // Restart
    m_restart = ip_ptr->m_restart;

    // Initial guess from file
    m_initGuess = ip_ptr->m_initGuess;
    m_initGuessFile = ip_ptr->m_initGuessFile;

    m_NeqTotal = ip_ptr->m_NeqTotal;

    m_nParts = ip_ptr->nParts;
    resizeGeometryArrays(m_nParts);

    nEq_lumen = ip_ptr->m_Neq[0];
    nv_lumen = nEq_lumen * m_lumenPoints;
    for(int i = 0; i<m_nParts; i++)
    {
        m_area2Vol[i] = ip_ptr->m_As[i];
        m_microPorosity[i] = ip_ptr->m_phi_micro[i];
        m_porosity[i] = ip_ptr->m_phi_gas[i];
        m_tortuosity[i] = ip_ptr->m_tau[i];
        m_Dp[i] = ip_ptr->m_Dp[i];
        m_PoreRad[i] = ip_ptr->m_rPore[i];
        m_microPoreRad[i] = ip_ptr->m_rMicroPore[i];

        //offset values        
        c_offset_U[i] = ip_ptr->offset_U[i];      // Axial velocity
        c_offset_T[i] = ip_ptr->offset_T[i];      // Temperature
        c_offset_Y[i] = ip_ptr->offset_Ygas[i];   // mass fractions of gas-phase species
        c_offset_Ts[i] = ip_ptr->offset_Ts[i];    // Solid phase temperature
        if(m_solveRadialVel)
        {
            c_offset_V[i] = ip_ptr->offset_V[i];
        }
        c_offset_theta[i] = ip_ptr->offset_theta[i]; // surface coverages
    }

    //Refine grid
    int refine = ip_ptr->m_refine;
    m_refineGrid = refine;
    m_maxGridPoints = ip_ptr->m_maxGridPoints;
    m_maxTimeSteps = ip_ptr->m_maxTimeSteps;
    m_minTimeSteps = ip_ptr->m_minTimeSteps;
    m_maxTimeStepSize = ip_ptr->m_maxTimeStepSize;
    m_minTimeStepSize = ip_ptr->m_minTimeStepSize;
    m_max_grid_ratio = ip_ptr->max_grid_ratio;
    m_max_delta = ip_ptr->max_delta;
    m_max_delta_slope = ip_ptr->max_delta_slope;
    m_prune = ip_ptr->prune;
    
    // Cantera Phases
    m_thermo = dynamic_cast<IdealGasPhase*> (ip_ptr->m_gas[0]);
    m_kin = (ip_ptr->m_kin[0]);
    m_trans = ip_ptr->m_tran[0];
    m_surfkin = ip_ptr->m_kin_surf[0];
	m_surf = ip_ptr->m_surf[0];

    // Thermo inputs
    T_in = ip_ptr->m_inletTemp;
    p_out = ip_ptr->m_outletPres;
    m_yInlet = ip_ptr->m_yInlet;
    m_xInlet = ip_ptr->m_xInlet;
    m_thetaInlet = ip_ptr->m_thetaInlet;
    m_inletVel = ip_ptr->m_inletVel;
    is_sccm = ip_ptr->is_sccm;
    if(is_sccm)
    {
        m_SCCM = ip_ptr->flowRate_sccm;
    }
    if(m_do_energy)
    {
        T_wall = ip_ptr->m_Twall;
        m_hcoeff = ip_ptr->m_hcoeff;
    }

    // Transport model
    m_tranModel = ip_ptr->m_tranModel;
    m_randPoreModel = ip_ptr->randPoreModel;
    
    if(m_tranModel == m_DGM)
    {
        m_tranDGM = ip_ptr->m_tranDGM[0];
        cout<<"\n Solving Dusty gas transport...";
    }
    else{
        cout<<"\n Solving Fickian diffusion transport...";
    }

    // Membrane parameters
    m_solveMem = ip_ptr->m_solveMem;
    if(m_solveMem)
    {
        m_memSpeciesIndex = ip_ptr->m_memSpeciesIndex;
        //m_perm = ip_ptr->m_perm;
        m_memThickness = ip_ptr->m_memThickness;
        m_presMem = ip_ptr->m_presMem;
        m_perm0 = ip_ptr->m_perm0;
        m_Ea_R = ip_ptr->m_Ea_R;
    }

    m_Rin = ip_ptr->m_Rin;
    m_length = ip_ptr->m_axialLength;

    if(m_solveSupport)
    {
        nEq_support = ip_ptr->m_Neq[1];
        m_Rout = ip_ptr->m_Rout;
        m_thickness_2 = ip_ptr->m_thickness_2;

        // Cantera Phases
        m_thermo_2 = dynamic_cast<IdealGasPhase*> (ip_ptr->m_gas[1]);
        m_kin_2 = ip_ptr->m_kin[1];
        m_trans_2 = ip_ptr->m_tran[1];
        m_surfkin_2 = ip_ptr->m_kin_surf[1];
        m_surf_2 = ip_ptr->m_surf[1];

        m_nsp_2 = ip_ptr->m_nsp[1];
        if(m_nsp_2 != m_nsp)
        {
            throw CanteraError("bvpQ2D::bvpQ2D", "Number of gas-phase species in the Lumen and the support must be the same");
        }
        m_nspSurf_2 = ip_ptr->m_nspSurf[1];
        m_gamma[1] = m_surf_2->siteDensity();
    }   

    if(m_do_channel)
    {
        // rho, u, T
        c_offset_ch_rho = 0;
        c_offset_ch_U = c_offset_ch_rho + 1;
        c_offset_ch_T = c_offset_ch_U + 1;
        nEq_channel = c_offset_ch_T + 1;
    }

    // Calculate inlet flux
    m_inletMassFlux = ip_ptr->m_inletFlux;

    // Set transport
    m_do_multicomponent = (m_trans->transportType() == "Multi");
    m_multidiff.resize(m_nsp*m_nsp*m_points, 0.0);
    
    // Calculate permeability
    for(int m = 0; m < m_nParts; m++)
    {
        KozenyCarmanPermeability(m);
    }

    // Inert gas index
    int Ar_index = m_thermo->speciesIndex("AR");

    //-------------- default solution bounds --------------------
    for(size_t i = 0; i<m_lumenPoints; i++)
    {
        int ind = i * nEq_lumen;
        setBounds(ind + c_offset_U[0], -1e20, 1e20); // no bounds on u    
        setBounds(ind + c_offset_T[0], 200.0, 2*m_thermo->maxTemp()); // temperature bounds
        setBounds(ind + c_offset_Ts[0], 200.0, 2*m_thermo->maxTemp()); // temperature bounds
        // mass fraction bounds
        for (size_t k = 0; k < m_nsp; k++) {
            setBounds(ind + c_offset_Y[0]+k, 1.0e-20, 1.0e5);
        }
        // surface coverages bounds
        for (size_t k = 0; k < m_nspSurf; k++) {
            setBounds(ind + c_offset_theta[0]+k, -1.0e-7, 1.0e5);
        }
        if(m_solveRadialVel)
        {
            setBounds(ind + c_offset_V[0], -1e20, 1e20); // no bounds on v
        }

        //-------------------- grid refinement -------------------------
        if(m_do_energy)
        {
            m_refiner->setActive(ind + c_offset_T[0], true);
            m_refiner->setActive(ind + c_offset_Ts[0], true);
        }
    }
    
    if(m_solveSupport)
    {
        for(size_t i = m_lumenPoints; i<m_radialPoints; i++)
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
        
            //-------------- default solution bounds --------------------
            setBounds(c_offset_U[1] + ind, -1e-20, 1e20); // Density
            setBounds(c_offset_T[1] + ind, 200.0, 2*m_thermo_2->maxTemp()); // temperature bounds
            setBounds(c_offset_Ts[1] + ind, 200.0, 2*m_thermo_2->maxTemp()); // temperature bounds
            m_refiner->setActive(c_offset_T[1] + ind, true);
            m_refiner->setActive(c_offset_Ts[1] + ind, true);
            for (size_t k = 0; k < m_nsp_2; k++) {
                setBounds(c_offset_Y[1] + k + ind, -1.0e-7, 1.0e5);
            }
            // surface coverages bounds
            for(size_t k = 0; k < m_nspSurf_2; k++) {
                setBounds(c_offset_theta[1] + k + ind, -1.0e-7, 1.0e5);
            }
            if(m_solveRadialVel)
            {
                setBounds(ind + c_offset_V[1], -1e20, 1e20); // no bounds on v
            }
            m_refiner->setActive(ind + c_offset_T[1], true);
            m_refiner->setActive(ind + c_offset_Ts[1], true);
        }
    }

    if(m_do_channel)
    {
        int ind = nEq_lumen*m_lumenPoints + nEq_support*m_supportPoints;
        setBounds(ind + c_offset_ch_U, -1e20, 1e20); // no bounds on u    
        setBounds(ind + c_offset_ch_T, 200.0, 2*m_thermo->maxTemp()); // temperature bounds
        setBounds(ind + c_offset_ch_rho, -1e-20, 1e5); // Density bounds
    }

    vector_double gr;
    for (size_t ng = 0; ng < m_axialPoints; ng++) {
        gr.push_back(L*ng/(m_axialPoints-1));
    }
    setupGrid(m_axialPoints, gr.data());
    setupRadialGrid();
    m_gamma[0] = m_surf->siteDensity();
    //cout<<"\n gamma = "<< m_surfkin->nReactions();

    // make a local copy of the species molecular weight vector
    m_wt.resize(m_nsp, 0.0);
    m_wt = m_thermo->molecularWeights();

    if(m_tranModel == m_DGM)
    {
        setDGMProperties();
    }
}   

/* 
    Setup the grid based on the finitevolume discretization
*/
void bvpQ2D::setupGrid(size_t n, const doublereal* z)
{
    resize(m_nv, n);
    double addConstant = 0;
    if(m_setupFirst)
    {
        addConstant = (z[1] - z[0])/2; // Shift the first point so that it becomes the cell center
        m_setupFirst = 0;
    }
    m_z[0] = addConstant + z[0];
    for (size_t j = 1; j < m_axialPoints; j++) {
        if (z[j] <= z[j-1]) {
            throw CanteraError("bvpQ2D::setupGrid", "grid points must be monotonically increasing");
        }
        m_z[j] = addConstant + z[j];
        m_dz[j-1] = m_z[j] - m_z[j-1];
        m_zWall[j - 1] = 0.5*(m_z[j] + m_z[j-1]); // right wall for cell (j-1)
    }
    // Half cell 
    m_zWall[m_axialPoints-1] = 2 * m_z[m_axialPoints-1] - m_zWall[m_axialPoints - 2];
}

/* 
    Resize all vectorsand arrays after the mesh refinement
*/
void bvpQ2D::resize(size_t ncomponents, size_t points)
{
    //points = points - 1;
    Domain1D::resize(ncomponents, points);
    m_axialPoints = points;
    m_points = m_axialPoints * m_totRadialPoints;
    m_rho.resize(m_axialPoints, m_totRadialPoints, 0.0);
    m_wtm.resize(m_axialPoints, m_totRadialPoints, 0.0);
    m_cp.resize(m_axialPoints, m_totRadialPoints, 0.0);
    m_press.resize(m_axialPoints, m_totRadialPoints, 0.0);
    m_ybar.resize(m_nsp, 0.0);
    
    m_flux.resize(m_nsp, m_points, 0.0);
    m_interfaceFlux.resize(m_nsp, m_points, 0.0);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_sdot.resize(m_nsp + m_nspSurf,m_points, 0.0);
    m_xMole.resize(m_nsp, m_points, 0.0);
    
    m_diff.resize(m_nsp*m_points, 0.0);
    m_diff_radial.resize(m_nsp*m_points, 0.0);
    m_tcon.resize(m_axialPoints, m_totRadialPoints, 0.0);
    m_visc.resize(m_axialPoints, m_totRadialPoints, 0.0);
    hk_interface.resize(m_nsp, m_points, 0.0);
    m_kExcess.resize(m_points,0.0);

    if (m_do_multicomponent || (m_tranModel == m_DGM) ) {
        m_multidiff.resize(m_nsp*m_nsp*m_points, 0.0);
        m_multidiff_radial.resize(m_nsp*m_nsp*m_points, 0.0);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    
    m_dz.resize(m_axialPoints-1, 0.0);
    m_zWall.resize(m_axialPoints, 0.0);
    m_z.resize(m_axialPoints, 0.0);
}

/* 
    Set thermodynamic state of the gas-phase at the cell center in the lumen
*/
void bvpQ2D::setGas(const doublereal* x, size_t j, size_t i)
{
    m_thermo->setTemperature(T(x,j,i));
    int ind = i * nEq_lumen;
    const doublereal* yy = x + m_nv*j + ind + c_offset_Y[0];
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press(j,i));
    m_rho(j,i) = m_thermo->density();
}

/* 
    Set thermodynamic state of the gas-phase at the cell center in the support
*/
void bvpQ2D::setGasSupport(const doublereal* x, size_t j, size_t i)
{
    m_thermo_2->setTemperature(T(x,j,i));
    int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
    const doublereal* yy = x + m_nv*j + ind + c_offset_Y[1];
    m_thermo_2->setMassFractions_NoNorm(yy);
    m_thermo_2->setPressure(m_press(j,i));
    m_rho(j,i) = m_thermo_2->density();
}

/* 
    Set thermodynamic state of the surface-phase at the cell center in the lumen
*/
void bvpQ2D::setSurf(const doublereal* x, size_t j, size_t i)
{
    m_surf->setTemperature(Ts(x,j,i));
    m_surf->setPressure(m_press(j,i));
    int ind = i * nEq_lumen;
    const doublereal* yy = x + m_nv*j + ind + c_offset_theta[0];
    m_surf->setCoveragesNoNorm(yy);
}

/* 
    Set thermodynamic state of the gas-phase at the cell center in the support
*/
void bvpQ2D::setSurfSupport(const doublereal* x, size_t j, size_t i)
{
    m_surf_2->setTemperature(Ts(x,j,i));
    m_surf_2->setPressure(m_press(j,i));
    int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
    const doublereal* yy = x + m_nv*j + ind + c_offset_theta[1];
    m_surf_2->setCoveragesNoNorm(yy); 
}

/* 
    Set thermodynamic state of the gas-phase at the cell boundary in the lumen
*/
void bvpQ2D::setGasAtMidpoint(const doublereal* x, size_t j, size_t i)
{
    m_thermo->setTemperature(0.5*(T(x,j,i)+T(x,j+1,i)));
    double pAvg = 0.5*(m_press(j,i) + m_press(j + 1,i));
    int ind = i * nEq_lumen;
    const doublereal* yyj = x + m_nv*j + ind + c_offset_Y[0];
    const doublereal* yyjp = x + m_nv*(j+1) + ind + c_offset_Y[0];
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(pAvg);
}

/* 
    Set thermodynamic state of the gas-phase at the cell boundary in the support
*/
void bvpQ2D::setGasAtMidpointSupport(const doublereal* x, size_t j, size_t i)
{
    m_thermo_2->setTemperature(0.5*(T(x,j,i)+T(x,j + 1,i)));
    double pAvg = 0.5*(m_press(j,i) + m_press(j + 1,i));
    int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
    const doublereal* yyj = x + m_nv*j + ind + c_offset_Y[1];
    const doublereal* yyjp = x + m_nv*(j + 1) + ind + c_offset_Y[1];
    for (size_t k = 0; k < m_nsp_2; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo_2->setMassFractions_NoNorm(m_ybar.data());
    m_thermo_2->setPressure(pAvg);
}

/* 
    Update transport and thermodynamic properties in the reactor
*/
void bvpQ2D::updateProperties(size_t jg, double* x, size_t jmin, size_t jmax)
{
    // properties are computed for grid points from j0 to j1
    size_t j0 = std::max<size_t>(jmin, 1) - 1;
    size_t j1 = std::min(jmax+1,m_axialPoints-1);
    updateThermo(x, j0, j1);
    if (jg == npos || m_force_full_update) {
        // update transport properties only if a Jacobian is not being
        // evaluated, or if specifically requested
        updateTransport(x, j0, j1);
        if(m_solveSupport || (m_lumenPoints > 1))
        {
            updateTransportRadial(x, j0, j1);
        }
    }

    // update the species diffusive mass fluxes whether or not a Jacobian is being evaluated
    updateDiffFluxes(x, j0, j1);
}

/* 
    Calculate axial diffusion fluxes in the reactor
*/
void bvpQ2D::updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
    for (size_t i = 0; i < m_radialPoints; i++) {
        if (m_tranModel == m_FICK)
        {
            if (m_do_multicomponent) {
                for (size_t j = j0; j <= (j1-1); j++) {
                    for (size_t k = 0; k < m_nsp; k++) {
                        double sum = 0.0;
                        for (size_t m = 0; m < m_nsp; m++) {
                            sum += m_wt[m] * m_multidiff[mindex(k,m,cellIndex(j,i))] * (X(x,m,j+1,i)-X(x,m,j,i));
                        }
                        m_flux(k,cellIndex(j,i)) = sum * m_diff[k+cellIndex(j,i)*m_nsp] / m_dz[j];
                        // Add convective flux
                        double rhou;
                        if(i>= m_lumenPoints)   //i.e. Support catalyst
                        {
                            setGasAtMidpointSupport(x,j,i);
                            rhou = (m_thermo_2->density())*0.5*(u(x,j,i) + u(x,j+1,i));
                        } else {
                            setGasAtMidpoint(x,j,i);
                            rhou = (m_thermo->density())*0.5*(u(x,j,i) + u(x,j+1,i));
                        }
                        m_flux(k,cellIndex(j,i)) += rhou * m_ybar[k];
                    }
                }
            } else 
            { //Fickian diffusion
                for (size_t j = j0; j <= (j1-1); j++) {
                    double sum = 0.0;
                    double wtm = m_wtm(j,i);
                    double rho = m_rho(j,i);
                    double dz = m_z[j+1] - m_z[j];
                    for (size_t k = 0; k < m_nsp; k++) {
                        m_flux(k,cellIndex(j,i)) = m_wt[k]*(rho*m_diff[k+m_nsp*cellIndex(j,i)]/wtm);
                        m_flux(k,cellIndex(j,i)) *= (X(x,k,j,i) - X(x,k,j+1,i))/m_dz[j];
                        sum -= m_flux(k,cellIndex(j,i));
                    }
                    // correction flux to ensure that \sum_k Y_k V_k = 0.
                    for (size_t k = 0; k < m_nsp; k++) {
                        //m_flux(k,cellIndex(j,i)) += sum*Y(x,k,j,i);
                        // Add convective flux
                        double rhou;
                        if(i>= m_lumenPoints) //i.e. Support catalyst
                        {
                            setGasAtMidpointSupport(x,j,i);
                            rhou = (m_thermo_2->density())*0.5*(u(x,j,i) + u(x,j+1,i));
                        }
                        else {
                            setGasAtMidpoint(x,j,i);
                            rhou = (m_thermo->density())*0.5*(u(x,j,i) + u(x,j+1,i));
                        }
                        m_flux(k,cellIndex(j,i)) += rhou * m_ybar[k];
                    }
                }
            }
        }
        /*else if (m_tranModel == m_DGM)
        {
            // Dusty gas model
            for (size_t j = j0; j < (j1-1); j++)
            {
                setGas(x,j,i); 
                m_tranDGM->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]); // unit: m2/s   
                m_thermo->saveState(state_prev);

                setGas(x,j+1,i); 
                m_thermo->saveState(state_curr);

                // Solve DGM transport to get mass fluxes (fluxes saved to m_flux)
                m_tranDGM->getMolarFluxes(&state_prev[0], &state_curr[0], (m_z[j+1] - m_z[j]), &m_flux(0,cellIndex(j,i))); 
                // Convert to mass fluxes
                for (size_t k = 0; k < m_nsp; k++)
                {
                    m_flux(k, cellIndex(j,i)) = m_wt[k] * m_flux(k, cellIndex(j,i));                  // unit: kg/m2/s
                }
            }
        }*/

        // Zero diffusion flux at the boundary
        for (size_t k = 0; k < m_nsp; k++) {
            // Add convective flux
            m_flux(k,cellIndex(m_axialPoints-1,i)) = m_rho(m_axialPoints-1,i)* u(x,m_axialPoints-1,i)*Y(x,k,m_axialPoints-1,i);
        }
    }
}

/* 
    Update transport properties in the axial direction
*/
void bvpQ2D::updateTransport(doublereal* x, size_t j0, size_t j1)
{
    for (size_t i = 0; i < m_radialPoints; i++) {
        if(m_tranModel == m_FICK)
        {
            if(m_do_multicomponent)
            {
                for (size_t j = j0; j < j1; j++) {
                    double rho, wtm;
                    if(i>= m_lumenPoints)
                    {
                        setGasAtMidpointSupport(x,j,i);
                        wtm = m_thermo_2->meanMolecularWeight();
                        rho = m_thermo_2->density();
                        m_visc(j,i) = m_trans_2->viscosity();
                        m_tcon(j,i) = m_trans_2->thermalConductivity();
                        m_trans_2->getMultiDiffCoeffs(m_nsp_2, &m_multidiff[mindex(0,0,cellIndex(j,i))]); // unit: m2/s
                    } else {
                        setGasAtMidpoint(x,j,i);
                        wtm = m_thermo->meanMolecularWeight();
                        rho = m_thermo->density();
                        m_visc(j,i) = m_trans->viscosity();
                        m_tcon(j,i) = m_trans->thermalConductivity();
                        m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,cellIndex(j,i))]); // unit: m2/s
                    }

                    // Use m_diff as storage for the factor outside the summation
                    for (size_t k = 0; k < m_nsp; k++)
                    {
                        m_diff[k+cellIndex(j,i)*m_nsp] = m_wt[k] * rho / (wtm*wtm);  //unit: kmol/m3
                    }
                }
            } else { // mixture averaged transport
                for (size_t j = j0; j < j1; j++) 
                {
                    if(i>= m_lumenPoints)
                    {
                        setGasAtMidpointSupport(x,j,i);

                        // Mixture-average coefficients
                        m_trans_2->getMixDiffCoeffs(&m_diff[cellIndex(j,i)*m_nsp_2]);

                        // Effective diffusion coefficients
                        calculateDiffCoeff(T(x,j,i), j, i, 0);
                        m_visc(j,i) = m_trans->viscosity();
                        m_tcon(j,i) = m_trans->thermalConductivity();
                    } else {
                        setGasAtMidpoint(x,j,i);

                        // Mixture-average coefficients
                        m_trans->getMixDiffCoeffs(&m_diff[cellIndex(j,i)*m_nsp]);

                        // Effective diffusion coefficients
                        calculateDiffCoeff(T(x,j,i), j, i, 0);
                        m_visc(j,i) = m_trans->viscosity();
                        m_tcon(j,i) = m_trans->thermalConductivity();
                    }
                }
            }
        }
        else if(m_tranModel == m_DGM)
        {
            /*for (size_t j = j0; j < j1; j++) 
            {
                setGasAtMidpoint(x,j,i);
                m_visc(j,i) = m_tranDGM->viscosity();
                m_tcon(j,i) = m_tranDGM->thermalConductivity();
            }*/
        }
    }
}

/* 
    Update transport properties in the radial direction
*/
void bvpQ2D::updateTransportRadial(doublereal* x, size_t j0, size_t j1)
{
    for (size_t i = 0; i < (m_radialPoints-1); i++) {
        if(m_tranModel == m_FICK)
        {
            if(m_do_multicomponent)
            {
                for (size_t j = j0; j <= j1; j++) {
                    double rho, wtm;
                    if(i>= m_lumenPoints)
                    {
                        setGasAtInterface(x,j,i);
                        wtm = m_thermo_2->meanMolecularWeight();
                        rho = m_thermo_2->density();
                        m_trans_2->getMultiDiffCoeffs(m_nsp_2, &m_multidiff_radial[mindex(0,0,cellIndex(j,i))]); // unit: m2/s
                    } else {
                        setGasAtInterface(x,j,i);
                        wtm = m_thermo->meanMolecularWeight();
                        rho = m_thermo->density();
                        m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff_radial[mindex(0,0,cellIndex(j,i))]); // unit: m2/s
                    }

                    // Use m_diff_radial as storage for the factor outside the summation
                    for (size_t k = 0; k < m_nsp; k++)
                    {
                        m_diff_radial[k+cellIndex(j,i)*m_nsp] = m_wt[k] * rho / (wtm*wtm);  //unit: kmol/m3
                    }
                }
            } else { // mixture averaged transport
                for (size_t j = j0; j <= j1; j++) 
                {
                    if(i>= m_lumenPoints)
                    {
                        setGasAtInterface(x,j,i);

                        // Mixture-average coefficients
                        m_trans_2->getMixDiffCoeffs(&m_diff_radial[cellIndex(j,i)*m_nsp_2]);

                        // Effective diffusion coefficients
                        calculateDiffCoeff(T(x,j,i), j, i, 1);
                    } else {
                        setGasAtInterface(x,j,i);

                        // Mixture-average coefficients
                        m_trans->getMixDiffCoeffs(&m_diff_radial[cellIndex(j,i)*m_nsp]);

                        // Effective diffusion coefficients
                        calculateDiffCoeff(T(x,j,i), j, i, 0);
                    }
                }
            }
        }
        else if(m_tranModel == m_DGM)
        {
            /*for (size_t j = j0; j < j1; j++) 
            {
                setGasAtMidpoint(x,j,i);
                m_visc(j,i) = m_tranDGM->viscosity();
                m_tcon(j,i) = m_tranDGM->thermalConductivity();
            }*/
        }
    }
}

/* 
    Calculate Knudsen coefficient
*/
Array2D bvpQ2D::calculateKnudsenDiffCoeff(double T, int m)
{
    // If the temperature is constant, Knudsen coefficients are calculated only once.
    
    Array2D D_Kn;
    double const1, const2 = 0.0;
    D_Kn.resize(m_nsp, 2, 0.0);

    const1 = (2.00 / 3.00) * (m_PoreRad[m] * m_porosity[m] / m_tortuosity[m]) * sqrt(8 * GasConstant * T / M_PI);
    if(m_randPoreModel)
    {
        const2 = (2.00 / 3.00) * (m_microPoreRad[m] * m_microPorosity[m] / m_tortuosity[m]) * sqrt(8 * GasConstant * T / M_PI);
    }
    
    for (int k = 0; k < m_nsp; k++)
    {
        D_Kn(k,0) = const1 * sqrt(1 / m_wt[k]);
        D_Kn(k,1) = const2 * sqrt(1 / m_wt[k]);
    }
    return D_Kn;
}

/* 
    Calculate effective diffusion coefficient
*/
void bvpQ2D::calculateDiffCoeff(double T, int j, int i, int n)
{
    //Get Knudsen coefficients    
    int m = 0;
    if(i>=m_lumenPoints)
    {        
        m = 1;
    }

    double phi_tau = m_porosity[m]/m_tortuosity[m];
    double microPhi_tau = m_microPorosity[m]/m_tortuosity[m];
    Array2D D_Kn = calculateKnudsenDiffCoeff(T,m);
    vector_double diffCoeff(m_nsp, 0.0), microDiffCoeff(m_nsp,0.0);

    for (int k = 0; k < m_nsp; k++)
    {
        diffCoeff[k] = 1 / (1 / (phi_tau*m_diff[k+m_nsp*cellIndex(j,i)]) + 1 / D_Kn(k,0));
        if(m_randPoreModel)
        {
            microDiffCoeff[k] = 1 / (1 / (microPhi_tau*m_diff[k+m_nsp*cellIndex(j,i)]) + 1 / D_Kn(k,1));
            
            // Calculate effective diffusion coefficient using Random-Pore model
            double const1 = m_microPorosity[m]*m_microPorosity[m]*(1 + 3*m_porosity[m])/(1 - m_porosity[m]);
            diffCoeff[k] = m_porosity[m]*m_porosity[m]*diffCoeff[k] + const1* microDiffCoeff[k];
        }
        // Update m_diff and m_diff_radial arrays
        if(n == 0)
        {
            m_diff[k+m_nsp*cellIndex(j,i)] = diffCoeff[k];
        } else {
            m_diff_radial[k+m_nsp*cellIndex(j,i)] = diffCoeff[k];
        }
    }
}

/* 
    Solver function
*/
void bvpQ2D::eval(size_t jg, doublereal* xg,
                  doublereal* rg, integer* diagg, doublereal rdt)
{
    // if evaluating a Jacobian, and the global point is outside the domain of
    // influence for this domain, then skip evaluating the residual
    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }
    resetBadValues(xg);
    //needJacUpdate();

    //Save pressure and density to p_prev and rho_prev
    p_prev.resize(m_axialPoints, m_totRadialPoints);
    rho_prev.resize(m_axialPoints, m_totRadialPoints);
    for(size_t j = 0; j < m_axialPoints; j++)
    {
        for(size_t i = 0; i < m_totRadialPoints; i++)
        {
            p_prev(j,i) = m_press(j,i);
            rho_prev(j,i) = m_rho(j,i);
        }
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = diagg + loc();

    size_t jmin, jmax;
    if (jg == npos) { // evaluate all points
        jmin = 0;
        jmax = m_axialPoints - 1;
    } else { // evaluate points for Jacobian
        size_t jpt = (jg == 0) ? 0 : jg - firstPoint();
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_axialPoints-1);
    }
    updateProperties(jg, x, jmin, jmax);
    interfaceFlux(x,jmin,jmax);
    if(m_solveSupport)
    {
        flux_LumenSupport(x, jmin, jmax);
    }
    for (size_t i = 0; i < m_totRadialPoints; i++) {
        if(m_do_channel && i ==(m_totRadialPoints-1))
        {
            evalAnnularChannel(x, rsd, diag, rdt, jmin, jmax);
        } else {
            evalResidual(x, rsd, diag, rdt, jmin, jmax, i);
        }
    }
    //getNH3Conversion(x);
    //showSolution(x);
} 

/* 
    Main residual function
*/
void bvpQ2D::evalResidual(double* x, double* rsd, int* diag,
                          double rdt, size_t jmin, size_t jmax, size_t i)
{
    //------------------------------------------------------------
    // evaluate the residual equations at all required grid points
    //------------------------------------------------------------

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         Inlet boundary
        //----------------------------------------------

        if (j == 0) {
            evalInletBoundary(x, rsd, diag, rdt, i);
        } else if (j == m_axialPoints - 1) {
            evalOuterBoundary(x, rsd, diag, rdt, i);
		} 
        else { // interior points
            
            //--------------------------------------------
            //    Continuity equation
            //
            //    d(\rho u)/dz = 0
            //---------------------------------------------
                        
            getWdot(x,j,i);
            getSdot(x,j,i);

            int m, ind, surfnsp;
            if(i < m_lumenPoints)
            {
                m = 0;
                ind = i * nEq_lumen;
                surfnsp = m_nspSurf;
            } else {
                m = 1;
                ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
                surfnsp = m_nspSurf_2;
            }

            // Calculate chemical and diffusion terms
            vector_double diffus(m_nsp), chem(m_nsp);
            double sumChem = 0, sumDiffus = 0, sumY = 0;

            for (size_t k = 0; k < m_nsp; k++) {
                diffus[k] = 2.0*(m_flux(k,cellIndex(j,i)) - m_flux(k,cellIndex(j-1,i))) / (m_z[j+1] - m_z[j-1]);
                chem[k] = m_area2Vol[m] * sdot(k+surfnsp,j,i) * m_wt[k];
                sumChem += chem[k];
                sumDiffus += diffus[k];
                sumY += Y(x,k,j,i);
            }
            // Flux through the membrane
            radialFluxTerm.resize(m_nsp,0.0);
            double memFluxTerm = 0.0, sumRadial = 0.0;
            //Radial flux
            for (size_t k = 0; k < m_nsp; k++)
            {
                if(i == 0 && m_radialPoints != 1)
                {
                    radialFluxTerm[k] = m_rWall[0]*m_interfaceFlux(k,cellIndex(j,0));
                    radialFluxTerm[k] *= 2/(m_rWall[0]*m_rWall[0] - 0);
                    sumRadial += radialFluxTerm[k];
                } else if (i == (m_radialPoints-1)) {
                    if(i != 0)
                    {
                        radialFluxTerm[k] = - m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
                        radialFluxTerm[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                        sumRadial += radialFluxTerm[k];
                    }
                    
                    if(k == m_memSpeciesIndex && m_solveMem)
                    {
                        memFluxTerm = getMembraneFlux(x, j) * m_rWall[i];
                        if(i != 0)
                        {
                            memFluxTerm *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                        } else{
                            memFluxTerm *= 2/(m_rWall[i]*m_rWall[i]);
                        }
                    }
                } else{
                    radialFluxTerm[k] = m_rWall[i]*m_interfaceFlux(k,cellIndex(j,i)) - m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
                    radialFluxTerm[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                    sumRadial += radialFluxTerm[k];
                }
            }
            double drho_dt = 1/m_porosity[m]*(sumChem - sumDiffus - memFluxTerm - sumRadial);
            rsd[index(c_offset_U[m] + ind, j)] = drho_dt;
            diag[index(c_offset_U[m] + ind, j)] = 0;

            //-----------------------------------------------------------------
            //    Species equations
            //
            //    \rho \phi dY_k/dt + \rho u dY_k/dz + d(J_k)/dz
            //       = A_s*sdot_k*W_k
            //-----------------------------------------------------------------

            for (size_t k = 0; k < m_nsp; k++) {
                rsd[index(c_offset_Y[m] + k + ind, j)] 
                = (chem[k] - diffus[k] - radialFluxTerm[k])/(m_porosity[m])
                   - rdt*(m_rho(j,i)*Y(x,k,j,i) - rho_prev(j,i)*Y_prev(k,j,i));
                if(k == m_memSpeciesIndex && m_solveMem)
                {
                    rsd[index(c_offset_Y[m] + k + ind, j)] -= memFluxTerm/m_porosity[m];
                }
                diag[index(c_offset_Y[m] + k + ind, j)] = 1;
            }

            //-----------------------------------------------------------------------
            //    Surface Species equations
            //
            //    d(\theta)/dt = sdot_k/gamma 
            //-----------------------------------------------------------------------

            double sumYsurf = 0;
            for (size_t k = 0; k < surfnsp; k++) {
                rsd[index(c_offset_theta[m] + k + ind, j)] 
                = sdot(k,j,i) / m_gamma[m] - rdt*(surfY(x,k,j,i) - surfY_prev(k,j,i)); 
                diag[index(c_offset_theta[m] + k + ind, j)] = 1;
                sumYsurf += surfY(x,k,j,i); 
            } 

            // Last species
            size_t k = surfnsp -1;
            rsd[index(c_offset_theta[m] + k + ind, j)] = 1 - sumYsurf;
            diag[index(c_offset_theta[m] + k + ind, j)] = 0;

            // residual equations if the energy equation is disabled
            if(m_do_energy)
            {
                double hk = 0;
                
                // Get enthalpies required to solve energy equation with membrane
                vector_double hk_RT(m_nsp, 0.0);

                //Radial flux
                vector_double radial_hkjk(m_nsp, 0.0);
                double sum_hkjk = 0.0, Temp;
                double qWall = 0.0;
                for (size_t k = 0; k < m_nsp; k++)
                {
                    if(i == 0 && m_radialPoints != 1)
                    {
                        setGasAtInterface(x,j,i);
                        radial_hkjk[k] = hk_interface(k,cellIndex(j,0))*m_rWall[0]*m_interfaceFlux(k,cellIndex(j,0));
                        radial_hkjk[k] *= 2/(m_rWall[0]*m_rWall[0] - 0);
                    } else if (i == (m_radialPoints-1)) {
                        if(i != 0)
                        {
                            setGasAtInterface(x,j,i-1);
                            radial_hkjk[k] = - hk_interface(k,cellIndex(j,i-1))*m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
                            radial_hkjk[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                        }
                        
                        if(i < m_lumenPoints)
                        {
                            setGas(x,j,i);
                            hk_RT = m_thermo->enthalpy_RT_ref();
                            Temp = m_thermo->temperature();
                        } else{
                            setGasSupport(x,j,i);
                            hk_RT = m_thermo_2->enthalpy_RT_ref();
                            Temp = m_thermo_2->temperature();
                        }
                        // Membrane heat flux
                        if(k == m_memSpeciesIndex && m_solveMem)
                        {
                            hk = hk_RT[m_memSpeciesIndex]*GasConstant*Temp;
                        }

                        // Heat transfer between the wall and the reactor
                        qWall = m_hcoeff*(T_wall - Temp)* m_rWall[i];
                        if(i != 0)
                        {
                            qWall *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                        } else{
                            qWall *= 2/(m_rWall[i]*m_rWall[i]);
                        }
                    } else{
                        setGasAtInterface(x,j,i);
                        radial_hkjk[k] = hk_interface(k,cellIndex(j,i))*m_rWall[i]*m_interfaceFlux(k,cellIndex(j,i));

                        setGasAtInterface(x,j,i-1);
                        radial_hkjk[k] -= hk_interface(k,cellIndex(j,i-1))*m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
                        radial_hkjk[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                    }
                    sum_hkjk += radial_hkjk[k];
                }

                // Calculate gradient
                rsd[index(c_offset_T[m] + ind, j)] = - grad_qGas(x,j,i) - qdotConv(x,j,i) - qdotSurf(x,j,i) - memFluxTerm*hk - sum_hkjk + qWall;
                rsd[index(c_offset_T[m] + ind, j)] /= (m_rho(j,i)*m_cp(j,i)*m_porosity[m]);
                rsd[index(c_offset_T[m]+ ind, j)] -= rdt*(T(x,j,i) - T_prev(j,i));
                diag[index(c_offset_T[m] + ind, j)] = 1;

                // Solid phase temperature
                rsd[index(c_offset_Ts[m] + ind, j)] =  - grad_qSolid(x,j,i) + qdotConv(x,j,i) + qdotSurf(x,j,i) - qdotEnv(x,j,i);
                rsd[index(c_offset_Ts[m] + ind, j)] /= (m_rhoSolid*m_CpSolid*(1 - m_porosity[m]));
                rsd[index(c_offset_Ts[m] + ind, j)] -= rdt*(Ts(x,j,i) - Ts_prev(j,i));
                diag[index(c_offset_Ts[m] + ind, j)] = 1;
            }
            else
            {
                rsd[index(c_offset_T[m] + ind, j)] = T(x,j,i) - T_in;
                diag[index(c_offset_T[m] + ind, j)] = 0;

                rsd[index(c_offset_Ts[m] + ind, j)] = Ts(x,j,i) - T_in;
                diag[index(c_offset_Ts[m] + ind, j)] = 0;
            }
        }
    }
}

/* 
    Residual function at the inlet cell
*/
void bvpQ2D::evalInletBoundary(double* x, double* rsd, int* diag, double rdt, size_t i)
{        
    getWdot(x,0,i);
    getSdot(x,0,i);

    int m, ind, surfnsp;
    if(i < m_lumenPoints)
    {
        m = 0;
        ind = i * nEq_lumen;
        surfnsp = m_nspSurf;
    } else {
        m = 1;
        ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
        surfnsp = m_nspSurf_2;
    } 

    // Calculate chemical and diffusion terms
    vector_double diffus(m_nsp), chem(m_nsp);
    double sumChem = 0, sumDiffus = 0, sumY = 0;

    for (size_t k = 0; k < m_nsp; k++) {
        if(i < m_radialPoints)
        {
            diffus[k] = (m_flux(k,cellIndex(0,i)) - m_inletMassFlux*m_yInlet[k]) / m_dz[0];
        } 
        else
        {
            diffus[k] = (m_flux(k,cellIndex(0,i))) / m_dz[0];
        }
        chem[k] = m_area2Vol[m] * sdot(k+surfnsp,0,i) * m_wt[k];
        sumChem += chem[k];
        sumDiffus += diffus[k];
        sumY += Y(x,k,0,i);
    }
    // Flux through the membrane
    radialFluxTerm.resize(m_nsp,0.0);
    double memFluxTerm = 0.0, sumRadial = 0.0;
    //Radial flux
    for (size_t k = 0; k < m_nsp; k++)
    {
        if(i == 0 && m_radialPoints != 1)
        {
            radialFluxTerm[k] = m_rWall[0]*m_interfaceFlux(k,cellIndex(0,0));
            radialFluxTerm[k] *= 2/(m_rWall[0]*m_rWall[0] - 0);
            sumRadial += radialFluxTerm[k];
        } else if (i == (m_radialPoints-1)) {
            if(i != 0)
            {
                radialFluxTerm[k] = - m_rWall[i-1]*m_interfaceFlux(k,cellIndex(0,i-1));
                radialFluxTerm[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                sumRadial += radialFluxTerm[k];
            }

            if(k == m_memSpeciesIndex && m_solveMem)
            {
                memFluxTerm = getMembraneFlux(x, 0) * m_rWall[i];
                if(i == 0)
                {
                    memFluxTerm *= 2/(m_rWall[i]*m_rWall[i]);
                } else {
                    memFluxTerm *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                }
            }
        } else{
            radialFluxTerm[k] = m_rWall[i]*m_interfaceFlux(k,cellIndex(0,i)) - m_rWall[i-1]*m_interfaceFlux(k,cellIndex(0,i-1));
            radialFluxTerm[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
            sumRadial += radialFluxTerm[k];
        }
    }

    double drho_dt = 1/m_porosity[m]*(sumChem - sumDiffus - memFluxTerm - sumRadial);
    rsd[index(c_offset_U[m] + ind, 0)] = drho_dt;
    diag[index(c_offset_U[m] + ind, 0)] = 0;

    //-----------------------------------------------------------------
    //    Species equations
    //
    //    \rho \phi dY_k/dt + \rho u dY_k/dz + d(J_k)/dz
    //       = A_s*sdot_k*W_k
    //-----------------------------------------------------------------

    for (size_t k = 0; k < m_nsp; k++) {
        rsd[index(c_offset_Y[m] + k + ind, 0)] 
        = (chem[k] - diffus[k] - radialFluxTerm[k])/(m_porosity[m])
            - rdt*(m_rho(0,i)*Y(x,k,0,i) - rho_prev(0,i)*Y_prev(k,0,i));
        if(k == m_memSpeciesIndex && m_solveMem)
        {
            rsd[index(c_offset_Y[m] + k + ind, 0)] -= memFluxTerm/m_porosity[m];
        }
        diag[index(c_offset_Y[m] + k + ind, 0)] = 1;
    }

    double sumYsurf = 0;
    for (size_t k = 0; k < surfnsp; k++)
    {
        diag[index(c_offset_theta[m] + k + ind, 0)] = 1;
        sumYsurf += surfY(x,k,0,i);
        rsd[index(c_offset_theta[m] + k + ind, 0)] = sdot(k,0,i) / m_gamma[m] - rdt*(surfY(x,k,0,i) - surfY_prev(k,0,i)); 
    }  
    size_t k = surfnsp -1;
    rsd[index(c_offset_theta[m] + k + ind, 0)] = 1 - sumYsurf;
    diag[index(c_offset_theta[m] + k + ind, 0)] = 0;

    if(m_do_energy)
    {
        solveEnergyatInlet(x, rsd, diag, rdt, memFluxTerm, i);
    } else {
        rsd[index(c_offset_T[0] + ind, 0)] = T(x,0,i) - T_in;
        diag[index(c_offset_T[0] + ind, 0)] = 0;

        rsd[index(c_offset_Ts[0] + ind, 0)] = Ts(x,0,i) - T_in;
        diag[index(c_offset_Ts[0] + ind, 0)] = 0;
    }
}

/* 
    Residual function at the outlet cell
*/
void bvpQ2D::evalOuterBoundary(double* x, double* rsd, int* diag, double rdt, size_t i)
{
    //------------------------------------------------
    //    Continuity equation
    //
    //    d(\rho u)/dz = 0
    //-------------------------------------------------
    int j = m_axialPoints - 1;
    int m, ind, surfnsp;
    if(i < m_lumenPoints)
    {
        m = 0;
        ind = i * nEq_lumen;
        surfnsp = m_nspSurf;
    } else {
        m = 1;
        ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
        surfnsp = m_nspSurf_2;
    }

    getWdot(x,j,i);
    getSdot(x,j,i);

    // Calculate chemical and diffusion terms
    vector_double diffus(m_nsp), chem(m_nsp);
    double sumChem = 0, sumDiffus = 0, sumFlux = 0;
    double z_out = 2*m_z[j] - m_z[j-1];

    for (size_t k = 0; k < m_nsp; k++) {
        if(i < m_radialPoints)
        {
            diffus[k] = (m_flux(k,cellIndex(j,i)) - m_flux(k,cellIndex(j-1,i))) / (m_zWall[j] - m_zWall[j-1]);
        } 
        else
        {
            diffus[k] = (0 - m_flux(k,cellIndex(j-1,i))) / (m_zWall[j] - m_zWall[j-1]);
        }
        chem[k] = m_area2Vol[m] * sdot(k+surfnsp,j,i) * m_wt[k];
        sumChem += chem[k];
        sumDiffus += diffus[k];
        sumFlux += m_flux(k,cellIndex(j-1,i));
    }
    // Flux through the membrane
    radialFluxTerm.resize(m_nsp,0.0);
    double memFluxTerm = 0.0, sumRadial = 0.0;
    //Radial flux
    for (size_t k = 0; k < m_nsp; k++)
    {
        if(i == 0 && m_radialPoints != 1)
        {
            setGasAtInterface(x,j,i);
            radialFluxTerm[k] = m_rWall[0]*m_interfaceFlux(k,cellIndex(j,0));
            radialFluxTerm[k] *= 2/(m_rWall[0]*m_rWall[0] - 0);
        } else if (i == (m_radialPoints-1)) {
            if(i != 0)
            {
                radialFluxTerm[k] = - m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
                radialFluxTerm[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
            }

            if(k == m_memSpeciesIndex && m_solveMem)
            {
                memFluxTerm = getMembraneFlux(x, j) * m_rWall[i];
                if(i == 0)
                {
                    memFluxTerm *= 2/(m_rWall[i]*m_rWall[i]);
                } else {;
                    memFluxTerm *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
                }
            }
        } else{
            radialFluxTerm[k] = m_rWall[i]*m_interfaceFlux(k,cellIndex(j,i)) - m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
            radialFluxTerm[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);            
        }
        sumRadial += radialFluxTerm[k];
    }

    double drho_dt = 1/m_porosity[m]*(sumChem - sumDiffus - memFluxTerm - sumRadial);

    // Free flow boundary
    diag[index(c_offset_U[m] + ind, j)] = 0;
    rsd[index(c_offset_U[m] + ind, j)] = drho_dt;

    //-----------------------------------------------------------------------
    //    Species equations
    //
    //    \rho \phi dY_k/dt + \rho u dY_k/dz + d(J_k)/dz
    //       = A_s*sdot_k*W_k
    //-----------------------------------------------------------------------
    double sum = 0;
    for (size_t k = 0; k < m_nsp; k++) {
        rsd[index(c_offset_Y[m] + k + ind, j)] 
            = (chem[k] - diffus[k] - radialFluxTerm[k])/(m_porosity[m])
            - rdt*(m_rho(j,i)*Y(x,k,j,i) - rho_prev(j,i)*Y_prev(k,j,i));
        if(k == m_memSpeciesIndex && m_solveMem)
        {
            rsd[index(c_offset_Y[m] + k + ind, j)] -= memFluxTerm/m_porosity[m];
        }
        diag[index(c_offset_Y[m] + k + ind, j)] = 1;
        sum += Y(x,k,j,i);
    }
    
    //-----------------------------------------------------------------------
    //    Surface Species equations
    //
    //    d(\theta)/dt = sdot_k/gamma 
    //-----------------------------------------------------------------------
    double sumYsurf = 0;
    for (size_t k = 0; k < surfnsp; k++) {
        rsd[index(c_offset_theta[m] + k + ind, j)] 
            = sdot(k,j,i) / m_gamma[m] - rdt*(surfY(x,k,j,i) - surfY_prev(k,j,i)); 
        diag[index(c_offset_theta[m] + k + ind, j)] = 1;
        sumYsurf += surfY(x,k,j,i);
    }  
    size_t k = surfnsp - 1;
    rsd[index(c_offset_theta[m] + k + ind, j)] = 1 - sumYsurf;
    diag[index(c_offset_theta[m] + k + ind, j)] = 0;         

    // residual equations if the energy equation is disabled
    if(m_do_energy)
    {
        solveEnergyatOutlet(x, rsd, diag, rdt, memFluxTerm, i);
    } else {
        rsd[index(c_offset_T[m] + ind, j)] = T(x,j,i) - T_in;
        diag[index(c_offset_T[m] + ind, j)] = 0;

        rsd[index(c_offset_Ts[m] + ind, j)] = Ts(x,j,i) - T_in;
        diag[index(c_offset_Ts[m] + ind, j)] = 0;
    }
}

/* 
    Calculate  a flux through a membrane at mesh point j
*/
double bvpQ2D::getMembraneFlux(double* x, size_t j)
{
    // Calculate permeability at current temperature
    double Tgas;
    int i;
    if(m_solveSupport && m_solveMem)
    {
        i = m_lumenPoints + m_supportPoints - 1;
        setGasSupport(x,j,i);
    } else if(!m_solveSupport && m_solveMem) {
        i = m_lumenPoints-1;
        setGas(x,j,i);
    }
    
    Tgas = T(x, j, i);
    m_perm = m_perm0 * exp(-(m_Ea_R) / Tgas);

    //Calculate partial pressure for membrane-permeable species
    double wt = m_wt[m_memSpeciesIndex];
    double pSpecies = m_press(j,i)*X(x, m_memSpeciesIndex, j, i);

    // Pressure inside the channel
    if(m_do_channel)
    {
        m_presMem = m_press(j,i+1);
    }

    double memFlux = (m_perm/m_memThickness) * (pSpecies - m_presMem) * wt;   // unit: kg/m2/s
    return memFlux;
}

/* 
    Initialize the solution vector
*/
void bvpQ2D::_getInitialSoln(double* x)
{
    if(!m_restart)
    {
        for (size_t j = 0; j < m_axialPoints; j++) {
            for (size_t i = 0; i < m_lumenPoints; i++)
            {
                vector_double cov(m_nspSurf, 0.0);
                int ind = i*nEq_lumen;
                // Set velocity
                x[index(c_offset_U[0] + ind, j)] = m_inletVel;

                //Set mass fractions, Temperature and pressure
                m_thermo->setMassFractions(m_yInlet.data());
                m_thermo->setState_TP(T_in, p_out);
                
                //Surface
                m_surf->setTemperature(T_in);
                m_surf->setCoveragesNoNorm(m_thetaInlet.data());
                m_surfkin->solvePseudoSteadyStateProblem(1,1);
                m_surf->getCoverages(cov.data());
                m_rho(j,i) = m_thermo->density();
                x[index(c_offset_T[0] + ind, j)] = m_thermo->temperature();
                x[index(c_offset_Ts[0] + ind, j)] = m_thermo->temperature();
                m_press(j,i) = m_thermo->pressure();
                m_thermo->getMassFractions(&x[index(c_offset_Y[0] + ind, j)]);
                m_surf->getCoverages(&x[index(c_offset_theta[0] + ind, j)]);
            }

            for (size_t i = m_lumenPoints; i < (m_lumenPoints + m_supportPoints); i++)
            {
                vector_double cov(m_nspSurf_2, 0.0);
                int ind = nv_lumen + (i - m_lumenPoints)*nEq_support;
                //Set velocity
                x[index(c_offset_U[1] + ind, j)] = m_inletVel;

                //Set mass fractions, Temperature and pressure
                m_thermo_2->setMassFractions(m_yInlet.data());
                m_thermo_2->setState_TP(T_in, p_out);
                
                //Surface
                m_surf_2->setTemperature(T_in);
                m_surfkin_2->advanceCoverages(1);
                m_surf_2->getCoverages(cov.data());
                m_rho(j,i) = m_thermo->density();
                x[index(c_offset_T[1] + ind, j)] = m_thermo_2->temperature();
                x[index(c_offset_Ts[1] + ind, j)] = m_thermo_2->temperature();
                m_press(j,i) = m_thermo_2->pressure();
                m_thermo_2->getMassFractions(&x[index(c_offset_Y[1] + ind, j)]);
                m_surf_2->getCoverages(&x[index(c_offset_theta[1] + ind, j)]);
            }

            // Annular channel
            if(m_do_channel)
            {
                m_thermo->setState_TPX(T_in, m_presMem, "H2:1.00");
                int ind = m_NeqTotal - nEq_channel;
                rho0_ch = m_thermo->density();
                x[index(c_offset_ch_U + ind, j)] = u0_ch;
                x[index(c_offset_ch_T + ind, j)] = T_in;
                x[index(c_offset_ch_rho + ind, j)] = rho0_ch;
            }
        }
    }
    if(m_initGuess)
    {
        readCSV(m_initGuessFile);
        interpolate(x);
    }
}

/* 
    Function to give component name based on the index in the solution vector
*/
string bvpQ2D::componentName(size_t m) const
{
    if(m < nv_lumen)        // Points inside the Lumen
    {
        size_t n = (m % nEq_lumen);
        switch (n) {
        case 0:
            return ("Velocity");
        case 1:
            return ("GasTemperature");
        case 2:
            return ("SolidTemperature");
        default:
            if (n >= c_offset_Y[0] && n < (c_offset_Y[0] + m_nsp)) {
                return (m_thermo->speciesName(n - c_offset_Y[0]));
            } else if (n >= c_offset_theta[0] && n < (c_offset_theta[0] + m_nspSurf))
            {
                return (m_surf->speciesName(n - c_offset_theta[0]));
            } else {
                return "<unknown>";
            }
        }
    } else if ((m - nv_lumen) < nEq_support*m_supportPoints) {
        size_t n = (m - nv_lumen) % nEq_support;
        size_t i = m_lumenPoints + (m - nv_lumen)/nEq_support;
        switch (n) {
        case 0:
            return ("Velocity_" + to_string(i));
        case 1:
            return ("GasTemperature_"+ to_string(i));
        case 2:
            return ("SolidTemperature_"+ to_string(i));
        default:
            if (n >= c_offset_Y[1] && n < (c_offset_Y[1] + m_nsp_2)) {
                return (m_thermo_2->speciesName(n - c_offset_Y[1])+ "_" + to_string(i));
            } else if (n >= c_offset_theta[1] && n < (c_offset_theta[1] + m_nspSurf_2))
            {
                return (m_surf_2->speciesName(n - c_offset_theta[1]) + "_" + to_string(i));
            } else {
                return "<unknown>";
            }
        }
    } else if ((m < m_NeqTotal) && m_do_channel)
    {
        // Annular channel
        size_t n = m - (m_NeqTotal - nEq_channel);
        switch(n) {
        case 0: 
            return ("Annular_Density");
        case 1:
            return ("Annular_Velocity");
        case 2:
            return ("Annular_Temperature");
        }
    } else {
        throw CanteraError("bvpQ2D::componentName",
                           "Wrong index" + m);
    }
}

/* 
    Function to give component index in the solution vector based on the name
*/
size_t bvpQ2D::componentIndex(const std::string& name) const
{
    if (name=="Velocity") {
        return 0;
    } else if (name=="GasTemperature") {
        return 1;
    } else if (name=="SolidTemperature") {
        return 2;
    } else {
        for (size_t n=c_offset_Y[0]; n<(m_nsp+c_offset_Y[0]); n++) {
            std::string name1 = componentName(n);
            boost::trim_left(name1);
            if (name1.compare(name) == 0) {
                return n;
            }
        }
        for (size_t n=c_offset_theta[0]; n<(m_nspSurf+c_offset_theta[0]); n++) {
            std::string name1 = componentName(n);
            boost::trim_left(name1);
            if (name1.compare(name) == 0) {
                return n;
            }
        }
        throw CanteraError("bvpQ2D::componentIndex",
                           "no component named " + name);
    }
}

/* 
    Function to print current solution
*/
void bvpQ2D::showSolution(const doublereal* x)
{
    Domain1D::showSolution(x);
}

/* 
    Function to reset the dependant variables if they are out of range
*/
void bvpQ2D::resetBadValues(doublereal* xg)
{
    double* x = xg + loc();
    for (size_t i = 0; i < m_totRadialPoints; i++) {
        int ind;
        if(i < m_lumenPoints)
        {
            ind = i * nEq_lumen;
            for (size_t j = 0; j < m_axialPoints; j++) {
                double* Y = x + m_nv*j + ind + c_offset_Y[0];
                m_thermo->setMassFractions_NoNorm(Y);
                m_thermo->getMassFractions(Y);
            }
            for (size_t j = 0; j < m_axialPoints; j++) {
                double* Z = x + m_nv*j + ind + c_offset_theta[0];
                m_surf->setCoveragesNoNorm(Z);
                m_surf->getCoverages(Z);
            }
        } else if (i < (m_lumenPoints + m_supportPoints)) {
            ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            for (size_t j = 0; j < m_axialPoints; j++) {
                double* Y = x + m_nv*j + ind + c_offset_Y[1];
                m_thermo->setMassFractions_NoNorm(Y);
                m_thermo->getMassFractions(Y);
            }
            for (size_t j = 0; j < m_axialPoints; j++) {
                double* Z = x + m_nv*j + ind + c_offset_theta[1];
                m_surf->setCoveragesNoNorm(Z);
                m_surf->getCoverages(Z);
            }
        }
    }
}

/* 
    Set the DGM properties if the dusty gas model is used.
*/
int bvpQ2D::setDGMProperties()
{
    m_thermo->setState_TPY(T_in, p_out, m_yInlet.data());
    m_tranDGM->setPorosity(m_porosity[0]);
    m_tranDGM->setTortuosity(m_tortuosity[0]);
    m_tranDGM->setMeanPoreRadius(m_PoreRad[0]);
    m_tranDGM->setMeanParticleDiameter(m_Dp[0]);
    cout<<"\n Finished setting up properties for the DGM transport \n";
    
    return 0;
}

/* 
    Calculate Kozeny-Carman permeability
*/
void bvpQ2D::KozenyCarmanPermeability(int m)
{
    /* Subroutine to calculate permeability */

    //Effective factor for the permeability
    double effFactor = m_porosity[m] / m_tortuosity[m];
    double volRatio = m_porosity[m] / (1 - m_porosity[m]);

    // Kozeny-Carman Permeability
    m_Bg[m] = pow(volRatio, 2) * pow(m_Dp[m], 2) * effFactor / 72;
}

/* Heat fluxes used in energy equation*/

/*
    This function calculates the gradient of gas-phase heat flux q_g
*/
double bvpQ2D::grad_qGas(const doublereal* x, size_t j, size_t i)
{
	vector_double h_right(m_nsp, 0.0), h_left(m_nsp,0.0);
    double RTright, RTleft;
    double phi;

    // Calculate \sum h_kj_k at the right and left interfaces
    if(i < m_lumenPoints)
    {
        phi = m_porosity[0];
        setGasAtMidpoint(x,j,i);
        h_right = m_thermo->enthalpy_RT_ref(); 
        RTright = GasConstant*(m_thermo->temperature());

        setGasAtMidpoint(x,j - 1,i);
        h_left = m_thermo->enthalpy_RT_ref();
        RTleft = GasConstant*(m_thermo->temperature());
    } else {
        phi = m_porosity[1];
        setGasAtMidpointSupport(x,j,i);
        h_right = m_thermo_2->enthalpy_RT_ref(); 
        RTright = GasConstant*(m_thermo_2->temperature());

        setGasAtMidpoint(x,j - 1,i);
        h_left = m_thermo_2->enthalpy_RT_ref();
        RTleft = GasConstant*(m_thermo_2->temperature());
    }

    // Calculate \sum grad(hkjk)
    double sum_hkjk_r = 0.0, sum_hkjk_l = 0.0;
    for (size_t k = 0; k < m_nsp; k++) 
    {
        h_right[k] *= RTright;
        h_left[k] *= RTleft;
        sum_hkjk_l += h_left[k]*m_flux(k,cellIndex(j-1,i));
        sum_hkjk_r += h_right[k]*m_flux(k,cellIndex(j,i));
    }
    double grad_hkjk;
    if(j == m_axialPoints - 1)
    {
        grad_hkjk = (sum_hkjk_r - sum_hkjk_l)/(m_z[j] - m_z[j-1]);
    } else {
        grad_hkjk = 2*(sum_hkjk_r - sum_hkjk_l)/(m_z[j+1] - m_z[j-1]);
    }
    double qGas_j = grad_hkjk + divHeatFlux(x,j,i)*phi;
    return qGas_j;
}

/*
    This function calculates the gradient of solid-phase heat flux q_g
*/
double bvpQ2D::grad_qSolid(const doublereal* x, size_t j, size_t i)
{
    double qSolid_p, qSolid_m;
    double T = 0.5*(Ts(x,j,i) + Ts(x,j-1,i));
    qSolid_m = -getConductivitySolid(T,i)*(Ts(x,j,i) - Ts(x,j-1,i))/(z(j) - z(j - 1));
    if(j == m_axialPoints-1)
    {
        qSolid_p = 0;
        return (0 - qSolid_m/(z(j) - z(j-1)));
    }
    else {
        T = 0.5*(Ts(x,j,i) + Ts(x,j+1,i));
        qSolid_p = -getConductivitySolid(T,i)*(Ts(x,j+1,i) - Ts(x,j,i))/(z(j+1) - z(j));
        return 2.0*(qSolid_p - qSolid_m)/(z(j+1) - z(j-1));
    }
}

/*
    This function calculates effective conductivity of the solid phase as follows:
	lambda_eff = phi_s*m_condSolid + cond_r
*/
double bvpQ2D::getConductivitySolid(double Ts, size_t i)
{
	double m_tconSolid, lambda_cond, const1;
    m_condSolid = solidConductivity(Ts);
    if(i<m_lumenPoints)
    {
        lambda_cond = (1 - m_porosity[0]) * m_condSolid;
        const1 = 4 * m_Dp[0] * Boltzmann * pow(Ts, 3);
    }
    else {
        lambda_cond = (1 - m_porosity[1]) * m_condSolid;
        const1 = 4 * m_Dp[1] * Boltzmann * pow(Ts, 3);
    }
	double lambda_star = m_condSolid / const1;
	double lambda_r = const1*(0.5756* m_emissivity * atan(1.5353*pow(lambda_star, 0.8011) / m_emissivity) + 0.1843);
	m_tconSolid =  (lambda_cond + lambda_r);
    return m_tconSolid;
}

/*
    This function calculates the net energy transfer between solid and gaseous phases due to convection
	qdotConv = hv*(Tg- Ts)
*/
double bvpQ2D::qdotConv(const doublereal* x, size_t j, size_t i)
{
	double qdotConv, hv;

    // calculates volumetric heat tranfer coefficient between the porous flow and the solid phase :
    // hv = Nu *lambda_gas/Dp^2
    if(i<m_lumenPoints)
    {
        hv = calculateNu(x,j,i)*m_tcon(j,i)/(m_Dp[0]*m_Dp[0]);
    }
    else{
        hv = calculateNu(x,j,i)*m_tcon(j,i)/(m_Dp[1]*m_Dp[1]);
    }
    
    qdotConv = hv*(T(x,j,i) - Ts(x,j,i));
	return qdotConv;
}

/*
    Calculate Nusselt number
*/
double bvpQ2D::calculateNu(const doublereal* x, size_t j, size_t i)
{
	// Calculate Nusselt number
	double Re, Pr, Nu, mu;

    // Calculate viscosity at the cell center
    if(i<m_lumenPoints)
    {
        setGas(x,j,i);
        mu = m_trans->viscosity();
        Re = abs(u(x,j,i) * m_rho(j,i) * m_Dp[0] / mu);
    } else{
        setGasSupport(x,j,i);
        mu = m_trans_2->viscosity();
        Re = abs(u(x,j,i) * m_rho(j,i) * m_Dp[1] / mu);
    }	
    Pr = mu * m_cp(j,i) / m_tcon(j,i);
	Nu = 2 + 1.1 * pow(Re, 0.6) * pow(Pr, (1 / 3));
    return Nu;
}

/*
    This function calculates radiative heat contribution qdotEnv:
    qdotEnv = sigma*epsilon*A_env*(Ts^4 - Tenv^4)
*/
double bvpQ2D::qdotEnv(const doublereal* x, size_t j, size_t i)
{
	double const1 = Boltzmann * m_emissivity *m_Aenv;
	double qdotEnv = const1 * (pow(Ts(x,j,i), 4) - pow(T_env, 4));
	return 0;//qdotEnv;
}

/*
    This function calculates the net energy transfer from heterogeneous surface reactions
*/
double bvpQ2D::qdotSurf(double* x, size_t j, size_t i)
{
	// Get sdot and enthalpy values at T = Tgas
    getSdot(x,j,i);
    vector_double h_RT_gas(m_nsp, 0.0), h_RT_surf(m_nspSurf, 0.0);
    int ind, m;
    int surfnsp;
    if(i< m_lumenPoints)
    {
        m = 0;
        surfnsp = m_nspSurf;
        h_RT_gas = m_thermo->enthalpy_RT_ref();
        ind = i * nEq_lumen;
        // Get sdot and enthalpy values at T = Tsurf
        m_thermo->setTemperature(Ts(x,j,i)); 
        const doublereal* yy = x + m_nv*j + ind + c_offset_Y[m];
        m_thermo->setMassFractions_NoNorm(yy);
        m_thermo->setPressure(m_press(j,i));
        h_RT_surf = m_thermo->enthalpy_RT_ref();
    } else{
        m = 1;
        surfnsp = m_nspSurf_2;
        h_RT_gas = m_thermo_2->enthalpy_RT_ref();
        ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
        // Get sdot and enthalpy values at T = Tsurf
        m_thermo_2->setTemperature(Ts(x,j,i)); 
        const doublereal* yy = x + m_nv*j + ind + c_offset_Y[m];
        m_thermo_2->setMassFractions_NoNorm(yy);
        m_thermo_2->setPressure(m_press(j,i));
        h_RT_surf.resize(m_nspSurf_2, 0.0);
        h_RT_surf = m_thermo_2->enthalpy_RT_ref();

    }
	double sum = 0;
    for (int k = 0; k < m_nsp; k++)
	{
		if (sdot(k+surfnsp,j,i) >= 0) {
			sum -= sdot(k+surfnsp,j,i) * h_RT_surf[k] * Ts(x,j,i);
		}
		else {
			sum -= sdot(k+surfnsp,j,i) * h_RT_gas[k] * T(x,j,i);
		}
	}
	double qdotSurf = GasConstant * sum * m_area2Vol[m];	
	return qdotSurf;
}

/*
    This function calculates conductivity of the solid phase
*/
double bvpQ2D::solidConductivity(double T)
{
    double cond_solid = 5.5 + 34.5*exp(-0.0033*(T - 273.15));
    return cond_solid;
}

/*
    Calculate pressure using Darcy's law
*/
void bvpQ2D::getPressure(const double* x)
{
    double mu, uAvg, dpdz;
    for (size_t i = 0; i < m_lumenPoints; i++) {
        m_press(m_axialPoints - 1,i) = p_out;
        for (int j = m_axialPoints -2; j >= 0; j--)
        {
            setGasAtMidpoint(x,j,i);
            mu = m_trans->viscosity();
            uAvg = u(x,j,i) + u(x,j + 1,i);
            dpdz = - (mu/m_Bg[0]) * uAvg;
            m_press(j,i) = m_press(j + 1, i) - dpdz *(m_z[j + 1] - m_z[j]);
        }
    }
    if(m_solveSupport)
    {
        for (size_t i = m_lumenPoints; i < m_radialPoints; i++) {
            m_press(m_axialPoints - 1,i) = p_out;
            for (int j = m_axialPoints -2; j >= 0; j--)
            {
                setGasAtMidpointSupport(x,j,i);
                mu = m_trans_2->viscosity();
                uAvg = u(x,j,i) + u(x,j + 1,i);
                dpdz = - (mu/m_Bg[1]) * uAvg;
                m_press(j,i) = m_press(j + 1, i) - dpdz *(m_z[j + 1] - m_z[j]);
            }
        }
    }
}

/*
    Solve energy equation at the inlet cell
*/
void bvpQ2D::solveEnergyatInlet(double* x, double* rsd, int* diag, double rdt, double flux, size_t i)
{              
    int m, ind;
    int j = 0;
    if(i < m_lumenPoints)
    {
        m = 0;
        ind = i * nEq_lumen;
    } else {
        m = 1;
        ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
    }
    
    // Get enthalpies required to solve energy equation with membrane
    vector_double hk_RT(m_nsp, 0.0);
    vector_double h_right(m_nsp, 0.0), h_left(m_nsp,0.0);

    //Radial flux
    vector_double radial_hkjk(m_nsp, 0.0);
    double sum_hkjk = 0.0;
    double phi, hk, Temp;
    double qWall = 0.0;
    for (size_t k = 0; k < m_nsp; k++)
    {
        if(i == 0 && m_radialPoints!= 1)
        {
            setGasAtInterface(x,0,i);
            radial_hkjk[k] = hk_interface(k,cellIndex(0,0))*m_rWall[0]*m_interfaceFlux(k,cellIndex(0,0));
            radial_hkjk[k] *= 2/(m_rWall[0]*m_rWall[0] - 0);
        } else if (i == (m_radialPoints-1)) {
            if(i != 0)
            {
                setGasAtInterface(x,j,i-1);
                radial_hkjk[k] = - hk_interface(k,cellIndex(0,i-1))*m_rWall[i-1]*m_interfaceFlux(k,cellIndex(0,i-1));
                radial_hkjk[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
            }

            // Membrane heat flux
            if(i < m_lumenPoints)
            {
                setGas(x,0,i);
                hk_RT = m_thermo->enthalpy_RT_ref();
                Temp = m_thermo->temperature();
            } else{
                setGasSupport(x,0,i);
                hk_RT = m_thermo_2->enthalpy_RT_ref();
                Temp = m_thermo_2->temperature();
            }
            if(k == m_memSpeciesIndex && m_solveMem)
            {
                hk = hk_RT[m_memSpeciesIndex]*GasConstant*Temp;
            }
            radial_hkjk[k] = 0;

            // Heat transfer between the wall and the reactor
            qWall = m_hcoeff*(T_wall - Temp)* m_rWall[i];
            if(i != 0)
            {
                qWall *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
            } else{
                qWall *= 2/(m_rWall[i]*m_rWall[i]);
            }
        } else {
            setGasAtInterface(x,0,i);
            radial_hkjk[k] = hk_interface(k,cellIndex(j,i))*m_rWall[i]*m_interfaceFlux(k,cellIndex(j,i));

            setGasAtInterface(x,0,i-1);
            radial_hkjk[k] -= hk_interface(k,cellIndex(j,i-1))*m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
            radial_hkjk[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]); 
        }
        sum_hkjk += radial_hkjk[k];
    }

    // Calculate \sum h_kj_k at the right interface
    double RTright, RTleft;
    if(i < m_lumenPoints)
    {
        phi = m_porosity[0];
        setGasAtMidpoint(x,0,i);
        h_right = m_thermo->enthalpy_RT_ref(); 
        RTright = GasConstant*(m_thermo->temperature());

        // Calculate \sum h_kj_k at the left interface
        m_thermo->setTemperature(T_in);
        m_thermo->setMassFractions_NoNorm(m_yInlet.data());
        m_thermo->setPressure(m_press(0,i));
        h_left = m_thermo->enthalpy_RT_ref();
        RTleft = GasConstant*T_in;
    } else {
        phi = m_porosity[1];
        setGasAtMidpointSupport(x,0,i);
        h_right = m_thermo_2->enthalpy_RT_ref(); 
        RTright = GasConstant*(m_thermo_2->temperature());

        // Calculate \sum h_kj_k at the left interface
        m_thermo_2->setTemperature(T_in);
        m_thermo_2->setMassFractions_NoNorm(m_yInlet.data());
        m_thermo_2->setPressure(m_press(0,i));
        h_left = m_thermo_2->enthalpy_RT_ref();
        RTleft = GasConstant*T_in;
    }

    // Calculate \sum grad(hkjk)
    double sum_hkjk_r = 0.0, sum_hkjk_l = 0.0;
    for (size_t k = 0; k < m_nsp; k++) 
    {
        h_right[k] *= RTright;
        h_left[k] *= RTleft;
        sum_hkjk_l += h_left[k]*m_inletMassFlux*m_yInlet[k];
        sum_hkjk_r += h_right[k]*m_flux(k,cellIndex(0,i));
    }
    double qGas_j = (sum_hkjk_r - sum_hkjk_l)/(m_z[1] - m_z[0]);

    // Calculate gradT
    double qGas_p = - m_tcon(0,i)*(T(x,1,i) - T(x,0,i))/(m_z[1] - m_z[0]);
    double qGas_m = 0;
    qGas_j += phi * (qGas_p - qGas_m)/(m_z[1] - m_z[0]);

    double Tavg = 0.5*(Ts(x,0,i) + Ts(x,1,i));
    double qSolid_p = - getConductivitySolid(Tavg,i)*(Ts(x,1,i) - Ts(x,0,i))/(m_z[1] - m_z[0]);
    double qSolid_m = 0;
    double qSolid_j = (qSolid_p - qSolid_m)/(m_z[1] - m_z[0]);

    // Calculate convective and radiative heat transfer across the inlet
    double q_env_gas = - m_hcoeff*(T(x,0,i) - T_env);
    double q_env_solid = - m_hcoeff*(Ts(x,0,i) - T_env);
    
    // Gas phase temperature
    rsd[index(c_offset_T[m] + ind, 0)] = - qGas_j - qdotConv(x,0,i) - qdotSurf(x,0,i) - flux*hk - sum_hkjk - q_env_gas/(m_z[1] - m_z[0]) + qWall;
    rsd[index(c_offset_T[m] + ind, 0)] /= (m_rho(0,i)*m_cp(0,i)*phi);
    rsd[index(c_offset_T[m] + ind, 0)] -= rdt*(T(x,0,i) - T_prev(0,i));
    diag[index(c_offset_T[m] + ind, 0)] = 1;
                
    // Solid phase temperature
    rsd[index(c_offset_Ts[m] + ind, 0)] =  - qSolid_j + qdotConv(x,0,i) + qdotSurf(x,0,i) - qdotEnv(x,0,i) + q_env_solid/(m_z[1] - m_z[0]);
    rsd[index(c_offset_Ts[m] + ind, 0)] /= (m_rhoSolid*m_CpSolid*(1 - phi));
    rsd[index(c_offset_Ts[m] + ind, 0)] -= rdt*(Ts(x,0,i) - Ts_prev(0,i));
    diag[index(c_offset_Ts[m] + ind, 0)] = 1;
}

/*
    Solve energy equation at the outlet cell
*/
void bvpQ2D::solveEnergyatOutlet(double* x, double* rsd, int* diag, double rdt, double flux, size_t i)
{
    int m, ind;
    int j = m_axialPoints-1;
    if(i < m_lumenPoints)
    {
        m = 0;
        ind = i * nEq_lumen;
    } else {
        m = 1;
        ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
    }
    
    // Get enthalpies required to solve energy equation with membrane
    vector_double hk_RT(m_nsp, 0.0);
    vector_double h_right(m_nsp, 0.0), h_left(m_nsp,0.0);

    //Radial flux
    vector_double radial_hkjk(m_nsp, 0.0);
    double sum_hkjk = 0.0;
    double phi, hk, Temp;
    double qWall = 0.0;
    for (size_t k = 0; k < m_nsp; k++)
    {
        if(i == 0 && m_radialPoints!= 1)
        {
            setGasAtInterface(x,j,0);
            radial_hkjk[k] = hk_interface(k,cellIndex(j,0))*m_rWall[0]*m_interfaceFlux(k,cellIndex(j,0));
            radial_hkjk[k] *= 2/(m_rWall[0]*m_rWall[0] - 0);
        } else if (i == (m_radialPoints-1)) {
            if(i != 0)
            {
                setGasAtInterface(x,j,i-1);
                radial_hkjk[k] = - hk_interface(k,cellIndex(j,i-1))*m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
                radial_hkjk[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
            }

            // Membrane heat flux
            if(i < m_lumenPoints)
            {
                setGas(x,j,i);
                hk_RT = m_thermo->enthalpy_RT_ref();
                Temp = m_thermo->temperature();
            } else{
                setGasSupport(x,j,i);
                hk_RT = m_thermo_2->enthalpy_RT_ref();
                Temp = m_thermo_2->temperature();
            }
            if(k == m_memSpeciesIndex && m_solveMem)
            {
                hk = hk_RT[m_memSpeciesIndex]*GasConstant*Temp;
            }

            // Heat transfer between the wall and the reactor
            qWall = m_hcoeff*(T_wall - Temp)* m_rWall[i];
            if(i != 0)
            {
                qWall *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
            } else{
                qWall *= 2/(m_rWall[i]*m_rWall[i]);
            }
        } else {
            setGasAtInterface(x,j,i);
            radial_hkjk[k] = hk_interface(k,cellIndex(j,i))*m_rWall[i]*m_interfaceFlux(k,cellIndex(j,i));

            setGasAtInterface(x,j,i-1);
            radial_hkjk[k] -= hk_interface(k,cellIndex(j,i-1))*m_rWall[i-1]*m_interfaceFlux(k,cellIndex(j,i-1));
            radial_hkjk[k] *= 2/(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]); 
        }
        sum_hkjk += radial_hkjk[k];
    }

    // Calculate \sum h_kj_k at the right and left interfaces

    double RTright, RTleft;
    if(i < m_lumenPoints)
    {
        // Calculate \sum h_kj_k at the right interface
        setGas(x,j,i);
        h_right = m_thermo->enthalpy_RT_ref(); 
        RTright = GasConstant*(m_thermo->temperature());

        // Calculate \sum h_kj_k at the left interface
        setGasAtMidpoint(x,j - 1,i);
        h_left = m_thermo->enthalpy_RT_ref();
        RTleft = GasConstant*(m_thermo->temperature());
    } else {
        // Calculate \sum h_kj_k at the right interface
        setGasSupport(x,j,i);
        h_right = m_thermo_2->enthalpy_RT_ref(); 
        RTright = GasConstant*(m_thermo_2->temperature());

        // Calculate \sum h_kj_k at the left interface
        setGasAtMidpointSupport(x,j - 1,i);
        h_left = m_thermo_2->enthalpy_RT_ref();
        RTleft = GasConstant*(m_thermo_2->temperature());
    }

    // Calculate \sum grad(hkjk)
    double sum_hkjk_r = 0.0, sum_hkjk_l = 0.0;
    for (size_t k = 0; k < m_nsp; k++) 
    {
        h_right[k] *= RTright;
        h_left[k] *= RTleft;
        sum_hkjk_l += h_left[k]*m_flux(k,cellIndex(j-1,i));
        sum_hkjk_r += h_right[k]*m_flux(k,cellIndex(j,i));
    }
    double qGas_j = (sum_hkjk_r - sum_hkjk_l)/(m_z[j] - m_z[j-1]);

    // Calculate gradT
    qGas_j += m_porosity[m] * divHeatFlux(x,j,i);

    // Calculate convective and radiative heat transfer across the outlet
    double q_env_gas = m_hcoeff*(T(x,j,i) - T_env);
    double q_env_solid = m_hcoeff*(Ts(x,j,i) - T_env) ;//+ Boltzmann * m_emissivity*(pow(Ts(x,j,i), 4) - pow(T_env, 4)));

    // Gas phase temperature
    rsd[index(c_offset_T[m] + ind, j)] = - qGas_j - qdotConv(x,j,i) - qdotSurf(x,j,i) - flux*hk - sum_hkjk - q_env_gas/(m_z[j] - m_z[j-1]) + qWall;
    rsd[index(c_offset_T[m] + ind, j)] /= (m_rho(j,i)*m_cp(j,i)*m_porosity[m]);
    rsd[index(c_offset_T[m] + ind, j)] -= rdt*(T(x,j,i) - T_prev(j,i));
    diag[index(c_offset_T[m] + ind, j)] = 1;

    // Solid phase temperature
    rsd[index(c_offset_Ts[m] + ind, j)] =  - grad_qSolid(x,j,i) + qdotConv(x,j,i) + qdotSurf(x,j,i) - qdotEnv(x,j,i) - q_env_solid/(m_z[j] - m_z[j-1]);
    rsd[index(c_offset_Ts[m] + ind, j)] /= (m_rhoSolid*m_CpSolid*(1 - m_porosity[m]));
    rsd[index(c_offset_Ts[m] + ind, j)] -= rdt*(Ts(x,j,i) - Ts_prev(j,i));
    diag[index(c_offset_Ts[m] + ind, j)] = 1;
}

/*
    Get ammonia conversion
*/
void bvpQ2D::getNH3Conversion(double* x)
{
	int j = m_axialPoints - 1;

    //get index of NH3, AR species
    int nh3_index = m_thermo->speciesIndex("NH3");
    int Ar_index = m_thermo->speciesIndex("AR");

    double mfr_in, mfr_out;
    double sumDenom = 0.0;
    double sumNum = 0.0;
    double area_i;
    for(int i = 0; i < m_radialPoints; i++)
    {
        if(i == 0)
        {
            area_i = M_PI*m_rWall[i]*m_rWall[i];
        } else{
            area_i = M_PI*(m_rWall[i]*m_rWall[i] - m_rWall[i-1]*m_rWall[i-1]);
        }
        mfr_in = area_i*m_rho(0,0)*u(x,0,0)*Y(x,nh3_index,0,0)/m_wt[nh3_index];
        mfr_out = area_i*m_rho(j,0)*u(x,j,0)*Y(x,nh3_index,j,0)/m_wt[nh3_index];
        sumDenom += mfr_in;
        sumNum += (mfr_in - mfr_out);
    }
    m_conversion = 100*sumNum/sumDenom;
}

/*
    Write ammonia conversion to a file
*/
void bvpQ2D::writeConversion()
{
    //Write conversion to the file
    std::ofstream outfile;
    outfile.open("ammonia_conversion.dat", ios::app);
    outfile<<"\n";
    outfile<<T_in<<"\t"<<p_out<<"\t"<<m_SCCM<<"\t"<<m_conversion;
    outfile.close();								// mol/s/m^3
}

/*
    Resize geometrical arrays based on whether support is present or not.
*/
void bvpQ2D::resizeGeometryArrays(int n)
{
    m_porosity.resize(n, 0.0);
    m_area2Vol.resize(n, 0.0);
    m_tortuosity.resize(n, 0.0);
    m_Dp.resize(n, 0.0);
    m_PoreRad.resize(n, 0.0);

    // Micro-pores
    m_microPorosity.resize(n, 0.0);
    m_microPoreRad.resize(n, 0.0);

    m_Bg.resize(n, 0.0);
    m_gamma.resize(n, 0.0);

    c_offset_U.resize(n, 0.0);
    c_offset_T.resize(n, 0.0);
    c_offset_Ts.resize(n, 0.0);
    c_offset_Y.resize(n, 0.0);
    c_offset_theta.resize(n, 0.0);
}

/*
    Calculate the flux at the lumen and the support interface.
*/
void bvpQ2D::interfaceFlux(double* x, size_t j0, size_t j1)
{
    double avgD = 0;
    for (size_t i = 0; i < m_radialPoints; i++) {
        if (m_tranModel == m_FICK)
        {
            if (m_do_multicomponent) {
                for (size_t j = j0; j <= j1; j++) {
                    if(i == m_radialPoints - 1)
                    {
                        for (size_t k = 0; k < m_nsp; k++) 
                        {
                            m_interfaceFlux(k,cellIndex(j,i)) = 0;
                        }
                    } 
                    else 
                    {
                        for (size_t k = 0; k < m_nsp; k++) 
                        {
                            double sum = 0.0;
                            for (size_t m = 0; m < m_nsp; m++) {
                                sum += m_wt[m] * m_multidiff_radial[mindex(k,m,cellIndex(j,i))] * (X(x,m,j,i + 1)-X(x,m,j,i));
                            }
                            m_interfaceFlux(k,cellIndex(j,i)) = sum * m_diff_radial[k+cellIndex(j,i)*m_nsp] / (m_r[i+1] - m_r[i]);
                        }
                    }
                }
            } else 
            { //Fickian diffusion
                for (size_t j = j0; j <= j1; j++) {
                    if(i == m_radialPoints - 1)
                    {
                        for (size_t k = 0; k < m_nsp; k++) 
                        {
                            m_interfaceFlux(k,cellIndex(j,i)) = 0;
                        }
                    } else {
                        double sum = 0.0;
                        double wtm = m_wtm(j,i);
                        double rho = m_rho(j,i);
                        double dr = m_r[i+1] - m_r[i];
                        for (size_t k = 0; k < m_nsp; k++) {
                            m_interfaceFlux(k,cellIndex(j,i)) = m_wt[k]*(rho*m_diff_radial[k+m_nsp*cellIndex(j,i)]/wtm);
                            m_interfaceFlux(k,cellIndex(j,i)) *= (X(x,k,j,i) - X(x,k,j,i + 1))/dr;
                            sum -= m_interfaceFlux(k,cellIndex(j,i));
                        }
                    }
                }
            }
        }
        /*else if (m_tranModel == m_DGM)
        {
            // Dusty gas model
            for (size_t j = j0; j <= j1; j++)
            {
                setGas(x,j); 
                m_tranDGM->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]); // unit: m2/s   
                m_thermo->saveState(state_prev);

                setGas(x,j+1); 
                m_thermo->saveState(state_curr);

                // Solve DGM transport to get mass fluxes (fluxes saved to m_flux)
                m_tranDGM->getMolarFluxes(&state_prev[0], &state_curr[0], (m_z[j+1] - m_z[j]), &m_flux(0,j)); 
                // Convert to mass fluxes
                for (size_t k = 0; k < m_nsp; k++)
                {
                    m_flux(k, j) = m_wt[k] * m_flux(k, j);                  // unit: kg/m2/s
                }
            }
        }*/
    }
}

/*
    Set thermodynamic state of the gas-phaseat the lumen-support interface
*/
void bvpQ2D::setGasAtInterface(const doublereal* x, size_t j, size_t i)
{
    double TAvg = 0.5*(T(x,j,i) + T(x,j,i+1));
    double pAvg = 0.5*(m_press(j,i) + m_press(j,i+1));
    int ind1, ind2;
    vector_double hk_RT(m_nsp,0.0);
    if(i < m_lumenPoints)
    {
        ind1 = i*nEq_lumen + c_offset_Y[0];
        ind2 = (i+1)*nEq_lumen + c_offset_Y[0];
    }
    else
    {
        ind1 = nv_lumen + (i - m_lumenPoints) * nEq_support + c_offset_Y[1];
        ind2 = nv_lumen + (i + 1 - m_lumenPoints) * nEq_support + c_offset_Y[1];
    }
    const doublereal* yyj_1 = x + m_nv*j + ind1;
    const doublereal* yyj_2 = x + m_nv*j + ind2;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj_1[k] + yyj_2[k]);
    }
    // Update m_thermo
    if(i < m_lumenPoints)
    {
        m_thermo->setTemperature(TAvg);
        m_thermo->setMassFractions_NoNorm(m_ybar.data());
        m_thermo->setPressure(pAvg);
        hk_RT = m_thermo->enthalpy_RT_ref();
        if(m_do_multicomponent)
        {
            m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff_radial[mindex(0,0,cellIndex(j,i))]);
        }
        
    } else {
        m_thermo_2->setTemperature(TAvg);
        m_thermo_2->setMassFractions_NoNorm(m_ybar.data());
        m_thermo_2->setPressure(pAvg);
        hk_RT = m_thermo_2->enthalpy_RT_ref();
        if(m_do_multicomponent)
        {
            m_trans_2->getMultiDiffCoeffs(m_nsp_2, &m_multidiff_radial[mindex(0,0,cellIndex(j,i))]);
        }
    }
    for(int k = 0;  k < m_nsp; k++)
    {
        hk_interface(k,cellIndex(j,i)) = hk_RT[k] * GasConstant * TAvg;
    }
}

/*
    Setup the grid in the radial direction
*/
void bvpQ2D::setupRadialGrid()
{
    m_r.resize(m_totRadialPoints, 0.0);
    m_rWall.resize(m_totRadialPoints, 0.0);

    // Lumen side
    double deltaR = m_Rin/m_lumenPoints;    
    m_r[0] = deltaR/2;
    m_rWall[0] = deltaR;
    for (size_t j = 1; j < m_lumenPoints; j++) {
        m_r[j] = m_r[j-1] + deltaR;
        m_rWall[j] = m_rWall[j - 1] + deltaR; // Upper wall for cell (j-1) or Lumen-support interface
    }

    if(m_solveSupport)
    {
        deltaR = (m_Rout - m_Rin)/m_supportPoints;
        m_r[m_lumenPoints] = m_Rin + deltaR/2;
        m_rWall[m_lumenPoints] = m_Rin + deltaR;
        for (size_t j = 1; j < m_supportPoints; j++){
            m_r[j + m_lumenPoints] = m_r[j-1 + m_lumenPoints] + deltaR;
            m_rWall[j + m_lumenPoints] = m_rWall[j - 1 + m_lumenPoints] + deltaR; // Upper wall for cell (j-1)
        }
    }
    if(m_do_channel)
    {
        int ind = m_lumenPoints + m_supportPoints;
        m_r[ind] = (m_rChannel - m_Rout)/2;
        m_rWall[ind] = m_rChannel;
    }
}

/*
    Calculate the diffusion flux at the interface boundary
*/
void bvpQ2D::flux_LumenSupport(double* x, size_t j0, size_t j1)
{
    int i = m_lumenPoints-1;
    double low = m_rWall[i] - m_r[i];
    double high = m_r[i+1] -  m_rWall[i];
    double avgD = 0;
    if (m_tranModel == m_FICK)
    {
        if (m_do_multicomponent) {
            for (size_t j = j0; j <= j1; j++) {
                setGas(x,j,i);
                m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff_radial[mindex(0,0,cellIndex(j,i))]);
                setGasSupport(x,j,i+1);
                m_trans_2->getMultiDiffCoeffs(m_nsp_2, &m_multidiff_radial[mindex(0,0,cellIndex(j,i+1))]);
                double rho = 0.5 * (m_rho(j,i) + m_rho(j,i+1));
                double wtm = 0.5 * (m_thermo->meanMolecularWeight() + m_thermo_2->meanMolecularWeight());
                double m_diff_avg;
                for (size_t k = 0; k < m_nsp; k++)
                {
                    double sum = 0.0;
                    for (size_t m = 0; m < m_nsp; m++) {
                        double m_multidiff_avg = (low*m_multidiff[mindex(k,m,cellIndex(j,i))] + high*m_multidiff[mindex(k,m,cellIndex(j,i+1))])/(low+high);
                        sum += m_wt[m] * m_multidiff_avg * (X(x,m,j,i + 1)-X(x,m,j,i));
                    }

                    // Use m_diff as storage for the factor outside the summation
                    m_diff_avg = m_wt[k] * rho / (wtm*wtm);  //unit: kmol/m3

                    m_interfaceFlux(k,cellIndex(j,i)) = sum * m_diff_avg / (m_r[i+1] - m_r[i]);
                }
            }
        } else 
        {
            //Fickian diffusion
            for (size_t j = j0; j <= j1; j++) {
                setGas(x,j,i);
                setGasSupport(x,j,i+1);

                double sum = 0.0;
                double wtm = (low*m_wtm(j,i) + high*m_wtm(j,i+1))/ (low+high);
                double rho = (low*m_rho(j,i) + high*m_rho(j,i+1))/ (low+high);
                double dr = m_r[i+1] - m_r[i];
                for (size_t k = 0; k < m_nsp; k++) {
                    m_interfaceFlux(k,cellIndex(j,i)) = m_wt[k]*(rho*m_diff_radial[k+m_nsp*cellIndex(j,i)]/wtm);
                    m_interfaceFlux(k,cellIndex(j,i)) *= (X(x,k,j,i) - X(x,k,j,i + 1))/dr;
                    sum -= m_interfaceFlux(k,cellIndex(j,i));
                }
                // correction flux to ensure that \sum_k Y_k V_k = 0.
                for (size_t k = 0; k < m_nsp; k++) {
                    //m_interfaceFlux(k,cellIndex(j,i)) += sum*Y(x,k,j,i);
                }
            }
        }
    }
        /*else if (m_tranModel == m_DGM)
        {
            // Dusty gas model
            for (size_t j = j0; j <= j1; j++)
            {
                setGas(x,j); 
                m_tranDGM->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]); // unit: m2/s   
                m_thermo->saveState(state_prev);

                setGas(x,j+1); 
                m_thermo->saveState(state_curr);

                // Solve DGM transport to get mass fluxes (fluxes saved to m_flux)
                m_tranDGM->getMolarFluxes(&state_prev[0], &state_curr[0], (m_z[j+1] - m_z[j]), &m_flux(0,j)); 
                // Convert to mass fluxes
                for (size_t k = 0; k < m_nsp; k++)
                {
                    m_flux(k, j) = m_wt[k] * m_flux(k, j);                  // unit: kg/m2/s
                }
            }
        }*/
}

/*
    Get the index of the excess species at each mesh point
*/
void bvpQ2D::getkExcess(const double* x)
{
    int ind, i, j;
    for(j = 0; j < m_axialPoints; j++)
    {
        for(i = 0; i < m_radialPoints; i++)
        {
            if(i < m_lumenPoints)
            {
                ind = i * nEq_lumen;
                const double* Y = x + m_nv*j + ind + c_offset_Y[0];
                m_kExcess[cellIndex(i,j)] = distance(Y, max_element(Y, Y + m_nsp));
            } else{
                ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
                const double* Y = x + m_nv*j + ind + c_offset_Y[1];
                m_kExcess[cellIndex(i,j)] = distance(Y, max_element(Y, Y + m_nsp));
            }
        }
    }
}

XML_Node& bvpQ2D::save(XML_Node& o, const doublereal* const sol)
{
    Array2D soln(m_nv, m_axialPoints, sol + loc());
    XML_Node& flow = Domain1D::save(o, sol);
    
    for(size_t iGrid =0; iGrid<m_totRadialPoints; iGrid++)
    {
        int ind, nsp, nsps, part;
        if(iGrid < m_lumenPoints)
        {
            ind = iGrid*nEq_lumen;
            nsp = m_nsp;
            nsps = m_nspSurf;
            part = 0;
        } else {
            ind = nv_lumen + (iGrid - m_lumenPoints) * nEq_support;
            nsp = m_nsp_2;
            nsps = m_nspSurf_2;
            part = 1;
        }
    
        doublereal *pres = new doublereal[m_axialPoints];
        for (int k = 0; k < m_axialPoints; k++)
        {
            pres[k] = m_press(iGrid, k);
        }

        XML_Node& gv = flow.addChild("grid_data");
        addFloat(flow, "pressure", *pres, "Pa", "pressure");

        addFloatArray(gv,"z",m_z.size(), m_z.data(), "m","length");
        vector_double x(soln.nColumns());
        soln.getRow(c_offset_U[part] + ind, x.data());
        addFloatArray(gv,"u",x.size(),x.data(),"m/s","velocity");

        if(m_solveRadialVel)
        {
            soln.getRow(c_offset_V[part], x.data());
            addFloatArray(gv,"V",x.size(),x.data(),"m/s","radialVel");
        }

        soln.getRow(c_offset_T[part], x.data());
        addFloatArray(gv,"T",x.size(),x.data(),"K","gas temperature");

        soln.getRow(c_offset_Ts[part], x.data());
        addFloatArray(gv,"Ts",x.size(),x.data(),"K","surf temperature");

        for (size_t k = 0; k < nsp; k++) {
            soln.getRow(c_offset_Y[part]+k, x.data());
            if(part == 0) {
                addFloatArray(gv,m_thermo->speciesName(k), x.size(),x.data(),"","massFraction");
            } else {
                addFloatArray(gv,m_thermo_2->speciesName(k), x.size(),x.data(),"","massFraction");
            }
        }
        for (size_t k = 0; k < nsps; k++) {
            soln.getRow(c_offset_theta[part]+k, x.data());
            if(part == 0) {
                addFloatArray(gv,m_surf->speciesName(k), x.size(),x.data(),"","surfCoverages");
            } else {
                addFloatArray(gv,m_surf_2->speciesName(k), x.size(),x.data(),"","surfCoverages");
            }
        }
    }

    XML_Node& ref = flow.addChild("refine_criteria");
    addFloat(ref, "ratio", refiner().maxRatio());
    addFloat(ref, "slope", refiner().maxDelta());
    addFloat(ref, "curve", refiner().maxSlope());
    addFloat(ref, "prune", refiner().prune());
    addFloat(ref, "grid_min", refiner().gridMin());
    return flow;
}

void bvpQ2D::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    int ind, iGrid;
    cout<<"\n Reading the saved solution...";
    Domain1D::restore(dom, soln, loglevel);

    std::vector<XML_Node*> str = dom.getChildren("string");
    for (size_t istr = 0; istr < str.size(); istr++) {
        const XML_Node& nd = *str[istr];
        writelog(nd["title"]+": "+nd.value()+"\n");
    }

    double pp = getFloat(dom, "pressure", "pressure");
    //setPressure(pp);
    std::vector<XML_Node*> d = dom.child("grid_data").getChildren("floatArray");
    vector_double x;
    size_t np = 0;
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        if (nm == "z") {
            getFloatArray(fa,x,false);
            np = x.size();
            if (loglevel >= 2) {
                writelog("Grid contains {} points.\n", np);
            }
            setupGrid(np, x.data());
        }
    }
    for (size_t n = 0; n < d.size(); n++) {
        const XML_Node& fa = *d[n];
        string nm = fa["title"];
        getFloatArray(fa,x,false);
        if (nm == "u") {
            debuglog("axial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("bvpQ2D::restore",
                                   "axial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_U[0],j)] = x[j];
            }
        } else if (nm == "z") {
            ; // already read grid
        } else if (nm == "V" && m_solveRadialVel) {
            debuglog("radial velocity   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("bvpQ2D::restore",
                                   "radial velocity array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_V[0],j)] = x[j];
            }
        } else if (nm == "T") {
            debuglog("gas temperature   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("bvpQ2D::restore",
                                   "temperature array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_T[0],j)] = x[j];
            }
        } else if (nm == "Ts") {
            debuglog("surf temperature   ", loglevel >= 2);
            if (x.size() != np) {
                throw CanteraError("bvpQ2D::restore",
                                   "temperature array size error");
            }
            for (size_t j = 0; j < np; j++) {
                soln[index(c_offset_Ts[0],j)] = x[j];
            }
        } else if (m_thermo->speciesIndex(nm) != npos) {
            debuglog(nm+"   ", loglevel >= 2);
            if (x.size() == np) {
                size_t k = m_thermo->speciesIndex(nm);
                for (size_t j = 0; j < np; j++) {
                    soln[index(k+c_offset_Y[0],j)] = x[j];
                }
            }
        } else if (m_surf->speciesIndex(nm) != npos) {
            debuglog(nm+"   ", loglevel >= 2);
            if (x.size() == np) {
                size_t k = m_surf->speciesIndex(nm);
                for (size_t j = 0; j < np; j++) {
                    soln[index(k+c_offset_theta[0],j)] = x[j];
                }
            }
        } else {
            // error
        }
    }

    if (dom.hasChild("refine_criteria")) {
        XML_Node& ref = dom.child("refine_criteria");
        refiner().setCriteria(getFloat(ref, "ratio"), getFloat(ref, "slope"),
                              getFloat(ref, "curve"), getFloat(ref, "prune"));
        refiner().setGridMin(getFloat(ref, "grid_min"));
    }
}

void bvpQ2D::interpolate(double* x)
{
    int i, ind;
    double z_current, zm, zp, slope;
    // Get the length of the z vector from the initial data
    nSteady = m_initData.at(0).second.size(); // m_initData[0].size();  
    cout<<"\n Number of grid points in the initial guess = "<<nSteady;
    for(int i = 0; i < m_axialPoints; i++)
    {
        // First calculate the range (i, i + 1) in which m_z[iGrid] lies
		z_current = m_z[i];
		ind = getMeshIndex(z_current);

        //Interpolate the initial solution vector m_y0
        zm = m_initData.at(0).second.at(ind);
        if(ind != (nSteady-1))
        {
            zp = m_initData.at(0).second.at(ind+1);
            slope = (z_current - zm) / (zp - zm);
        }

        // 
        int jVar = m_initData.size();
        int jInd;
        double ym, yp;
        std::string name = {};
        for(int j = 1; j < jVar; j++)
        {
            name = m_initData.at(j).first;
            jInd = componentIndex(name);
            ym = m_initData.at(j).second.at(ind);
            if(ind != (nSteady-1))
            {
                yp = m_initData.at(j).second.at(ind+1);
                x[index(jInd, i)] = ym + slope * (yp - ym);
            } else {
                x[index(jInd, i)] = ym;
            }
        }

        // Assign the same values along the radial direction
        for (size_t j = 0; j < m_lumenPoints; j++)
        {
            for(int k = 0; k < nEq_lumen; k++)
            {
                x[index(j*nEq_lumen + k, i)] = x[index(k, i)];
            }
        }
    }
}

int bvpQ2D::getMeshIndex(double z)
{
	int i, ind;
	if (z == 0)
	{
		ind = 0;
	}
	else {
		for (i = 0; i < (nSteady - 1); i++)
		{
			if (m_initData.at(0).second.at(i+1) >= z)
			{
				ind = i;
				break;
			}
            if (z >= m_initData.at(0).second.at(nSteady-1))
            {
                ind = nSteady-1;
            }
		}
	}
	return ind;
}

void bvpQ2D::evalAnnularChannel(double* x, double* rsd, int* diag, double rdt, size_t jmin,size_t jmax)
{
    int ind = m_NeqTotal - nEq_channel;
    double memFluxTerm = 0, convec = 0;
    double rhou, rho, dz;
    double drhou_dt, drhou2dz, dpdz, stressTerm;
    int i = m_totRadialPoints - 1;  // Last radial point which is inside the annular channel
    int r_in = m_rWall[i-1];
    
    // Inlet of the annular channel
    rsd[index(c_offset_ch_rho + ind, 0)] = m_rho(0,i) - rho0_ch;
    diag[index(c_offset_ch_rho + ind, 0)] = 0;

    rsd[index(c_offset_ch_U + ind, 0)] = u(x,0,i) - u0_ch;
    diag[index(c_offset_ch_U + ind, 0)] = 0;

    rsd[index(c_offset_ch_T + ind, 0)] = T(x,0,i) - T_in;
    diag[index(c_offset_ch_T + ind, 0)] = 0;

    for(int j = 1; j < m_axialPoints; j++)
    {
        rho = m_rho(j,i);
        rhou = rho * u(x,j,i);
        dz = (m_z[j] - m_z[j-1]);
        
        // Continuity equation
        convec = (rhou - m_rho(j-1,i)*u(x,j-1,i))/dz;
        memFluxTerm = getMembraneFlux(x, j) * m_rWall[i-1];
        memFluxTerm *= 2/(m_rChannel*m_rChannel - r_in*r_in);
        double drho_dt = memFluxTerm - convec;
        rsd[index(c_offset_ch_rho + ind, j)] = drho_dt - rdt * (rho - rho_prev(j,i));
        diag[index(c_offset_ch_rho + ind, j)] = 1;

        // Momentum equation
        drhou2dz = (rhou*u(x,j,i) - m_rho(j-1,i)*u(x,j-1,i)*u(x,j-1,i))/dz;
        dpdz = (m_press(j,i) - m_press(j-1,i))/dz;
        stressTerm = getFrictionFactor(x,j,i)* 0.5 * rhou * u(x,j,i);
        
        rsd[index(c_offset_ch_U + ind, j)] = (- dpdz - stressTerm - drhou2dz - u(x,j,i)*drho_dt)/m_rho(j,i);
        rsd[index(c_offset_ch_U + ind, j)] -= rdt * (u(x,j,i) - u_prev(j,i));
        diag[index(c_offset_ch_U + ind, j)] = 1;

        // Energy equation
        if(m_do_energy)
        {
            double dTdz = (T(x,j,i) - T(x,j-1,i))/dz;
            double q_conv = m_hcoeff * (T(x,j,i) - T_wall);
            setChannelState(x,j,i);
            double cp = m_thermo->cp_mole();
            double q_mem = memFluxTerm * m_thermo->enthalpy_mass();
            rsd[index(c_offset_ch_T + ind, j)] = (q_conv - q_mem - rhou*cp*dTdz)/(rho*cp);
            rsd[index(c_offset_ch_T + ind, j)] -= rdt * (T(x,j,i) - T_prev(j,i));
            diag[index(c_offset_ch_T + ind, j)] = 1;
        } else {
            rsd[index(c_offset_ch_T + ind, j)] = T(x,j,i) - T(x,j,i-1); //dTdt;
            diag[index(c_offset_ch_T + ind, j)] = 0;
        }

        if(isnan(rsd[index(c_offset_ch_T + ind, j)]))
        {
            throw CanteraError("evalAnnularChannel::T", "NaN error");
        }
        if(isnan(rsd[index(c_offset_ch_rho + ind, j)]))
        {
            throw CanteraError("evalAnnularChannel::rho", "NaN error");
        }
        if(isnan(rsd[index(c_offset_ch_U + ind, j)]))
        {
            throw CanteraError("evalAnnularChannel::u", "NaN error");
        }
    }
}

void bvpQ2D::setChannelState(const doublereal* x, size_t j, size_t i)
{
    double Temp_ch, rho_ch; 
    int ind = m_NeqTotal - nEq_channel;
    Temp_ch = x[index(c_offset_ch_T + ind, j)];
    rho_ch = x[index(c_offset_ch_rho + ind, j)];
    vector_double xx(m_nsp, 0.0);
    xx[m_memSpeciesIndex] = 1.00;
    m_thermo->setState_TRX(Temp_ch, rho_ch, xx.data());
}

double bvpQ2D::getFrictionFactor(double* x, int j, int i)
{
    double xi, Re, Ref, Dh, var1, var2;
    xi = m_rWall[i-1]/m_rChannel;
    Dh = 2 * m_rChannel * (1-xi);
    Re = m_rho(j,i)*u(x,j,i)*Dh/m_visc(j,i);

    var1 = (1 - xi)*(1 + xi)*log(xi);
    var2 = (1 + pow(xi,4)) - 2*xi*xi + var1* (1 + xi*xi);
    Ref = 16*pow((1-xi),2)*var1/var2;
    return (Ref/Re);
}