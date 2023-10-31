#ifndef BVP_DGM_H
#define BVP_DGM_H

#include <fstream>
#include "BoundaryValueProblem.h"
#include "canteraInput.h"

#include <cantera/oneD/Domain1D.h>
#include <cantera/base/Array.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/thermo/SurfPhase.h>
#include <cantera/kinetics/InterfaceKinetics.h>
#include <cantera/numerics/FuncEval.h>
#include <cantera/transport/DustyGasTransport.h>
#include <cantera/transport.h>
#include <cantera/numerics/Func1.h>

class Cantera::Transport;
using namespace Cantera;
using namespace std;
using namespace BVP;
class Cantera::Func1;

class bvpQ2D : public BVP::BoundaryValueProblem
{
public:
    //--------------------------------
    // construction and destruction
    //--------------------------------

    //! Create a new flow domain.
    //! @param filename Input filename
    //! @param nsp Number of species.
    //! @param points Initial number of grid points
    //friend class user_func;

    // Constructor and destructor functions
	bvpQ2D(std::string filename, inputClass *ptr, int nv, int points, double Length, int energyFlag);
	//virtual ~bvpQ2D();

    virtual void setupGrid(size_t n, const doublereal* z);
    void setupRadialGrid();
    //! Change the grid size. Called after grid refinement.
    virtual void resize(size_t components, size_t points);
    void resizeGeometryArrays(int n);

    //! Set the gas and surface object state to be consistent with the solution at point j.
    void setGas(const doublereal* x, size_t j, size_t i);
    void setSurf(const doublereal* x, size_t j, size_t i);
    void setGasSupport(const doublereal* x, size_t j, size_t i);
    void setSurfSupport(const doublereal* x, size_t j, size_t i);

    //! Set the gas state to be consistent with the solution at the midpoint
    //! between j and j + 1.
    void setGasAtMidpoint(const doublereal* x, size_t j, size_t i);
    void setGasAtMidpointSupport(const doublereal* x, size_t j,size_t i);

    //! Update the properties (thermo, transport, and diffusion flux).
    //! This function is called in eval after the points which need
    //! to be updated are defined.
    virtual void updateProperties(size_t jg, double* x, size_t jmin, size_t jmax);

    //! Update the diffusive mass fluxes.
    virtual void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1);

    // Calculate Knudsen diffusion and effective diffusion coefficient
    Array2D calculateKnudsenDiffCoeff(double T, int m);
    void calculateDiffCoeff(double T, int j, int i, int n);    // n = 0 for axial, n = 1 for radial

    //! Update the transport properties at grid points in the range from `j0`
    //! to `j1`, based on solution `x`.
    virtual void updateTransport(doublereal* x, size_t j0, size_t j1);
    virtual void updateTransportRadial(doublereal* x, size_t j0, size_t j1);

    /**
     * Update the thermodynamic properties from point j0 to point j1
     * (inclusive), based on solution x.
     */
    void updateThermo(doublereal* x, size_t j0, size_t j1) {
        getPressure(x);
        for (size_t i = 0; i < m_lumenPoints; i++) {
            for (size_t j = j0; j <= j1; j++) {
                setGas(x,j,i);
                setSurf(x,j,i);
                m_rho(j,i) = m_thermo->density();
                m_wtm(j,i) = m_thermo->meanMolecularWeight();
                m_press(j,i) = m_thermo->pressure();
                m_cp(j,i) = m_thermo->cp_mass();
                vector_double xmole(m_nsp, 0.0);
                m_thermo->getMoleFractions(&xmole[0]);
                for(int k = 0; k< m_nsp; k ++)
                {
                    m_xMole(k, cellIndex(j,i)) = xmole[k];
                }
            }
        }

        if(m_solveSupport)
        {
            for (size_t i = m_lumenPoints; i < (m_lumenPoints + m_supportPoints); i++) 
            {
                getPressure(x);
                for (size_t j = j0; j <= j1; j++) {
                    setGasSupport(x,j,i);
                    setSurfSupport(x,j,i);
                    m_rho(j,i) = m_thermo_2->density();
                    m_wtm(j,i) = m_thermo_2->meanMolecularWeight();
                    m_press(j,i) = m_thermo_2->pressure();
                    m_cp(j,i) = m_thermo_2->cp_mass();
                    vector_double xmole(m_nsp_2, 0.0);
                    m_thermo_2->getMoleFractions(&xmole[0]);
                    for(int k = 0; k< m_nsp_2; k++)
                    {
                        m_xMole(k, cellIndex(j,i)) = xmole[k];
                    }
                }
            }
        }

        // Annular region
        if(m_do_channel)
        {
            size_t i = m_totRadialPoints - 1;
            for (size_t j = j0; j <= j1; j++) {
                setChannelState(x,j,i);
                m_rho(j,i) = m_thermo->density();
                m_press(j,i) = m_thermo->pressure();
                //m_cp(j,i) = m_thermo->cp_mass();
                //m_visc(j,i) = m_trans->viscosity();
                //m_tcon(j,i) = m_trans->thermalConductivity();
            }
        }
    }

    //DGM
    int setDGMProperties();

    //Heat fluxes used in energy equation
    double grad_qGas(const doublereal* x, size_t j, size_t i);
    double grad_qSolid(const doublereal* x, size_t j, size_t i);
    double qdotConv(const doublereal* x, size_t j, size_t i);
    double qdotEnv(const doublereal* x, size_t j, size_t i);
    double qdotSurf(double* x, size_t j,size_t i);
    double calculateNu(const doublereal* x, size_t j, size_t i);

    // Calculate solid phase conductivity
    double getConductivitySolid(double, size_t);
    double solidConductivity(double);

    //Calculate pressure using Darcy flow
    void getPressure(const double* x);

    // NH3 conversion
    void getNH3Conversion(double* x);
    void writeConversion();

    // Read restart file
    XML_Node& save(XML_Node& o, const doublereal* const sol);
    void restore(const XML_Node& dom, double* soln, int loglevel);

    // Inital guess from the diffusion-free steady solution
    void interpolate(double* x);
    int getMeshIndex(double);

    // Annular region
    virtual void evalAnnularChannel(double* x, double* rsd, int* diag, double rdt, size_t jmin,size_t jmax);
    void setChannelState(const doublereal* x, size_t j, size_t i);
    double getFrictionFactor(double*, int, int); // Friction factor calculation for annular region

    //! @name Solution components
    //! @{

    doublereal u(const doublereal* x, size_t j, size_t i) const {
        int ind;
        if(i < m_lumenPoints)
        {
            return x[index(c_offset_U[0] + i*nEq_lumen, j)];
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return x[index(c_offset_U[1] + ind, j)];
        } else if((i == m_totRadialPoints-1)&& m_do_channel)
        {
            // Annular channel
            ind = m_NeqTotal - nEq_channel;
            return x[index(c_offset_ch_U + ind, j)];
        } else {
            throw CanteraError("bvpQ2D::u", "Index out of bounds for U");
        }
    }

    doublereal u_prev(size_t j, size_t i) const {
        int ind;
        if(i < m_lumenPoints)
        {
            return prevSoln(c_offset_U[0] + i*nEq_lumen, j);
        } else if (i < (m_lumenPoints + m_supportPoints))
        {
            ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return prevSoln(c_offset_U[1] + ind, j);
        } else if((i == m_totRadialPoints-1)&& m_do_channel)
        {
            // Annular channel
            ind = m_NeqTotal - nEq_channel;
            return prevSoln(c_offset_ch_U + ind, j);
        }
        else {
            throw CanteraError("bvpQ2D::u_prev", "Index out of bounds for u_prev");
        }
    }

    /*doublereal rho_u(const doublereal* x, size_t j) {
        return m_rho[j]*x[index(c_offset_U[0], j)];
    }*/
    
    doublereal T(const doublereal* x, size_t j, size_t i) const {
        int ind;
        if(i < m_lumenPoints)
        {
            return x[index(c_offset_T[0] + i*nEq_lumen, j)];
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return x[index(c_offset_T[1] + ind, j)];
        } else if((i == m_totRadialPoints-1)&& m_do_channel)
        {
            // Annular channel
            ind = m_NeqTotal - nEq_channel;
            return x[index(c_offset_ch_T + ind, j)];
        } else {
            throw CanteraError("bvpQ2D::T", "Index out of bounds for T");
        }
    }

    
    doublereal Ts(const doublereal* x, size_t j, size_t i) const {
        if(i < m_lumenPoints)
        {
            return x[index(c_offset_Ts[0] + i*nEq_lumen, j)];
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return x[index(c_offset_Ts[1] + ind, j)];
        } else {
            throw CanteraError("bvpQ2D::Ts", "Index out of bounds for Ts");
        }
    }

    doublereal T_prev(size_t j, size_t i) const {
        int ind;
        if(i < m_lumenPoints)
        {
            return prevSoln(c_offset_T[0] + i*nEq_lumen, j);
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return prevSoln(c_offset_T[1] + ind, j);
        } else if((i == m_totRadialPoints-1)&& m_do_channel)
        {
            // Annular channel
            ind = m_NeqTotal - nEq_channel;
            return prevSoln(c_offset_ch_T + ind, j);
        } else {
            throw CanteraError("bvpQ2D::T_prev", "Index out of bounds for T_prev");
        }
    }

    doublereal Ts_prev(size_t j, size_t i) const {
        if(i < m_lumenPoints)
        {
            return prevSoln(c_offset_Ts[0] + i*nEq_lumen, j);
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return prevSoln(c_offset_Ts[1] + ind, j);
        } else {
            throw CanteraError("bvpQ2D::Ts_prev", "Index out of bounds for Ts_prev");
        }
    }

    doublereal p(const doublereal* x, size_t j, size_t i) const {
        return m_press(j,i);
    }
    
    doublereal Y(const doublereal* x, size_t k, size_t j, size_t i) const {
        if(i < m_lumenPoints)
        {
            return x[index(c_offset_Y[0] + i*nEq_lumen + k, j)];
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return x[index(c_offset_Y[1] + ind + k, j)];
        } else{
            throw CanteraError("bvpQ2D::Y", "Index out of bounds for Y");
        }
    }

    doublereal& Y(doublereal* x, size_t k, size_t j, size_t i) {
        if(i < m_lumenPoints)
        {
            return x[index(c_offset_Y[0] + i*nEq_lumen + k, j)];
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return x[index(c_offset_Y[1] + ind + k, j)];
        } else{
            throw CanteraError("bvpQ2D::Y", "Index out of bounds for Y");
        }
    }

    doublereal Y_prev(size_t k, size_t j, size_t i) const {
        if(i < m_lumenPoints)
        {
            return prevSoln(c_offset_Y[0] + k + i*nEq_lumen, j);
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return prevSoln(c_offset_Y[1] + k + ind, j);
        } else {
            throw CanteraError("bvpQ2D::Y_prev", "Index out of bounds for Y_prev");
        }
    }

    doublereal& surfY(doublereal* x, size_t k, size_t j, size_t i) {
        if(i < m_lumenPoints)
        {
            return x[index(c_offset_theta[0] + i*nEq_lumen + k, j)];
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return x[index(c_offset_theta[1] + ind + k, j)];
        } else {
            throw CanteraError("bvpQ2D::surfY", "Index out of bounds for surfY");
        }
    }

    doublereal surfY(const doublereal* x, size_t k, size_t j, size_t i) const {
        if(i < m_lumenPoints)
        {
            return x[index(c_offset_theta[0] + i*nEq_lumen + k, j)];
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return x[index(c_offset_theta[1] + ind + k, j)];
        } else {
            throw CanteraError("bvpQ2D::surfY", "Index out of bounds for surfY");
        }
    }

    doublereal surfY_prev(size_t k, size_t j, size_t i) const {
        if(i < m_lumenPoints)
        {
            return prevSoln(c_offset_theta[0] + k + i*nEq_lumen, j);
        } else if(i < (m_lumenPoints+m_supportPoints))
        {
            int ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
            return prevSoln(c_offset_theta[1] + k + ind, j);
        } else {
            throw CanteraError("bvpQ2D::surfY_prev", "Index out of bounds for surfY_prev");
        }
    }

    doublereal X(const doublereal* x, size_t k, size_t j, size_t i) const {
        return m_wtm(j,i)*Y(x,k,j,i)/m_wt[k];
    }

    /*doublereal flux(size_t k, size_t j) const {
        return m_flux(k, j);
    }*/

    /*doublereal dYdz(const doublereal* x, size_t k, size_t j, size_t n) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (Y(x,k,jloc,n) - Y(x,k,jloc-1,n))/m_dz[jloc-1];
    }*/

    /*doublereal dpdz(const doublereal* x, size_t j, size_t i) const {
        // Upwind differences
        return (p(x,j + 1,i) - p(x,j,i))/m_dz[j];
    }*/

    /*doublereal dTdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (T(x,jloc) - T(x,jloc-1))/m_dz[jloc-1];
    }*/

    /*doublereal dTsdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (Ts(x,jloc) - Ts(x,jloc-1))/m_dz[jloc-1];
    }*/

    doublereal divHeatFlux(const doublereal* x, size_t j, size_t i) const {
        double qGas_p, qGas_m, cond_p, cond_m;
        //setGasAtMidpoint(x,j);
        cond_p = m_tcon(j,i);
        //setGasAtMidpoint(x,j-1);
        cond_m = m_tcon(j-1,i);

        qGas_m = -cond_m*(T(x,j,i) - T(x,j-1,i))/(z(j) - z(j - 1));
        if(j == m_axialPoints-1)
        {
            qGas_p = 0;
            return (-qGas_m/(z(j) - z(j-1)));
        }
        else {
            qGas_p = -cond_p*(T(x,j+1,i) - T(x,j,i))/(z(j+1) - z(j));
            return 2.0*(qGas_p - qGas_m)/(z(j+1) - z(j-1));
        }
    }

    size_t cellIndex(size_t j, size_t i) const
    {
        // jth cell in axial direction and ith cell in radial directon
        return (j*m_totRadialPoints + i);
    }

    //! Write the initial solution estimate into array x.
    virtual void _getInitialSoln(double* x);

    // Membrane flux as the boundary condition
    double getMembraneFlux(double* x, size_t j);
    void interfaceFlux(double* x, size_t j0, size_t j1);
    void setGasAtInterface(const doublereal* x, size_t j, size_t i);
    void flux_LumenSupport(double* x, size_t j0, size_t j1);

    // Kozeny-Carman permeability
    void KozenyCarmanPermeability(int);

    // Get coverages of point j
    vector_fp getCov(int j)
    {
        vector_fp cov(m_nspSurf, 0.0);
        for(int i = 0; i<m_nspSurf; i++)
        {
            cov[i] = m_cov(i,j);
        }
        return cov;
    }

    // Set coverages of point j
    void setCov(vector_fp cov, int j)
    {
        for(int i = 0; i<m_nspSurf; i++)
        {
            m_cov(i,j) = cov[i];
        }
    }

    /*//! Index of the species on the left boundary with the largest mass fraction
    size_t leftExcessSpecies() const {
        return m_kExcessLeft;
    }

    //! Index of the species on the right boundary with the largest mass fraction
    size_t rightExcessSpecies() const {
        return m_kExcessRight;
    }*/
    
    /*!
     *  Evaluate the residual function for axisymmetric stagnation flow. If
     *  j == npos, the residual function is evaluated at all grid points.
     *  Otherwise, the residual function is only evaluated at grid points
     *  j-1, j, and j+1. This option is used to efficiently evaluate the
     *  Jacobian numerically.
     */
    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt);

    //! Evaluate the residual function. This function is called in eval
    //! after updateProperties is called.
    virtual void evalResidual(double* x, double* rsd, int* diag, double rdt, size_t jmin, size_t jmax, size_t i);
    virtual void evalInletBoundary(double* x, double* rsd, int* diag, double rdt, size_t i);
    void evalOuterBoundary(double* x, double* rsd, int* diag, double rdt, size_t i);
    virtual void solveEnergyatInlet(double* x, double* rsd, int* diag, double rdt, double memflux, size_t i);
    virtual void solveEnergyatOutlet(double* x, double* rsd, int* diag, double rdt, double memflux, size_t i);

    virtual void showSolution(const doublereal* x);
    virtual void resetBadValues(double* xg);
    //virtual void updateCoverages(double* xg);

    /**
     * Solve the boundary value problem.
     * @param loglevel controls amount of diagnostic output.
     */
    void solve(int loglevel=0, Func1* f = nullptr) {
        if (m_sim == 0) {
            start();
        }
        bool refine = m_refineGrid!= 0;
        m_sim->setMaxGridPoints(-1, m_maxGridPoints);
        m_sim->setMaxTimeStepCount(m_maxTimeSteps);
        const int ts[1] = {m_minTimeSteps};
        m_sim->setTimeStep(m_minTimeStepSize,1,ts);
        m_sim->setMaxTimeStep(m_maxTimeStepSize);
        m_sim->setMinTimeStep(m_minTimeStepSize);
        m_sim->setRefineCriteria(1, m_max_grid_ratio, m_max_delta, m_max_delta_slope, m_prune);
        if(m_restart)
        {
            m_sim->restore("savedSoln.xml", {}, loglevel);
        }
        m_sim->setSteadyCallback(f);
        m_sim->solve(loglevel, refine);
        m_sim->save("savedSoln.xml", {}, {}, loglevel);
    }

    void restoreTimeSteppingSolution()
    {
        m_sim->restoreTimeSteppingSolution();
    }

    void restoreSteadySolution()
    {
        m_sim->restoreSteadySolution();
    }

    /*void setSteadyCallback(user_func* f)
    {
        f->eval(0);
    }*/

    void writeCSV(int err, bool dotitles = true, std::string ztitle = "z") const {
        std::string filename;
        if (err == (-1))
        {
            if(is_sccm)
            {
                filename = "Err_solution_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(int(m_SCCM)) + "sccm.csv";
            }
            else {
                filename = "Err_solution_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(double(m_inletVel)) + "vel.csv";
            }
        } else if (err == 0) {
            if(is_sccm)
            {
                filename = "solution_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(int(m_SCCM)) + "sccm.csv";
            }
            else {
                filename = "solution_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(double(m_inletVel)) + "vel.csv";
            }
        } else {
            if(is_sccm)
            {
                filename = "steadysolve_sol_" + to_string(int(err)) + "_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(int(m_SCCM)) + "sccm.csv";
            }
            else {
                filename = "steadysolve_sol_" + to_string(int(err)) + "_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(double(m_inletVel)) + "vel.csv";
            }
        }      
        std::ofstream f(filename);
        int np = nPoints();
        int nc = nComponents();
        int n, m, ind, neq;
        if (dotitles) {
            f << ztitle << ", ";
            for(int i = 0; i < m_radialPoints; i++)
            {
                f << "Pressure_"<<i << ", ";
            }
            for(int i = 0; i < m_radialPoints; i++)
            {
                if(i < m_lumenPoints)
                {
                    ind = i* nEq_lumen;
                    neq = nEq_lumen;
                } else {
                    neq = nEq_support;
                    ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
                }
                for (m = 0; m < neq; m++)
                {
                    f << componentName(ind+m) << "_"<<i;
                    f << ", ";
                }
            }
            if(m_do_channel)
            {
                ind = m_NeqTotal - nEq_channel;
                for (m = 0; m < nEq_channel; m++) 
                {
                    f << componentName(ind + m) << "_ch";
                    f << ", ";
                }
            }
            f << std::endl;
        }
        for (n = 0; n < np; n++) {
            f << z(n) << ", " ;
            for(int i = 0; i < m_radialPoints; i++)
            {
                f << m_press(n,i) << ", ";
            }
            for (m = 0; m < nc; m++) {
                f << m_sim->value(1, m, n);
                if (m != nc - 1) {
                    f << ", ";
                }
            }
            f << std::endl;
        }
    }

    void writeCSVMole(int err, bool dotitles = true, std::string ztitle = "z") const {
        std::string filename;
        if (err == (-1))
        {
            if(is_sccm)
            {
                filename = "Err_moleFractions_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(int(m_SCCM)) + "sccm.csv";
            }
            else {
                filename = "Err_moleFractions_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(double(m_inletVel)) + "vel.csv";
            }
        } else if (err == 0) {
            if(is_sccm)
            {
                filename = "moleFractions_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(int(m_SCCM)) + "sccm.csv";
            }
            else {
                filename = "moleFractions_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(double(m_inletVel)) + "vel.csv";
            }
        } else {
            if(is_sccm)
            {
                filename = "steadysolve_mole_" + to_string(int(err)) + "_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(int(m_SCCM)) + "sccm.csv";
            }
            else {
                filename = "steadysolve_mole_" + to_string(int(err)) + "_" + to_string(int(T_in)) + "K_" + to_string(int(p_out/1e5)) + "Bar_" + to_string(double(m_inletVel)) + "vel.csv";
            }
        }
        std::ofstream f(filename);
        int np = nPoints();
        int nc = nComponents();
        int n, m, ind, neq;
        if (dotitles) {
            f << ztitle << ", ";
            for(int i = 0; i < m_radialPoints; i++)
            {
                f << "Pressure_"<<i << ", ";
            }
            for(int i = 0; i < m_radialPoints; i++)
            {
                if(i < m_lumenPoints)
                {
                    ind = i* nEq_lumen;
                    neq = nEq_lumen;
                } else {
                    neq = nEq_support;
                    ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
                }
                for (m = 0; m < neq; m++)
                {
                    f << componentName(ind+m) << "_"<<i;
                    f << ", ";
                }
            }
            if(m_do_channel)
            {
                ind = m_NeqTotal - nEq_channel;
                for (m = 0; m < nEq_channel; m++) 
                {
                    f << componentName(ind + m) << "_ch";
                    f << ", ";
                }
            }
            f << std::endl;
        }
        for (n = 0; n < np; n++) {
            f << z(n)<< ", ";
            for(int i = 0; i < m_radialPoints; i++)
            {
                f << m_press(n,i) << ", ";
            }
            for(int i = 0; i < m_radialPoints; i++)
            {
                int ind;
                if(i < m_lumenPoints)
                {
                    ind = i* nEq_lumen;
                } else {
                    ind = nv_lumen + (i - m_lumenPoints) * nEq_support;
                }
                
                for (m = 0; m < c_offset_Y[0]; m++)
                {
                    f << m_sim->value(1, m + ind, n) << ", ";
                }
                for (m = 0; m < m_nsp; m++) {
                    f << m_xMole(m,cellIndex(n,i));
                    f << ", ";
                    //f << std::endl;
                }    
                for (m = c_offset_theta[0]; m < nEq_lumen; m++) {
                    f << m_sim->value(1, m + ind, n);
                    f << ", ";
                    //f << std::endl;
                }
            }
            if(m_do_channel)
            {
                int ind = m_NeqTotal - nEq_channel;
                for (m = 0; m < nEq_channel; m++) 
                {
                    f << m_sim->value(1, m + ind, n);
                    f << ", ";
                }
            } 
            f << std::endl;
        }
    }

    void readCSV(std::string filename){
        // The data is stored in m_initData array
        std::ifstream f(filename);
        // Check if file is open
        if(!f.is_open()) throw std::runtime_error("Could not open file");

        std::string line, col;
        double val;

        // Read the first line with column names
        if(f.good())
        {
            // Read the first line in the file
            std::getline(f, line);
            std::stringstream ss(line);

            // Extract column name
            while(std::getline(ss, col, ',')) {
            
                // Initialize and add <colname, double vector> pairs to result
                m_initData.push_back({col, std::vector<double> {}});
            }
        }

        // Read remaining lines
        while(std::getline(f, line)) {
        
            int ind = 0;
            std::stringstream ss(line);
            
            // Extract data
            while(ss >> val){
            
                // Add the current value to the 'ind' column's values vector
                m_initData.at(ind).second.push_back(val);
            
                // If the next token is a comma, ignore it and move on
                if(ss.peek() == ',') ss.ignore();
            
                 // Increment the column index
                ind++;
            }
        }

    // Close file
    f.close();
}

    //Component
    string componentName(size_t n) const;
    size_t componentIndex(const std::string& name) const;
    void getkExcess(const double* x);

protected:
    doublereal wdot(size_t k, size_t j, size_t i) const {
        return m_wdot(k,cellIndex(j,i));
    }

    //! Write the net production rates at point `j` into array `m_wdot`
    void getWdot(doublereal* x, size_t j, size_t i) {
        if(i< m_lumenPoints)
        {
            setGas(x,j,i);
            m_kin->getNetProductionRates(&m_wdot(0,cellIndex(j,i)));
        } else {
            setGasSupport(x,j,i);
            m_kin_2->getNetProductionRates(&m_wdot(0,cellIndex(j,i)));
        }
    }

    doublereal sdot(size_t k, size_t j, size_t i) const {
        return m_sdot(k,cellIndex(j,i));
    }

    //! Write the net surface production rates at point `j` into array `m_sdot`
    void getSdot(doublereal* x, size_t j, size_t i) {
        if(i< m_lumenPoints)
        {
            m_sdot.resize(m_nsp + m_nspSurf,m_points, 0.0);
            setGas(x,j,i);
            setSurf(x,j,i);
            m_surfkin->getNetProductionRates(&m_sdot(0,cellIndex(j,i)));
        } else{
            m_sdot.resize(m_nsp_2 + m_nspSurf_2,m_points, 0.0);
            setGasSupport(x,j,i);
            setSurfSupport(x,j,i);
            m_surfkin_2->getNetProductionRates(&m_sdot(0,cellIndex(j,i)));
        }
        
    }

    size_t mindex(size_t k, size_t j, size_t m) {
        return m*m_nsp*m_nsp + m_nsp*j + k;
    }

    //---------------------------------------------------------
    //             member data
    //---------------------------------------------------------

    // Offsets of solution components in the solution array.
    std::vector<int> c_offset_U;        // Axial velocity
    std::vector<int> c_offset_V;        // Radial velocity
    std::vector<int> c_offset_T;        // Temperature
    std::vector<int> c_offset_Ts;       // Solid phase temperature
    std::vector<int> c_offset_Y;        // mass fractions of gas-phase species
    std::vector<int> c_offset_theta;    // surface coverages
    size_t c_offset_ch_U;               // Velocity inside annular channel
    size_t c_offset_ch_rho;             // Density inside annular channel
    size_t c_offset_ch_T;               // Temperature inside annular channel

    size_t m_axialPoints = 0;
    size_t m_radialPoints = 1;
    size_t m_totRadialPoints = 1;
    size_t m_lumenPoints = 1;
    size_t m_supportPoints = 1;
    size_t m_points = 0;
    int m_nParts = 1;
    int nv_lumen = 0;

    // grid parameters
    vector_fp m_dz, m_zWall;
    vector_fp m_r, m_rWall;

    size_t m_nsp = 0;
    size_t m_nspSurf = 0;
    size_t m_nv = 0;
    int m_do_energy = 0;

    IdealGasPhase* m_thermo;
    GasKinetics* m_kin;
    Transport* m_trans;
    InterfaceKinetics* m_surfkin;
	SurfPhase* m_surf;
    DustyGasTransport* m_tranDGM;

    // Support
    IdealGasPhase* m_thermo_2;
    GasKinetics* m_kin_2;
    Transport* m_trans_2;
    InterfaceKinetics* m_surfkin_2;
	SurfPhase* m_surf_2;

    size_t m_nsp_2 = 0;
    size_t m_nspSurf_2 = 0;
    size_t nEq_lumen = 0;
    size_t nEq_support = 0;
    size_t m_NeqTotal = 0;
    size_t nEq_channel = 0;

    // mixture thermo properties
    Array2D m_rho;
    Array2D m_wtm;
    Array2D m_press;
    Array2D p_prev;
    Array2D rho_prev;
    vector_fp rho_prev_2;

    // species thermo properties
    vector_fp m_wt;
    vector_fp m_wt_2;
    Array2D m_cp;

    //Surface site-density
    vector_double m_gamma;

    // transport properties
    vector_fp m_diff;
    vector_fp m_diff_radial;
    Array2D m_visc;
    Array2D m_tcon;
    vector_fp m_multidiff;
    vector_fp m_multidiff_radial;
    Array2D m_dthermal;
    Array2D m_flux;
    Array2D m_interfaceFlux;
    Array2D m_xMole;
    vector_double radialFluxTerm;

    // production rates
    Array2D m_wdot, m_sdot;

    // Membrane parameters
    int m_memSpeciesIndex = -1, m_solveMem;
    double m_perm, m_memThickness, m_presMem, m_perm0, m_Ea_R;
    double m_Rin, m_length, m_memArea2Vol;

    // Restart
    int m_restart;

    // Initial guess from solution of diffusion-free code
    int m_initGuess = 0, nSteady = 0;
    std::string m_initGuessFile;
    std::vector<std::pair<std::string, std::vector<double>>> m_initData;

    // Surface solve
    Array2D m_cov;

    //! Index of gas and surface phase species with a large mass fraction at each boundary, for which
    //! the mass fraction may be calculated as 1 minus the sum of the other mass
    //! fractions
    vector_fp m_kExcess;
    size_t m_kExcessLeft;
    size_t m_kExcessRight;
    //size_t m_surfExcessLeft;
    //size_t m_surfExcessRight;

    vector_double m_porosity, m_area2Vol, m_tortuosity, m_Dp, m_PoreRad, m_Bg, m_microPorosity, m_microPoreRad;
    double m_Rout, m_thickness_2;
    int m_solveSupport = 0;
    int m_solveRadialVel = 0;

    int m_refineGrid = 0;
    int m_maxGridPoints = 0;
    int m_maxTimeSteps = 0;
    int m_minTimeSteps = 0;
    double m_minTimeStepSize = 0;
    int m_maxTimeStepSize = 0;

    //Inlet values
    double rho_fixed, T_in, T_wall, p_out, m_inletVel, m_inletMassFlux;
    double m_hcoeff;
    vector_fp m_yInlet, m_xInlet;
    vector_fp m_thetaInlet;
    double m_conversion;

    //Multi/Mix transport
    int m_do_multicomponent = 0;

    //DGM/Fickian model
    int m_tranModel;

    // Random pore model
    int m_randPoreModel = 0;

    // Solve solid phase energy equation
    int m_solveSolid = 0;
    double m_condSolid, m_rhoSolid, m_CpSolid, m_emissivity, T_env, m_Aenv;

    double m_max_grid_ratio, m_prune, m_max_delta_slope, m_max_delta;
    Array2D hk_interface;

    // Annular channel flows
    int m_do_channel;
    double m_rChannel;
    double rho0_ch, u0_ch;

    // Dusty gas model
    vector_double state_prev = {};
	vector_double state_curr = {};

    private:
    vector_fp m_ybar;
    double const m_FICK = 0;
    double const m_DGM = 1;

    int m_setupFirst = 1;
    int is_sccm = 0;
    double m_SCCM;

    public:
    int refine_count = 0;
};
#endif

class user_func : public Cantera::Func1
{
public:
    friend class bvpQ2D;
    user_func(bvpQ2D* ptr)
    {
        Func1();
        //cout<<"\n refinegrid == "<< ptr->m_refineGrid;
        main_ptr = ptr;
        //m_func = &write_instance;
    }

    virtual doublereal eval(doublereal a) const
    {
        main_ptr->refine_count += 1;
        main_ptr->writeCSV(main_ptr->refine_count);
        main_ptr->writeCSVMole(main_ptr->refine_count);
        return 0;
    }

    //protected:
    bvpQ2D* main_ptr;  
};