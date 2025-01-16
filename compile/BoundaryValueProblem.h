/* 
    This file is taken from Cantera sourse code
    (https://github.com/Cantera/cantera/tree/main/samples/cxx/bvp)
*/

/// @file BoundaryValueProblem.h
/// Simplified interface to the capabilities provided by Cantera to
/// solve boundary value problems.

#ifndef BVP_H
#define BVP_H

#include "cantera/onedim.h"
#include <cantera/numerics/Func1.h>
#include <fstream>

/// Namespace for the boundary value problem package.
namespace BVP
{

// default grid refinement parameters
const double max_grid_ratio = 5.0; ///< max ratio of neighboring grid intervals
const double max_delta = 0.1; ///< max difference in function values
const double max_delta_slope = 0.1; ///< max difference in slopes
const double prune = 0.01; ///< don't remove grid points

/**
 * Used to specify component-specific options for method
 * setComponent of method BoundaryValueProblem. An instance of
 * class Component should be created for each solution component,
 * and its values set appropriately.
 */
class Component
{
public:
    double lower; ///< lower bound
    double upper; ///< upper bound
    double rtol; ///< relative error tolerance
    double atol; ///< absolute error tolerance
    bool refine; ///< make this component active for grid refinement
    std::string name; ///< component name

    /**
     * Constructor. Sets default values.
     */
    Component() : lower(0.0), upper(1.0), rtol(1.0e-9), atol(1.0e-12),
        refine(true) {}
};

/**
 * Base class for boundary value problems. This class is designed
 * to provide a simplified interface to the capabilities Cantera
 * provides to solve boundary value problems. Classes for specific
 * boundary value problems should be derived from this one.
 *
 * Class BoundaryValueProblem derives from Cantera's Domain1D
 * class.
 */
class BoundaryValueProblem :
    public Cantera::Domain1D,
    public std::enable_shared_from_this<BoundaryValueProblem>
{

public:

    /**
     * Constructor. This constructor begins with a uniform grid of
     * np points starting at zmin, and ending at zmax.
     *
     * @param nv Number of solution components
     * @param np Number of grid points in initial grid
     * @param zmin Location of left-hand side of domain
     * @param zmax Location of right-hand side of domain
     */
    BoundaryValueProblem(int nv, int np, double zmin, double zmax) :
        m_left(0), m_right(0), m_sim(0)
    {
        // Create the initial uniform grid
        Cantera::vector<double> z(np);
        int iz;
        for (iz = 0; iz < np; iz++) {
            z[iz] = zmin + iz*(zmax - zmin)/(np-1);
        }
        setupGrid(np, z.data());
        resize(nv, np);

        // Add dummy terminator domains on either side of this one.
        m_left = std::make_shared<Cantera::Empty1D>();
        m_right = std::make_shared<Cantera::Empty1D>();
    }

    /**
     * Constructor. This alternate constructor starts with a
     * specified grid, unlike the above that uses a uniform grid
     * to start. The array z must contain the z coordinates of np
     * grid points.
     */
    BoundaryValueProblem(int nv, int np, double* z) :
        m_left(0), m_right(0), m_sim(0)
    {
        setupGrid(np, z);
        resize(nv, np);

        // Add dummy terminator domains on either side of this one.
        m_left = std::make_shared<Cantera::Empty1D>();
        m_right = std::make_shared<Cantera::Empty1D>();
    }

    /**
     *  Set parameters and options for solution component @e n.
     *  This method should be invoked for each solution component
     *  before calling 'solve'. The parameter values should first
     *  be set by creating an instance of class Component, and
     *  setting its member data appropriately.
     *
     *  @param n Component number.
     *  @param c Component parameter values
     */
    void setComponent(size_t n, Component& c) {
        if (n >= m_nv) {
            throw Cantera::CanteraError("BoundaryValueProblem::setComponent",
                                        "Illegal solution component number");
        }
        // set the upper and lower bounds for this component
        setBounds(n, c.lower, c.upper);
        // set the error tolerances
        setSteadyTolerances(c.rtol, c.atol, n);
        setTransientTolerances(c.rtol, c.atol, n);
        // specify whether this component should be considered in
        // refining the grid
        m_refiner->setActive(n, c.refine);
        // set a default name if one has not been entered
        if (c.name == "") {
            c.name = fmt::format("Component {}", n);
        }
        setComponentName(n, c.name);
    }

    /**
     * Solve the boundary value problem.
     * @param loglevel controls amount of diagnostic output.
     */
    void solve(int loglevel=0) {
        if (!m_sim) {
            start();
        }
        bool refine = true;
        m_sim->solve(loglevel, refine);
    }

    /**
     * Write the solution to a CSV file.
     * @param filename CSV file name.
     * @param ztitle Title for 'z' column.
     * @param dotitles If true, begin with a row of column titles.
     */
    void writeCSV(std::string filename = "output.csv",
                  bool dotitles = true, std::string ztitle = "z") const {
        std::ofstream f(filename);
        int np = nPoints();
        int nc = nComponents();
        int n, m;
        if (dotitles) {
            f << ztitle << ", ";
            for (m = 0; m < nc; m++) {
                f << componentName(m);
                if (m != nc - 1) {
                    f << ", ";
                }
            }
            f << std::endl;
        }
        for (n = 0; n < np; n++) {
            f << z(n) << ", ";
            for (m = 0; m < nc; m++) {
                f << m_sim->value(1, m, n);
                if (m != nc - 1) {
                    f << ", ";
                }
            }
            f << std::endl;
        }
    }

    /**
     * Initial value of solution component \a n at initial grid
     * point \a j. The default is zero for all components at all
     * grid points. Overload in derived classes to specify other
     * choices for initial values.
     */
    double initialValue(size_t n, size_t j) override {
        return 0.0;
    }

protected:
    std::shared_ptr<Cantera::Domain1D> m_left; ///< dummy terminator
    std::shared_ptr<Cantera::Domain1D> m_right; ///< dummy terminator
    std::shared_ptr<Cantera::Sim1D> m_sim; ///< controller for solution
    
    /**
     * Set up the problem. Creates the solver instance, and sets
     * default grid refinement parameters. This method is called
     * internally, and does not need to be invoked explicitly in
     * derived classes.
     */
    void start() {
        std::vector<std::shared_ptr<Cantera::Domain1D>> domains {
            m_left, shared_from_this(), m_right
        };

        // create the Sim1D instance that will control the
        // solution process
        m_sim = std::make_shared<Cantera::Sim1D>(domains);

        // set default grid refinement parameters
        m_sim->setRefineCriteria(1, max_grid_ratio, max_delta,
                                 max_delta_slope, prune);
    }

};
}
#endif
