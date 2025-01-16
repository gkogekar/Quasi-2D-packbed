#include "main.h"
using namespace std;

int main()
{
	// This is the main subroutine which initilizes, solves and 
	// then post-processes the boundary value problem

	PrintInfo();

	// Get input file
	string filename = "input.dat";

    // Create instance for Input class
    inputClass input_instance;
    input_instance.readInputFile(filename);
    static inputClass *ip_ptr = &input_instance;

    // Required parameters
    int nv = (ip_ptr->m_NeqTotal);
    int points = ip_ptr->m_meshPoints;
    double L = ip_ptr->m_axialLength;
    int energyFlag = ip_ptr->m_do_energy;
    double ATOL_Steady = ip_ptr->m_Atol_ss;
    double RTOL_Steady = ip_ptr->m_Rtol_ss;
    double ATOL_Transient = ip_ptr->m_Atol_ts;
    double RTOL_Transient = ip_ptr->m_Rtol_ts;
    int loglevel = ip_ptr->m_loglevel;
    int refine = ip_ptr->m_refine;
    
    // Create BVP_instance instance
    auto bvp_ptr = std::make_shared<bvpQ2D>(filename, ip_ptr, nv, points, L, energyFlag);
    
    // Create user_func class
    user_func write_instance(bvp_ptr);
    static user_func *write_ptr = &write_instance;
    //write_ptr->write_eval(0);
    	
	// Solve function
    try {
        // Solve the equations, refining the grid as needed        
        cout<<"\n number of components = " << bvp_ptr->nComponents();
        bvp_ptr->setSteadyTolerances(RTOL_Steady, ATOL_Steady, -1);
        bvp_ptr->setTransientTolerances(RTOL_Transient, ATOL_Transient, -1);
        bvp_ptr->solve(loglevel, write_ptr);
        
        // Write the solution to a CSV file.
        bvp_ptr->writeCSV(0);
        bvp_ptr->writeCSVMole(0);
        bvp_ptr->writeConversion();
        cout<<"\n Steady state is reached. Program exits successfully."<<endl;
        return 0;
    } catch (Cantera::CanteraError& err) {
        std::cerr << err.what() << std::endl;
        // Restore the last successful solution
        bvp_ptr->restoreSteadySolution();
        //bvp_ptr->restoreTimeSteppingSolution();
        bvp_ptr->writeCSV(-1);
        bvp_ptr->writeCSVMole(-1);
        cout<<"\n Error"<<endl;
        return -1;
    }
}

void PrintInfo(void)
{
    cout<<"\nThis code simulates a steady, quasi-2D boundary-value problem across a catalytic membrane reactor \n \n";
};
