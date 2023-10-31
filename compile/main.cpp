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
    bvpQ2D BVP_instance(filename, ip_ptr, nv, points, L, energyFlag);
    static bvpQ2D *bvp_ptr = &BVP_instance;

    // Create user_func class
    user_func write_instance(bvp_ptr);
    static user_func *write_ptr = &write_instance;
    //write_ptr->eval(0);
    	
	// Solve function
    try {
        // Solve the equations, refining the grid as needed        
        cout<<"\n number of components = " << BVP_instance.nComponents();
        BVP_instance.setSteadyTolerances(RTOL_Steady, ATOL_Steady, -1);
        BVP_instance.setTransientTolerances(RTOL_Transient, ATOL_Transient, -1);
        BVP_instance.solve(loglevel, write_ptr);
        
        // Write the solution to a CSV file.
        BVP_instance.writeStats();
        BVP_instance.writeCSV(0);
        BVP_instance.writeCSVMole(0);
        //BVP_instance.writeConversion();
        cout<<"\n Steady state is reached. Program exits successfully."<<endl;
        return 0;
    } catch (Cantera::CanteraError& err) {
        std::cerr << err.what() << std::endl;
        // Restore the last successful solution
        BVP_instance.restoreSteadySolution();
        //BVP_instance.restoreTimeSteppingSolution();
        BVP_instance.writeStats();
        BVP_instance.writeCSV(-1);
        BVP_instance.writeCSVMole(-1);
        cout<<"\n Error"<<endl;
        return -1;
    }
}

void PrintInfo(void)
{
    cout<<"\nThis code simulates a steady, one dimensional boundary-value problem across a catalyst support \n \n";
};
