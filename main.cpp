#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quad = false;
    const bool t_elemap = false;
    const bool t_shapefunction = false;
    const bool t_Ke = false;
    const bool t_F= false;

    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if (t_quad) Tests::test_quadrature(2);
    if (t_elemap) Tests :: test_elementmapping(false, 4);
    if (t_shapefunction) Tests :: test_ShapeFunctions( 1,  1);
    if (t_Ke) Tests :: test_Ke();
    if (t_F) Tests :: test_F();
}

void run_simu()
{

    const bool simu_pure_dirichlet = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) {
        //Simu::pure_dirichlet_pb("data/square.mesh", verbose);
        //Simu :: dirichlet_ts_pb("data/square.mesh",verbose);
        //Simu :: pb_sinus_bump("data/square.mesh",verbose);
        //Simu :: sol_exacte("data/square.mesh",verbose);
        Simu :: ecart_sinus_bump("data/square.mesh",verbose);
    }
}

int main( int argc, const char * argv[] )
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
