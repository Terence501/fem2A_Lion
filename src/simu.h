#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh;
            
            mesh.load (mesh_filename);
            
            std::vector<double> F(mesh.nb_vertices(),0);
            std::cout << "F calculée" << std::endl;
            
            SparseMatrix K(mesh.nb_vertices());
            
            for (int triangle = 0; triangle < mesh.nb_triangles(); ++ triangle){
            	ElementMapping elt_mapping = ElementMapping(mesh,false,triangle);
            	ShapeFunctions sh_functions = ShapeFunctions(2,1);
            	Quadrature quad = Quadrature :: get_quadrature(2);
            	DenseMatrix Ke;
            	assemble_elementary_matrix(elt_mapping,sh_functions,quad,unit_fct,Ke);
            	local_to_global_matrix(mesh, triangle , Ke, K);
            }
            std::cout << "K calculée" << std::endl;
            
            std::vector <double> values(mesh.nb_vertices()) ;
            for (int i= 0; i < mesh.nb_vertices(); ++ i){
            	values[i]= xy_fct(mesh.get_vertex(i));
            }
            
            std::vector <bool> attribute_is_dirichlet(2,false) ;
            attribute_is_dirichlet[1] = true ;
            mesh.set_attribute(unit_fct,1,true);
            std::vector<double> u(mesh.nb_vertices());
            
            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet,values,K,F );
            solve(K,F,u);
            
            std::cout << "Système résolu" << std::endl;
            std :: string export_name = "pureDirichletsquare";
            mesh.save(export_name+=".mesh");
            save_solution(u,export_name+".bb");
        }

    }

}
