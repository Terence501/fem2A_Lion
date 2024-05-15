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
        
        double sinus_bump_fct(vertex v)
        {
        	const double pi = std::acos(-1);
        	return 2* (pi*pi) * std::sin(pi*v.x) * std::sin(pi*v.y);
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
            mesh.save(export_name + ".mesh");
            save_solution(u,export_name + ".bb");
        }
        
        
        
        
        
	void dirichlet_ts_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh;
            
            mesh.load (mesh_filename);
            
            std::vector<double> F(mesh.nb_vertices());
            SparseMatrix K(mesh.nb_vertices());
            
            for (int triangle = 0; triangle < mesh.nb_triangles(); ++ triangle){
            	ElementMapping elt_mapping = ElementMapping(mesh,false,triangle);
            	ShapeFunctions sh_functions = ShapeFunctions(2,1);
            	Quadrature quad = Quadrature :: get_quadrature(2);
            	DenseMatrix Ke;
            	std :: vector <double> Fe;
            	assemble_elementary_matrix(elt_mapping,sh_functions,quad,unit_fct,Ke);
            	assemble_elementary_vector(elt_mapping,sh_functions,quad,unit_fct,Fe);
            	local_to_global_matrix(mesh, triangle , Ke, K);
            	local_to_global_vector(mesh, false, triangle , Fe, F);
            }
            
            std::cout << "K calculée" << std::endl;
            
            std::vector <double> values(mesh.nb_vertices(),0) ;

            
            std::vector <bool> attribute_is_dirichlet(2,false) ;
            attribute_is_dirichlet[1] = true ;
            mesh.set_attribute(unit_fct,1,true);
            std::vector<double> u(mesh.nb_vertices());
            
            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet,values,K,F );
            solve(K,F,u);
            
            std::cout << "Système résolu" << std::endl;
            std :: string export_name = "dirichlet_ts_square";
            mesh.save(export_name + ".mesh");
            save_solution(u,export_name + ".bb");
        }
	
	
	
	void pb_sinus_bump( const std::string& mesh_filename, bool verbose )
	{
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh;
            
            mesh.load (mesh_filename);
            
            std::vector<double> F(mesh.nb_vertices());
            SparseMatrix K(mesh.nb_vertices());
            
            for (int triangle = 0; triangle < mesh.nb_triangles(); ++ triangle){
            	ElementMapping elt_mapping = ElementMapping(mesh,false,triangle);
            	ShapeFunctions sh_functions = ShapeFunctions(2,1);
            	Quadrature quad = Quadrature :: get_quadrature(2);
            	DenseMatrix Ke;
            	std :: vector <double> Fe;
            	assemble_elementary_matrix(elt_mapping,sh_functions,quad,unit_fct,Ke);
            	assemble_elementary_vector(elt_mapping,sh_functions,quad,sinus_bump_fct,Fe);
            	local_to_global_matrix(mesh, triangle , Ke, K);
            	local_to_global_vector(mesh, false, triangle , Fe, F);
            }
            
            std::cout << "K calculée" << std::endl;
            
            std::vector <double> values(mesh.nb_vertices(),0) ;

            
            std::vector <bool> attribute_is_dirichlet(2,false) ;
            attribute_is_dirichlet[1] = true ;
            mesh.set_attribute(unit_fct,1,true);
            std::vector<double> u(mesh.nb_vertices());
            
            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet,values,K,F );
            solve(K,F,u);
            
            std::cout << "Système résolu" << std::endl;
            std :: string export_name = "test_sinus_bump_square";
            mesh.save(export_name + ".mesh");
            save_solution(u,export_name + ".bb");
        }
        
                                              
	void sol_exacte(const std::string& mesh_filename,bool verbose){
   		Mesh mesh;
            	mesh.load(mesh_filename);
           
            	const double pi = std::acos(-1);
            	std::vector<double> u(mesh.nb_vertices(),0);
            	for (int i =0 ; i < mesh.nb_vertices() ; ++ i ){
            	u[i] = std :: sin(pi*mesh.get_vertex(i).x)*std::sin(pi*mesh.get_vertex(i).y);
            	}	
            	std :: string export_name = "sinus_exact_square";
           	mesh.save (export_name + ".mesh");
           	save_solution(u, export_name+".bb");
           	}


	void ecart_sinus_bump(const std::string& mesh_filename, bool verbose){
	
            std::cout << "Solving a pure Dirichlet problem" << std::endl;
            if ( verbose ) {
                std::cout << " with lots of printed details..." << std::endl;
            }
            
            Mesh mesh;
            
            mesh.load (mesh_filename);
            
            std::vector<double> F(mesh.nb_vertices());
            SparseMatrix K(mesh.nb_vertices());
            
            for (int triangle = 0; triangle < mesh.nb_triangles(); ++ triangle){
            	ElementMapping elt_mapping = ElementMapping(mesh,false,triangle);
            	ShapeFunctions sh_functions = ShapeFunctions(2,1);
            	Quadrature quad = Quadrature :: get_quadrature(2);
            	DenseMatrix Ke;
            	std :: vector <double> Fe;
            	assemble_elementary_matrix(elt_mapping,sh_functions,quad,unit_fct,Ke);
            	assemble_elementary_vector(elt_mapping,sh_functions,quad,sinus_bump_fct,Fe);
            	local_to_global_matrix(mesh, triangle , Ke, K);
            	local_to_global_vector(mesh, false, triangle , Fe, F);
            }
            
            std::cout << "K calculée" << std::endl;
            
            std::vector <double> values(mesh.nb_vertices(),0) ;

            
            std::vector <bool> attribute_is_dirichlet(2,false) ;
            attribute_is_dirichlet[1] = true ;
            mesh.set_attribute(unit_fct,1,true);
            std::vector<double> u(mesh.nb_vertices());
            
            apply_dirichlet_boundary_conditions(mesh,attribute_is_dirichlet,values,K,F );
            solve(K,F,u);
	    
	    std::vector< double > ecart(mesh.nb_vertices(),0);
	    const double pi = std::acos(-1);
	    double somme;
	    for (int i =0 ; i < mesh.nb_vertices() ; ++i ){
		ecart[i] = std::abs((std :: sin(pi*mesh.get_vertex(i).x)*std::sin(pi*mesh.get_vertex(i).y))-u[i]);
		somme+=ecart[i];
		   }
	   std :: string export_name = "ecart_sinus_square";
           mesh.save (export_name + ".mesh");
	   save_solution(ecart, export_name+".bb");
	   somme =somme/mesh.nb_vertices();
	   std::cout << "Erreur totale divisée par le nombre de triangles : " << somme << std::endl;

	}
}
}

