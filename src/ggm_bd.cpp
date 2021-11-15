// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C)  Laura Codazzi, Alessandro Colombi, Matteo Gianella                           |
//                                                                                                 |
//     Part of this code is adapted from BDgraph package (C Reza Mohammadi)                        |
//                                                                                                 |
//     FGM is free software: you can redistribute it and/or modify it under                        |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Laura Codazzi (laura.codazzi@tuhh.de),                                          |
//				   Alessandro Colombi (a.colombi10@campus.unimib.it),                              |
//				   Matteo Gianella (matteo.gianella@polimi.it)                                     |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

#include "matrix.h"
#include "rgwish.h"

extern "C" {
// birth-death MCMC for Gaussian Graphical models  
// for case D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_bdmcmc_map( int *iter, int *burnin, int G[], double g_prior[], int *FGM_ptr_gprior_length, double Ts[], double K[], 
                    int *p, double *threshold, int all_graphs[], double all_weights[], double K_hat[], 
                    char *sample_graphs[], double graph_weights[], int *size_sample_g,
                    int *b, int *b_star, double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, count_all_g = 0;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int i, j, ij, one = 1, dim = *p, pxp = dim * dim;
	double Dsij, sum_weights = 0.0, weight_C, sum_rates;
	bool this_one;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );

	// - - allocation for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// - - - - - - - - - - - - - - 

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	// FGM modification,
	int FGM_gprior_length = *FGM_ptr_gprior_length;
	int PriorType(-1); //integer that defines the type of prior. 1 is for BetaBernoulli, 2 for Bernoulli
	(FGM_gprior_length==2) ? PriorType=1 : PriorType=2; //1 is BetaBernoulli, 2 is Bernoulli

	// Counting size of nodes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}
	
	// For finding the index of rates 
	vector<int> index_row( qp );
	vector<int> index_col( qp );
	int counter = 0;
	vector<double> Dsijj( pxp ); 
	int FGM_E(0);
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
		    ij = j * dim + i;
		    
		    if(G[ij]==1){
		    	FGM_E++;
		    }
		    if(PriorType==2){

		    	if( ( g_prior[ ij ] != 0.0 ) or ( g_prior[ ij ] != 1.0 ) )
		    	{
		        	index_row[ counter ] = i;
		        	index_col[ counter ] = j;
		        	counter++;
		    	}	
		    }
		    else{
		    	index_row[ counter ] = i;
		        index_col[ counter ] = j;
		    	counter++;
		    }

		    // for calculating the birth/death rates
		    Dsij        = Ds[ ij ];
		    Dsijj[ ij ] = Dsij * Dsij / Ds[ j * dim + j ]; 
		}


	int sub_qp = counter;
	vector<double> rates( sub_qp );

	vector<double> log_ratio_g_prior;
	if(PriorType==2)
	{
		log_ratio_g_prior.resize(pxp);
			for( j = 1; j < dim; j++ )
			{
				for( i = 0; i < j; i++ )
				{
					ij = j * dim + i;
					log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
				}	
			}
		
	}
	else{
		log_ratio_g_prior.resize(2);
		log_ratio_g_prior[0] = g_prior[0];
		log_ratio_g_prior[1] = g_prior[1];
	}
		
	

// - - Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		//if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
// - - - STEP 1: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		
		if(PriorType==2){

			rates_bdmcmc_parallel( &rates[0], &log_ratio_g_prior[0], &PriorType, G, &FGM_E, &index_row[0], &index_col[0], &sub_qp, Ds, &Dsijj[0], &sigma[0], &K[0], b, &dim );
		}
		else{
			rates_bdmcmc_parallel( &rates[0], &log_ratio_g_prior[0], &PriorType, G, &FGM_E, &index_row[0], &index_col[0], &sub_qp, Ds, &Dsijj[0], &sigma[0], &K[0], b, &dim );
		}
		
		
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &sub_qp );
		selected_edge_i = index_row[ index_selected_edge ];
		selected_edge_j = index_col[ index_selected_edge ];

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			counter = 0;	
			for( j = 1; j < dim; j++ )
				for( i = 0; i < j; i++ )
					char_g[ counter++ ] = G[ j * dim + i ] + '0'; 

			weight_C = 1.0 / sum_rates;
			
			//for( i = 0; i < pxp; i++ ) K_hat[i] += K[i] * weight_C;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat[0], &one );			

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[ count_all_g ] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[ i ] == string_g )
				{
					graph_weights[ i ] += all_weights[ count_all_g ];
					all_graphs[ count_all_g ] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[ size_sample_graph ] = string_g;
				graph_weights[ size_sample_graph ]   = all_weights[ count_all_g ];
				all_graphs[ count_all_g ]            = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
    			
		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		(G[ selected_edge_ij ] == 1) ? FGM_E-- : FGM_E++; //update number of links currently in the graph
		G[ selected_edge_ij ] = 1 - G[ selected_edge_ij ];
		G[ selected_edge_i * dim + selected_edge_j ] = G[ selected_edge_ij ];
		if( G[ selected_edge_ij ] )
		{ 
			++size_node[ selected_edge_i ]; 
			++size_node[ selected_edge_j ]; 
		}else{ 
			--size_node[ selected_edge_i ]; 
			--size_node[ selected_edge_j ]; 
		}

// - - - STEP 2: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
// - - End of main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - | 

	#pragma omp parallel for
	for( i = 0; i < size_sample_graph; i++ ) 
	{
		sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	
	#pragma omp parallel for
	for( i = 0; i < pxp; i++ )
		K_hat[ i ] /= sum_weights;
}
        

              
} // End of exturn "C"
