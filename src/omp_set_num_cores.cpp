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

#include "util.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

extern "C" {
	void omp_set_num_cores( int *cores ) 
	{
	    #ifdef _OPENMP
	        omp_set_num_threads( *cores );
	    #else
	        Rprintf( "  This OS does not support multi-threading for the BDgraph package  \n" ); 
	    #endif
	}
}
