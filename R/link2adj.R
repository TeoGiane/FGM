# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C)  Laura Codazzi, Alessandro Colombi, Matteo Gianella                           |
#                                                                                                 |
#     Part of this code is adapted from BDgraph package (C Reza Mohammadi)                        |
#                                                                                                 |
#     FGM is free software: you can redistribute it and/or modify it under                        |
#     the terms of the GNU General Public License as published by the Free                        |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
#                                                                                                 |
#     Maintainer: Laura Codazzi (laura.codazzi@tuhh.de),                                          |
#                 Alessandro Colombi (a.colombi10@campus.unimib.it),                              |
#                 Matteo Gianella (matteo.gianella@polimi.it)                                     |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Creating an adjacency matrix based on links                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

link2adj = function( link, p = NULL )
{
    if( !is.matrix( link ) & !is.data.frame( link ) ) stop( " Input 'link' must be a matrix or dataframe" )
    if( is.data.frame( link ) ) link <- data.matrix( link )
    
    if( ncol( link ) != 2 ) stop( " Input 'link' must have only 2 columns" )
    if( nrow( link ) < 1 ) stop( " Input 'link' must have at least one row" )
    
    if( !is.null( p ) ) 
        if( max( link ) > p ) stop( " Value of 'p' is not matched with input 'link'" )
    
    if( is.null( p ) ) p = max( link )
    
    adj = matrix( 0, p, p )
    
    for( i in 1:nrow( link ) )
        adj[ link[ i, 1 ], link[ i, 2 ] ] = 1 
    
    return( adj ) 
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
