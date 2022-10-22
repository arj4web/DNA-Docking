/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.cuh"

__global__ void zero1_interaction_grid(cufftReal *grid,int grid_size)
{
    int x=threadIdx.x+(blockDim.x*blockIdx.x);
    int y=threadIdx.y+(blockDim.y*blockIdx.y);
    int z=threadIdx.z+(blockDim.z*blockIdx.z);

    if(z<grid_size&&x<grid_size&&y<grid_size)grid[gaddress(x,y,z,grid_size)] = (cufftReal)0;
}

__global__ void interaction_grid(cufftReal *grid, Amino_Acid *Residue,float grid_span , int grid_size ,int steps,int ydim)
{
  int residue=threadIdx.y+(blockDim.y*blockIdx.y);
  int atom=threadIdx.x+(blockDim.x*blockIdx.x);
    int x_step , y_step , z_step ;

     float		x_centre , y_centre , z_centre ;

  /* Variables */

     float         distance,one_span;
     one_span = grid_span / (float)grid_size ;

     distance = 1.8 ;

if(residue<ydim){

    if((residue>0)&&(atom>0)&&(atom<=Residue[residue].size))
    {

        
        int x = gord(Residue[residue].Atom[atom].coord[1] , grid_span , grid_size );
        int y = gord(Residue[residue].Atom[atom].coord[2] , grid_span , grid_size );
        int z = gord(Residue[residue].Atom[atom].coord[3] , grid_span , grid_size );

        for( x_step = max( ( x - steps ) , 0 ) ; x_step <= min( ( x + steps ) , ( grid_size - 1 ) ) ; x_step ++ ) {

            x_centre  = gcentre( x_step , grid_span , grid_size ) ;

        for( y_step = max( ( y - steps ) , 0 ) ; y_step <= min( ( y + steps ) , ( grid_size - 1 ) ) ; y_step ++ ) {

          y_centre  = gcentre( y_step , grid_span , grid_size ) ;

          for( z_step = max( ( z - steps ) , 0 ) ; z_step <= min( ( z + steps ) , ( grid_size - 1 ) ) ; z_step ++ ) {

            z_centre  = gcentre( z_step , grid_span , grid_size ) ;

            if( pythagoras(Residue[residue].Atom[atom].coord[1] ,Residue[residue].Atom[atom].coord[2] ,Residue[residue].Atom[atom].coord[3] , x_centre , y_centre , z_centre ) < distance ) grid[gaddress(x_step,y_step,z_step,grid_size)] = (cufftReal)1 ;

          }
        }
     


    }

    }
  }
}

void discretise_structure( struct Structure This_Structure , float grid_span , int grid_size , cufftReal *grid, int size1 ) {

/************/

  /* Counters */

 

  /* Co-ordinates */


  int	steps;


  /* Variables */

  float         distance , one_span ;

/************/

  one_span = grid_span / (float)grid_size ;

  distance = 1.8 ;

/************/
dim3 numblocks(((grid_size-1)/threadperblock3D.x)+1,((grid_size-1)/threadperblock3D.y)+1,((grid_size-1)/threadperblock3D.z)+1);



zero1_interaction_grid<<<numblocks,threadperblock3D>>>(grid,grid_size);
cudaDeviceSynchronize();


/************/
struct Amino_Acid *Residue,*d_Residue;
Residue = (struct Amino_Acid*)malloc((This_Structure.length+1)*sizeof(Amino_Acid));
int a=0;
for (int i = 1; i <=This_Structure.length; i++)
{
  Residue[i]=This_Structure.Residue[i];
  cudaMalloc(&Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom));
  cudaMemcpy(Residue[i].Atom,This_Structure.Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom),cudaMemcpyHostToDevice);
  a=max(a,This_Structure.Residue[i].size);
  
}
cudaMalloc((void**)&d_Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid));
cudaMemcpy(d_Residue,Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid),cudaMemcpyHostToDevice);

  dim3 numblock1((a/threadperblock2D.x)+1,(This_Structure.length/threadperblock2D.y)+1);
  steps = (int)( ( distance / one_span ) + 1.5 ) ;
  interaction_grid<<<numblock1,threadperblock2D>>>(grid, d_Residue, grid_span,grid_size,steps,This_Structure.length+1);
  cudaDeviceSynchronize();
  cudaFree(d_Residue);
  
  free(Residue);
  /************/

  return ;

}



/************************/




__global__ void surface_grid( float grid_span , int grid_size , cufftReal *grid , float surface , float internal_value ) {


/************/

  /* Counters */

  int	x=threadIdx.x+(blockIdx.x*blockDim.x) , y=threadIdx.y+(blockIdx.y*blockDim.y), z=threadIdx.z +(blockIdx.z*blockDim.z);
  int	steps , x_step , y_step , z_step ;

  /* Variables */

  float		one_span ;

  int	at_surface ;

/************/
if(z<grid_size&&x<grid_size&&y<grid_size){

  one_span = grid_span / (float)grid_size ;

/************/

  /* Surface grid atoms */

  steps = (int)( ( surface / one_span ) + 1.5 ) ;

  
        if( (int)grid[gaddress(x,y,z,grid_size)] == 1 ) {

          at_surface = 0 ;

          for( x_step = max( x - steps , 0 ) ; x_step <= min( x + steps , grid_size - 1 ) ; x_step ++ ) {
            for( y_step = max( y - steps , 0 ) ; y_step <= min( y + steps , grid_size - 1 ) ; y_step ++ ) {
              for( z_step = max( z - steps , 0 ) ; z_step <= min( z + steps , grid_size - 1 ) ; z_step ++ ) {

                if( (int)grid[gaddress(x_step,y_step,z_step,grid_size)] == 0 ) {

                  if( ( (float)( ( ( x_step - x ) * ( x_step - x ) ) + ( ( y_step - y ) * ( y_step - y ) ) + ( ( z_step - z ) * ( z_step - z ) ) ) * one_span * one_span ) < ( surface * surface ) ) at_surface = 1 ;

                }

              }
            }
          }

          if( at_surface == 0 ) grid[gaddress(x,y,z,grid_size)] = (cufftReal)internal_value ;

        }
}
/************/

  return ;

}
