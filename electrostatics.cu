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
__device__ int my_strcmp(const char *str_a, const char *str_b, unsigned len = 256){
  int match = 0;
  unsigned i = 0;
  unsigned done = 0;
  while ((i < len) && (match == 0) && !done){
    if ((str_a[i] == 0) || (str_b[i] == 0)) done = 1;
    else if (str_a[i] != str_b[i]){
      match = i+1;
      if ((int)str_a[i] - (int)str_b[i] < 0) match = 0 - (i + 1);}
    i++;}
  return match;
  }

  __device__ int my_strncmp(const char *s1, const char *s2, size_t n)
  {
      unsigned char c1 = '\0';
  unsigned char c2 = '\0';
  if (n >= 4)
    {
      size_t n4 = n >> 2;
      do
        {
          c1 = (unsigned char) *s1++;
          c2 = (unsigned char) *s2++;
          if (c1 == '\0' || c1 != c2)
            return c1 - c2;
          c1 = (unsigned char) *s1++;
          c2 = (unsigned char) *s2++;
          if (c1 == '\0' || c1 != c2)
            return c1 - c2;
          c1 = (unsigned char) *s1++;
          c2 = (unsigned char) *s2++;
          if (c1 == '\0' || c1 != c2)
            return c1 - c2;
          c1 = (unsigned char) *s1++;
          c2 = (unsigned char) *s2++;
          if (c1 == '\0' || c1 != c2)
            return c1 - c2;
        } while (--n4 > 0);
      n &= 3;
    }
  while (n > 0)
    {
      c1 = (unsigned char) *s1++;
      c2 = (unsigned char) *s2++;
      if (c1 == '\0' || c1 != c2)
        return c1 - c2;
      n--;
    }
  return c1 - c2;
  }
__global__ void assign_charges_on_GPU(struct Amino_Acid *Residue,int ydim)
{
  int residue=threadIdx.y+(blockDim.y*blockIdx.y);
  int atom=threadIdx.x+(blockDim.x*blockIdx.x);


if(residue<ydim){
  if((residue>0)&&(atom>0)&&(atom<=Residue[residue].size)){

    Residue[residue].Atom[atom].charge = 0.0;
    /* peptide backbone */

    if( my_strcmp(Residue[residue].Atom[atom].atom_name , " N  ",3 ) == 0 ) {
        if( my_strcmp(Residue[residue].res_name , "PRO",3 ) == 0 ) {
          Residue[residue].Atom[atom].charge = -0.10 ;
        } else {
          Residue[residue].Atom[atom].charge =  0.55 ;
          if( residue == 1 )Residue[residue].Atom[atom].charge = 1.00 ;
        }
      }


    if( my_strcmp( Residue[residue].Atom[atom].atom_name , " O  ",3 ) == 0 ) {
        Residue[residue].Atom[atom].charge = -0.55 ;
        if( residue == ydim-1)Residue[residue].Atom[atom].charge = -1.00 ;
      }
     /* charged residues */

      if( ( my_strcmp( Residue[residue].res_name , "ARG",3 ) == 0 ) && ( my_strncmp(Residue[residue].Atom[atom].atom_name , " NH" , 3 ) == 0 ) ) Residue[residue].Atom[atom].charge =  0.50 ;
      if( ( my_strcmp( Residue[residue].res_name , "ASP",3 ) == 0 ) && ( my_strncmp(Residue[residue].Atom[atom].atom_name , " OD" , 3 ) == 0 ) ) Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( my_strcmp( Residue[residue].res_name , "GLU",3 ) == 0 ) && ( my_strncmp(Residue[residue].Atom[atom].atom_name , " OE" , 3 ) == 0 ) ) Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( my_strcmp( Residue[residue].res_name , "LYS",3 ) == 0 ) && ( my_strcmp( Residue[residue].Atom[atom].atom_name , " NZ ",3 ) == 0 ) )Residue[residue].Atom[atom].charge =  1.00 ;

  }
  }
}
void assign_charges( struct Structure This_Structure ) {

/************/

  /* Counters */

  int a=0 ;

/************/
struct Amino_Acid *Residue,*d_Residue;
Residue = (struct Amino_Acid*)malloc((This_Structure.length+1)*sizeof(Amino_Acid));
for (int i = 1; i <= This_Structure.length; i++)
{
  Residue[i]=This_Structure.Residue[i];
  cudaMalloc((void**)&Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom));
  cudaMemcpy(Residue[i].Atom,This_Structure.Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom),cudaMemcpyHostToDevice);
  a=max(a,This_Structure.Residue[i].size);
  
}
cudaMalloc((void**)&d_Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid));
cudaMemcpy(d_Residue,Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid),cudaMemcpyHostToDevice);

dim3 numblocks((a/threadperblock2D.x)+1,(This_Structure.length/threadperblock2D.y)+1);
assign_charges_on_GPU<<<numblocks,threadperblock2D>>>(d_Residue,This_Structure.length+1);
cudaDeviceSynchronize();
cudaMemcpy(Residue,d_Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid),cudaMemcpyDeviceToHost);
for (int i = 1; i <= This_Structure.length; i++)
{
 
  cudaMemcpy(This_Structure.Residue[i].Atom,Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom),cudaMemcpyDeviceToHost);
  cudaFree(Residue[i].Atom);
  
}
cudaFree(d_Residue);
free(Residue);
/************/

}



/************************/
__global__ void zero_interaction_grid(cufftReal *grid,int grid_size)
{
    int x=threadIdx.x +(blockDim.x*blockIdx.x);
    int y=threadIdx.y+ (blockDim.y*blockIdx.y);
    int z=threadIdx.z + (blockDim.z*blockIdx.z);


    if(z<grid_size&&x<grid_size&&y<grid_size)
    {
      grid[gaddress(x,y,z,grid_size)] = (cufftReal)0;
    }
}

__global__ void electric_fieldonGPU(cufftReal *grid,int grid_size,float grid_span,int size1, Amino_Acid *Residue)
{
    int x=threadIdx.x +(blockDim.x*blockIdx.x);
    int y=threadIdx.y+ (blockDim.y*blockIdx.y);
    int z=threadIdx.z+(blockDim.z*blockIdx.z);
    float		distance ;
   float epsilon ;
  
    
    if(z<grid_size&&x<grid_size&&y<grid_size){
    if (y==0&&z==0)
    {
      printf( "." );
    }
    float x_centre  = gcentre( x , grid_span , grid_size ) ;
    float y_centre  = gcentre( y , grid_span , grid_size ) ;
    float z_centre  = gcentre( z , grid_span , grid_size ) ;
    float phi=0;
        for( int residue = 1 ; residue <= size1 ; residue ++ ) {
          for( int atom = 1 ; atom <=Residue[residue].size ; atom ++ ) {

            if(Residue[residue].Atom[atom].charge != 0 ) {

              distance = pythagoras( Residue[residue].Atom[atom].coord[1] , Residue[residue].Atom[atom].coord[2] , Residue[residue].Atom[atom].coord[3] , x_centre , y_centre , z_centre ) ;
         
              if( distance < 2.0 ) distance = 2.0 ;

              if( distance >= 2.0 ) {

                if( distance >= 8.0 ) {

                  epsilon = 80 ;

                } else { 

                  if( distance <= 6.0 ) { 

                    epsilon = 4 ;
             
                  } else {

                    epsilon = ( 38 * distance ) - 224 ;

                  }

                }
  
                phi += (Residue[residue].Atom[atom].charge / ( epsilon * distance ) ) ;

              }

            }

          }
        }

    grid[gaddress(x,y,z,grid_size)] = (cufftReal)phi ;
  }
}



void electric_field( struct Structure This_Structure , float grid_span , int grid_size , cufftReal *grid ) {

dim3 numblocks(((grid_size-1)/threadperblock3D.x)+1,((grid_size-1)/threadperblock3D.y)+1,((grid_size-1)/threadperblock3D.z)+1);

zero_interaction_grid<<<numblocks,threadperblock3D>>>(grid,grid_size);
cudaDeviceSynchronize();
struct Amino_Acid *Residue,*d_Residue;
Residue = (struct Amino_Acid*)malloc((This_Structure.length+1)*sizeof(Amino_Acid));
for (int i = 1; i <= This_Structure.length; i++)
{
  Residue[i]=This_Structure.Residue[i];
  cudaMalloc((void**)&Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom));
  cudaMemcpy(Residue[i].Atom,This_Structure.Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom),cudaMemcpyHostToDevice);

}
cudaMalloc((void**)&d_Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid));
cudaMemcpy(d_Residue,Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid),cudaMemcpyHostToDevice);
  



/************/

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  printf( "  electric field calculations ( one dot / grid sheet ) " ) ;

  electric_fieldonGPU<<<numblocks,threadperblock3D>>>(grid,grid_size,grid_span,This_Structure.length,d_Residue);
  cudaDeviceSynchronize();

  printf( "\n" ) ;
  cudaFree(d_Residue);
  free(Residue);
 


/************/

  return ;

}



/************************/
__global__ void helper_point_charge_GPU(Amino_Acid *Residue,float one_span,int x_low,int y_low,int z_low,float a,float b,float c,int grid_size,cufftReal *grid,int residue,int atom )
{
  float		x_corner , y_corner , z_corner ;
  int x=threadIdx.x + x_low;
  int y=threadIdx.y + y_low;
  int z=threadIdx.z + z_low;
  int x_high = x_low + 1 ;
  int y_high = y_low + 1 ;
  int z_high = z_low + 1 ;
  

  x_corner = one_span * ( (float)( x - x_high ) + .5 ) ;
  y_corner = one_span * ( (float)( y - y_high ) + .5 ) ;
  z_corner = one_span * ( (float)( z - z_high ) + .5 ) ;
  float w = ( ( x_corner + a ) * ( y_corner + b ) * ( z_corner + c ) ) / ( 8.0 * x_corner * y_corner * z_corner ) ;
  grid[gaddress(x,y,z,grid_size)] += (cufftReal)( w*Residue[residue].Atom[atom].charge ) ;

}

__global__ void point_charge_GPU(Amino_Acid *Residue,float one_span,float grid_span,int grid_size,cufftReal *grid,int ydim )
{
    int residue=threadIdx.y+(blockIdx.y*blockDim.y);
    int atom=threadIdx.x+(blockIdx.x*blockDim.x);
    int x_low, y_low, z_low;
    if(residue<ydim){
    if((residue>0)&&(atom>0)&&(atom<=Residue[residue].size))
    {
        if(Residue[residue].Atom[atom].charge != 0 ) {

        x_low = gord( Residue[residue].Atom[atom].coord[1] - ( one_span / 2 ) , grid_span , grid_size ) ;
        y_low = gord( Residue[residue].Atom[atom].coord[2] - ( one_span / 2 ) , grid_span , grid_size ) ;
        z_low = gord( Residue[residue].Atom[atom].coord[3] - ( one_span / 2 ) , grid_span , grid_size ) ;

        float a = Residue[residue].Atom[atom].coord[1] - gcentre( x_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        float b = Residue[residue].Atom[atom].coord[2] - gcentre( y_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        float c = Residue[residue].Atom[atom].coord[3] - gcentre( z_low , grid_span , grid_size ) - ( one_span / 2 ) ;

        dim3 threadPerblock(2,2,2);
        helper_point_charge_GPU<<<1,threadPerblock>>>(Residue,one_span,x_low,y_low,z_low,a,b,c,grid_size,grid,residue,atom);
        cudaDeviceSynchronize();
 
        }
    }
  }
}

void electric_point_charge( struct Structure This_Structure , float grid_span , int grid_size , cufftReal *grid ) {

/************/

  /* Variables */

  float		one_span ;

/************/
int a=0;
dim3 numblocks(((grid_size-1)/threadperblock3D.x)+1,((grid_size-1)/threadperblock3D.y)+1,((grid_size-1)/threadperblock3D.z)+1);





zero_interaction_grid<<<numblocks,threadperblock3D>>>(grid,grid_size);
cudaDeviceSynchronize();

/************/
struct Amino_Acid *Residue,*d_Residue;
Residue = (struct Amino_Acid*)malloc((This_Structure.length+1)*sizeof(Amino_Acid));
for (int i = 1; i <= This_Structure.length; i++)
{
  Residue[i]=This_Structure.Residue[i];
  cudaMalloc((void**)&Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom));
  cudaMemcpy(Residue[i].Atom,This_Structure.Residue[i].Atom,(This_Structure.Residue[i].size+1)*sizeof(struct Atom),cudaMemcpyHostToDevice);
  a=max(a,This_Structure.Residue[i].size);
  
}
cudaMalloc((void**)&d_Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid));
cudaMemcpy(d_Residue,Residue,(This_Structure.length+1)*sizeof(struct Amino_Acid),cudaMemcpyHostToDevice);
dim3 numblock1((a/threadperblock2D.x)+1,(This_Structure.length/threadperblock2D.y)+1);
  one_span = grid_span / (float)grid_size ;
  point_charge_GPU<<<numblock1,threadperblock2D>>>(d_Residue,one_span,grid_span,grid_size,grid,This_Structure.length+1);
  cudaDeviceSynchronize();
  free(Residue);
  cudaFree(d_Residue);

/************/

  return ;

}



/************************/



__global__ void electric_field_zero_core( int grid_size , cufftReal *elec_grid , cufftReal *surface_grid , float internal_value ) {

/************/

  /* Co-ordinates */

  int	x=threadIdx.x+(blockIdx.x*blockDim.x),y=threadIdx.y+(blockIdx.y*blockDim.y),z=threadIdx.z+(blockIdx.z*blockDim.z);

/************/

  if(z<grid_size)if( surface_grid[gaddress(x,y,z,grid_size)] == (cufftReal)internal_value ) elec_grid[gaddress(x,y,z,grid_size)] = (cufftReal)0 ;



/************/


}
