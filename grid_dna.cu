#include "dna_structures.cuh"
__global__ void zero1_interaction_grid(cufftReal *grid,int grid_size)
{
    int x=threadIdx.x+(blockDim.x*blockIdx.x);
    int y=threadIdx.y+(blockDim.y*blockIdx.y);
    int z=threadIdx.z+(blockDim.z*blockIdx.z);

    if(z<grid_size&&x<grid_size&&y<grid_size)grid[gaddress(x,y,z,grid_size)] = (cufftReal)0;
}
__global__ void interaction_grid(cufftReal *grid, Nucleic_Acid *Nucleotide,float grid_span , int grid_size ,int steps,int ydim)
{
  int nucleotide=threadIdx.y+(blockDim.y*blockIdx.y);
  int atom=threadIdx.x+(blockDim.x*blockIdx.x);
    int x_step , y_step , z_step ;

     float		x_centre , y_centre , z_centre ;

  /* Variables */

     float         distance,one_span;
     one_span = grid_span / (float)grid_size ;

     distance = 1.8 ;

if(nucleotide<ydim){

    if((nucleotide>0)&&(atom>0)&&(atom<=Nucleotide[nucleotide].size))
    {

        
        int x = gord(Nucleotide[nucleotide].atom[atom].coord[1] , grid_span , grid_size );
        int y = gord(Nucleotide[nucleotide].atom[atom].coord[2] , grid_span , grid_size );
        int z = gord(Nucleotide[nucleotide].atom[atom].coord[3] , grid_span , grid_size );

        for( x_step = max( ( x - steps ) , 0 ) ; x_step <= min( ( x + steps ) , ( grid_size - 1 ) ) ; x_step ++ ) {

            x_centre  = gcentre( x_step , grid_span , grid_size ) ;

        for( y_step = max( ( y - steps ) , 0 ) ; y_step <= min( ( y + steps ) , ( grid_size - 1 ) ) ; y_step ++ ) {

          y_centre  = gcentre( y_step , grid_span , grid_size ) ;

          for( z_step = max( ( z - steps ) , 0 ) ; z_step <= min( ( z + steps ) , ( grid_size - 1 ) ) ; z_step ++ ) {

            z_centre  = gcentre( z_step , grid_span , grid_size ) ;

            if( pythagoras(Nucleotide[nucleotide].atom[atom].coord[1] ,Nucleotide[nucleotide].atom[atom].coord[2] ,Nucleotide[nucleotide].atom[atom].coord[3] , x_centre , y_centre , z_centre ) < distance ) grid[gaddress(x_step,y_step,z_step,grid_size)] = (cufftReal)1 ;

          }
        }
     


    }

    }
  }
}

void discretise_dna(struct DNA_structure dna, float grid_span, int grid_size, cufftReal *grid){
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
struct Nucleic_Acid *nucleotide,*d_nucleotide;
nucleotide = (struct Nucleic_Acid*)malloc((dna.length+1)*sizeof(Nucleic_Acid));
int a=0;
for (int i = 1; i <=dna.length; i++)
{
  nucleotide[i]=dna.nucleotide[i];
  cudaMalloc(&nucleotide[i].atom,(dna.nucleotide[i].size+1)*sizeof(struct Atom));
  cudaMemcpy(nucleotide[i].atom,dna.nucleotide[i].atom,(dna.nucleotide[i].size+1)*sizeof(struct Atom),cudaMemcpyHostToDevice);
  a=max(a,dna.nucleotide[i].size);
  
}
cudaMalloc((void**)&d_nucleotide,(dna.length+1)*sizeof(struct Nucleic_Acid));
cudaMemcpy(d_nucleotide,nucleotide,(dna.length+1)*sizeof(struct Nucleic_Acid),cudaMemcpyHostToDevice);

  dim3 numblock1((a/threadperblock2D.x)+1,(dna.length/threadperblock2D.y)+1);
  steps = (int)( ( distance / one_span ) + 1.5 ) ;
  interaction_grid<<<numblock1,threadperblock2D>>>(grid, d_nucleotide, grid_span,grid_size,steps,dna.length+1);
  cudaDeviceSynchronize();
  cudaFree(d_nucleotide);
  
  free(nucleotide);
  /************/

  return ;

   
}

__global__ void surface_grid( float grid_span , int grid_size , cufftReal *grid , float surface , float internal_value ){
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