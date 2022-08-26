#include "dna_structures.h"

void discretise_dna(struct DNA_structure dna, float grid_span, int grid_size, fftw_real *grid){
    int x, y, z;
    int steps, x_step, y_step, z_step;
    float x_centre, y_centre, z_centre;
    float distance = 1.8, one_span;
    one_span = grid_span/(float)grid_size;

    for(x=0; x<grid_size; x++){
        for(y=0; y<grid_size; y++){
            for(z=0; z<grid_size; z++){
                grid[gaddress(x, y, z, grid_size)] = (fftw_real)0;
            }
        }
    }    
    steps = (int)( ( distance / one_span ) + 1.5 );

    for(int nucleotide=1; nucleotide<=dna.length; nucleotide++){
        for(int atom=1; atom<=dna.nucleotide[nucleotide].size; atom++){

            x = gord(dna.nucleotide[nucleotide].atom[atom].coord[1], grid_span, grid_size);
            y = gord(dna.nucleotide[nucleotide].atom[atom].coord[2], grid_span, grid_size);
            z = gord(dna.nucleotide[nucleotide].atom[atom].coord[3], grid_span, grid_size);

            for( x_step = max( ( x - steps ) , 0 ) ; x_step <= min( ( x + steps ) , ( grid_size - 1 ) ) ; x_step ++ ){
                x_centre  = gcentre( x_step , grid_span , grid_size ) ;

                   for( y_step = max( ( y - steps ) , 0 ) ; y_step <= min( ( y + steps ) , ( grid_size - 1 ) ) ; y_step ++ ) {
                        y_centre  = gcentre( y_step , grid_span , grid_size ) ;

                        for( z_step = max( ( z - steps ) , 0 ) ; z_step <= min( ( z + steps ) , ( grid_size - 1 ) ) ; z_step ++ ) {
                            z_centre  = gcentre( z_step , grid_span , grid_size ) ;

                            if(pythagoras(dna.nucleotide[nucleotide].atom[atom].coord[1], dna.nucleotide[nucleotide].atom[atom].coord[2], dna.nucleotide[nucleotide].atom[atom].coord[3], x_centre, y_centre, z_centre)<distance)
                                grid[gaddress(x_step,y_step,z_step,grid_size)] = (fftw_real)1;
                        }
                   }
            }
        }
    }
    return;
}

void surface_grid( float grid_span , int grid_size , fftw_real *grid , float surface , float internal_value ){
  int	x, y, z;
  int	steps, x_step, y_step, z_step;
  float		one_span;
  int	at_surface;
  one_span = grid_span / (float)grid_size;
  steps = (int)( ( surface / one_span ) + 1.5 );

  for(x = 0; x < grid_size; x ++){
    for(y = 0; y < grid_size; y ++){
      for(z = 0; z < grid_size; z ++){

        if((int)grid[gaddress(x,y,z,grid_size)] == 1){

          at_surface = 0 ;

          for(x_step = max( x - steps, 0) ; x_step <= min(x + steps, grid_size - 1); x_step ++){
            for(y_step = max( y - steps, 0) ; y_step <= min(y + steps, grid_size - 1); y_step ++){
              for(z_step = max( z - steps, 0) ; z_step <= min(z + steps, grid_size - 1); z_step ++){

                if((int)grid[gaddress(x_step,y_step,z_step,grid_size)] == 0){

                  if(((float)( ( ( x_step - x ) * ( x_step - x ) ) + ( ( y_step - y ) * ( y_step - y ) ) + ( ( z_step - z ) * ( z_step - z ) ) ) * one_span * one_span ) < ( surface * surface ) ) at_surface = 1 ;

                }

              }
            }
          }

          if(at_surface == 0) grid[gaddress(x,y,z,grid_size)] = (fftw_real)internal_value ;

        }

      }
    }
  }
  return ;

}