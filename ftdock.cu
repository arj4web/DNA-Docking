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
#include <cuda_runtime.h>




__global__ void convolution(cufftComplex *static_fsg,cufftComplex *multiple_fsg,cufftComplex *mobile_fsg,cufftComplex *static_elec_fsg,cufftComplex *mobile_elec_fsg,cufftComplex *multiple_elec_fsg,int electrostatics,int global_grid_size)
{
  int fx=threadIdx.x+(blockDim.x*blockIdx.x);
  int fy=threadIdx.y+(blockDim.y*blockIdx.y);
  int fz=threadIdx.z+(blockDim.z*blockIdx.z);
  if(fx<global_grid_size&&fy<global_grid_size&&fz<(((global_grid_size)/2)+1)){
  int fxyz = fz + ( global_grid_size/2 + 1 ) * ( fy + global_grid_size * fx ) ;

          multiple_fsg[fxyz].x =
           static_fsg[fxyz].x * mobile_fsg[fxyz].x + static_fsg[fxyz].y * mobile_fsg[fxyz].y;
          multiple_fsg[fxyz].y =
           static_fsg[fxyz].y * mobile_fsg[fxyz].x - static_fsg[fxyz].x * mobile_fsg[fxyz].y;
           
          if( electrostatics == 1 ) {
            multiple_elec_fsg[fxyz].x =
             static_elec_fsg[fxyz].x * mobile_elec_fsg[fxyz].x + static_elec_fsg[fxyz].y * mobile_elec_fsg[fxyz].y ;
            multiple_elec_fsg[fxyz].y =
             static_elec_fsg[fxyz].y * mobile_elec_fsg[fxyz].x - static_elec_fsg[fxyz].x * mobile_elec_fsg[fxyz].y ;
          } 
}
}
__global__ void init_score(Score *d_Scores,int dim_x)
{
      int i=threadIdx.x+(blockDim.x*blockIdx.x);
      if(i<dim_x)
      { 
        d_Scores[i].score = 0 ;
        d_Scores[i].rpscore = 0.0 ;
        d_Scores[i].coord[1] = 0 ;
        d_Scores[i].coord[2] = 0 ;
        d_Scores[i].coord[3] = 0 ;
      }

}
__global__ void get_score(Score *d_Scores,cufftReal *convoluted_grid,cufftReal *convoluted_elec_grid, int electrostatics,int keep_per_rotation,int global_grid_size)
{
  int x=threadIdx.x+(blockDim.x*blockIdx.x);
  int y=threadIdx.y+(blockDim.y*blockIdx.y);
  int z=threadIdx.z+(blockDim.z*blockIdx.z);
  int i=0;
if(x<global_grid_size&&y<global_grid_size&&z<global_grid_size){
    int fx = x ;
    if( fx > ( global_grid_size / 2 ) ) fx -= global_grid_size ;
    int fy = y ;
    if( fy > ( global_grid_size / 2 ) ) fy -= global_grid_size ;
    int fz = z ;
    if( fz > ( global_grid_size / 2 ) ) fz -= global_grid_size ;

    int xyz = z + ( 2 * ( global_grid_size / 2 + 1 ) ) * ( y + global_grid_size * x ) ;

    if( ( electrostatics == 0 ) || ( convoluted_elec_grid[xyz] < 0 ) ) {

         /* Scale factor from FFTs */
        if( (int)convoluted_grid[xyz] != 0 ) {
            convoluted_grid[xyz] /= ( global_grid_size * global_grid_size * global_grid_size ) ;
        }

        if( (int)convoluted_grid[xyz] > d_Scores[keep_per_rotation-1].score ) {

          i = keep_per_rotation - 2 ;

          while( ( (int)convoluted_grid[xyz] > d_Scores[i].score ) && ( i >= 0 ) ) {
                d_Scores[i+1].score    = d_Scores[i].score ;
                d_Scores[i+1].rpscore  = d_Scores[i].rpscore ;
                d_Scores[i+1].coord[1] = d_Scores[i].coord[1] ;
                d_Scores[i+1].coord[2] = d_Scores[i].coord[2] ;
                d_Scores[i+1].coord[3] = d_Scores[i].coord[3] ;
                i -- ;
              }

              d_Scores[i+1].score    = (int)convoluted_grid[xyz] ;
              if( ( electrostatics != 0 ) && ( convoluted_elec_grid[xyz] < 0.1 ) ) {
                d_Scores[i+1].rpscore  = (float)convoluted_elec_grid[xyz] ;
              } else {
                d_Scores[i+1].rpscore  = (float)0 ;
              }
              d_Scores[i+1].coord[1] = fx ;
              d_Scores[i+1].coord[2] = fy ;
              d_Scores[i+1].coord[3] = fz ;

            }

          }

  }
}
int main( int argc , char *argv[] ) {

  /* index counters */

  int	i ;

  /* Command line options */

  char		*output_file_name ;
  char		*static_file_name ;
  char		*mobile_file_name ;
  int		global_grid_size ;
  int		angle_step ;
  float		surface ;
  float		internal_value ;
  int		electrostatics ;
  int		keep_per_rotation ;
  int 		kept_scores ;
  int		rescue ;
  int		calculate ;
  float		reverse_calculated_one_span ;

  char		*default_global_grid_size ;
  char		*default_angle_step ;
  char		*default_surface ;
  char		*default_internal_value ;
  char		*default_electrostatics ;
  char		*default_keep_per_rotation ;

  /* File stuff */

  FILE		*ftdock_file ;
  char		line_buffer[100] ;
  int		id , id2 , SCscore ;
  float		RPscore ;
  int		x , y , z , z_twist , theta , phi ;

  /* Angles stuff */

  struct Angle	Angles ;
  int		first_rotation , rotation ;

  /* Structures */

  struct Structure	Static_Structure , Mobile_Structure ;
  struct Structure	Origin_Static_Structure , Origin_Mobile_Structure ;
  struct Structure	Rotated_at_Origin_Mobile_Structure ;



  /* Grid stuff */

  float		grid_span , one_span ;

  cufftReal	*static_grid;
  cufftReal	*mobile_grid;
  cufftReal	*convoluted_grid;

  cufftReal	*static_elec_grid;
  cufftReal	*mobile_elec_grid;
  cufftReal	*convoluted_elec_grid;

  /* FFTW stuff */

  cufftHandle	p , pinv ;

  cufftComplex  *static_fsg ;
  cufftComplex  *mobile_fsg ;
  cufftComplex  *multiple_fsg ;

  cufftComplex  *static_elec_fsg;
  cufftComplex  *mobile_elec_fsg;
  cufftComplex  *multiple_elec_fsg;
  cufftResult result;
  cudaError_t error;

  /* Scores */

  struct Score	*Scores, *d_Scores ;
  float		max_es_value ;

/************/

  /* Its nice to tell people what going on straight away */

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;


  printf( "\n          3D-Dock Suite (March 2001)\n" ) ;
  printf( "          Copyright (C) 1997-2000 Gidon Moont\n" ) ;
  printf( "          This program comes with ABSOLUTELY NO WARRANTY\n" ) ;
  printf( "          for details see license. This program is free software,\n"); 
  printf( "          and you may redistribute it under certain conditions.\n\n"); 

  printf( "          Biomolecular Modelling Laboratory\n" ) ;
  printf( "          Imperial Cancer Research Fund\n" ) ;
  printf( "          44 Lincoln's Inn Fields\n" ) ;
  printf( "          London WC2A 3PX\n" ) ;
  printf( "          +44 (0)20 7269 3348\n" ) ;
  printf( "          http://www.bmm.icnet.uk/\n\n" ) ;


  printf( "Starting FTDock (v2.0) global search program\n" ) ;


/************/

  /* Memory allocation */

  if( ( ( output_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( static_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( mobile_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ) {
    GENERAL_MEMORY_PROBLEM 
  }

/************/

  /* Command Line defaults */

  strcpy( output_file_name , "ftdock_global.dat" ) ;
  strcpy( static_file_name , " --static file name was not provided--" ) ;
  strcpy( mobile_file_name , " --mobile file name was not provided--" ) ;
  global_grid_size = 128 ;
  angle_step = 12 ;
  surface = 1.3 ;
  internal_value = -15 ;
  electrostatics = 1 ;
  keep_per_rotation = 3 ;
  rescue = 0 ;
  calculate = 1 ;
  reverse_calculated_one_span = 0.7 ;

  default_global_grid_size = "(default calculated)" ;
  default_angle_step = "(default)" ;
  default_surface = "(default)" ;
  default_internal_value = "(default)" ;
  default_electrostatics = "(default)" ;
  default_keep_per_rotation = "(default)" ;

  /* Command Line parse */

  for( i = 1 ; i < argc ; i ++ ) {

    if( strcmp( argv[i] , "-out" ) == 0 ) {
      i ++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
        printf( "Bad command line\n" ) ;
        exit( EXIT_FAILURE ) ;
      }
      strcpy( output_file_name , argv[i] ) ;
    } else {
      if( strcmp( argv[i] , "-static" ) == 0 ) {
        i ++ ;
        if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
          printf( "Bad command line\n" ) ;
          exit( EXIT_FAILURE ) ;
        }
        strcpy( static_file_name , argv[i] ) ;
      } else {
        if( strcmp( argv[i] , "-mobile" ) == 0 ) {
          i ++ ;
          if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
            printf( "Bad command line\n" ) ;
            exit( EXIT_FAILURE ) ;
          }
          strcpy( mobile_file_name , argv[i] ) ;
        } else {
          if( strcmp( argv[i] , "-grid" ) == 0 ) {
            i ++ ;
            if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
              printf( "Bad command line\n" ) ;
              exit( EXIT_FAILURE ) ;
            }
            sscanf( argv[i] , "%d" , &global_grid_size ) ;
            if( ( global_grid_size % 2 ) != 0 ) {
              printf( "Grid size must be even\n" ) ;
              exit( EXIT_FAILURE ) ;
            }
            default_global_grid_size = "(user defined)" ;
            calculate = 0 ;
          } else {
            if( strcmp( argv[i] , "-angle_step" ) == 0 ) {
              i ++ ;
              if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
                printf( "Bad command line\n" ) ;
                exit( EXIT_FAILURE ) ;
              }
              sscanf( argv[i] , "%d" , &angle_step ) ;
              default_angle_step = "(user defined)" ;
            } else {
              if( strcmp( argv[i] , "-surface" ) == 0 ) {
                i ++ ;
                if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
                  printf( "Bad command line\n" ) ;
                  exit( EXIT_FAILURE ) ;
                }
                sscanf( argv[i] , "%f" , &surface ) ;
                default_surface = "(user defined)" ;
              } else {
                if( strcmp( argv[i] , "-internal" ) == 0 ) {
                  i ++ ;
                  if( i == argc ) {
                    printf( "Bad command line\n" ) ;
                    exit( EXIT_FAILURE ) ;
                  }
                  sscanf( argv[i] , "%f" , &internal_value ) ;
                  default_internal_value = "(user defined)" ;
                } else {
                  if( strcmp( argv[i] , "-noelec" ) == 0 ) {
                    electrostatics = 0 ;
                    default_electrostatics = "(user defined)" ;
                  } else {
                    if( strcmp( argv[i] , "-keep" ) == 0 ) {
                      i ++ ;
                      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
                        printf( "Bad command line\n" ) ;
                        exit( EXIT_FAILURE ) ;
                      }
                      sscanf( argv[i] , "%d" , &keep_per_rotation ) ;
                      default_keep_per_rotation = "(user defined)" ;
                    } else {
                      if( strcmp( argv[i] , "-rescue" ) == 0 ) {
                        rescue = 1 ;
                      } else {
                        if( strcmp( argv[i] , "-calculate_grid" ) == 0 ) {
                          i ++ ;
                          if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
                            printf( "Bad command line\n" ) ;
                            exit( EXIT_FAILURE ) ;
                          }
                          calculate = 1 ;
                          default_global_grid_size = "(user defined calculated)" ;
                          sscanf( argv[i] , "%f" , &reverse_calculated_one_span ) ;
                        } else {
                          printf( "Bad command line\n" ) ;
                          exit( EXIT_FAILURE ) ;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

  }

/************/

  /* Rescue option */

  if( rescue == 1 ) {

    printf( "RESCUE mode\n" ) ;

    if( ( ftdock_file = fopen( "scratch_parameters.dat" , "r" ) ) == NULL ) {
      printf( "Could not open scratch_parameters.dat for reading.\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    }

    calculate = 0 ;

    default_global_grid_size = "(read from rescue file)" ;
    default_angle_step = "(read from rescue file)" ;
    default_surface = "(read from rescue file)" ;
    default_internal_value = "(read from rescue file)" ;
    default_electrostatics = "(read from rescue file)" ;

    while( fgets( line_buffer , 99 , ftdock_file ) ) {

      if( strncmp( line_buffer , "Static molecule" , 15 ) == 0 ) sscanf( line_buffer , "Static molecule :: %s" , static_file_name ) ;
      if( strncmp( line_buffer , "Mobile molecule" , 15 ) == 0 ) sscanf( line_buffer , "Mobile molecule :: %s" , mobile_file_name ) ;
      if( strncmp( line_buffer , "Output file name" , 16 ) == 0 ) sscanf( line_buffer , "Output file name :: %s" , output_file_name ) ;
      if( strncmp( line_buffer , "Global grid size" , 16 ) == 0 ) sscanf( line_buffer , "Global grid size :: %d" , &global_grid_size ) ;
      if( strncmp( line_buffer , "Global search angle step" , 24 ) == 0 ) sscanf( line_buffer , "Global search angle step :: %d" , &angle_step ) ;
      if( strncmp( line_buffer , "Global surface thickness" , 24 ) == 0 ) sscanf( line_buffer , "Global surface thickness :: %f" , &surface ) ;
      if( strncmp( line_buffer , "Global internal deterrent value" , 31 ) == 0 ) sscanf( line_buffer , "Global internal deterrent value :: %f" , &internal_value ) ;
      if( strncmp( line_buffer , "Electrostatics                     ::     on" , 44 ) == 0 ) electrostatics = 1 ;    
      if( strncmp( line_buffer , "Electrostatics                     ::    off" , 44 ) == 0 ) electrostatics = 0 ;    
      if( strncmp( line_buffer , "Global keep per rotation" , 25 ) == 0 ) sscanf( line_buffer , "Global keep per rotation :: %d" , &keep_per_rotation ) ;

    }

    fclose( ftdock_file ) ;

    if( ( ftdock_file = fopen( "scratch_scores.dat" , "r" ) ) == NULL ) {
      printf( "Could not open scratch_scores.dat for reading.\nDying\n\n" ) ;
      exit( EXIT_FAILURE ) ;
    }

    fgets( line_buffer , 99 , ftdock_file ) ;

    while( fgets( line_buffer , 99 , ftdock_file ) ) {

      sscanf( line_buffer , "G_DATA %d " , &first_rotation ) ;

    }

    fclose( ftdock_file ) ;

    first_rotation ++ ;

    printf( "Will be starting from rotation %d\n" , first_rotation ) ;

/************/

  } else {

    first_rotation = 1 ;

  }

/************/

  /* Do these things first so that bad inputs will be caught soonest */

  /* Read in Structures from pdb files */
  Static_Structure = read_pdb_to_structure( static_file_name ) ;
  Mobile_Structure = read_pdb_to_structure( mobile_file_name ) ;

  if( Mobile_Structure.length > Static_Structure.length ) {
    printf( "WARNING\n" ) ;
    printf( "The mobile molecule has more residues than the static\n" ) ;
    printf( "Are you sure you have the correct molecules?\n" ) ;
    printf( "Continuing anyway\n" ) ;
  }
  
/************/

  /* Get angles */
  Angles = generate_global_angles( angle_step ) ;

  printf( "Total number of rotations is %d\n" , Angles.n ) ;

/************/

  /* Assign charges */

  if( electrostatics == 1 ) {
    printf( "Assigning charges\n" ) ;
    assign_charges( Static_Structure );
    assign_charges( Mobile_Structure ) ;
    

  }

/************/

  /* Store new structures centered on Origin */

  Origin_Static_Structure = translate_structure_onto_origin( Static_Structure ) ;
  Origin_Mobile_Structure = translate_structure_onto_origin( Mobile_Structure ) ;
  /* Free some memory */

  for( i = 1 ; i <= Static_Structure.length ; i ++ ) {
    free( Static_Structure.Residue[i].Atom ) ;
  }
  free( Static_Structure.Residue ) ;

  for( i = 1 ; i <= Mobile_Structure.length ; i ++ ) {
    free( Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Mobile_Structure.Residue ) ;

/************/

  /* Calculate Grid stuff */

  grid_span = total_span_of_structures( Origin_Static_Structure , Origin_Mobile_Structure ) ;

  if( calculate == 1 ) {
    printf( "Using automatic calculation for grid size\n" ) ;
    global_grid_size = (int)( grid_span / reverse_calculated_one_span ) ;
    if( ( global_grid_size % 2 ) != 0 ) global_grid_size ++ ;
  }

  one_span = grid_span / (float)global_grid_size ;

  printf( "Span = %.3f angstroms\n" , grid_span ) ;
  printf( "Grid size = %d\n" , global_grid_size ) ;
  printf( "Each Grid cube = %.5f angstroms\n" , one_span ) ;

/************/

  /* Memory Allocation */
  int size1 =global_grid_size * global_grid_size * ( 2 * ( global_grid_size / 2 + 1 ) );

  if( ( Scores = ( struct Score * ) malloc ( ( keep_per_rotation + 2 ) * sizeof( struct Score ) ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
     cudaMalloc((void**)&d_Scores,( keep_per_rotation + 2 ) * sizeof( struct Score ));

  if(
    (cudaMalloc((void**)&static_grid,size1* sizeof( cufftReal )) == cudaErrorMemoryAllocation)
    ||
    (cudaMalloc((void**)&mobile_grid,size1* sizeof( cufftReal )) == cudaErrorMemoryAllocation)
    ||
    (cudaMalloc((void**)&convoluted_grid,size1* sizeof( cufftReal ))== cudaErrorMemoryAllocation)
    ) {
    printf( "Not enough memory for surface grids\nUse (sensible) smaller grid size\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }
  

  static_fsg = ( cufftComplex * ) static_grid ;
  mobile_fsg = ( cufftComplex * ) mobile_grid ;
  multiple_fsg = ( cufftComplex * ) convoluted_grid ;

  if( electrostatics == 1 ) {

  if(
    (cudaMalloc((void**)&static_elec_grid,size1* sizeof( cufftReal )) == cudaErrorMemoryAllocation)
    ||
    (cudaMalloc((void**)&mobile_elec_grid,size1* sizeof( cufftReal )) == cudaErrorMemoryAllocation)
    ||
    (cudaMalloc((void**)&convoluted_elec_grid,size1* sizeof( cufftReal ))== cudaErrorMemoryAllocation)
    ) {
    printf( "Not enough memory for surface grids\nUse (sensible) smaller grid size\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  } 
  else {
      /* all ok */
      printf( "Electrostatics are on\n" ) ;
    }

    static_elec_fsg = ( cufftComplex * ) static_elec_grid ;
    mobile_elec_fsg = ( cufftComplex * ) mobile_elec_grid ;
    multiple_elec_fsg = ( cufftComplex * ) convoluted_elec_grid ;

  }

/************/

  /* Create FFTW plans */

  printf( "Creating plans\n" ) ;
  result = cufftPlan3d(&p, global_grid_size , global_grid_size , global_grid_size ,
                               CUFFT_R2C) ;

  result = cufftPlan3d(&pinv, global_grid_size , global_grid_size , global_grid_size ,
                               CUFFT_C2R) ;

/************/

  printf( "Setting up Static Structure\n" ) ;

  /* Discretise and surface the Static Structure (need do only once) */
  
  discretise_structure( Origin_Static_Structure , grid_span , global_grid_size , static_grid,size1);
  printf( "  surfacing grid\n") ;
  dim3 numblocks(((global_grid_size-1)/threadperblock3D.x)+1,((global_grid_size-1)/threadperblock3D.y)+1,((global_grid_size-1)/threadperblock3D.z)+1);
  surface_grid<<<numblocks,threadperblock3D>>>( grid_span , global_grid_size , static_grid , surface , internal_value) ;
  cudaDeviceSynchronize();
  /* Calculate electic field at all grid nodes (need do only once) */
  if( electrostatics == 1 ) {
    electric_field( Origin_Static_Structure , grid_span , global_grid_size , static_elec_grid ) ;
    electric_field_zero_core<<<numblocks,threadperblock3D>>>( global_grid_size , static_elec_grid , static_grid , internal_value) ;
    cudaDeviceSynchronize();
  }

  /* Fourier Transform the static grids (need do only once) */
  printf( "  one time forward FFT calculations\n" ) ;
  result = cufftExecR2C( p , static_grid , static_fsg ) ;
  cudaDeviceSynchronize();
  if( electrostatics == 1 ) {
    result =cufftExecR2C( p , static_elec_grid , static_elec_fsg ) ;
    cudaDeviceSynchronize();
  }

  printf( "  done\n" ) ;

/************/

  /* Store paramaters in case of rescue */

  if( ( ftdock_file = fopen( "scratch_parameters.dat" , "w" ) ) == NULL ) {
    printf( "Could not open scratch_parameters.dat for writing.\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

  fprintf( ftdock_file, "\nGlobal Scan\n" ) ;

  fprintf( ftdock_file, "\nCommand line controllable values\n" ) ;
  fprintf( ftdock_file, "Static molecule                    :: %s\n" , static_file_name ) ;
  fprintf( ftdock_file, "Mobile molecule                    :: %s\n" , mobile_file_name ) ;
  fprintf( ftdock_file, "Output file name                   :: %s\n" , output_file_name ) ;
  fprintf( ftdock_file, "\n" ) ;
  fprintf( ftdock_file, "Global grid size                   :: %6d      %s\n" , global_grid_size , default_global_grid_size ) ;
  fprintf( ftdock_file, "Global search angle step           :: %6d      %s\n" , angle_step , default_angle_step ) ;
  fprintf( ftdock_file, "Global surface thickness           :: %9.2f   %s\n" , surface , default_surface ) ;
  fprintf( ftdock_file, "Global internal deterrent value    :: %9.2f   %s\n" , internal_value , default_internal_value ) ;
  if( electrostatics == 1 ) {
    fprintf( ftdock_file, "Electrostatics                     ::     on      %s\n" , default_electrostatics ) ;
  } else {
    fprintf( ftdock_file, "Electrostatics                     ::    off      %s\n" , default_electrostatics ) ;
  }
  fprintf( ftdock_file, "Global keep per rotation           :: %6d      %s\n" , keep_per_rotation , default_keep_per_rotation ) ;

  fprintf( ftdock_file, "\nCalculated values\n" ) ;
  fprintf( ftdock_file, "Global rotations                   :: %6d\n" , Angles.n ) ;
  fprintf( ftdock_file, "Global total span (angstroms)      :: %10.3f\n" , grid_span ) ;
  fprintf( ftdock_file, "Global grid cell span (angstroms)  :: %10.3f\n" , one_span ) ;

  fclose( ftdock_file ) ;

/************/

  /* Main program loop */

  max_es_value = 0 ;

  printf( "Starting main loop through the rotations\n" ) ;

  for( rotation = first_rotation ; rotation <= Angles.n ; rotation ++ ) {

    printf( "." ) ; 
    

    if( ( rotation % 50 ) == 0 ) printf( "\nRotation number %5d\n" , rotation ) ;

    /* Rotate Mobile Structure */
    Rotated_at_Origin_Mobile_Structure = rotate_structure( Origin_Mobile_Structure , (int)Angles.z_twist[rotation] , (int)Angles.theta[rotation] , (int)Angles.phi[rotation] ) ;

    /* Discretise the rotated Mobile Structure */
    discretise_structure( Rotated_at_Origin_Mobile_Structure , grid_span , global_grid_size , mobile_grid,size1) ;
  
    /* Electic point charge approximation onto grid calculations ( quicker than filed calculations by a long way! ) */
    if( electrostatics == 1 ) {
     
      electric_point_charge( Rotated_at_Origin_Mobile_Structure , grid_span , global_grid_size , mobile_elec_grid) ;

    }
      

    /* Forward Fourier Transforms */
    result = cufftExecR2C( p , mobile_grid , mobile_fsg ) ;
    cudaDeviceSynchronize();
    if( electrostatics == 1 ) {
          result = cufftExecR2C( p , mobile_elec_grid , mobile_elec_fsg ) ;
          cudaDeviceSynchronize();
    }

/************/

    /* Do convolution of the two sets of grids
       convolution is equivalent to multiplication of the complex conjugate of one
       fourier grid with other (raw) one
       hence the sign changes from a normal complex number multiplication
    */
   dim3 numblockconvo(((global_grid_size-1)/threadperblock3D.x)+1,((global_grid_size-1)/threadperblock3D.y)+1,((global_grid_size/2)/threadperblock3D.z)+1);
   convolution<<<numblockconvo,threadperblock3D>>>(static_fsg,multiple_fsg,mobile_fsg,static_elec_fsg,mobile_elec_fsg,multiple_elec_fsg,electrostatics,global_grid_size);
   cudaDeviceSynchronize();
  
  
    /* Reverse Fourier Transform */
    result = cufftExecC2R( pinv , multiple_fsg , convoluted_grid ) ;
    cudaDeviceSynchronize();
    if( electrostatics == 1 ) {
      result = cufftExecC2R( pinv , multiple_elec_fsg , convoluted_elec_grid ) ;
      cudaDeviceSynchronize();
    }
   

/************/

    /* Get best scores */
 
    init_score<<<((keep_per_rotation-1)/1024)+1,1024>>>(d_Scores,keep_per_rotation);
    cudaDeviceSynchronize();
      
    get_score<<<numblocks,threadperblock3D>>>(d_Scores,convoluted_grid,convoluted_elec_grid,electrostatics,keep_per_rotation,global_grid_size);
    cudaDeviceSynchronize();
  
    cudaMemcpy(Scores,d_Scores,( keep_per_rotation + 2 ) * sizeof( struct Score ),cudaMemcpyDeviceToHost);

    
    if( rotation == 1 ) {
      if( ( ftdock_file = fopen( "scratch_scores.dat" , "w" ) ) == NULL ) {
        printf( "Could not open scratch_scores.dat for writing.\nDying\n\n" ) ;
        exit( EXIT_FAILURE ) ;
      }
    } else {
      if( ( ftdock_file = fopen( "scratch_scores.dat" , "a" ) ) == NULL ) {
        printf( "Could not open scratch_scores.dat for writing.\nDying\n\n" ) ;
        exit( EXIT_FAILURE ) ;
      }
    }
  

    for( i = 0 ; i < keep_per_rotation ; i ++ ) {

      max_es_value = min( max_es_value , Scores[i].rpscore ) ;
      fprintf( ftdock_file, "G_DATA %6d   %6d    %7d       %.0f      %4d %4d %4d      %4d%4d%4d\n" ,
                rotation , 0 , Scores[i].score , (double)Scores[i].rpscore , Scores[i].coord[1] , Scores[i].coord[2] , Scores[i].coord[3] ,
                 Angles.z_twist[rotation] , Angles.theta[rotation]  , Angles.phi[rotation] ) ;

    }

    fclose( ftdock_file ) ;

    /* Free some memory */
    for( i = 1 ; i <= Rotated_at_Origin_Mobile_Structure.length ; i ++ ) {
      free( Rotated_at_Origin_Mobile_Structure.Residue[i].Atom ) ;
    }
    free( Rotated_at_Origin_Mobile_Structure.Residue ) ;

  }
   
  /* Finished main loop */

/************/

  /* Free the memory */
   cudaFree(d_Scores);

  result = cufftDestroy(p) ;
  result = cufftDestroy( pinv ) ;

  error = cudaFree( static_grid ) ;
  error =cudaFree( mobile_grid ) ;
  error =cudaFree( convoluted_grid ) ;

  if( electrostatics == 1 ) {
    error =cudaFree( static_elec_grid ) ;
    error =cudaFree( mobile_elec_grid ) ;
    error =cudaFree( convoluted_elec_grid ) ;
  }

  for( i = 1 ; i <= Origin_Static_Structure.length ; i ++ ) {
    free( Origin_Static_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Static_Structure.Residue ) ;

  for( i = 1 ; i <= Origin_Mobile_Structure.length ; i ++ ) {
    free( Origin_Mobile_Structure.Residue[i].Atom ) ;
  }
  free( Origin_Mobile_Structure.Residue ) ;

/************/

  /* Read in all the scores */

  printf( "\nReloading all the scores\n" ) ;

  if( ( ftdock_file = fopen( "scratch_scores.dat" , "r" ) ) == NULL ) {
    printf( "Could not open scratch_scores.dat for reading.\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

  if( ( Scores = ( struct Score * ) realloc ( Scores , ( 1 + keep_per_rotation ) * Angles.n * sizeof( struct Score ) ) ) == NULL ) {
    printf( "Not enough memory left for storing scores\nProbably keeping too many per rotation\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

  kept_scores = 0 ;

  while( fgets( line_buffer , 99 , ftdock_file ) ) {

    sscanf( line_buffer , "G_DATA %d %d %d %f  %d %d %d  %d %d %d" , &id , &id2 , &SCscore , &RPscore ,
                                                                     &x , &y , &z , &z_twist , &theta , &phi ) ;

    Scores[kept_scores].score    = SCscore ;
    Scores[kept_scores].rpscore  = RPscore ;
    Scores[kept_scores].coord[1] = x ;
    Scores[kept_scores].coord[2] = y ;
    Scores[kept_scores].coord[3] = z ;
    Scores[kept_scores].angle[1] = z_twist ;
    Scores[kept_scores].angle[2] = theta ;
    Scores[kept_scores].angle[3] = phi ;

    kept_scores ++ ;

  }

  fclose( ftdock_file ) ;

  kept_scores -- ;

  qsort_scores( Scores , 0 , kept_scores ) ;

/************/

  /* Writing data file */

  if( ( ftdock_file = fopen( output_file_name , "w" ) ) == NULL ) {
    printf( "Could not open %s for writing.\nDying\n\n" , output_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }

  fprintf( ftdock_file, "FTDOCK data file\n" ) ;

  fprintf( ftdock_file, "\nGlobal Scan\n" ) ;

  fprintf( ftdock_file, "\nCommand line controllable values\n" ) ;
  fprintf( ftdock_file, "Static molecule                    :: %s\n" , static_file_name ) ;
  fprintf( ftdock_file, "Mobile molecule                    :: %s\n" , mobile_file_name ) ;
  fprintf( ftdock_file, "\n" ) ;
  fprintf( ftdock_file, "Global grid size                   :: %6d      %s\n" , global_grid_size , default_global_grid_size ) ;
  fprintf( ftdock_file, "Global search angle step           :: %6d      %s\n" , angle_step , default_angle_step ) ;
  fprintf( ftdock_file, "Global surface thickness           :: %9.2f   %s\n" , surface , default_surface ) ;
  fprintf( ftdock_file, "Global internal deterrent value    :: %9.2f   %s\n" , internal_value , default_internal_value ) ;
  if( electrostatics == 1 ) {
    fprintf( ftdock_file, "Electrostatics                     ::     on      %s\n" , default_electrostatics ) ;
  } else {
    fprintf( ftdock_file, "Electrostatics                     ::    off      %s\n" , default_electrostatics ) ;
  }
  fprintf( ftdock_file, "Global keep per rotation           :: %6d      %s\n" , keep_per_rotation , default_keep_per_rotation ) ;

  fprintf( ftdock_file, "\nCalculated values\n" ) ;
  fprintf( ftdock_file, "Global rotations                   :: %6d\n" , Angles.n ) ;
  fprintf( ftdock_file, "Global total span (angstroms)      :: %10.3f\n" , grid_span ) ;
  fprintf( ftdock_file, "Global grid cell span (angstroms)  :: %10.3f\n" , one_span ) ;

  fprintf( ftdock_file, "\nData\n" ) ;
  fprintf( ftdock_file , "Type       ID    prvID    SCscore        ESratio         Coordinates            Angles\n\n" ) ;

  if( electrostatics == 1 ) {

    for( i = 0 ; i <= min( kept_scores , ( NUMBER_TO_KEEP - 1 ) ) ; i ++ ) {

      fprintf( ftdock_file, "G_DATA %6d   %6d    %7d       %8.3f      %4d %4d %4d      %4d%4d%4d\n" ,
               i + 1 , 0 , Scores[i].score , 100 * ( Scores[i].rpscore / max_es_value ) ,
               Scores[i].coord[1] , Scores[i].coord[2] , Scores[i].coord[3] ,
               Scores[i].angle[1] , Scores[i].angle[2] , Scores[i].angle[3] ) ;

    }

  } else {

    for( i = 0 ; i <= min( kept_scores , ( NUMBER_TO_KEEP - 1 ) ) ; i ++ ) {

      fprintf( ftdock_file, "G_DATA %6d   %6d    %7d       %8.3f      %4d %4d %4d      %4d%4d%4d\n" ,
               i + 1 , 0 , Scores[i].score , 0.0 ,
               Scores[i].coord[1] , Scores[i].coord[2] , Scores[i].coord[3] ,
               Scores[i].angle[1] , Scores[i].angle[2] , Scores[i].angle[3] ) ;

    }

  }

  fclose( ftdock_file ) ;
    
/************/

  printf( "\n\nFinished\n\n" ) ;

  return( 0 ) ;

} /* end main */
