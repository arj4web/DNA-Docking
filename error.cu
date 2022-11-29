#include "structures.cuh"

__device__ float pythagoras( float x1 , float y1 , float z1 , float x2 , float y2 , float z2 ) {

  return sqrt( ( ( x1 - x2 ) * ( x1 - x2 ) ) + ( ( y1 - y2 ) * ( y1 - y2 ) ) + ( ( z1 - z2 ) * ( z1 - z2 ) ) ) ;

}
struct Structure read_pdb_to_structure( char *pdb_file_name ) {

/************/

  /* Variables */

  /* Counters */
  int	n_residues ;	/* number of residues */
  int	res_size ;	/* number of atoms in single residue */

  /* File stuff */
  FILE	*pdb_file ;
  char	line_buffer[100] ;

  /* What the data is going into */
  struct Structure		This_Structure ;

  /* Variables from the PDB file */
  int	serial ;
  char		atom_name[5] ;
  char		res_name[4] ;
  char		chainID[2] ;
  char		res_seq_plus_iCode[6] ;
  float		coord_x , coord_y , coord_z ;
  float		occupancy, temp_factor ;
  char		olc[2] ;

  /* Comparison values */
  char	present_res_seq_plus_iCode[6] ;

/************/

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  /* File handling */

  /* Open file */
  printf( "  reading parsed pdb file: %s\n", pdb_file_name ) ;
  if( ( pdb_file = fopen( pdb_file_name, "r" ) ) == NULL ) {
    printf( "This file does not exist here, or is unreadable.\nDying\n\n" ) ;
    exit( EXIT_FAILURE ) ;
  }

/************/

  /* Initialisations */

  /* Counters */
  n_residues = 0 ;
  res_size = 0 ;

  /* Comparison values */
  strcpy( present_res_seq_plus_iCode , ">" ) ;

  /* Memory allocation */
  if( ( This_Structure.Residue = ( struct Amino_Acid * ) malloc ( sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

/************/

  /* Read PDB file */

  /* The Atoms */

  while( fgets( line_buffer, 85, pdb_file ) ) {

    if( strncmp( line_buffer, "ATOM", 4 ) == 0 ) {

      /* Have an ATOM */

      /* Get Values */

      /* the following may seem silly, but sscanf convention means that two
         float fields with no white space between them, where the first is
         less than the maximum field width, mucks up everything.
      */

      sscanf( line_buffer +  6 , "%5d" , &serial ) ;
      sscanf( line_buffer + 30 , "%8f" , &coord_x ) ;
      sscanf( line_buffer + 38 , "%8f" , &coord_y ) ;
      sscanf( line_buffer + 46 , "%8f" , &coord_z ) ;
      sscanf( line_buffer + 54 , "%6f" , &occupancy ) ;
      sscanf( line_buffer + 60 , "%6f" , &temp_factor ) ;
      

      strncpy( atom_name,		line_buffer+12,	4 ) ;
      strncpy( res_name,		line_buffer+17,	3 ) ;
      strncpy( chainID,			line_buffer+21,	1 ) ;
      strncpy( res_seq_plus_iCode,	line_buffer+22,	5 ) ;
      strncpy( olc,			line_buffer+77,	1 ) ;

      strncpy( atom_name + 4,		"\0", 1 ) ;
      strncpy( res_name + 3,		"\0", 1 ) ;
      strncpy( chainID + 1,		"\0", 1 ) ;
      strncpy( res_seq_plus_iCode + 5,	"\0", 1 ) ;
      strncpy( olc + 1,			"\0", 1 ) ;

/************/

      /* New Residue */

      if( strcmp( res_seq_plus_iCode , present_res_seq_plus_iCode ) != 0 ) {

        /* have next residue */

        /* Store old info */
        This_Structure.Residue[n_residues].size = res_size ;

        /* Increment, Reset numbers */
        n_residues ++ ;
        res_size = 0 ;

        /* Memory management */
        if( ( This_Structure.Residue = (struct Amino_Acid * ) realloc ( This_Structure.Residue, ( n_residues + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
          GENERAL_MEMORY_PROBLEM
        }
        if( ( This_Structure.Residue[n_residues].Atom = ( struct Atom * ) malloc ( sizeof_Atom ) ) == NULL ) {
          GENERAL_MEMORY_PROBLEM
        }

        /* Store new info */
        strcpy( This_Structure.Residue[n_residues].res_seq_plus_iCode , res_seq_plus_iCode );
        strcpy( This_Structure.Residue[n_residues].res_name ,           res_name ) ;
        strcpy( This_Structure.Residue[n_residues].chainID ,            chainID ) ;
        strcpy( This_Structure.Residue[n_residues].olc,                 olc ) ;
        
      }

      strcpy( present_res_seq_plus_iCode , res_seq_plus_iCode ) ;

/************/

      /* Put Atoms into Structure */

      res_size ++ ;

      if( ( This_Structure.Residue[n_residues].Atom = ( struct Atom * ) realloc ( This_Structure.Residue[n_residues].Atom, ( res_size + 1 ) * sizeof_Atom ) ) == NULL ) {
        GENERAL_MEMORY_PROBLEM
      }

      This_Structure.Residue[n_residues].Atom[res_size].serial = serial ;
      strcpy( This_Structure.Residue[n_residues].Atom[res_size].atom_name, atom_name ) ;
      This_Structure.Residue[n_residues].Atom[res_size].coord[1] = coord_x ;
      This_Structure.Residue[n_residues].Atom[res_size].coord[2] = coord_y ;
      This_Structure.Residue[n_residues].Atom[res_size].coord[3] = coord_z ;
      This_Structure.Residue[n_residues].Atom[res_size].occupancy = occupancy ;
      This_Structure.Residue[n_residues].Atom[res_size].temp_factor = temp_factor ;

/************/

    }

  } /* got to end of pdb file */

/************/

  /* Clean up */

  This_Structure.Residue[n_residues].size = res_size ;
  This_Structure.length = n_residues ;
  strcpy( This_Structure.ident , pdb_file_name  );

  /* Finish off */

  fclose( pdb_file ) ;

  return This_Structure ;

}
struct Structure duplicate_structure( struct Structure This_Structure ) {

/************/

  /* Variables */
  struct Structure	New_Structure ;

  /* Counters */
  int		residue , atom ;

/************/

  if( ( New_Structure.Residue = ( struct Amino_Acid * ) malloc ( ( This_Structure.length + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

  strcpy( New_Structure.ident , This_Structure.ident ) ;
  New_Structure.length = This_Structure.length ;

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {

    New_Structure.Residue[residue] = This_Structure.Residue[residue] ;

    if( ( New_Structure.Residue[residue].Atom = ( struct Atom * ) malloc ( ( This_Structure.Residue[residue].size + 1 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      New_Structure.Residue[residue].Atom[atom] = This_Structure.Residue[residue].Atom[atom] ;

    }

  }

  return New_Structure ;

/************/

}
__global__ void calculate_error(Amino_Acid *residue1,Amino_Acid *residue2,int ydim,float *RMS)
{
        int residue=threadIdx.y+(blockDim.y*blockIdx.y);
        int atom=threadIdx.x+(blockDim.x*blockIdx.x);
        float error;
        if(residue< ydim){
            if((residue>0)&&(atom>0)&&(atom<=residue1[residue].size)){
            error=pythagoras(residue1[residue].Atom[atom].coord[1],residue1[residue].Atom[atom].coord[2],residue1[residue].Atom[atom].coord[3],residue2[residue].Atom[atom].coord[1],residue2[residue].Atom[atom].coord[2],residue2[residue].Atom[atom].coord[3]);
            error=error/pythagoras(residue1[residue].Atom[atom].coord[1],residue1[residue].Atom[atom].coord[2],residue1[residue].Atom[atom].coord[3],0.0,0.0,0.0);
            int counter=0;
            for(int i=1;i<residue;i++)
            {
                counter+=residue1[i].size;
            }
            RMS[counter+atom]=error;
            
            }
      }
}
int main(int argc, char *argv[])
{
  
    char		*output_file_name ;
    char		*orignal_file ;
    char		*translated_file;
    struct Structure	Orignal_Structure , Translated_Structure ;

    int i=0;
    
  if( ( ( output_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( orignal_file  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ||
      ( ( translated_file  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL ) ) {
    GENERAL_MEMORY_PROBLEM 
  }

    for( i = 1 ; i < argc ; i ++ ) 
    {
        printf("thdfghhjyjjnygere there\n");
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
                
        strcpy( orignal_file , argv[i] ) ;
        printf("thdfghhjyjjnygere there\n"); 
      } else {
        if( strcmp( argv[i] , "-mobile" ) == 0 ) {
       
          i ++ ;
          if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
            printf( "Bad command line\n" ) ;
            exit( EXIT_FAILURE ) ;
          }
          strcpy( translated_file , argv[i] ) ;
          
        } 
        else{ 
            printf( "Bad command line\n" ) ;
            exit( EXIT_FAILURE ) ;
        }
        }
        }
   
}

    Orignal_Structure = read_pdb_to_structure( orignal_file) ;
    Translated_Structure = read_pdb_to_structure(translated_file) ;
printf("there there\n");
    int a=0;
  Orignal_Structure = duplicate_structure( Orignal_Structure ) ;
  Translated_Structure = duplicate_structure( Translated_Structure ) ;
  struct Amino_Acid *Residue1,*d_Residue1, *Residue2,*d_Residue2;
  Residue1 = (struct Amino_Acid*)malloc((Orignal_Structure.length+1)*sizeof(Amino_Acid));
  Residue2 = (struct Amino_Acid*)malloc((Translated_Structure.length+1)*sizeof(Amino_Acid));
  int t=0;
  for (int i = 1; i <=Orignal_Structure.length; i++)
  {
     
    Residue1[i]=Orignal_Structure.Residue[i];
    Residue2[i]=Translated_Structure.Residue[i];
    cudaMalloc((void**)&Residue1[i].Atom,(Orignal_Structure.Residue[i].size+1)*sizeof(struct Atom));
    cudaMalloc((void**)&Residue2[i].Atom,(Translated_Structure.Residue[i].size+1)*sizeof(struct Atom));
    cudaMemcpy(Residue1[i].Atom,Orignal_Structure.Residue[i].Atom,(Orignal_Structure.Residue[i].size+1)*sizeof(struct Atom),cudaMemcpyHostToDevice);
    cudaMemcpy(Residue2[i].Atom,Translated_Structure.Residue[i].Atom,(Translated_Structure.Residue[i].size+1)*sizeof(struct Atom),cudaMemcpyHostToDevice);
    a=max(a,Orignal_Structure.Residue[i].size);
    t+=Orignal_Structure.Residue[i].size;
    
  }
  printf("The value of t is : %d\n",t);
  cudaMalloc((void**)&d_Residue1,(Orignal_Structure.length+1)*sizeof(struct Amino_Acid));
  cudaMalloc((void**)&d_Residue2,(Translated_Structure.length+1)*sizeof(struct Amino_Acid));
  cudaMemcpy(d_Residue1,Residue1,(Orignal_Structure.length+1)*sizeof(struct Amino_Acid),cudaMemcpyHostToDevice);
  cudaMemcpy(d_Residue2,Residue2,(Translated_Structure.length+1)*sizeof(struct Amino_Acid),cudaMemcpyHostToDevice);

  dim3 numblocks((a/threadperblock2D.x)+1,(Orignal_Structure.length/threadperblock2D.y)+1);
  float *d_RMS,*RMS;
  cudaMalloc((void**)&d_RMS,(t+1)*sizeof(float));
  calculate_error<<<numblocks,threadperblock2D>>>(d_Residue1,d_Residue2,Orignal_Structure.length+1,d_RMS);
  cudaDeviceSynchronize();
  RMS=(float *)malloc((t+1)*sizeof(float));
  cudaMemcpy(RMS,d_RMS,(t+1)*sizeof(float),cudaMemcpyDeviceToHost);

  for ( i = 1; i < t+1; i++)
  {
    /* code */
    if(RMS[i]>0)printf("%f\n",RMS[i]);
  }
  
  cudaFree(d_Residue1);
  cudaFree(d_Residue2);
  cudaFree(d_RMS);


}


             