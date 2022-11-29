#include "structures.cuh"

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
            RMS[counter+atom]=error*100;
            
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
        if( strcmp( argv[i] , "-out" ) == 0 ) {
            i ++ ;
            if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
            printf( "Bad command line\n" ) ;
            exit( EXIT_FAILURE ) ;
            }
            strcpy( output_file_name , argv[i] ) ;
    } else {
      if( strcmp( argv[i] , "-primary" ) == 0 ) {
        i ++ ;
        if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
          printf( "Bad command line\n" ) ;
          exit( EXIT_FAILURE ) ;
        }
                
        strcpy( orignal_file , argv[i] ) ;
      } else {
        if( strcmp( argv[i] , "-secondary" ) == 0 ) {
       
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

  printf("Following are the notable errors in the 2 PDB's\n");
  for ( i = 1; i < t+1; i++)
  {
    if(RMS[i]>0)printf("%f\n",RMS[i]);
  }
  
  cudaFree(d_Residue1);
  cudaFree(d_Residue2);
  cudaFree(d_RMS);


}


             