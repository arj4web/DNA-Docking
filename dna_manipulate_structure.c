#include "dna_structures.h"

struct DNA_structure read_dna_file(char *filename){
    FILE *fp;
    fp = fopen(filename, "r");
    char line_buffer[85];


    ///// variables for Atom
    int id;
    char name_of_atom[5];
    float x, y, z;
    float occupancy;
    float temp_factor;
    float charge;
    char mol[2];

    //// variables for Nucleotide

    char nucleotide_name[4];
    char chain_no[2];
    char nuCode[6];
    int nc;
    int size;

    // variables for whole DNA Structure
    char dna_name[256];
    int size_of_dna;

    // Counter Variables
    int no_of_atoms = 0;
    int no_of_nucleotide = 0;

    // variable to keep count of nucleotides
    char present_nuCode[6];
    strcpy(present_nuCode, "#");

    struct DNA_structure dna;
    dna.nucleotide = (struct Nucleic_Acid *)malloc(sizeof(struct Nucleic_Acid));

    while (fgets(line_buffer, 85, fp)!=NULL)
    {
        if(strncmp(line_buffer, "ATOM", 4) == 0){
            sscanf( line_buffer +  6 , "%5d" , &id ) ;
            sscanf( line_buffer + 30 , "%8f" , &x ) ;
            sscanf( line_buffer + 38 , "%8f" , &y ) ;
            sscanf( line_buffer + 46 , "%8f" , &z ) ;
            sscanf( line_buffer + 54 , "%6f" , &occupancy ) ;
            sscanf( line_buffer + 60 , "%6f" , &temp_factor ) ;
            sscanf( line_buffer + 82 , "%2d" , &nc ) ;

            strncpy( name_of_atom,		line_buffer+12,	4 ) ;
            strncpy( nucleotide_name,		line_buffer+17,	3 ) ;
            strncpy( chain_no,			line_buffer+21,	1 ) ;
            strncpy( nuCode,	line_buffer+22,	5 ) ;
            strncpy( mol,			line_buffer+77,	1 ) ;

            strncpy( name_of_atom + 4,		"\0", 1 ) ;
            strncpy( nucleotide_name + 3,		"\0", 1 ) ;
            strncpy( chain_no + 1,		"\0", 1 ) ;
            strncpy( nuCode + 5,	"\0", 1 ) ;
            strncpy( mol + 1,			"\0", 1 ) ;

            if(strcmp(nuCode, present_nuCode)!=0){ //New Nucleotide
                
                dna.nucleotide[no_of_nucleotide].size = no_of_atoms;

                //reinitialize
                no_of_atoms = 0;
                no_of_nucleotide += 1 ;

                dna.nucleotide = (struct Nucleic_Acid*)realloc(dna.nucleotide, (no_of_nucleotide +1)*sizeof(struct Nucleic_Acid));

                dna.nucleotide[no_of_nucleotide].atom = (struct Atom*)malloc(sizeof(struct Atom));

                strcpy(dna.nucleotide[no_of_nucleotide].chainID, chain_no);
                dna.nucleotide[no_of_nucleotide].nc = nc;
                
                strcpy(dna.nucleotide[no_of_nucleotide].nucleicAcid_name, nucleotide_name);
                strcpy(dna.nucleotide[no_of_nucleotide].nuCode, nuCode);       
            }
            // put the present no of nucleotide
            strcpy(present_nuCode, nuCode);
            //Create New Atom
            no_of_atoms += 1;

            dna.nucleotide[no_of_nucleotide].atom = (struct Atom*)realloc(dna.nucleotide[no_of_nucleotide].atom, (no_of_atoms+1)*sizeof(struct Atom));

            strcpy(dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].atom_name, name_of_atom);
            dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].charge = charge;
            dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].coord[1] = x;
            dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].coord[2] = y;
            dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].coord[3] = z;
            dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].occupancy = occupancy;
            dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].serial = id;
            dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].temp_factor = temp_factor;
            strcpy(dna.nucleotide[no_of_nucleotide].atom[no_of_atoms].molecule, mol);
        }
    }
    
    dna.nucleotide[no_of_nucleotide].size = no_of_atoms;
    dna.length = no_of_nucleotide;
    strcpy(dna.name, filename);
    fclose(fp);
    //printf("%s", dna.nucleotide[1].atom[1].olc);
    return dna;

}

void write_dna_file(struct DNA_structure dna, char *filename){
    FILE *fp;
    fp = fopen(filename, "w");
    for(int nucleotide = 1; nucleotide<=dna.length; nucleotide++){
        for(int atom = 1; atom<=dna.nucleotide[nucleotide].size; atom++){
            fprintf(fp, "ATOM  %5d %4s %3s %1s %5s   %8.3f %8.3f %8.3f %6.2f %6.2f            %1s\n", 
            dna.nucleotide[nucleotide].atom[atom].serial,
            dna.nucleotide[nucleotide].atom[atom].atom_name,
            dna.nucleotide[nucleotide].nucleicAcid_name,
            dna.nucleotide[nucleotide].chainID,
            dna.nucleotide[nucleotide].nuCode,
            dna.nucleotide[nucleotide].atom[atom].coord[1],
            dna.nucleotide[nucleotide].atom[atom].coord[2],
            dna.nucleotide[nucleotide].atom[atom].coord[3],
            dna.nucleotide[nucleotide].atom[atom].occupancy,
            dna.nucleotide[nucleotide].atom[atom].temp_factor,
            dna.nucleotide[nucleotide].atom[atom].molecule
            );

        }
    }
}

struct DNA_structure duplicate_dna_structure(struct DNA_structure dna){
    struct DNA_structure new_dna;
    struct Nucleic_Acid d_nucleotide,nucleotide;

    new_dna.length = dna.length;
    nucleotide = (struct Nucleic_Acid*)malloc(sizeof(struct Nucleic_Acid));
    new_dna.nucleotide = (struct Nucleic_Acid*)malloc(sizeof(struct Nucleic_Acid));

    for(int nucleotide = 1; nucleotide <= dna.length; nucleotide++){
        new_dna.nucleotide = dna.nucleotide;
        for(int atom=1; atom<=dna.nucleotide[nucleotide].size; atom++){
            new_dna.nucleotide[nucleotide].atom[atom]=dna.nucleotide[nucleotide].atom[atom];
        }
    }
    return new_dna;
}

struct DNA_structure translate_dna_structure(struct DNA_structure dna, float x_shift, float y_shift, float z_shift){
    struct DNA_structure new_dna;
    new_dna = duplicate_dna_structure(dna);
    for(int nucleotide=1; nucleotide<=new_dna.length; nucleotide++){
        for(int atom=1; atom<=dna.nucleotide[nucleotide].size; atom++){
            dna.nucleotide[nucleotide].atom[atom].coord[1] += x_shift;
            dna.nucleotide[nucleotide].atom[atom].coord[2] += y_shift;
            dna.nucleotide[nucleotide].atom[atom].coord[3] += z_shift;
        }
    }    
    return new_dna;
}

struct DNA_structure translate_dna_to_origin(struct DNA_structure dna){
    int total_no_of_atoms=0;
    float average_x=0.0, average_y=0.0, average_z=0.0;
    struct DNA_structure new_dna;
    new_dna = duplicate_dna_structure(dna);

    for(int nucleotide=1; nucleotide<=new_dna.length; nucleotide++){
        for(int atom=1; atom<=dna.nucleotide[nucleotide].size; atom++){
            total_no_of_atoms += 1;
            average_x += dna.nucleotide[nucleotide].atom[atom].coord[1];
            average_y += dna.nucleotide[nucleotide].atom[atom].coord[2];
            average_z += dna.nucleotide[nucleotide].atom[atom].coord[3];
        }
    }
    average_x = average_x/(float)total_no_of_atoms;
    average_y = average_y/(float)total_no_of_atoms;
    average_z = average_z/(float)total_no_of_atoms;

    for(int nucleotide=1; nucleotide<=new_dna.length; nucleotide++){
        for(int atom=1; atom<=dna.nucleotide[nucleotide].size; atom++){
            dna.nucleotide[nucleotide].atom[atom].coord[1] -= average_x;
            dna.nucleotide[nucleotide].atom[atom].coord[2] -= average_y;
            dna.nucleotide[nucleotide].atom[atom].coord[3] -= average_z;
        }
    }
    return new_dna;
}

struct DNA_structure rotate_dna(struct DNA_structure dna, int z_twist, int theta, int phi){
    struct DNA_structure new_dna;
    new_dna = duplicate_dna_structure(dna);
    float post_z_twist_x, post_z_twist_y, post_z_twist_z;
    float post_theta_x, post_theta_y, post_theta_z;

    for(int nucleotide=1; nucleotide<=dna.length; nucleotide++){
        for(int atom=1; atom<=dna.nucleotide[nucleotide].size; atom++){
         // Z axis twist
         post_z_twist_x = new_dna.nucleotide[nucleotide].atom[atom].coord[1] * cos(0.017453293*z_twist) - new_dna.nucleotide[nucleotide].atom[atom].coord[1] * sin(0.017453293*z_twist);
         post_z_twist_y = new_dna.nucleotide[nucleotide].atom[atom].coord[2] * sin(0.017453293*z_twist) + new_dna.nucleotide[nucleotide].atom[atom].coord[2] * cos(0.017453293*z_twist);
         post_z_twist_z = new_dna.nucleotide[nucleotide].atom[atom].coord[3];

        // X-Z plane twist
         post_theta_x = post_z_twist_x  * cos(0.017453293*theta) + post_z_twist_x * sin(0.017453293*theta);
         post_theta_y = post_z_twist_y;
         post_theta_z = - post_z_twist_z * sin(0.017453293*theta) + post_z_twist_z * cos(0.017453293*theta);

         // X axis twist
         new_dna.nucleotide[nucleotide].atom[atom].coord[1] = post_theta_x * cos(0.017453293*z_twist) - post_theta_x * sin(0.017453293*z_twist);
         new_dna.nucleotide[nucleotide].atom[atom].coord[2] = post_theta_y * sin(0.017453293*z_twist) + post_theta_y * cos(0.017453293*z_twist);
         new_dna.nucleotide[nucleotide].atom[atom].coord[3] = post_theta_z;

        }
    } 
    return new_dna;
}

struct DNA_structure merge_dnas(struct DNA_structure dna1, struct DNA_structure dna2){
    int new_nucleotide;
    struct DNA_structure new_dna;
    strcpy(new_dna.name, "Complex.pdb");
    new_dna.length = dna1.length + dna2.length;

    new_dna.nucleotide = (struct Nucleic_Acid*)malloc((dna1.length + dna2.length + 1)* sizeof(struct Nucleic_Acid));

    for(int nucleotide=1; nucleotide<=dna1.length; nucleotide++){
        new_dna.nucleotide[nucleotide].atom = (struct Atom*)malloc((dna1.nucleotide[nucleotide].size+1)*sizeof(struct Atom));
        strcpy(new_dna.nucleotide[nucleotide].nucleicAcid_name, dna1.nucleotide[nucleotide].nucleicAcid_name);
        strcpy(new_dna.nucleotide[nucleotide].chainID, dna1.nucleotide[nucleotide].chainID);
        strcpy(new_dna.nucleotide[nucleotide].nuCode, dna1.nucleotide[nucleotide].nuCode);
        new_dna.nucleotide[nucleotide].nc = dna1.nucleotide[nucleotide].nc;
        new_dna.nucleotide[nucleotide].size = dna1.nucleotide[nucleotide].size;

        new_dna.nucleotide[nucleotide].atom = (struct Atom*)malloc((dna1.nucleotide[nucleotide].size+1)*sizeof(struct Atom));

        for(int atom=1; atom<=dna1.nucleotide[nucleotide].size; atom++){
            new_dna.nucleotide[nucleotide].atom[atom] = dna1.nucleotide[nucleotide].atom[atom];
        }
    }

    for(int nucleotide=1; nucleotide<=dna2.length; nucleotide++){
        new_nucleotide = nucleotide + dna1.length;
        strcpy(new_dna.nucleotide[new_nucleotide].nucleicAcid_name, dna2.nucleotide[nucleotide].nucleicAcid_name);
        strcpy(new_dna.nucleotide[new_nucleotide].chainID, dna2.nucleotide[nucleotide].chainID);
        strcpy(new_dna.nucleotide[new_nucleotide].nuCode, dna2.nucleotide[nucleotide].nuCode);
        new_dna.nucleotide[new_nucleotide].nc = dna2.nucleotide[nucleotide].nc;
        new_dna.nucleotide[new_nucleotide].size = dna2.nucleotide[nucleotide].size;

        new_dna.nucleotide[new_nucleotide].atom = (struct Atom*)malloc((dna2.nucleotide[nucleotide].size+1)*sizeof(struct Atom));
        
        for(int atom=1; atom<=dna2.nucleotide[nucleotide].size; atom++){
            new_dna.nucleotide[new_nucleotide].atom[atom] = dna2.nucleotide[nucleotide].atom[atom];
        }
    }
    return new_dna;
}

float radius_of_dna(struct DNA_structure dna){
    float present, largest=0;

    for(int nucleotide=1; nucleotide<=dna.length; nucleotide++){
        for(int atom=1; atom<=dna.nucleotide[nucleotide].size; atom++){
            present = dna.nucleotide[nucleotide].atom[atom].coord[1]*dna.nucleotide[nucleotide].atom[atom].coord[1] + dna.nucleotide[nucleotide].atom[atom].coord[2]*dna.nucleotide[nucleotide].atom[atom].coord[2] + dna.nucleotide[nucleotide].atom[atom].coord[3]*dna.nucleotide[nucleotide].atom[atom].coord[3];
            if(present>largest) largest = present;
        }
    }
    return largest;
}

float span_of_dna(struct DNA_structure dna1, struct DNA_structure dna2){
    return 1+((radius_of_dna(dna1)+radius_of_dna(dna2))*2);
}


int main(){
    struct DNA_structure dna1, dna2, new_dna;
   dna2 = read_dna_file("1.pdb");
   dna1 = read_dna_file("2.pdb");
   //new_dna = merge_dnas(dna2, dna1); 
   //write_dna_file(new_dna, "see.pdb");
   printf("%f", span_of_dna(dna1, dna2));
   return 0;   
}

