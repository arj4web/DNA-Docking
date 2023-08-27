import numpy as np
import os
import pandas as pd
import sys
import subprocess

def centroid(p):
    
    """ Returns the centroid of the array """
    
    return p.mean(axis = 0)

def rmsd(p, q):
    
    """ Calculates the RMSD of the metrices P and Q. P and Q have to of same dimension. """
    
    diff = p - q
    return np.sqrt((diff * diff).sum() / p.shape[0])

def kabsch(p, q):
    
    """ Rotates p onmtyo q generating covariance and opptinmal matrix """
    
    c = np.dot(np.transpose(p), q)
    v, s, w = np.linalg.svd(c)
    d = (np.linalg.det(v) * np.linalg.det(w)) < 0.0
    if d:
        s[-1] = -s[-1]
        v[:-1] = -v[:-1]
    u = np.dot(v, w)
    return u

def kabsch_rotation(p, q):
    
    """ Rotates matrix p onto q """
    u = kabsch(p, q)
    
    p = np.dot(p, u)
    return p

def kabsch_rmsd(p, q)->float:
    
    """ Calculates the RMSD using Kabsch Algorithm """
    
    # Translate the matrix onto origin
    
    q = q - centroid(q)
    p = p - centroid(p)
    
    p = kabsch_rotation(p, q)
    val = rmsd(p, q)
    return val

def read_pdb(filename:str):
    
    """ Reads the specified PDB's backbone atoms and returns their coordinate as a matrix """
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    atoms = [line for line in lines if ((line[0:4] == 'ATOM'))]
    backbone_atoms = []
    for atom in atoms:
        backbone_atoms.append(atom)
    no_of_atoms = len(backbone_atoms)
    chains = []
    coordinate = []

    for line in backbone_atoms:
        chain = line[20:22].strip()
        coordinate.append([ float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ])
        if(chain not in chains):
            chains.append(chain)
            
    return coordinate, no_of_atoms, chains

def cal_rmsd(file1:str, file2:str):
    
    """ The main programme """
    
    coord_1, no_atoms_1, chains_1 = read_pdb(file1)
    coord_2, no_atoms_2, chains_2 = read_pdb(file2)
    
    if( (no_atoms_1 == no_atoms_2) and  (len(chains_1) == len(chains_2)) ):
        print('Calculating...')
        
    else:
        print(f'Given two PDBs are not equivalent\n Programme Terminating...\nFile1: {file1},File2: {file2}')
        return(999)
    
    p = np.array(coord_1)
    q = np.array(coord_2)
    
    val = kabsch_rmsd(p, q)
    return (round(val, 2))

###############################################################
def similar_pdb(file_1, file_2):
    atoms_1 = []
    atoms_2 = []
    chains_1 = []
    chains_2 = []
    modified_pdb = []
    with open(file_1, "r") as f1:
        lines_1 = f1.readlines()
        for line in lines_1:
            if(line[:6] == 'ENDMDL'):
                break
            if(line[0:4] == 'ATOM'):
                atoms_1.append(line)
    for atom in atoms_1:
        chain = atom[21:22].strip()
        if(chain not in chains_1):
            chains_1.append(chain)

    with open(file_2) as f2:
        lines_2 = f2.readlines()
        for line in lines_2:
            if(line[:6] == 'ENDMDL'):
                break
            if(line[0:4] == 'ATOM'):
                atoms_2.append(line)
    for atom in atoms_2:
        chain = atom[21:22].strip()
        if(chain not in chains_2):
            chains_2.append(chain)

    if(len(atoms_1) == len(atoms_2)):
        if( chains_1[0] == chains_2[1] and chains_1[1] == chains_2[0]):
            different_chains = [[]  for _ in range(len(chains_1))]
            for i in range(len(chains_1)):
                for atom in atoms_1:
                    if atom[20:22].strip() == chains_1[i]:
                        different_chains[i].append(atom) 
            for line in different_chains[1]:
                modified_pdb.append(line)
            for line in different_chains[0]:
                modified_pdb.append(line)
        else:
            modified_pdb = atoms_1
    with open(file_1.split('.')[0]+'_modified.pdb', 'w') as f:
        for atom in modified_pdb:
            f.write(atom)

################### Main Program ##############################
def main_prog(file_path:str, filename:str):
    try:
        os.chdir(file_path)
        files = os.listdir()
        pdbfiles = [file for file in files if file[-4:] == '.pdb']
        result = pd.DataFrame()
        similar_pdb(filename, 'Complex_1g.pdb')
        file_name = filename.split('.')[0]+'_modified.pdb'
        for file in files:
            val = cal_rmsd(file_name, file)
            result = pd.concat([result, pd.DataFrame({'Complex Name': [file], 'RMSD': [val]})])

        os.chdir(file_path)
        r = result.sort_values(['RMSD'], ascending = True)[0:5]
        with open('rmsdGPU.txt', 'w')as f:
            f.write(str(r))
        new_name = filename.split('.')[0]+'_rmsd.csv'
        result.to_csv(new_name)
        print(f'File Written: {new_name}')
    except:
        print(f'Keep "{filename}" and "Complex_1g.pdb" present in: {file_path}')

if __name__ == "__main__":
    main_prog(sys.argv[1], sys.argv[2])
