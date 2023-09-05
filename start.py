# Open the file in read mode
import subprocess
import os

# Replace 'your_command' with the actual command you want to execute

# Run the command in the terminal


with open('names.txt', 'r') as file:
    # Read lines and strip newline characters
    names = [line.strip() for line in file.readlines()]

# Print the names
base_dir = "/mnt/c/Users/adiso/OneDrive/Desktop/DNA-Docking/"

for name in names:
    name_dir = os.path.join(base_dir, "Test_Data_CUDA", name)
    os.chdir(name_dir)
    command="cp " + name + "* ../../"
    subprocess.run(command, shell=True)
    os.chdir(base_dir)
    command="./ftdock -noelec -mode 1 -static " + name + "_Protein_rotated.pdb" + " -mobile " + name + "_DNA_rotated.pdb"
    subprocess.run(command, shell=True)
    command = "mkdir test"
    subprocess.run(command, shell=True)
    command = "cp " + name + "* ftdock_global.dat build test/"
    subprocess.run(command, shell=True)
    os.chdir(os.path.join(base_dir, "test"))
    command = "./build -b1 1 -in ftdock_global.dat"
    subprocess.run(command, shell=True)
    command = "rm build ftdock_global.dat " + name + "_Protein_rotated.pdb " + name + "_DNA_rotated.pdb"
    subprocess.run(command, shell=True)
    os.chdir(base_dir)
    command = "python rmsd_new.py /mnt/c/Users/adiso/OneDrive/Desktop/DNA-Docking/test/ " +name+".pdb"
    subprocess.run(command, shell=True)
    os.chdir(os.path.join(base_dir, "test"))
    command = "cp rmsdGPU.txt " +name+"_rmsd.csv ../ftdock_global.dat ../Test_Data_CUDA/"+name+"/"
    subprocess.run(command, shell=True)
    os.chdir(base_dir)
    command="rm -rf test"
    subprocess.run(command, shell=True)
    command="rm "+name+"*"
    subprocess.run(command, shell=True)

