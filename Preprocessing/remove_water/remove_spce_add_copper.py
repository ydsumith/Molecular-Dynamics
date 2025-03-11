import numpy as np

# Global variables to store LAMMPS data
num_atoms = 0  # Total number of atoms
num_atom_types = 0  # Total number of atom types
num_bonds = 0  # Total number of bonds
num_bond_types = 0  # Total number of bond types
num_angles = 0  # Total number of angles
num_angle_types = 0  # Total number of angle types
box = np.zeros((3, 2))  # Box dimensions (xlo, xhi, ylo, yhi, zlo, zhi)
masses = {}  # Dictionary to store atom masses
atoms = []  # List to store atom coordinates
atom_types = []  # List to store atom types
bonds = []  # List to store bond information
angles = []  # List to store angle information
molecule_ids = []  # List to store molecule IDs
charges = []  # List to store atom charges

def write_lammps_data(filename):
    """
    Writes the filtered LAMMPS data to a new file.
    
    Args:
        filename (str): The name of the output file.
    """
    with open(filename, 'w') as file:
        # Write header and counts
        file.write(f"\n{num_atoms} atoms\n")
        file.write(f"{num_atom_types} atom types\n")
        file.write(f"{num_bonds} bonds\n")
        file.write(f"{num_bond_types} bond types\n")
        file.write(f"{num_angles} angles\n")
        file.write(f"{num_angle_types} angle types\n\n")
        
        # Write box dimensions
        file.write(f"{box[0, 0]} {box[0, 1]} xlo xhi\n")
        file.write(f"{box[1, 0]} {box[1, 1]} ylo yhi\n")
        file.write(f"{box[2, 0]} {box[2, 1]} zlo zhi\n\n")
        
        # Write masses
        file.write("Masses\n\n")
        for atom_type, mass in masses.items():
            file.write(f"{atom_type} {mass}\n")
        
        # Write atoms
        file.write("\n Atoms\n\n")
        for i, (x, y, z) in enumerate(atoms):
            file.write(f"{i+1} {molecule_ids[i]} {atom_types[i]} {charges[i]} {x} {y} {z}\n")

        # Write bonds
        file.write("\n Bonds\n\n")
        for i, (atom1, atom2) in enumerate(bonds):
            file.write(f"{i+1} 1 {atom1} {atom2}\n")
        
        # Write angles
        file.write("\n Angles\n\n")
        for i, (atom1, atom2, atom3) in enumerate(angles):
            file.write(f"{i+1} 1 {atom1} {atom2} {atom3}\n")

def filter_molecules_by_region(z_min, z_max, cylinder_center, cylinder_radius):
    """
    Removes entire molecules if any atom is outside the given z-range or outside the specified cylindrical region.
    """
    global atoms, atom_types, bonds, angles, num_atoms, num_bonds, num_angles, molecule_ids, charges

    print(num_atoms)

    new_atoms = []  # New list to store filtered atoms
    new_atom_types = []  # New list to store filtered atom types
    new_molecule_ids = []  # New list to store filtered molecule IDs
    new_bonds = []  # New list to store filtered bonds
    new_angles = []  # New list to store filtered angles
    new_charges = []  # New list to store filtered charges
    new_atom_id = 1  # Counter for new atom IDs
    new_mol_id = 1  # Counter for new molecule IDs
    
    # Iterate over atoms in groups of 3 (assuming water molecules)
    for i in range(0, num_atoms, 3):
        x = atoms[i][0]
        y = atoms[i][1]
        z = atoms[i][2]
        
        # Check if the atom is within the z-range and cylindrical region
        if (z_min <= z <= z_max and np.sqrt((x - cylinder_center[0])**2 + (y - cylinder_center[1])**2) <= cylinder_radius) or (z > z_max) or (z < z_min):
            # Add the water molecule to the new lists
            new_atoms.append(atoms[i])
            new_atoms.append(atoms[i+1])
            new_atoms.append(atoms[i+2])

            # Assign atom types (2 for oxygen, 1 for hydrogen)
            new_atom_types.append(2)
            new_atom_types.append(1)
            new_atom_types.append(1)

            # Assign charges (-0.8476 for oxygen, 0.4238 for hydrogen)
            new_charges.append(-0.8476)
            new_charges.append(0.4238)
            new_charges.append(0.4238)
            
            # Add bonds (O-H bonds)
            new_bonds.append([new_atom_id, new_atom_id+1])
            new_bonds.append([new_atom_id, new_atom_id+2])
            
            # Add angle (H-O-H angle)
            new_angles.append([new_atom_id+1, new_atom_id, new_atom_id+2])

            # Assign molecule IDs
            new_molecule_ids.append(new_mol_id)
            new_molecule_ids.append(new_mol_id)
            new_molecule_ids.append(new_mol_id)
            
            # Increment counters
            new_atom_id += 3
            new_mol_id += 1
            
        if i < 5:
            print(new_charges[i])

    # Update global variables with filtered data
    atoms = np.array(new_atoms)
    atom_types = np.array(new_atom_types)
    molecule_ids = np.array(new_molecule_ids)
    bonds = np.array(new_bonds)
    angles = np.array(new_angles)
    charges = np.array(new_charges)
    num_atoms = len(atoms)
    num_bonds = len(bonds)
    num_angles = len(angles)

    print("natoms: ", num_atoms)

def add_copper(filename):
    """
    Adds copper atoms from a LAMMPS data file to the existing system.
    """
    global num_atoms, num_atom_types, num_bonds, num_bond_types, num_angles, num_angle_types
    global box, masses, atoms, atom_types, bonds, angles, molecule_ids, charges

    with open(filename, 'r') as file:
        lines = file.readlines()

    in_atoms = False
    new_atoms = []  # New list to store copper atoms
    new_atom_types = []  # New list to store copper atom types
    new_molecule_ids = []  # New list to store copper molecule IDs
    new_charges = []  # New list to store copper charges
    
    # Find the highest current atom ID to offset new atoms correctly
    max_atom_id = num_atoms

    for line in lines:
        if "atom types" in line:
            num_atom_types = max(num_atom_types, int(line.split()[0]))  # Update number of atom types
        elif "Atoms" in line:
            in_atoms = True
            continue  # Skip the "Atoms" header line
        elif line.strip() == "" or line.startswith("#"):
            continue
        
        if in_atoms:
            parts = line.split()
            if len(parts) < 4:
                continue  # Skip malformed lines
            
            atom_id = max_atom_id + int(parts[0])  # Offset new atom IDs
            atom_type = 3  # Copper atom type
            x, y, z = map(float, parts[2:5])
            z = z + 20.0  # Shift copper atoms along the z-axis
            
            # Append copper atom data
            new_atoms.append([x, y, z])
            new_atom_types.append(atom_type)
            new_molecule_ids.append(0)  # No molecular ID for copper
            new_charges.append(0.0)  # Assume neutral copper atoms

    
    # Append the copper atoms to the existing atoms
    atoms = np.append(atoms, new_atoms, axis=0)
    atom_types = np.append(atom_types, new_atom_types)
    molecule_ids = np.append(molecule_ids, new_molecule_ids)
    charges = np.append(charges, new_charges)

    # Update num_atoms after adding copper atoms
    num_atoms = len(atoms)

def read_lammps_data(filename):
    """
    Reads LAMMPS data from a file and stores it in global variables.
    """
    global num_atoms, num_atom_types, num_bonds, num_bond_types, num_angles, num_angle_types
    global box, masses, atoms, atom_types, bonds, angles, molecule_ids, charges

    with open(filename, 'r') as file:
        lines = file.readlines()

    # Flags to identify sections
    in_masses = False
    in_atoms = False
    in_bonds = False
    in_angles = False

    for line in lines:
        if "atoms" in line:
            num_atoms = int(line.split()[0])
        elif "atom types" in line:
            num_atom_types = int(line.split()[0])
        elif "bonds" in line:
            num_bonds = int(line.split()[0])
        elif "bond types" in line:
            num_bond_types = int(line.split()[0])
        elif "angles" in line:
            num_angles = int(line.split()[0])
        elif "angle types" in line:
            num_angle_types = int(line.split()[0])
        elif "xlo xhi" in line:
            box[0] = list(map(float, line.split()[:2]))
        elif "ylo yhi" in line:
            box[1] = list(map(float, line.split()[:2]))
        elif "zlo zhi" in line:
            box[2] = list(map(float, line.split()[:2]))
        elif "Masses" in line:
            in_masses = True
            in_atoms = False
            in_bonds = False
            in_angles = False
            continue
        elif "Atoms" in line:
            in_masses = False
            in_atoms = True
            in_bonds = False
            in_angles = False
            continue
        elif "Bonds" in line:
            in_masses = False
            in_atoms = False
            in_bonds = True
            in_angles = False
            continue
        elif "Angles" in line:
            in_masses = False
            in_atoms = False
            in_bonds = False
            in_angles = True
            continue

        if in_masses:
            if line.strip() == "" or line.startswith("#"):
                continue
            atom_type, mass = map(float, line.split()[:2])
            masses[int(atom_type)] = mass

        if in_atoms:
            if line.strip() == "" or line.startswith("#"):
                continue
            atom_id, mol_id, atom_type, charge, x, y, z = map(float, line.split()[:7])
            atoms.append([x, y, z])
            atom_types.append(int(atom_type))
            molecule_ids.append(int(mol_id))

        if in_bonds:
            if line.strip() == "" or line.startswith("#"):
                continue
            bond_id, bond_type, atom1, atom2 = map(int, line.split()[:4])
            bonds.append([atom1, atom2])

        if in_angles:
            if line.strip() == "" or line.startswith("#"):
                continue
            angle_id, angle_type, atom1, atom2, atom3 = map(int, line.split()[:5])
            angles.append([atom1, atom2, atom3])

    print("Natoms: ", num_atoms)
    atoms = np.array(atoms)
    atom_types = np.array(atom_types)
    molecule_ids = np.array(molecule_ids)
    bonds = np.array(bonds)
    angles = np.array(angles)

# Example usage
z_min, z_max = 17.0, 172.0
cylinder_center = (30.0, 30.0)
cylinder_radius = 8.0

read_lammps_data("spcewater_2nm.lammps")
print("Natoms: ", num_atoms)

filter_molecules_by_region(z_min, z_max, cylinder_center, cylinder_radius)

add_copper("coppertube_2nm.lmp")

write_lammps_data("filtered_output.lammps")
print("Filtered data written to filtered_output.lammps")
