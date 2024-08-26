# Interactions-search
Description:

This Python script identifies close interactions between a protein and a ligand from provided PDB files. Given the PDB files for the protein and the ligand, the script analyzes spatial proximity to detect potential interactions between the ligand and a specific chain of the protein.
Inputs:

    -l : The PDB file containing the ligand structure.
    -r : The PDB file containing the protein structure.
    -c : The chain identifier within the protein PDB that interacts with the ligand.

Functionality:

    The script reads the PDB files for both the protein and the ligand.
    It extracts the coordinates of atoms in the specified protein chain and the ligand.
    By calculating the distances between these atoms, the script identifies close interactions based on a predefined cutoff distance.
    The results include a list of atoms from both the protein and the ligand that are within the interaction range, along with the corresponding distances.

Output:

    A list of interactions detected between the specified protein chain and the ligand, including the atoms involved and their respective distances.

This tool is useful for analyzing protein-ligand interactions in structural biology, aiding in the study of binding sites and interaction dynamics.
