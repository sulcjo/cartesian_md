import pickle
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.DSSP import DSSP
import Bio

handle = '3AI5'
ref_path = '/run/media/sulcjo/sulcjo-data/IOCB/alphafold/alphafold2_out/python/3AI5_.pdb'
model_path = '/run/media/sulcjo/sulcjo-data/IOCB/alphafold/alphafold2_out/python/3AI5_af.pdb'

### This will get the AA sequence from PDB, but only from the first structure. If there's more, it will simply ignore any additional structures
### which is handled by the condition that subsequent residue index can't be lower (i.e. 1) than the last
with open(model_path) as pdb_file:
    residues = []
    pdb_lines = pdb_file.readlines()
    remember_last_resindex = 0
    for line in pdb_lines:
        split_line = line.split()

        if 'ATOM' in split_line[0] and remember_last_resindex != split_line[5] and int(remember_last_resindex) < int(split_line[5]):
            residues.append(split_line[4] + split_line[3] + split_line[5])
            remember_last_resindex = split_line[5]
    pdb_file.seek(0)
    pdb_file = pdb_file.read()
sequence = residues

# Transform into FASTA format
def three_to_one(sequence):
    assignment_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                       'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                       'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                       'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    modified_sequence = []
    for residue in sequence:
        modified_sequence.append(assignment_dict[residue[1:4]] + residue[4:])
    return (modified_sequence)
sequence = three_to_one(sequence)


p = PDBParser(QUIET=True)
ref_structure = p.get_structure(handle, ref_path)
model_structure = p.get_structure(f'{handle}_model', model_path)

ref_model = ref_structure[0]
sample_model = model_structure[0]




ref_atoms = []
sample_atoms = []


atoms_to_be_aligned = range(1, len(sequence)+1)

for sample_chain in sample_model:
    for sample_res in sample_chain:
        if sample_res.get_id()[1] in atoms_to_be_aligned:
            sample_atoms.append(sample_res['CA'])

for ref_chain in ref_model:
    for ref_res in ref_chain:
        # Check if residue number ( .get_id() ) is in the list
        if ref_res.get_id()[1] in atoms_to_be_aligned:
            # Append CA atom to list
            ref_atoms.append(ref_res['CA'])





super_imposer = Superimposer()
super_imposer.set_atoms(ref_atoms, sample_atoms)
super_imposer.apply(sample_model.get_atoms())
total_model_rms = (super_imposer.rms)
print(f'{handle} Ca rmsd to ref structure is {total_model_rms} with rottrans {super_imposer.rottran}')