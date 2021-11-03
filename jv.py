import Bio.PDB
from Bio import SeqIO
import numpy
from Bio.SeqUtils import seq1

pdb_code = '3AI5'
pdb_filename = "/run/media/sulcjo/sulcjo-data/IOCB/alphafold2_out/jv_send/3AI5.pdb"
pdb_out_filename = "/run/media/sulcjo/sulcjo-data/IOCB/alphafold2_out/jv_send/3AI5_Pyaligned.pdb"

structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)

ref_model = structure[0]
alt_model = structure[1]








"""
def DSSP_assign(pdb_path):
    p = PDBParser()
    # print(pdb_path)
    model_number_pdb = int(modelno) + 1
    structure = p.get_structure('pdb', file=pdb_path)

    model = structure[0]
    dssp = DSSP(model, pdb_path)

    ss = []
    for i, residue in enumerate(sequence):
        try:
            a_key = list(dssp.keys())[i]
            ss.append(dssp[a_key][2])
        except:
            ss.append('x')

    return (ss)

    # Assign colors to SS


def get_sec_str_colors(ss):
    colors = []
    color_assignments = {'-': 'grey', 'E': 'blue', 'T': 'purple', 'S': 'purple', 'H': 'orange'}
    for struc in ss:
        try:
            colors.append(color_assignments[struc])
        except:
            colors.append('white')
    return (colors)


ss_model = DSSP_assign(pdb_path)
ss_ref = DSSP_assign(solved_structure_path)
color_ss = get_sec_str_colors(ss_model)
color_ss_ref = get_sec_str_colors(ss_ref)
"""