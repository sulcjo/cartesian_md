import Bio.PDB
from Bio import SeqIO
import numpy
from Bio.SeqUtils import seq1


structure_AF2 = Bio.PDB.PDBParser().get_structure('1AMM', '/run/media/sulcjo/sulcjo-data/IOCB/alphafold/pdb_reviewers_paper/AF2/1AMM_1.pdb')
structure_ref = Bio.PDB.PDBParser().get_structure('ref', '/run/media/sulcjo/sulcjo-data/IOCB/alphafold/pdb_reviewers_paper/AF2/1amm.pdb1')


ref_model = structure_ref[0]
alt_model = structure_AF2[0]









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
