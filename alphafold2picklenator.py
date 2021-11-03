import pickle
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.DSSP import DSSP
"""
handle = 'run1_pdztrp'
models = [f'result_model_{i}.pkl' for i in [1,2,3,4,5]]
models_pdbs = [f'relaxed_model_{i}.pdb' for i in range(1,6)]
path = '/home/sulcjo/Desktop/IOCB/pdz_trpc/alphafold2_out/pdz_trp/outputs-8897441.meta-pbs.metacentrum.cz/pdztrp'
solved_structure_path = '/home/sulcjo/Desktop/IOCB/pdz_trpc/NMR_structures/pdztrp/pdztrp_nmr_0.1.pdb'
"""

handle = 'run1_trppdz'
models = [f'result_model_{i}.pkl' for i in [1,2,3,4,5]]
models_pdbs = [f'relaxed_model_{i}.pdb' for i in range(1,6)]
path = '/home/sulcjo/Desktop/IOCB/pdz_trpc/alphafold2_out/trp_pdz/outputs-8910776.meta-pbs.metacentrum.cz/trppdz'
solved_structure_path = '/home/sulcjo/Desktop/IOCB/pdz_trpc/NMR_structures/trppdz/trppdz_nmr_0.3.pdb'



plt.rcParams['font.size'] = 12




"""
Prediction confidence plots
"""

# Unpickle AlphaProt2 output and plot

nrows = 3
ncols = 2

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(40,20))
row = 0
col = 0

for modelno, model in enumerate(models):
    with open(path+'/'+model, 'rb') as model_file:
        pickletini = pickle.load(model_file)
    y_data = pickletini['plddt']

    pdb_path = path+f'/relaxed_model_{int(modelno)+1}.pdb'

    # Get pdb sequence
    with open(pdb_path) as pdb_file:
        residues = []
        pdb_lines = pdb_file.readlines()
        remember_last_resindex = 0
        for line in pdb_lines:
            split_line = line.split()

            if 'ATOM' in split_line[0] and remember_last_resindex != split_line[5]:
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


    # Assign DSSP

    def DSSP_assign(pdb_path):

        p = PDBParser()
        #print(pdb_path)
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



    x_data = sequence



    if modelno == nrows:
        row = 0
        col += 1

    axs[row][col].plot(sequence, y_data)
    axs[row][col].set_title(f'plddt score for {model}')
    #axs[row][col].set_xticklabels(sequence, rotation=90)
    #Top - Reference Sstructure
    axs[row][col].scatter(sequence, [106 for i in range(0, len(sequence))],
                color=color_ss_ref, s=150, marker='s', linewidths=0.25,
                edgecolors='black')

    #Bottom - Model Sstructure
    axs[row][col].scatter(sequence, [102 for i in range(0, len(sequence))],
                          color=color_ss, s=150, marker='v', linewidths=0.25,
                          edgecolors='black')



    axs[row][col].set(xlim=(0,len(sequence)+1) , ylim=(min(y_data)-10,112))
    axs[row][col].set_xticklabels([label for label in sequence[::2]])
    axs[row][col].set_xticks(ticks=[tick for tick in range(0, len(sequence), 2)])
    plt.setp(axs[row][col].get_xticklabels(), rotation=90)
    row += 1

for ss_color, ss_type in zip(['grey', 'blue', 'purple', 'orange'],
                                 ['Coil', r'$\beta$-sheet', 'Turn/Bend', r'$\alpha$-Helix']):
    axs[2][1].scatter(-100,-100, color=ss_color, marker='s', label=ss_type, s=15)
    axs[2][1].set(xlim=(-1,1),ylim=(-1,1))
    axs[2][1].axis('off')

axs[2][1].scatter(-100, -100, color='black', marker='s', label='Ref. str.')
axs[2][1].scatter(-100, -100, color='black', marker='v', label='Mod. str.')
axs[2][1].legend(loc='lower left')



plt.subplots_adjust(hspace=0.400)
plt.title('plddt for model')
plt.savefig(path+f'/{handle}_{model}_confidence.png')

"""
Align structures, plot RMSD as function of resi
https://gist.github.com/andersx/6354971
"""

p = PDBParser()
ref_structure = p.get_structure('ref_structure', solved_structure_path)
ref_atoms = []
atoms_to_be_aligned = range(1, len(sequence)+1)


"""
Plot distances
"""
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(40,20))

row = 0
col = 0


for modelpdbno, model_pdb in enumerate(models_pdbs):
    sample_structure = p.get_structure(model_pdb, path+'/'+model_pdb)
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]
    ref_atoms = []
    sample_atoms = []

    for ref_chain in ref_model:
        for ref_res in ref_chain:
            # Check if residue number ( .get_id() ) is in the list
            if ref_res.get_id()[1] in atoms_to_be_aligned:
            # Append CA atom to list
                ref_atoms.append(ref_res['CA'])
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned:
                sample_atoms.append(sample_res['CA'])

    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())
    total_model_rms = (super_imposer.rms)
    print(f'{model_pdb} Ca rmsd to ref structure is {total_model_rms}')

    """
    Get distance between ref and predicted model residues
    """
    distances = []
    for ref_chain, sample_chain in zip(ref_model, sample_model):
        for ref_resi, sample_resi in zip(ref_chain, sample_chain):
            distances.append(abs(ref_resi['CA']-sample_resi['CA']))


    """
    Plot it
    """
    #print(modelno)
    #print(row)
    #print(col)
    if modelpdbno == nrows:
        row = 0
        col += 1

    upper_lim_plot = max(distances)

    axs[row][col].plot(sequence, distances,color='red')
    axs[row][col].set_title(f'Abs CA distance ref/model for {model_pdb}')
    # axs[row][col].set_xticklabels(sequence, rotation=90)
    # Top - Reference Sstructure
    axs[row][col].scatter(sequence, [upper_lim_plot+2 for i in range(0, len(sequence))],
                          color=color_ss_ref, s=150, marker='s', linewidths=0.25,
                          edgecolors='black')

    # Bottom - Model Sstructure
    axs[row][col].scatter(sequence, [upper_lim_plot+1 for i in range(0, len(sequence))],
                          color=color_ss, s=150, marker='v', linewidths=0.25,
                          edgecolors='black')

    # Bottom - Model Sstructure


    axs[row][col].set(xlim=(0, len(sequence) + 1), ylim=(-0.1, upper_lim_plot+3))
    axs[row][col].set_ylabel('|CA-CA dist| / A')
    #axs[row][col].set_xlabel('Res no')
    plt.setp(axs[row][col].get_xticklabels(), rotation=90)
    axs[row][col].set_xticklabels([label for label in sequence[::2]])
    axs[row][col].set_xticks(ticks=[tick for tick in range(0, len(sequence), 2)])
    row += 1

for ss_color, ss_type in zip(['grey', 'blue', 'purple', 'orange'],
                                 ['Coil', r'$\beta$-sheet', 'Turn/Bend', r'$\alpha$-Helix']):
    axs[2][1].scatter(-100,-100, color=ss_color, marker='s', label=ss_type, s=15)
    axs[2][1].set(xlim=(-1,1),ylim=(-1,1))
    axs[2][1].axis('off')

axs[2][1].scatter(-100, -100, color='black', marker='s', label='Ref. str.')
axs[2][1].scatter(-100, -100, color='black', marker='v', label='Mod. str.')
axs[2][1].legend(loc='lower left')
plt.setp(axs[row][col].get_xticklabels(), rotation=90)



plt.subplots_adjust(hspace=0.400)
plt.title(f'distances for model {model_pdb}')
plt.savefig(path+f'/{handle}_{model}_dist.png')



plt.show()
