import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import operator
from matplotlib import gridspec

def get_pairwise_matrix(path):
    dataframe = pd.read_csv(path)
    return(dataframe)

def get_energy_matrix(path_to_folder, type):
    if type == 'total':
        filename = 'energies_intEnMeanTotal.dat'
    elif type == 'elec':
        filename = 'energies_intEnMeanElec.dat'
    elif type == 'vdw':
        filename = 'energies_intEnMeanVdW.dat'


    with open(f'{path_to_folder}{filename}') as file:
        lines = file.readlines()
        lines_split = [line.split() for line in lines ]
        df = pd.DataFrame(lines_split).astype(float)
    return(df)

def three_to_one(sequence):
    assignment_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    new_sequence = []
    for residue in sequence:
        new_sequence.append(residue[0] + assignment_dict[residue[1:4]] + residue[4:])
    return(new_sequence)

def obtain_sequence(path, aminoacids = 'triple'):
    with open(path) as pdb_file:
        residues = []
        pdb_lines = pdb_file.readlines()
        remember_last_resindex = 0
        for line in pdb_lines:
            split_line = line.split()

            if 'ATOM' in split_line[0] and remember_last_resindex != split_line[5]:
                residues.append(split_line[4] + split_line[3] + split_line[5])
                remember_last_resindex = split_line[5]
        if aminoacids == 'single':
            residues = three_to_one(residues)


        return(residues)

def get_best_pairs(dataframe, dataframe_elec = None, dataframe_vdw = None, pairs = 20):
    statistics_dataframe = dataframe.describe().T

    if dataframe_elec and dataframe_vdw:
        elec_stats = dataframe_elec.describe().T
        elec_stats.columns = ['elec' + column for column in elec_stats.columns]
        vdw_stats = dataframe_vdw.describe().T
        vdw_stats.columns = ['vdw' + column for column in vdw_stats.columns]
        statistics_dataframe = pd.concat([statistics_dataframe, elec_stats, vdw_stats], axis = 1)


    best = statistics_dataframe.nsmallest(n=pairs, columns='mean')
    best = best.drop_duplicates()
    return(best)

""" Calculate total energies for every residue and obtain best residues , from file:///home/sulcjo/Desktop/myomedin/grinn_linux_v110_hf1/24_10E8_model_0_last5ns_noPBCnojump/energies_intEnMeanTotalList.dat"""
def calculate_energy_per_residue(path):
    with open(path) as file:
        lines = file.readlines()
    #Obtain list of unique residues (first column from the file)
    unique_residues = []
    for line in lines:
        residue_candidate = line.split()[0]
        holder = residue_candidate if residue_candidate not in unique_residues else None
        if holder:
            unique_residues.append(holder)
    #Calculate total interaction energies for each unique residue
    residue_dict = {}
    for unique_residue in unique_residues:
        holder = []

        for line in lines:
            if unique_residue == line.split()[0]:
                holder.append(float(line.split()[2]))
        residue_dict[unique_residue] = sum(holder)

    return(residue_dict)

def sort_dictionary(dict):
    sorted_tuples = sorted(dict.items(), key=operator.itemgetter(1))
    sorted_dict = {k: v for k, v in sorted_tuples}
    return(sorted_dict)

def write_resis_into_file(dict, path):
    file = open(path, 'a')

    for key in dict.keys():
        file.write((key + ' : ' + str(dict[key]) + ' kcal/mol' + '\n'))
    return(True)

def plot_heatmap(path_to_folder, sequence, name):

    df = get_energy_matrix(path_to_folder)


    df.replace(0, np.nan, inplace=True)


    ax = sns.heatmap(df, cmap='RdBu', vmin=-10, vmax=10, xticklabels=sequence, yticklabels=sequence)


    ax.set_xticks(ticks = [i for i in range(0, len(sequence), 3)])
    ax.set_xticklabels(sequence[::3])

    ax.set_yticks(ticks = [i for i in range(0, len(sequence), 3)])
    ax.set_yticklabels(sequence[::3])

    plt.xticks(ticks = [i for i in range(0, len(sequence), 3)], labels=sequence[::3])
    plt.yticks(ticks = [i for i in range(0, len(sequence), 3)], labels=sequence[::3])


    plt.title(name)
    ax.set(xlim = (240, 350), ylim = (0, 239))

"""
x = df.index.values
y = df.columns.values
x, y = np.meshgrid(x, y)
z = df.values


fig = plt.figure(figsize=(15,15))
ax3d = fig.add_subplot(projection='3d')
surf = ax3d.plot_surface(x,y,z, cmap='viridis', vmin=-40, vmax=40)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax3d.set(xlim = (240, 350), ylim = (0, 239), zlim = (-40, 40))
#ax3d.set(zlim = (-40, 40))

"""

def obtain_secondary_structures(path_to_folder):
    from Bio.PDB import PDBParser
    from Bio.PDB.DSSP import DSSP
    p = PDBParser()
    structure = p.get_structure('Myo', f'{path_to_folder}system_dry.pdb')
    model = structure[0]
    dssp = DSSP(model, f'{path_to_folder}system_dry.pdb')

    ss = []
    for i, residue in enumerate(sequence):
        try:
            a_key = list(dssp.keys())[i]
            ss.append(dssp[a_key][2])
            #ss.append(dssp[a_key][2] + '\n' + dssp[a_key][1])
        except:
            ss.append('x')

    return(ss)

def get_sec_str_colors(ss):
    colors = []
    color_assignments = {'-' : 'grey', 'E' : 'blue', 'T' : 'purple', 'S' : 'purple', 'H' : 'orange'}
    for struc in ss:
        try:
            colors.append(color_assignments[struc])
        except:
            colors.append('white')
    return(colors)

def obtain_sum_interactions(dataframe):
    sum_interactions = []
    for row_name, row in dataframe.items():
        sum_interactions.append(sum(row))

    return(sum_interactions)

def plot_heatmap_with_resis(path_to_folder, sequence, name, protein_range = (0, 239), ligand_range = (240, 350), vh_lines = True, sec_str = True, mark_resis = None, type = 'total' ):
    # Change step of ticks in heatmap
    ticks_step = 5

    plt.rcParams.update({'font.sans-serif': 'Verdana'})
    # Obtain pairwise interaction dataset, calculate total energies per residue
    df = get_energy_matrix(path_to_folder, type=type)
    sum_interactions = obtain_sum_interactions(df)
    df.replace(0, np.nan, inplace=True)



    # Prepare figure as a grid of 4 subplots
    fig = plt.figure(figsize=(15, 15))
    #plt.tight_layout()
    spec = gridspec.GridSpec(ncols=3, nrows=2, width_ratios=(4, 40, 1), height_ratios=(20, 3))
    ax_cbar = fig.add_subplot(spec[2])
    colors = ['red' if x < 0 else 'grey' if x == 0 else 'blue' for x in sum_interactions]
    # Plot heatmap using Seaborn
    ax3 = fig.add_subplot(spec[1])
    ax3 = sns.heatmap(df, cmap='RdBu', center=0, robust=True,
                      cbar_kws={'label': 'kcal/mol'}, cbar=True, cbar_ax=ax_cbar, linewidths=0.25)
    # ax3 = plt.imshow(df, cmap='RdBu')



    # Plot ax1, upper-right plot showing ligand residues energies
    ax1 = fig.add_subplot(spec[4], sharex=ax3)

    ax1.scatter(sequence, sum_interactions, color=colors, marker='x')

        # Because heatmap is an image, it doesn't fill the whole subplot frame, therefore ax1 wouldn't correspond to heatmap axis
        # Add a little bit to the xlim (beyond 110 residues of MYO) so it compacts the frame
    #ax1.set(xlim=(0, 135))
        # Hide x-axis ticks (resis)
    plt.setp(ax1.get_xticklabels(), visible=False)
    #ax1.xaxis.set_major_locator(ticker.NullLocator())
    #ax1.tick_params(labelrotation=90)
        # Hide right, top and bottom parts of the frame
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_ylabel('kcal/mol')
    # Similar process for the ax2, showing Ab residues left of the heatmap
    ax2 = fig.add_subplot(spec[0], sharey=ax3)
    ax2.scatter(sum_interactions, sequence, color=colors, marker='x')
    #ax2.set(xlim=(-60, 60))
    #ax1.set(ylim=(-60,60))
    #ax2.yaxis.set_major_locator(ticker.NullLocator())
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_xlabel('kcal/mol')
    plt.setp(ax2.get_yticklabels(), visible=False)

    ax2.tick_params(axis='both', which='both', length=0)
    ax1.tick_params(axis='both', which='both', length=0)



    ax3.set_xticks(ticks=[i for i in range(0, len(sequence), ticks_step)])
    ax3.set_yticks(ticks=[i for i in range(0, len(sequence), ticks_step)])
    ax3.set(xlim=ligand_range, ylim=protein_range)

    # Additional stuff
    if vh_lines:
        ax1.vlines(sequence, ymin=[0 for i in sequence], ymax=sum_interactions, color=colors)
        ax2.hlines(sequence, xmin=[0 for i in sequence], xmax=sum_interactions, color=colors)

    if mark_resis:
        for mark_resi in mark_resis:
            if int(mark_resi) < ligand_range[0]:
                ax3.scatter(ligand_range[0] + 2, mark_resi, color='red', marker='s', s=150, linewidths=0.25,
                            edgecolors='black')

            elif int(mark_resi) >= ligand_range[0]:
                ax3.scatter(mark_resi, protein_range[0] + 2, color='red', marker='s', s=150, linewidths=0.25,
                            edgecolors='black')
            else:
                print('Mark residue values out of range')



    if sec_str:
        secondary_structure = obtain_secondary_structures(path_to_folder)
        """ax4 = ax3.twiny()
        ax4.plot(range(ligand_range[1]), [None for i in range(ligand_range[1])])
        ax4.cla()

        ax4.set_xticks(ticks=[i for i in range(0, len(secondary_structure))])
        ax4.set_xticklabels(secondary_structure)
        #ax3.tick_params(axis='x', which='minor', direction='out', length=30, bottom=False, top=True)
        #ax3.tick_params(axis='x', which='major', direction='in', bottom=True, top=False)
        """
        # Plot ligand SS (x-axis)
        ax3.scatter(df.index.values[ligand_range[0]:ligand_range[1]], [-2 for i in df.index.values[ligand_range[0]:ligand_range[1]]], color=get_sec_str_colors(secondary_structure[ligand_range[0]:ligand_range[1]]), s=150, marker='s', linewidths=0.25, edgecolors='black')

        # Plot protein  SS (y-axis)
        ax3.scatter([(ligand_range[0] - 0.5) for i in df.index.values], df.index.values, color=get_sec_str_colors(secondary_structure), s=150, marker='s', linewidths=0.25, edgecolors='black')
        ax3.set(ylim=(protein_range[0] - 2.5, protein_range[1]))
        ax3.set(xlim=(ligand_range[0] - 1, ligand_range[1]))




    if mark_resis or sec_str:
        ax_ss_leg = fig.add_subplot(spec[3])
    else:
        ax_ss_leg = None


    if sec_str:
        # Construct a legend - make a non-visible scatter plot to obtain labels, put it into lower left subplot and hide everything plot-related, then create a legend
        for ss_color, ss_type in zip(['grey', 'blue', 'purple', 'orange'], ['Coil', r'$\beta$-sheet', 'Turn/Bend', r'$\alpha$-Helix']):
            ax_ss_leg.scatter(-1000,-1000,color=ss_color,marker='s',label=ss_type, s=150)

    if mark_resis:
        ax_ss_leg.scatter(-1000,1000,color='red',marker='s',label='Marked residues', s=150)

    if ax_ss_leg:
        ax_ss_leg.legend(ncol=1, loc='upper center', fontsize=20,frameon=False)
        ax_ss_leg.set(xlim=(0,1), ylim=(0.1))
        ax_ss_leg.set_xticks([])
        ax_ss_leg.set_yticks([])
        for spine_pos in ax_ss_leg.spines:
            ax_ss_leg.spines[spine_pos].set_visible(False)





    # Bottom frame line
    ax3.axhline(y=protein_range[0], color='k', linewidth=2)

    # Top frame line
    ax3.axhline(y=protein_range[1], color='k', linewidth=2)

    # Left frame line
    ax3.axvline(x=ligand_range[0], color='k', linewidth=2)

    # Right frame line
    ax3.axvline(x=ligand_range[1], color='k', linewidth=2)





    # Obtain best pairs to put into text
    if type == 'total':
        filename = 'energies_intEnTotal.csv'
    elif type == 'elec':
        filename = 'energies_intEnElec.csv'
    elif type == 'vdw':
        filename = 'energies_intEnVdW.csv'

    pairs_df = get_pairwise_matrix(f'{path_to_folder}{filename}')
    best_pairs = get_best_pairs(pairs_df, pairs=7).T
    best_pairs_text = ''
    for row in best_pairs:
        best_pairs_text += f'{row} {best_pairs[row]["mean"].round(decimals=1)} kcal/mol \n'

    ax3.scatter(x=-100, y=-100, label=best_pairs_text, marker='x', s=0)
    ax3.legend(loc='upper right', fontsize=14, frameon=False)
    #ax3.text(s=best_pairs_text, x=0.5, y=0.5, fontsize=16)
    #ax_best_pairs = fig.add_subplot(spec[5])
    #ax_best_pairs.text(0, 0,s=best_pairs_text, fontsize=20, ha='left')
    #ax_best_pairs.set_xticks([])
    #ax_best_pairs.set_yticks([])
    #for spine_pos in ax_best_pairs.spines:
        #ax_best_pairs.spines[spine_pos].set_visible(False)



    fig.tight_layout()
    plt.subplots_adjust(wspace=0.300)
    plt.suptitle(name + '\n' + type)
    #plt.show()



#analyses = ['24_10E8_model_0_last5ns_noPBCnojump', '25_10E8_model_6', '92_10E8_model_4', '158_10E8_model_8', 'WT_10E8_model_6']
#analysis = '158_10E8_model_8'
#path_to_folder = f'/home/sulcjo/Desktop/myomedin/grinn_linux_v110_hf1/{analysis}/'


models = [6]
variant = 'WT'
for model in models:

    path_to_folder = f'/home/sulcjo/Desktop/myomedin/mmpbsa_gmxmmpbsa/{variant}/{model}/grinn_output/'
    analysis = f'{variant}_{model}'


    sequence = obtain_sequence(f'{path_to_folder}system_dry.pdb', aminoacids='triple')
    #plot_heatmap(path_to_folder, sequence, name='24_10E8_model_0 - 10E8 Fv')
    mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
    plot_heatmap_with_resis(path_to_folder, sequence,
                            name=f'{analysis}', protein_range=(0,239), ligand_range=(240,350),
                            vh_lines = True, sec_str = True, mark_resis=mark_resis, type='total')

    """
    plot_heatmap_with_resis(path_to_folder, sequence,
                            name=f'{analysis}', protein_range=(0,239), ligand_range=(240,350),
                            vh_lines = True, sec_str = True, mark_resis=mark_resis, type='elec')
    
    plot_heatmap_with_resis(path_to_folder, sequence,
                            name=f'{analysis}', protein_range=(0,239), ligand_range=(240,350),
                            vh_lines = True, sec_str = True, mark_resis=mark_resis, type='vdw')
    """
    plt.savefig(f'/home/sulcjo/Desktop/myomedin/mmpbsa_gmxmmpbsa/{variant}/{model}.png')
plt.show()
