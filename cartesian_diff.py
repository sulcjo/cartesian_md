import sys, argparse, os
# ujson is recommended but not needed, speeds up json vectors loading
try:
    import ujson as json
except ImportError:
    try:
        import simplejson as json
    except ImportError:
        import json

import numpy as np
import pandas as pd
from alive_progress import alive_bar
import scipy.spatial as ss
from cartesian import __prepare_matplotlib
import re
import matplotlib.pyplot as plt

__prepare_matplotlib()

# This script creates pickle files with histograms and a csv of spatial analysis of vectors
def analyse_space(vector):
    histograms = {}
    output_df = pd.DataFrame()
    atom_keys = list(vector.keys())
    with alive_bar(len(atom_keys)) as bar:
        print('Calculating vector explored volumes')
        for atom_key in atom_keys:
            # Histogramming for plots
            x = [i[0] for i in vector[atom_key]]
            y = [i[1] for i in vector[atom_key]]
            z = [i[2] for i in vector[atom_key]]

            # Create histograms
            #hist_atom, edges = np.histogramdd((x,y,z))
            #histograms[atom_key] = hist_atom

            # Calculating exploration volume
            arr = np.array([x,y,z])
            arr = np.transpose(arr)
            try:
                hull = ss.ConvexHull(arr)
            except ss.qhull.QhullError:
                print('QHullError, check your vectors! (are you sure you are not subtracting two identical ones?')
                exit()

            # Obtained hull volumes, simplices (points), vertices (lines) of the shape
            #print(f'Vol={hull.volume} nm^3')
            #print(f'Simplices={hull.simplices}')
            #print(f'Vertices={hull.vertices}')

            outvol = float(hull.volume)
            outvolstep = float(outvol / (len(x)))

            #which_vector = f'Vector {vectin}'
            #out_series = pd.Series({which_vector: atom_key,
            #                  f'Volume {vectin} / nm^3': outvol,
            #                  f'Vol {vectin}/step': outvolstep})

            out_series = pd.Series({'atom' : atom_key,
                                    'vol': outvol,
                                    'volstep': outvolstep
                                    })

            output_df = output_df.append(out_series, True)
            bar()

    return(output_df)

def write_to_pdb_beta(pdb, delta):
    # This assigns the value of deltaVolumeExplored (obtained by |VE(2)-VE(1)| where VE is volume explored)
    # to a pdb file, where it's put into the last column (usually occupied by temperature factors) to be
    # visualized at your leisure (PyMol, vmd, Chimera...)

    try:
        with open(pdb,'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"{pdb} doesn't exist, will exit")
        exit()

    atom_ind = 0

    new_lines = []

    print('Writing .pdb')
    for line in lines:

        if 'ENDMDL' in line:
            new_line = line
            new_line = new_line.replace('\n', '')
            new_lines.append(new_line + '\n')
            #atom_ind = 0
            break # exits the loop after writing the first frame
        elif 'ATOM' in line:
            factor = delta.iloc[atom_ind]
            factor = "{:10.4f}".format(factor)
            # This is an ugly hack, replace later
            new_line = line.replace(' 0.00 ', factor)
            atom_ind += 1
        else:
            new_line = line

        new_line = new_line.replace('\n','')
        new_lines.append(new_line+'\n')

    new_pdb = ''.join(new_lines)

    return(new_pdb)

def main(argv=sys.argv[1:]):

    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=str, help='Input vector.json for trajectory 1', required=True)
    parser.add_argument("--s", type=str, help='OPTIONAL Input vector.json for trajectory 2', default=False)
    parser.add_argument("--o", type=str, help='Output directory, defaults to names of trajectories used separated by _', required=False, default=False)
    #parser.add_argument("--ohf", type=str, help='Output histograms for the --f file', default='histogram1')
    #parser.add_argument("--ohs", type=str, help='OPTIONAL Output histograms for the --s file', default='histogram2')
    parser.add_argument("--plot", type=bool, help='Plot spatial stuff, defaults to False', default=False)
    parser.add_argument("--pdbs", type=str, help='OPTIONAL Input structure in .pdb format of the second trajectory', default=False)
    parser.add_argument("--resi", type=str, help='OPTIONAL atoms-to-residues assignment file.json. Will turn on residue mode', required=False, default=False)
    global args
    args = parser.parse_args(argv)

    # Set up default output directory
    if not args.o:
        traj1_name = args.f.replace('.json', '')

        if args.s:
            traj2_name = args.s.replace('.json', '')
            args.o = f'{traj1_name}_{traj2_name}'
        else:
            args.o = traj1_name

    if not os.path.exists(args.o):
        os.makedirs(args.o)
    ###

    try:
        with open(args.f) as file:
            global vectors1
            vectors1 = json.load(file)

    except FileNotFoundError:
        print(f'{args.f} not found, will quit.')
        exit()

    if args.s:
        try:
            with open(args.s) as file:
                global vectors2
                vectors2 = json.load(file)
        except:
            print(f'{args.s} not found, will proceed with only one vector file analysis, obtained from --f1')
            args.f2 = False


    traj1_name = args.f.replace('.json','')
    output_df1 = analyse_space(vectors1)
    output_df1 = output_df1.rename(columns={
        'atom': traj1_name,
        'vol': f'V({traj1_name})',
        'volstep': f'V({traj1_name})/step'
    })

    if args.s:

        traj2_name = args.s.replace('.json', '')
        output_df2 = analyse_space(vectors2)
        output_df2 = output_df2.rename(columns={
            'atom': traj2_name,
            'vol': f'V({traj2_name})',
            'volstep': f'V({traj2_name})/step'
        })

        output_df = output_df1.join(output_df2)

        tot_explored_volume_1 = output_df[f'V({traj1_name})'].sum()
        tot_explored_volume_2 = output_df[f'V({traj2_name})'].sum()

        delta = output_df[f'V({traj2_name})'] - output_df[f'V({traj1_name})']
        delta = pd.DataFrame(delta, columns=['V_delta'])

        output_df = output_df.join(delta)
        output_df.loc[output_df.index[0], f'SUMV({traj1_name})'] = tot_explored_volume_1
        output_df.loc[output_df.index[0], f'SUMV({traj1_name})/step'] = tot_explored_volume_1/len(list(vectors1.keys()))
        output_df.loc[output_df.index[0], f'SUMV({traj2_name})'] = tot_explored_volume_2
        output_df.loc[output_df.index[0], f'SUMV({traj2_name})/step'] = tot_explored_volume_2/len(list(vectors2.keys()))

    else:
        output_df = output_df1
        tot_explored_volume_1 = output_df[f'V({traj1_name})'].sum()
        output_df.loc[output_df.index[0], f'SUMV({traj1_name})'] = tot_explored_volume_1
        output_df.loc[output_df.index[0], f'SUMV({traj1_name})/step'] = tot_explored_volume_1/len(list(vectors1.keys()))

    output_df.to_csv(f'{args.o}/diff_atom.csv')

    output_residue_df = pd.DataFrame()

    if args.resi:
        output_residue_df = pd.DataFrame()
        try:
            with open(args.resi) as file:
                resi_assignment = json.load(file)
        except FileNotFoundError:
            print(f'{args.resi} not found, will quit.')
            exit()

        residue_keys = list(resi_assignment.keys())
        output_df = output_df.set_index(traj1_name)

        for i, key in enumerate(residue_keys):
            atoms_in_residue = resi_assignment[key]
            min_atom=atoms_in_residue[0]
            max_atom=atoms_in_residue[-1]
            sum_vol1=output_df.loc[f'atom {min_atom}':f'atom {max_atom}', f'V({traj1_name})'].sum()
            sum_vol1_step=output_df.loc[f'atom {min_atom}':f'atom {max_atom}', f'V({traj1_name})/step'].sum()

            if args.s:
                sum_vol2 = output_df.loc[f'atom {min_atom}':f'atom {max_atom}', f'V({traj2_name})'].sum()
                sum_vol2_step = output_df.loc[f'atom {min_atom}':f'atom {max_atom}', f'V({traj2_name})/step'].sum()
                residue_df = pd.DataFrame({
                    f'residue' : key,
                    f'V({traj1_name})' : sum_vol1,
                    f'V({traj1_name})/step' : sum_vol1_step,
                    f'V({traj2_name})' : sum_vol2,
                    f'V({traj2_name})/step' : sum_vol2_step
                }, index=[key])
            else:
                residue_df = pd.DataFrame({
                    f'residue': key,
                    f'V({traj1_name})': sum_vol1,
                    f'V({traj1_name})/step': sum_vol1_step,
                }, index=[key])

            output_residue_df = output_residue_df.append(residue_df)

        if args.s:
            delta_resi = output_residue_df[f'V({traj2_name})'] - output_residue_df[f'V({traj1_name})']
            delta_resi = pd.DataFrame(delta_resi, columns=['V_delta'])
            output_residue_df = output_residue_df.join(delta_resi)

        output_residue_df.to_csv(f'{args.o}/diff_resi.csv')

    if args.pdbs:
        if not args.s:
            print('(!!) Will print B-factors for a single trajectory only in nm^3')
            delta = output_df[f'V({traj1_name})']
            diff_pdb = 'REMARK Volume explored by atom in cartesian space writen in Beta-factors in nm^3 \n'
        else:
            print('(!!) Will print B-factors deltas between the two trajectories (s - f) in nm^3')
            delta = output_df[f'V({traj2_name})'] - output_df[f'V({traj1_name})']
            diff_pdb = 'REMARK Difference of volumes explored by atoms in two trajectories in cartesian space writen in Beta-factors in nm^3 \n'


        diff_pdb += write_to_pdb_beta(pdb=args.pdbs, delta=delta)
        with open(f'{args.o}/{args.o}_diff.pdb', 'w') as file:
            file.write(diff_pdb)

    if args.plot:

        # Plot atom-wise diff
        fig, ax = plt.subplots(figsize=(15,10))
        vol_atoms_1 = output_df[f'V({traj1_name})']
        name_atoms_1 = list(vectors1.keys())
        name_atoms_1 = [re.findall(r'\d+', i)[0] for i in name_atoms_1]
        name_atoms_1 = list(map(int, name_atoms_1))
        ax.plot(name_atoms_1, vol_atoms_1, color='blue', label=args.f)
        if args.s:
            vol_atoms_2 = output_df[f'V({traj2_name})']
            name_atoms_2 = list(vectors2.keys())
            name_atoms_2 = [re.findall(r'\d+', i)[0] for i in name_atoms_2]
            name_atoms_2 = list(map(int, name_atoms_2))
            ax.plot(name_atoms_2, vol_atoms_2, color='orange', label=args.s)
        ax.set_xlabel('Atom number')
        ax.set_ylabel('Explored volume / nm^3')
        ax.legend()
        plt.savefig(f'{args.o}/diff_plot_atom.png', dpi=300)

        if args.resi:
            fig_resi, ax_resi = plt.subplots(figsize=(15,10))
            vol_atoms_1 = output_residue_df[f'V({traj1_name})']
            name_atoms_1 = output_residue_df.index
            plt.xticks(rotation=90)
            ax_resi.plot(name_atoms_1, vol_atoms_1, color='blue', label=args.f)
            if args.s:
                vol_atoms_2 = output_residue_df[f'V({traj2_name})']
                name_atoms_2 = output_residue_df.index
                plt.xticks(rotation=90)
                ax_resi.plot(name_atoms_2, vol_atoms_2, color='orange', label=args.s)
            ax_resi.set_xlabel('Residue number')
            ax_resi.set_ylabel('(sum by atoms) Explored volume / nm^3')
            ax_resi.legend()
            plt.savefig(f'{args.o}/diff_plot_resi.png', dpi=300)
        plt.show()

if __name__ == "__main__":
    main()