import sys, argparse
import json
import numpy as np
import pandas as pd
from alive_progress import alive_bar
import scipy.spatial as ss
import pickle
import matplotlib.pyplot as plt
import re
import multiprocessing as mp

# This script creates pickle files with histograms and a csv of spatial analysis of vectors
def analyse_space(vector, vectin=1):
    histograms = {}
    output_df = pd.DataFrame()
    atom_keys = list(vector.keys())
    with alive_bar(len(atom_keys)) as bar:
        print('Histograming vectors and calculating volumes')
        for atom_key in atom_keys:
            # Histogramming for plots
            x = [i[0] for i in vector[atom_key]]
            y = [i[1] for i in vector[atom_key]]
            z = [i[2] for i in vector[atom_key]]




            hist_atom, edges = np.histogramdd((x,y,z))
            histograms[atom_key] = hist_atom

            # Calculating exploration volume
            arr = np.array([x,y,z])
            arr = np.transpose(arr)


            try:
                hull = ss.ConvexHull(arr)
            except ss.qhull.QhullError:
                print('QHullError, check your vectors! (are you sure you are not subtracting two identical ones?')
                exit()

            print(f'Vol={hull.volume} nm^3')
            #print(f'Simplices={hull.simplices}')
            #print(f'Vertices={hull.vertices}')

            outvol = float(hull.volume)
            outvolstep = float(outvol / (len(x)))

            which_vector = f'Vector {vectin}'
            out_series = pd.Series({which_vector: atom_key,
                              f'Volume {vectin} / nm^3': outvol,
                              f'Vol {vectin}/step': outvolstep})

            output_df = output_df.append(out_series, True)

            bar()


    return(output_df, histograms)

def color_structure(space_df, delta=False):
    # This assigns the value of deltaVolumeExplored (obtained by |VE(2)-VE(1)| where VE is volume explored)
    # to a pdb file, where it's put into the last column (usually occupied by temperature factors) to be
    # visualized at your leisure (PyMol, vmd, Chimera...)
    if not args.s:
        print('(!!) Will print B-factors for a single trajectory only in nm^3')
        delta = space_df['Volume 1 / nm^3']
        new_pdb = 'REMARK Volume explored by atom in cartesian space writen in Beta-factors in nm^3'
    else:
        print('(!!) Will print B-factors deltas between the two trajectories (s - f) in nm^3')
        delta = space_df['Volume 2 / nm^3'] - space_df['Volume 1 / nm^3']
        new_pdb = 'REMARK Difference of volumes explored by atoms in two trajectories in cartesian space writen in Beta-factors in nm^3'

    def listToString(s):

        # initialize an empty string
        str1 = ""

        # traverse in the string
        for ele in s:
            str1 += f'{ele} '

            # return string
        return str1

    try:
        with open(args.pdbs,'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"{args.pdbs} doesn't exist, will exit")
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
            break

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

    with open(f'{args.pdbs}_diff.pdb','w') as file:
        file.write(new_pdb)

def main(argv=sys.argv[1:]):

    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=str, help='Input vector.json for trajectory 1', required=True)
    parser.add_argument("--s", type=str, help='OPTIONAL Input vector.json for trajectory 2', default=False)
    parser.add_argument("--o", type=str, help='Analysis output file', default='cartesian_diff.out')
    parser.add_argument("--ohf", type=str, help='Output histograms for the --f file', default='histogram1')
    parser.add_argument("--ohs", type=str, help='OPTIONAL Output histograms for the --s file', default='histogram2')
    parser.add_argument("--plot", type=bool, help='Plot spatial stuff, defaults to False', default=False)
    parser.add_argument("--pdbs", type=str, help='OPTIONAL Input structure in .pdb format of the second trajectory', default=False)
    parser.add_argument("--resi", type=str, help='OPTIONAL atoms-to-residues assignment file.json. Will turn on residue mode', required=False)
    global args
    args = parser.parse_args(argv)

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

    output_df1, histograms1 = analyse_space(vectors1, vectin=1)

    if args.s:

        output_df2, histograms2 = analyse_space(vectors2, vectin=2)

        output_df = output_df1.join(output_df2)
        tot_explored_volume_1 = output_df['Volume 1 / nm^3'].sum()
        tot_explored_volume_2 = output_df['Volume 2 / nm^3'].sum()

        delta = output_df['Volume 2 / nm^3'] - output_df['Volume 1 / nm^3']
        delta = pd.DataFrame(delta, columns=['Vol2-Vol1'])
        output_df = output_df.join(delta)
        output_df_rename = output_df.rename(columns={'Vector 1': f'{args.f}', 'Volume 1 / nm^3': f'Volume {args.f} / nm^3',
                                  'Vol 1/step': f'{args.f} volume/trajectory_step', 'Vector 2': f'{args.s}', 'Volume 2 / nm^3': f'Volume {args.s} / nm^3',
                                  'Vol 2/step': f'{args.s} volume/trajectory_step', 'Vol2-Vol1':f'V({args.f})-V({args.s})'})
        output_df_rename.loc[output_df_rename.index[0], f'{args.f} total volume / nm^3'] = tot_explored_volume_1
        output_df_rename.loc[output_df_rename.index[0], f'{args.f} total volume/atom'] = tot_explored_volume_1/len(list(vectors1.keys()))
        output_df_rename.loc[output_df_rename.index[0], f'{args.s} total volume / nm^3'] = tot_explored_volume_2
        output_df_rename.loc[output_df_rename.index[0], f'{args.s} total volume/atom'] = tot_explored_volume_2/len(list(vectors2.keys()))

    else:
        output_df = output_df1
        tot_explored_volume_1 = output_df['Volume 1 / nm^3'].sum()
        output_df_rename = output_df.rename(columns={'Vector 1': f'{args.f}', 'Volume 1 / nm^3': f'Volume {args.f} / nm^3', 'Vol 1/step':f'{args.f} volume/trajectory_step'})
        output_df_rename.loc[output_df_rename.index[0], f'{args.f} total volume / nm^3'] = tot_explored_volume_1
        output_df_rename.loc[output_df_rename.index[0], f'{args.f} total volume/atom'] = tot_explored_volume_1/len(list(vectors1.keys()))

    output_df_rename.to_csv(args.o)

    output_residue_df = pd.DataFrame()

    if args.resi:
        try:
            with open(args.resi) as file:
                resi_assignment = json.load(file)
        except FileNotFoundError:
            print(f'{args.resi} not found, will quit.')
            exit()

        residue_keys = list(resi_assignment.keys())
        output_df = output_df.set_index('Vector 1')



        for i, key in enumerate(residue_keys):
            atoms_in_residue = resi_assignment[key]
            min_atom=atoms_in_residue[0]
            max_atom=atoms_in_residue[-1]
            sum_vol1=output_df.loc[f'atom {min_atom}':f'atom {max_atom}', 'Volume 1 / nm^3'].sum()
            sum_vol1_step=output_df.loc[f'atom {min_atom}':f'atom {max_atom}', 'Vol 1/step'].sum()

            if args.s:
                sum_vol2 = output_df.loc[f'atom {min_atom}':f'atom {max_atom}', 'Volume 2 / nm^3'].sum()
                sum_vol2_step = output_df.loc[f'atom {min_atom}':f'atom {max_atom}', 'Vol 2/step'].sum()
                residue_df = pd.DataFrame({
                    f'residue' : key,
                    f'{args.f} total volume / nm^3' : sum_vol1,
                    f'{args.f} total volume/atom' : sum_vol1_step,
                    f'{args.s} total volume / nm^3' : sum_vol2,
                    f'{args.s} total volume/atom' : sum_vol2_step
                }, index=[key])
            else:
                residue_df = pd.DataFrame({
                    f'residue': key,
                    f'{args.f} total volume / nm^3': sum_vol1,
                    f'{args.f} total volume/atom': sum_vol1_step,
                }, index=[key])

            output_residue_df = output_residue_df.append(residue_df)

        if args.s:
            delta_resi = output_residue_df[f'{args.s} total volume / nm^3'] - output_residue_df[f'{args.f} total volume / nm^3']
            delta_resi = pd.DataFrame(delta_resi, columns=[f'V({args.s})-V({args.f})'])
            output_residue_df = output_residue_df.join(delta_resi)

        output_residue_df.to_csv(f'{args.o}_resi.csv')





    with open(args.ohf, 'wb') as file:
        pickle.dump(histograms1, file)
    if args.s:
        with open(args.ohs, 'wb') as file:
            pickle.dump(histograms2, file)
    if args.pdbs:
        color_structure(output_df)


    if args.plot:

        fig, ax = plt.subplots()
        if args.resi:
            vol_atoms_1 = output_residue_df[f'{args.f} total volume / nm^3']
            name_atoms_1 = output_residue_df.index
            plt.xticks(rotation=90)
        else:
            vol_atoms_1 = output_df['Volume 1 / nm^3']
            name_atoms_1 = list(vectors1.keys())
            name_atoms_1 = [re.findall(r'\d+', i)[0] for i in name_atoms_1]
            name_atoms_1 = list(map(int, name_atoms_1))
        ax.plot(name_atoms_1, vol_atoms_1, color='blue', label=args.f)

        if args.s:
            if args.resi:
                vol_atoms_2 = output_residue_df[f'{args.s} total volume / nm^3']
                name_atoms_2 = output_residue_df.index
                plt.xticks(rotation=90)
            else:
                vol_atoms_2 = output_df['Volume 2 / nm^3']
                name_atoms_2 = list(vectors2.keys())
                name_atoms_2 = [re.findall(r'\d+', i)[0] for i in name_atoms_2]
                name_atoms_2 = list(map(int, name_atoms_2))

            ax.plot(name_atoms_2, vol_atoms_2, color='orange', label=args.s)

        ax.set_xlabel('Atom number')
        ax.set_ylabel('Explored volume / nm^3')
        ax.legend()
        plt.show()

if __name__ == "__main__":
    main()


# Analyses vector in cartesian space, basis for one-vector analysis


# Put to main after testing (remove global from vectors)

