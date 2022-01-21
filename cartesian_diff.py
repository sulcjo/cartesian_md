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

            hull = ss.ConvexHull(arr)
            #print(f'Vol={hull.volume} nm^3')
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

def color_structure(space_df):
    # This assigns the value of deltaVolumeExplored (obtained by |VE(2)-VE(1)| where VE is volume explored)
    # to a pdb file, where it's put into the last column (usually occupied by temperature factors) to be
    # visualized at your leisure (PyMol, vmd, Chimera...)
    if not args.s:
        print('Comparing structures requires two of them, please provide one more structure using the --s flag, will quit.')
        exit()

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

    delta = space_df['Volume 2 / nm^3'] - space_df['Volume 1 / nm^3']

    new_pdb = ''
    atom_ind = 0


    with alive_bar(len(lines)) as bar:
        print('Writing .pdb')
        for line in lines:
            """
            if 'ENDMDL' in line:
                new_line = line
                atom_ind = 0

            elif 'ATOM' in line:
                # All the columns until temp-factor
                split = line.split()
                new_line = split[0] + '      ' + split[1] + '  ' + split[2] + '   ' + split[3] + '    ' + split[4] + '      ' + split[5] + '  ' + split[6] + '  ' + split[7] + '  ' + split[8] + '  '
                # Temp factor / delta of explored volume
                factor = str(abs(round(delta.iloc[atom_ind]*100,2)))

                new_line += factor
                # Atom type (C, N, O...)
                new_line += '           ' + split[9]
                atom_ind += 1
            else:
                new_line = line
            """
            if 'ENDMDL' in line:
                new_line = line
                atom_ind = 0
            elif 'ATOM' in line:
                factor = str(abs(round(delta.iloc[atom_ind] * 100, 2)))
                # This is an ugly hack, replace later
                new_line = line.replace(' 0.00 ', factor)
                atom_ind += 1
            else:
                new_line = line

            new_line = new_line.replace('\n','')
            new_pdb += new_line
            new_pdb += '\n'

            bar()


    with open(f'{args.pdbs}_colored.pdb','w') as file:
        file.write(new_pdb)







def main(argv=sys.argv[1:]):

    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=str, help='Input vector.json for trajectory 1', required=True)
    parser.add_argument("--s", type=str, help='OPTIONAL Input vector.json for trajectory 2', default=False)
    parser.add_argument("--o", type=str, help='Analysis output file', default='analysis.out')
    parser.add_argument("--ohf", type=str, help='Output histograms for the --f file', default='histogram1')
    parser.add_argument("--ohs", type=str, help='OPTIONAL Output histograms for the --s file', default='histogram2')
    parser.add_argument("--plot", type=bool, help='Plot spatial stuff, defaults to False', default=False)
    parser.add_argument("--pdbs", type=str, help='Input structure in .pdb format of the second trajectory', default=False)
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
    else:
        output_df = output_df1

    output_df.to_csv(args.o)
    with open(args.ohf, 'wb') as file:
        pickle.dump(histograms1, file)
    if args.s:
        with open(args.ohs, 'wb') as file:
            pickle.dump(histograms2, file)

    if args.pdbs:
        #pool = mp.Pool(mp.cpu_count())
        #pool.map(color_structure, output_df)
        #pool.close()
        #pool.join()
        color_structure(output_df)

    if args.plot:
        #print(output_df['Volume 1 / nm^3'])
        vol_atoms_1 = output_df['Volume 1 / nm^3']
        name_atoms_1 = list(vectors1.keys())
        name_atoms_2 = list(vectors2.keys())
        name_atoms_1 = [re.findall(r'\d+', i)[0] for i in name_atoms_1]
        name_atoms_2 = [re.findall(r'\d+', i)[0] for i in name_atoms_2]
        name_atoms_1 = list(map(int, name_atoms_1))
        name_atoms_2 = list(map(int, name_atoms_2))
        vol_atoms_2 = output_df['Volume 2 / nm^3']

        fig, ax = plt.subplots()
        ax.plot(name_atoms_1, vol_atoms_1, color='blue', label=args.f)
        ax.plot(name_atoms_2, vol_atoms_2, color='orange', label=args.s)
        ax.set_xlabel('Atom number')
        ax.set_ylabel('Explored volume / nm^3')
        ax.legend()
        plt.show()

if __name__ == "__main__":
    main()


# Analyses vector in cartesian space, basis for one-vector analysis


# Put to main after testing (remove global from vectors)

