import argparse
import numpy as np
import pandas as pd
import pyarrow
import re
import sys

"""
python cartesian.py [--f file containing coordinates of atoms <cart.xvg>] [--s file containing coordinates of molecular COMs <mol.xvg>]

Cartesian requires an output file from gmx traj -f TRAJ.xtc -s TOP.tpr -ox CARTESIAN.xvg (vectors of atoms).
Optionally, it also requires an output from the same command, but with -com flag specified (vectors of COMs, beware - always 
compare comparable, that is COMs of the same part of the molecule for both systems!). Alternative to this
is just suplying a single coord.xvg for each of the trajectories, which were previously aligned.

Input can also be read from .pdb file (again fitted for both rotation and translation). This usage doesn't allow for using --s 
parameter for COM fitting. Frames of the trajectory are separated by ENDMDL keyword in pdb. Coordinates from .pdb are recalculated from
Angstroms to nm.

Please make sure that your trajectory is free of PBC and is fitted both for rotational as well as translational
movement of you molecule of interest. The molecule should also be (usually) exploring the equilibrium
part of the trajectory.

Consider omitting hydrogens out of the analysis, these are either constrained or too flexible and usually don't carry meaningful
information.
"""

def __prepare_matplotlib():
    import matplotlib.pyplot as plt
    # Prepare matplotlib parameters
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def parse_from_xvg(path, com=False):
    """
    parse_cartesian() parses a dataset of this type:
        # comments for .xvg reader
        @ notes about the plot for .xvg reader (names, legends etc.)
        25000	4.282	6.226	3.555	4.271	6.317	3.685	4.413	6.326	3.587	4.342	6.113	3.765	4.55 ...
        25010	...
        ...     ...

    The datalines (no @ or #) read as: at time 25000, the (X Y Z) of the first atom was (4.282 6.226 3.555), of the second
    atom (4.271 6.317 3.685), of the third atom (4.413 6.326 3.587) ...

    Prepares the input for assign_to_dict() which runs in parallel.
    1) Read the file, split comment, note and data lines
    2) Get atom names from the file, which are assigned as keys into x/y/z_cart dictionaries as x_cart={'atom 2 X' : [POS1, POS2 ...]}
    3) Get atom indices as such x_ind=[2,7,8,9,11], because the output atoms aren't necessarily a linear progressiin
    4) Compare lens of different pre-prepared dicts and lists, if they're not all of same length, exit
    5) Call assign_to_dict with a multiprocessing.Pool()
    6) Return assigned dictionaries with coords x_cart, y_cart, z_cart and times

    :return: x_cart, y_cart, z_cart, times
    """

    try:
        with open(path) as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f'{path} does not exist')
        exit()

    # Read only lines not containing '#'
    comments = tuple([line for line in lines if '@' in line])
    data = tuple([line for line in lines if '@' not in line and '#' not in line])

    # Get atom numbers from the file using regex, create 3 dictionaries containing X Y Z values for atoms
    pattern = '[a-z]+\s+[0-9]+\s+[XYZ]*'

    """
    This doesn't work, generated dictionaries can't have one of the lists appended without modifying all the other lists,
    this is normal in Python    
    x_cart = dict.fromkeys([re.search(pattern, line)[0] for line in comments if 'atom' in line and 'X' in line], [])
    y_cart = dict.fromkeys([re.search(pattern, line)[0] for line in comments if 'atom' in line and 'Y' in line], [])
    z_cart = dict.fromkeys([re.search(pattern, line)[0] for line in comments if 'atom' in line and 'Z' in line], [])
    """

    # Get keys
    # Creating a tuple from list comprehension is faster than other possibilities. https://stackoverflow.com/questions/16940293/why-is-there-no-tuple-comprehension-in-python
    if com:
        keys = ['com_x','com_y','com_z']
    else:
        x_keys = tuple([re.search(pattern, line)[0] for line in comments if 'atom' in line and 'X' in line])
        y_keys = tuple([i.replace('X','Y') for i in x_keys])
        z_keys = tuple([i.replace('X','Z') for i in x_keys])
        keys = [i for i in zip(x_keys, y_keys, z_keys)]
        keys = [x for xs in keys for x in xs]
    #####

    # Split lines and transform into pandas
    vectors_df = pd.DataFrame(data)
    vectors_df = vectors_df[0].str.split('\t', expand=True)
    vectors_df.iloc[:,-1] = vectors_df.iloc[:,-1].str.replace('\n','')
    times_df = vectors_df[0] # => simulation times
    vectors_df = vectors_df.drop(columns=vectors_df.columns[0], axis=1)
    vectors_df.columns = keys
    #####

    return(vectors_df.astype(np.float32))

def parse_from_pdb(path):

    #path='/run/timeshift/backup/IOCB/md/FDs/MSM/fitted_trajectories/VECTORS/pdz/run_1/short.pdb'
    try:
        from Bio.PDB import PDBParser
        import warnings
        with warnings.catch_warnings(): # Suppress these pesky BioPython warnings
            warnings.simplefilter("ignore")
            parser = PDBParser()
            structure = parser.get_structure("struct", path)
    except FileNotFoundError:
        print(f'{path} does not exist')
        exit()
    except ModuleNotFoundError:
        print('Missing optional dependency biopython (pip3 install biopython), cannot read .pdb, will quit.')
        exit()
    coordinates_df = pd.DataFrame()

    for i, model in enumerate(structure):

        single_row = {}
        atoms = model.get_atoms()
        for j, atom in enumerate(atoms):
            coords = list(atom.get_coord())
            single_row[f'Atom {j} X'] = coords[0] / 10 # From A to nm
            single_row[f'Atom {j} Y'] = coords[1] / 10
            single_row[f'Atom {j} Z'] = coords[2] / 10

        if i == 0:
            coordinates_df = pd.DataFrame(single_row, index=[i])
        else:
            single_row_series = pd.DataFrame(single_row, index=[i])
            coordinates_df = pd.concat([coordinates_df, single_row_series])

    print(coordinates_df)
    exit()

    # Separate .pdb into frames by ENDMDL



def recalculate_vectors_com(coordinates_df, com_df):

    print(coordinates_df)
    print(com_df)
    num_cols = len(coordinates_df.columns)
    num_atoms = int(num_cols/3)

    for column in coordinates_df.columns:
        if 'X' in column:
            i=0
        elif 'Y' in column:
            i=1
        else:
            i=2
        coordinates_df[column] = coordinates_df[column] - com_df.iloc[:, i]

    return(coordinates_df)

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=str, help='Input traj (.xvg from GROMACS gmx traj or .pdb of trajectory) (aligned for rotation+translation or for rotation in conjuction with --s COM file)', required=True)
    parser.add_argument("--s", type=str, help='OPTIONAL Input traj -com -ox filepath', required=False)
    parser.add_argument("--o", type=str, help='Output.json (json)', required=False, default='cartesian_outfile')
    parser.add_argument("--resi", type=str, help='OPTIONAL .pdb file for residue-to-atom assignment', required=False)
    global args
    args = parser.parse_args(argv)
    if not args.o:
        args.o = f'{args.f.strip(".xvg").strip(".pdb")}'

    if '.xvg' in args.f:
        traj_format = 'xvg'
    elif '.pdb' in args.f:
        traj_format = 'pdb'
        print('Usage of .pdb format does not allow for COM fitting, will ignore --s argument if present')
        args.s = False
    else:
        print('Unknown format of the --f input file (accepted: .xvg, .pdb), will quit.')
        exit()

    #x_cart, y_cart, z_cart, times = parse_cartesian(args.f)


    if args.s: # If we're using COM fitting
        print('######')
        print('(!!) Recalculating vectors with COM positions as (0,0,0)')
        print('######')
        com_df = parse_from_xvg(args.s, com=True)
        coordinates_df = parse_from_xvg(args.f)
        #vectors = get_vectors(x_cart, y_cart, z_cart, x_com, y_com, z_com) REMOVE
        vectors = recalculate_vectors_com(coordinates_df, com_df)
    elif traj_format == 'pdb':
        print('######')
        print('(!!) Using aligned trajectories instead of COM fitting')
        print('(!!) Reading from .pdb')
        print('######')
        vectors = parse_from_pdb(args.f)
    else:
        print('######')
        print('(!!) Using aligned trajectories instead of COM fitting')
        print('(!!) Reading from GROMACS .xvg')
        print('######')
        vectors = parse_from_xvg(args.f)

    with open(args.o,'wb') as fp:
        print(f'Outputting vectors to {args.o}')
        vectors.to_parquet(fp, compression='snappy')

    if args.resi:
        try:
            with open(args.resi) as file:
                # The whole trajectory
                lines = file.readlines()
                # Get the first frame
                frame_lines = []
                for line in lines:
                    if 'TER' in line or 'ENDMDL' in line:
                        break
                    elif 'ATOM' in line:
                        frame_lines.append(line.split()[1:5])
                print('######')
                print('Assigning atoms to residues')
                print('######')
        except FileNotFoundError:
            print(f'{args.resi} does not exist, no assignment will be done')
            exit()

        # Assign atoms to residues
        resis_df = pd.DataFrame(frame_lines, columns=['atom_id','atom_type','resi_type','resi_id'])
        with open(f'{args.o}_resis', 'wb') as fp:
            print(f'Outputting residue assignments to {args.o}_resis')
            resis_df.to_parquet(fp, compression='snappy')

if __name__ == '__main__':
    main()




