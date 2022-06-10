"""
Combined cartesian_diff and cartesian_grid modules into cartesian_ana
"""

import sys, argparse, os
import json
import numpy as np
import pandas as pd
from alive_progress import alive_bar
import scipy.spatial as ss
from cartesian import __prepare_matplotlib
import re
import matplotlib.pyplot as plt
import math
import matplotlib
matplotlib.use('qtagg')
import seaborn as sb
import matplotlib.ticker as ticker
try:
    import multiprocessing as mp
    multiprocess = True
except ImportError:
    multiprocess = False

from alive_progress import alive_bar
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from cartesian import __prepare_matplotlib
from cartesian_diff import write_to_pdb_beta
from matplotlib.widgets import TextBox, Slider, Button
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import math

__prepare_matplotlib()

def analyse_space(vector):
    """
    :param vector: Cartesian.py produced MD trajectory atomic vectors
    :return: Pandas DataFrame of volumes explored by atoms, approximated with Convex Hull algorithm in 3D-space
    """
    output_df = pd.DataFrame(columns=['atom', 'vol', 'volstep'])
    num_atoms = math.ceil(len(vector.columns) / 3)


    out_dict_list = []

    with alive_bar(num_atoms) as bar:
        print('Calculating vector explored volumes')

        # Atoms in vector dataframe have data in three columns (x,y,z,)

        for i in range(num_atoms):

            lower = i*3 # Three columns per atom, we're counting triples

            atom_coordinates = np.array(vector.iloc[:, lower:lower+3]).astype(np.float32)

            try:
                hull = ss.ConvexHull(atom_coordinates)
            except ss.qhull.QhullError:
                print('QHullError, check your vectors! (are you sure you are not subtracting two identical ones?')
                exit()


            # Obtained hull volumes, simplices (points), vertices (lines) of the shape
            #print(f'Vol={hull.volume} nm^3')
            #print(f'Simplices={hull.simplices}')
            #print(f'Vertices={hull.vertices}')

            outvol = np.float32(hull.volume)
            outvolstep = np.float32(outvol / (len(vector.iloc[:, lower]))) # steps = number of rows for the current atom X coordinate

            out_dict_list.append({'atom':int(i), 'vol':outvol, 'volstep':outvolstep})

            bar()

    output_df = pd.DataFrame(out_dict_list)

    return(output_df)

def grid(vectors):


    divided = (vectors / args.grid).astype(np.int16)  # Divide all cells by args.grid #

    import time
    start_time = time.time()

    """
    Atomic coordinates are assigned into 3D histogram bins with binning parameter set by args.bin.
    In the unlikely case that a certain coordinate is on the binning boundary (binning=1 ; coordinate=2.0 for example)
    the coordinate always goes into the "bin on the left". Returns two DataFrames, one containing unique grids for each atom XYZ value,
    the other one containing corresponding counts (how many times the vector
    visited these grid members during the trajectory)
    :param vectors:
    :return: Atomic coordinates as bin positions in XYZ, only unique members and their counts (two DFs)
    """
    # Contrary to the old solution, all values that are on the bin border (bin=1 ; values like 2.0, 3.0, 1.0) go "to the left" bin
    # i.e. Coordinate=2.0 / Bin=1.0 ends up in bin #2

    number_of_triples = int(len(divided.columns)/3)
    list_of_uniq_dfs = []
    list_of_count_dfs = []

    # np.unique slower than pd.unique, which uses hash-tables to lookup the unique values
    # but pandas can't work on three numerical values at once unlike numpy which can find unique values in 3,n arrays
    # join XYZ columns into strings (one string for one atom in one observation)
    # then find unique strings using pd.unique. Then split back into numerical DataFrame with str.split(',')

    for triple in range(number_of_triples):

        # Column triple iterator
        lowcol = triple*3

        # XYZ values from three columns into a single column. String of form 'x,y,z'
        triple_arr = divided.iloc[:, lowcol:lowcol+3].to_numpy().reshape(-1,3)

        if args.gridbackend == 0:
            # Variant np.unique solution
            uniq, count = np.unique(triple_arr, axis=0, return_counts=True)

            uniq = tuple(uniq.astype(np.int32)) # It works with tuples, with arrays pandas just expand everything
            count = tuple(count.astype(np.int32))

        else:
            # pd.unique solution
            triplize = lambda tr: f'{tr[0]},{tr[1]},{tr[2]}'
            #triplized = np.apply_along_axis(triplize, 1, triple_arr)
            triplized = tuple(map(triplize, triple_arr))

            uniq = pd.value_counts(triplized)

            count = uniq.values
            detriplize = lambda tr: tuple(np.array(tr.split(',')).astype(np.int32))
            uniq = tuple(map(detriplize, uniq.index.values))

        list_of_uniq_dfs.append(uniq)
        list_of_count_dfs.append(count)

    uniq_df = pd.DataFrame(list_of_uniq_dfs).T
    count_df = pd.DataFrame(list_of_count_dfs).T


    print("--- %s Uniques ided---" % (time.time() - start_time))

    # 11.8 seconds with np.unique(triple_arr)
    # 15.1 seconds with pd.unique and lambdas for modifying using maps

    return(uniq_df, count_df)

def grid_rms(grid_unq1, grid_unq2, grid_count1, grid_count2):
    #print(grid_unq1)
    #print(grid_count1)

    #grid_unq1 = pd.read_parquet('/run/timeshift/backup/IOCB/md/FDs/MSM/fitted_trajectories/VECTORS/pdz/run_1/old_method_unique_grids/grids1')
    #grid_unq2 = pd.read_parquet('/run/timeshift/backup/IOCB/md/FDs/MSM/fitted_trajectories/VECTORS/pdz/run_1/old_method_unique_grids/grids2')
    #grid_count1 = pd.read_parquet('/run/timeshift/backup/IOCB/md/FDs/MSM/fitted_trajectories/VECTORS/pdz/run_1/old_method_unique_grids/grids1_count')
    #grid_count2 = pd.read_parquet('/run/timeshift/backup/IOCB/md/FDs/MSM/fitted_trajectories/VECTORS/pdz/run_1/old_method_unique_grids/grids2_count')


    import time
    start_time = time.time()

    # How many triples
    colsno_1 = len(grid_unq1.columns)
    colsno_2 = len(grid_unq2.columns)

    # Prepare RMSmetric list (one value for one compared atom)
    rms_lst = []

    number_of_triples = min(colsno_1, colsno_2) # Handles cases where one of the datasets has less atoms

    if colsno_1 != colsno_2:
        print('WARNING grid-datasets used for RMS calculation are of different size, RMS value may be faulty. Will iterate through'
              'atoms until one of the datasets runs out.')

    # two lens, ceil and /2 because datasets can theoretically have different size
    # this way the bar really calculates progress in such a case

    print("--- %s A: initial check&prep ---" % (time.time() - start_time))


    print(f'Calculating inter-grid RMS using gridsize {args.grid} nm')
    for col1, col2 in zip(grid_unq1.columns, grid_unq2.columns): # This works even if one of the datasets has less atoms
        rmsd_sum = 0



        triple1 = pd.DataFrame(grid_unq1[col1].dropna().tolist())



        triple2 = pd.DataFrame(grid_unq2[col2].dropna().tolist())
        xs2, ys2, zs2 = tuple(triple2[0]), tuple(triple2[1]), tuple(triple2[2])


        triple1_counts = tuple(grid_count1[col1].dropna().tolist())


        triple2_counts = pd.DataFrame(grid_count2[col2].dropna().tolist())




        sum_counts1, sum_counts2 = sum(triple1_counts), sum(triple2_counts[0])

        samples = sum_counts1*sum_counts2


        # For i-th entry in first trajectory grids, for j-th entry in second trajectory grids
        i = 0

        matrix_x = []
        matrix_y = []
        matrix_z = []
        matrix_weight = []


        while i < len(triple1):

            # X-term (single i with all j), Y-term, Z-term
            x_term_i = triple1.iloc[i, 0] - xs2 # Create a series for i=1 minus all possible j
            y_term_i = triple1.iloc[i, 1] - ys2
            z_term_i = triple1.iloc[i, 2] - zs2

            # Weights for all pairs
            i_weight = np.array(triple1_counts[i] * triple2_counts, dtype=np.int32)


            matrix_x.append(x_term_i)
            matrix_y.append(y_term_i)
            matrix_z.append(z_term_i)
            matrix_weight.append(i_weight)

            i += 1


        for ix, iy, iz, iweight in zip(matrix_x, matrix_y, matrix_z, matrix_weight):
            #print(ix)
            #print(iy)
            #print(iz)


            ix = ix**2
            iy = iy**2
            iz = iz**2

            ixyz = ix+iy+iz
            ixyz_w = ixyz*iweight

            rmsd_sum += ixyz_w.sum()

        rmsd = math.sqrt((1 / samples) * rmsd_sum)

        rms_lst.append(rmsd)

    print("--- %s C: Sums calculated ---" % (time.time() - start_time))

    fig, axs = plt.subplots()
    import seaborn as sb
    plotdata = pd.DataFrame(rms_lst)



    #backend0 = pd.read_csv('/run/timeshift/backup/IOCB/md/FDs/MSM/fitted_trajectories/VECTORS/pdz/run_1/uniques_backend0.csv')
    #sb.lineplot(data=backend0.iloc[:, 1], ax=axs, color='blue')

    #df = pd.read_csv('//run/timeshift/backup/IOCB/md/FDs/MSM/fitted_trajectories/VECTORS/pdz/run_1/old_method_unique_grids/rms_grid01_skipno.csv')
    df = pd.read_csv('/run/timeshift/backup/IOCB/md/FDs/MSM/fitted_trajectories/VECTORS/pdz/run_1/old_method_unique_grids/correct_grids_old.csv')
    sb.lineplot(data=df.iloc[:, 1], ax=axs, color='red')


    plotdata = plotdata
    sb.lineplot(data=plotdata, ax=axs, color='green')
    print(rms_lst)






    plt.show()

    exit()

    return (rms_lst)

def internal_grid_rms(grid, grid_count):
    atom_keys = list(grid.keys())
    with alive_bar(len(atom_keys)) as bar:
        print('Calculating internal RMS')
        rmsd_dict = {}
        # For each atom
        for atom in atom_keys:
            rmsd_sum = 0
            atom_grid = grid[atom]
            atom_grid_count = grid_count[atom]

            atom_xs = [x[0] for x in atom_grid]
            atom_ys = [x[1] for x in atom_grid]
            atom_zs = [x[2] for x in atom_grid]
            current_atom_grid_count = [count for count in atom_grid_count]

            # print(atom_ys)
            # print(atom_zs)

            i = 0
            j = 0
            while i < len(atom_xs):
                # print(i,j)
                # print(rmsd)
                if j == len(atom_xs) - 1:
                    i += 1
                    j = 0
                    continue
                if i != j:
                    rmsd_sum += (current_atom_grid_count[i] * current_atom_grid_count[j]) * (
                                atom_xs[i] - atom_xs[j]) ** 2 + (atom_ys[i] - atom_ys[j]) ** 2 + (
                                            atom_zs[i] + atom_zs[j]) ** 2
                j += 1

            rmsd = math.sqrt((1 / 2) * (1 / len(atom_xs)) * rmsd_sum)
            rmsd_dict[atom] = rmsd
            bar()
    return (rmsd_dict)

def parse_violin(start, stop, vectors1, traj1_name, vectors2=False, traj2_name=False):
    violin_df = pd.DataFrame()
    vector1_keys = list(vectors1.keys())

    coords = []
    datapoint_id = []
    axis_id = []
    traj_id = []
    atom_id = []

    # Fix later to index by vector2_keys

    for atom in vector1_keys[start:stop]:
        coords.append([x[0] for x in vectors1[atom]])
        datapoint_id.append(np.repeat(f'{atom}_x', len(vectors1[atom])))
        axis_id.append(np.repeat('x', len(vectors1[atom])))
        traj_id.append(np.repeat(traj1_name, len(vectors1[atom])))
        atom_id.append(np.repeat(atom, len(vectors1[atom])))

        if args.s:
            coords.append([x[0] for x in vectors2[atom]])
            datapoint_id.append(np.repeat(f'{atom}_x', len(vectors2[atom])))
            axis_id.append(np.repeat('x', len(vectors2[atom])))
            traj_id.append(np.repeat(traj2_name, len(vectors2[atom])))
            atom_id.append(np.repeat(atom, len(vectors2[atom])))

        coords.append([y[1] for y in vectors1[atom]])
        datapoint_id.append(np.repeat(f'{atom}_y', len(vectors1[atom])))
        axis_id.append(np.repeat('y', len(vectors1[atom])))
        traj_id.append(np.repeat(traj1_name, len(vectors1[atom])))
        atom_id.append(np.repeat(atom, len(vectors1[atom])))
        if args.s:
            coords.append([y[1] for y in vectors2[atom]])
            datapoint_id.append(np.repeat(f'{atom}_y', len(vectors2[atom])))
            axis_id.append(np.repeat('y', len(vectors2[atom])))
            traj_id.append(np.repeat(traj2_name, len(vectors2[atom])))
            atom_id.append(np.repeat(atom, len(vectors2[atom])))

        coords.append([z[2] for z in vectors1[atom]])
        datapoint_id.append(np.repeat(f'{atom}_z', len(vectors1[atom])))
        axis_id.append(np.repeat('z', len(vectors1[atom])))
        traj_id.append(np.repeat(traj1_name, len(vectors1[atom])))
        atom_id.append(np.repeat(atom, len(vectors1[atom])))
        if args.s:
            coords.append([z[2] for z in vectors2[atom]])
            datapoint_id.append(np.repeat(f'{atom}_z', len(vectors2[atom])))
            axis_id.append(np.repeat('z', len(vectors2[atom])))
            traj_id.append(np.repeat(traj2_name, len(vectors2[atom])))
            atom_id.append(np.repeat(atom, len(vectors2[atom])))

    coords = [item for sublist in coords for item in sublist]  # flatten the coords list
    datapoint_id = [item for sublist in datapoint_id for item in sublist]
    axis_id = [item for sublist in axis_id for item in sublist]
    traj_id = [item for sublist in traj_id for item in sublist]
    atom_id = [item for sublist in atom_id for item in sublist]

    violin_df['coords'] = coords
    violin_df['datapoint_id'] = datapoint_id
    violin_df['axis_id'] = axis_id
    violin_df['traj_id'] = traj_id
    violin_df['atom_id'] = atom_id

    return (violin_df)

def plot_cubes(data, color, label, ax):
    x_data = [x[0] for x in data]
    y_data = [y[1] for y in data]
    z_data = [z[2] for z in data]

    def cuboid_data(o, size):
        X = [[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
             [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
             [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
             [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
             [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
             [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]]]
        X = np.array(X).astype(float)
        for i in range(3):
            X[:, :, i] *= size[i]
        X += np.array(o)
        return X

    def plotCubeAt(positions, colors=None, **kwargs):
        g = []
        for p, c in zip(positions, colors):
            g.append(cuboid_data(p, size=[args.grid, args.grid, args.grid]))
        return Poly3DCollection(np.concatenate(g),
                                facecolors=np.repeat(colors, 6, axis=0), **kwargs)

    ### POSITION OF THE CUBE
    positions = [(x * args.grid, y * args.grid, z * args.grid) for x, y, z in zip(x_data, y_data, z_data)]
    colors = [color for i in range(0, len(positions))]
    ###

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect([1, 1, 1])
    pc = plotCubeAt(positions, colors=colors, edgecolor="k")
    pc.set_zsort('min')  # not needed?
    ax.add_collection3d(pc)
    ax.scatter(-100, -100, -100, color=color, label=label, marker='s', s=50)

    # plt.show()

def set_ax_lims_3d(grids1, grids2, ax):
    """ Flatten the grids datasets, join them and find the lowest and highest gridpoint. Set ax limits and ticks """
    """
    Since plot_cubes() function determines the cube size using args.grid, we have to count with that (grids contain 
    "grid size" units, while the plot is in real nanometers (grid-size*args.grid)
    """

    """ Flatten the grids datasets, join them and find the lowest and highest gridpoint. Set ax limits and ticks """
    flat_list = [item for sublist in grids1 for item in sublist]
    flat_list_2 = [item for sublist in grids2 for item in sublist]
    flat_list = flat_list + flat_list_2

    grids_x = [x[0] for x in grids1] + [x[0] for x in grids2]
    grids_y = [y[1] for y in grids1] + [y[1] for y in grids2]
    grids_z = [z[2] for z in grids1] + [z[2] for z in grids2]

    lower, upper = (min(flat_list) * args.grid) - args.grid, (max(flat_list) * args.grid) + args.grid
    # print(lower, upper)

    ax.set_xlabel('Bin(x) / nm')
    ax.set_ylabel('Bin(y) / nm')
    ax.set_zlabel('Bin(z) / nm')
    # ax.autoscale_view()
    lower_int = int(math.floor(lower) - 1000 * args.grid)
    upper_int = int(math.ceil(upper) + 1000 * args.grid)

    ticks = [float(i * args.grid) for i in range(lower_int, upper_int)]
    ax.set_xticks(ticks=ticks)
    ax.set_yticks(ticks=ticks)
    ax.set_zticks(ticks=ticks)
    ax.set(xlim=(lower, upper), ylim=(lower, upper), zlim=(lower, upper))
    #ax.set_box_aspect((np.ptp(grids_x), np.ptp(grids_y), np.ptp(grids_z)))

    xmin = min(grids_x)*args.grid-args.grid
    xmax = max(grids_x)*args.grid+args.grid
    ymin = min(grids_y)*args.grid-args.grid
    ymax = max(grids_y)*args.grid+args.grid
    zmin = min(grids_z)*args.grid-args.grid
    zmax = max(grids_z)*args.grid+args.grid

    #ax.set(xlim=(min(grids_x)*args.grid, max(grids_x)*args.grid)+args.grid, ylim=(min(grids_y)*args.grid-args.grid, max(grids_y)*args.grid)+args.grid, zlim=(min(grids_z)*args.grid-args.grid, max(grids_z)*args.grid)+args.grid)
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), zlim=(zmin, zmax))

def write_to_pdb_beta(pdb, delta):
    """
    :param pdb: Reference .pdb file of the trajectory
    :param delta: Series of a parameter to be written into B-field of the .pdb file
    :return: Textlines of a modified .pdb file, with B-values being filled with the provided delta
    """
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
            #print(factor)


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

# Helper function for better boolean argument handling
def str2bool(v):
    if isinstance(v, bool):
        return(v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return(True)
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return(False)
    else:
        raise(argparse.ArgumentTypeError('Boolean value expected.'))

# Helper function
strlower = lambda v: v.lower()

def main(argv=sys.argv[1:]):

    ## General I/O
    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=str, help='Input vector.json for trajectory 1', required=True)
    parser.add_argument("--s", type=str, help='OPTIONAL Input vector.json for trajectory 2', default=False)
    parser.add_argument("--o", type=str, help='Output directory, defaults to names of trajectories used separated by _', required=False, default=False)
    parser.add_argument("--pdbs", type=str, help='OPTIONAL Input structure in .pdb format of the second trajectory', default=False)
    parser.add_argument("--resi", type=str, help='OPTIONAL atoms-to-residues assignment file.json. Will turn on residue mode', required=False, default=False)


    # cartesian grid unique arguments
    parser.add_argument('--method', type=strlower, help='Method - grid, grid_scan, violin, volume, all; defaults to all', required=False, default='grid',
                        choices=('grid','grid_scan','violin','volume','all'))

    ## Grid method specifics
    parser.add_argument('--grid', type=float, help='Grid size in nm, defaults to 0.1 nm', required=False, default=0.1)
    parser.add_argument('--pop_threshold', type=int,
                        help='Bins with populations lower than this will be disregarded for plotting, defaults to 10 (good for throwing away SSAP caused artefacts)',
                        required=False, default=10)
    parser.add_argument('--mp', type=int, help='Nthreads to use for grid calculations, defaults to 1', required=False, default=False)
    parser.add_argument('--gridbackend', type=int, help='Backend to use for unique grid assignment', required=False, default=1, choices=(0,1))

    ## Plotting
    parser.add_argument("--plot", type=str2bool, help='Plot spatial stuff, defaults to False', const=True, default=True, nargs='?')
    parser.add_argument('--plot_3d', type=str2bool, help='OPTIONAL plot results in a 3D-plot', const=True, required=False, default=True, nargs='?')
    parser.add_argument('--plot_positions', type=str2bool, help='OPTIONAL plot positional violin plots', const=True, required=False, default=True, nargs='?')
    parser.add_argument('--plot_diff', type=str2bool, help='OPTIONAL plot explored volume by atom', const=True, required=False, default=True, nargs='?')
    parser.add_argument('--plot_violin', type=str2bool, help='OPTIONAL plot violin plots of spatial positions', const=True, required=False, default=True, nargs='?')
    ##

    global args
    args = parser.parse_args(argv)

    # Set up trajectory names from file names
    expression = "[\w-]+?((?=\.)|(?=$))"
    traj1_name = re.search(expression, args.f).group(0)
    _ = lambda: re.search(expression, args.s).group(0) if args.s else False
    traj2_name = _()

    # Load vectors
    try:
        global vectors1
        vectors1 = pd.read_parquet(args.f).astype(np.float32)
    except FileNotFoundError:
        print(f'{args.f} not found, will quit.')
        exit()
    ##

    # Prepare helper file for cartesian_batch.py if comparing two different trajectories
    if args.s:
        try:
            global vectors2
            vectors2 = pd.read_parquet(args.s).astype(np.float32)

            """
            This part only works if you're using two trajectories, perhaps modify later, so it prints out an outputs.txt
            file for cartesian_batch.py even if we're only interested in obtaining distributions with a single batch of trajectories
            """

            with open('outputs.txt', 'at') as file: # This assumes that --method=all !

                diffname = 'diff_atom.csv'
                confname = f'{args.o}_grid_g{args.grid}_p{args.pop_threshold}.csv'
                grids1name = f'{traj1_name}_grid_g{args.grid}.cart'
                grids2name = f'{traj2_name}_grid_g{args.grid}.cart'
                grids1count = f'{traj1_name}_grid_g{args.grid}_count.cart'
                grids2count = f'{traj2_name}_grid_g{args.grid}_count.cart'

                file.write(f'{args.f},{args.s},{args.o},{diffname},{confname},{grids1name},{grids2name},{grids1count},{grids2count}\n')

        except:
            print(f'{args.s} not found, will proceed with only one vector file analysis, obtained from --f')
            args.s = False
    ##

    # Setup output directory name and create it if necessary
    _ = lambda: args.o if args.o else (f'{traj1_name}_{traj2_name}_out' if args.s else f'{traj1_name}_out')
    args.o = _()
    print(args.o)
    if not os.path.exists(args.o):  # Fix this, only take the name of the args.f and args.s, not the whole path!
        os.makedirs(args.o)
    ##

    if args.method == 'volume' or args.method == 'all':
        output_df1 = analyse_space(vectors1)

        output_df1 = output_df1.rename(columns={
            'atom': traj1_name,
            'vol': f'V({traj1_name})',
            'volstep': f'V({traj1_name})/step'
        })

        if args.s:

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
            output_df.loc[output_df.index[0], f'SUMV({traj1_name})/step'] = tot_explored_volume_1 / len(vectors1)
            output_df.loc[output_df.index[0], f'SUMV({traj2_name})'] = tot_explored_volume_2
            output_df.loc[output_df.index[0], f'SUMV({traj2_name})/step'] = tot_explored_volume_2 / len(vectors2)
        else:
            output_df = output_df1
            tot_explored_volume_1 = output_df[f'V({traj1_name})'].sum()
            output_df.loc[output_df.index[0], f'SUMV({traj1_name})'] = tot_explored_volume_1
            output_df.loc[output_df.index[0], f'SUMV({traj1_name})/step'] = tot_explored_volume_1 / len(vectors1)

        output_df.to_csv(f'{args.o}/diff_atom.csv')

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

        if args.plot_diff:
            # Plot atom-wise diff
            fig, ax = plt.subplots(figsize=(15, 10))
            vol_atoms_1 = output_df[f'V({traj1_name})']
            name_atoms_1 = output_df[traj1_name]
            ax.plot(name_atoms_1, vol_atoms_1, color='blue', label=traj1_name)
            # Rework plotting with seaborn instead of matplotlib
            if args.s:
                vol_atoms_2 = output_df[f'V({traj2_name})']
                name_atoms_2 = output_df[traj2_name]
                ax.plot(name_atoms_2, vol_atoms_2, color='orange', label=traj2_name)
            ax.set_xlabel('Atom number')
            ax.set_ylabel('Explored volume / nm^3')
            ax.legend()
            plt.savefig(f'{args.o}/diff_plot_atom.png', dpi=300)

            plt.show() # use only one plt.show() at the end
        global VOLDF
        VOLDF = output_df

    if args.method == 'grid' or args.method == 'all':

        atomic_grid_uniq, atomic_grid_count = grid(vectors1)
        nongridded_vectors1_num = vectors1.count(numeric_only=True, axis=0).sum() # how many vectors before
        gridded_vectors1_num = atomic_grid_uniq.count(numeric_only=True, axis=0).sum() # how many vectors after binning

        nongridded_vectors2_num = 0
        gridded_vectors2_num = 0

        if args.s:
            atomic_grid_uniq2, atomic_grid_count2 = grid(vectors2)
            nongridded_vectors2_num = vectors2.count(numeric_only=True, axis=0).sum()
            gridded_vectors2_num = atomic_grid_uniq2.count(numeric_only=True, axis=0).sum()

            rms = grid_rms(atomic_grid_uniq, atomic_grid_uniq2, atomic_grid_count, atomic_grid_count2)
            exit()

            rms_out = pd.Series(rms).T
            rms_out = pd.DataFrame(rms_out, columns=[f'{traj1_name}/{traj2_name}'])
            rms_out.to_csv(f'{args.o}/{args.o}_grid_g{args.grid}_p{args.pop_threshold}.csv')
        else:
            grid2_size = 0
            grid2_unq_size = 0

        print(f'Total amount of coordinates: {nongridded_vectors1_num + nongridded_vectors2_num}, '
              f'after binning {gridded_vectors1_num + gridded_vectors2_num}')

        # Save grids for later use with cartesian_batch
        with open(f'{args.o}/{traj1_name}_grid_g{args.grid}.json', 'w') as file:
            json.dump(atomic_grid_uniq, file, cls=NumpyEncoder)
        with open(f'{args.o}/{traj1_name}_grid_g{args.grid}_count.json', 'w') as file:
            json.dump(grid_count, file, cls=NumpyEncoder)

        if args.s:
            with open(f'{args.o}/{traj2_name}_grid_g{args.grid}.json', 'w') as file:
                json.dump(atomic_grid_2_uniq, file, cls=NumpyEncoder)
            with open(f'{args.o}/{traj2_name}_grid_g{args.grid}_count.json', 'w') as file:
                json.dump(grid_count_2, file, cls=NumpyEncoder)
        ###




        global GRIDF
        GRIDF = rms_out

        print(f'Total amount of vectors {vectors_num*2}')
        print(f'Total amount of gridpoints {grid1_size+grid2_size}')
        print(f'Total amount of unique gridpoints {grid1_unq_size+grid2_unq_size}')

        # Calculate internal RMSD of grid-set 1 and 2
        #int_rms_1 = internal_grid_rms(atomic_grid_uniq, grid_count)
        #int_rms_2 = internal_grid_rms(atomic_grid_2_uniq, grid_count_2)
        #print(int_rms_1)
        #print(int_rms_2)

        if args.pdbs:
            print('(!!) Will print B-factors of conformational deltas between the two trajectories in arbitrary units')
            delta = pd.Series(rms)

            grid_pdb = 'REMARK Differences in exploring cartesian space per atom writen in Beta-factors in arbitrary units \n'
            grid_pdb += write_to_pdb_beta(pdb=args.pdbs, delta=delta)

            with open(f'{args.o}/{args.o}_grid.pdb', 'w') as file:
                file.write(grid_pdb)

        if args.plot_3d and args.s:

            if args.resi:
                try:
                    with open(args.resi) as file:
                        resi_assignment = json.load(file)
                except FileNotFoundError:
                    print(f'{args.resi} not found, will not use.')
                    args.resi = False

                residue_keys = list(resi_assignment.keys())

            atomic_grid_keys = list(atomic_grid.keys())
            atomic_grid_keys_2 = list(atomic_grid_2.keys())
            count_dict_keys = list(grid_count.keys())
            count_dict_keys_2 = list(grid_count_2.keys())

            # Adjust bottom to make room for Buttons
            fig, ax = plt.subplots(figsize=(15, 12))
            plt.axis('off')
            ax = fig.add_subplot(projection='3d')
            plt.subplots_adjust(bottom=0.25)



            def submit(expression, minpop=0):
                global current_plot
                current_plot = expression
                # print(atomic_grid_uniq[f'atom {current_plot+1}'])

                """ Get (X,Y,Z) unique grid coordinates for TRAJ1"""
                grids1 = tuple(map(tuple, atomic_grid_uniq[atomic_grid_keys[current_plot]]))
                """ Get (X,Y,Z) unique grid coordinates for TRAJ2"""
                grids2 = tuple(map(tuple, atomic_grid_2_uniq[atomic_grid_keys_2[current_plot]]))
                """ Have to keep order, so can't use sets with intersections etc. """

                grids1_unq = [grid for grid in grids1 if grid not in grids2]  # Coords where only grids 1 reside
                grids2_unq = [grid for grid in grids2 if grid not in grids1]  # Coords where only grids 2 reside
                grids_intersect = [grid for grid in grids1 if grid in grids2]  # Intersection of the two

                ax.cla()  # clean ax

                if len(grids1_unq) > 0:
                    plot_cubes(data=grids1_unq, color='red',
                               label=f'{traj1_name} only', ax=ax)  # Plot cubes which are only present in dataset 1
                if len(grids2_unq) > 0:
                    plot_cubes(data=grids2_unq, color='blue',
                               label=f'{traj2_name} only', ax=ax)  # Plot cubes which are only present in dataset 2
                if len(grids_intersect) > 0:
                    plot_cubes(data=grids_intersect, color='purple',
                               label=f'intersection', ax=ax)  # Intersections (i.e. both 1 and 2 have a point here)

                ax.legend(loc='upper left')
                set_ax_lims_3d(grids1, grids2, ax=ax)
                ax.set_title(f'Atom {expression + 1}\ngridsize {args.grid}, pop. thresh. {args.pop_threshold}')

                # ax.autoscale_view()

                plt.draw()

            def save_all(prefix):
                global current_plot
                for plot in range(0, len(atomic_grid_keys)):
                    submit(plot)
                    plt.savefig(f'{args.o}/gridplot_atom{current_plot}.png')

            global current_plot
            current_plot = 0
            submit(current_plot)

            axbox = fig.add_axes([0.3, 0.05, 0.6, 0.075])
            save_axbox = fig.add_axes([0.8, 0.075, 0.2, 0.05])
            save_all_axbox = fig.add_axes([0.8, 0.025, 0.2, 0.05])
            # text_box = TextBox(axbox, "Atom # / Resi #", textalignment="center")
            # text_box.on_submit(submit)
            # text_box.set_val("0")  # Trigger `submit` with the initial string.

            slider = Slider(ax=axbox, label='Atom/Residue #', valmin=0, valmax=len(atomic_grid_keys) - 1, valinit=0,
                            valstep=1)
            slider.on_changed(submit)

            save_button = Button(ax=save_axbox, label='Save current')
            save_button.on_clicked(lambda x: plt.savefig(f'{args.o}/gridplot_atom{current_plot}.png'))

            save_all_button = Button(ax=save_all_axbox, label='Save all')
            save_all_button.on_clicked(save_all)

            # Show
            plt.show()

    if str.lower(args.method) == 'violin' or (str.lower(args.method) == 'all' and args.plot_violin):
        """
        Using non-unique grid data, plot a violin plot to show x,y,z and density positions in a grid using violin plot       
        """
        if args.resi:
            try:
                with open(args.resi) as file:
                    resi_assignment = json.load(file)
                    residue_keys = list(resi_assignment.keys())
            except FileNotFoundError:
                print(f'{args.resi} not found, will not use.')
                args.resi = False

        def plot_violin(indexer):
            indexer = int(indexer)

            global chunk_size
            global axs_list

            # Clean up after previous plot
            try:
                for i in axs_list:
                    i.remove()
            except:
                axs_list = []

            global current_atoms_violin

            # catplot is a figure-level function, doesn't accept target axes, can't use it for this usage
            if args.s:
                violin_df = parse_violin(chunk_size * indexer, ((chunk_size * indexer) + chunk_size), vectors1=vectors1,
                                         traj1_name=traj1_name, vectors2=vectors2,
                                         traj2_name=traj2_name)  # start from atom 0, end with atom 3 (excluding)
                current_atoms_violin_arr = violin_df['atom_id'].unique()
                current_atoms_violin = ''
                for atom in current_atoms_violin_arr:
                    current_atoms_violin += atom
                current_atoms_violin.replace(' ', '_')

                # order = violin_df['atom_id'].unique()
                # hue_order = violin_df['traj_id'].unique()
                # sb.catplot(data=violin_df, x='datapoint_id', y='coords', hue='traj_id', col='atom_id', kind='violin', split='traj_id', sharex=False, inner='quartiles')
            else:
                violin_df = parse_violin(chunk_size * indexer, ((chunk_size * indexer) + chunk_size), vectors1=vectors1,
                                         traj1_name=traj1_name)
                # sb.catplot(data=violin_df, x='datapoint_id', y='coords', hue='axis_id', cut=0, inner='quartiles')

            for i, (n, grp) in enumerate(violin_df.groupby("atom_id")):
                _ = fig.add_subplot(1, chunk_size, i + 1)
                axs_list.append(_)
                # 1-row, 3-cols, index
                # sb.countplot(x="sex", hue="smoker", data=grp,
                #              order=order, hue_order=hue_order, ax=ax)
                if args.s:
                    axs_list[i] = sb.violinplot(x='datapoint_id', y='coords', hue='traj_id', data=grp, split='True',
                                                cut=0,
                                                inner='quartile')
                else:
                    axs_list[i] = sb.violinplot(x='datapoint_id', y='coords', hue='axis_id', data=grp, cut=0,
                                                inner='quartile')
                axs_list[i].set_title(f"atom = {n}")
                axs_list[i].get_legend().remove()
                # axs_list[i].get_xaxis().set_visible(False)
            axs_list[-1].legend()

            plt.draw()

        def plot_violin_complex(indexer):
            indexer = int(indexer)
            global ax1
            global ax2
            global ax3
            global ax4
            try:
                ax1.remove()
                ax2.remove()
                ax3.remove()
                ax4.remove()
            except:
                pass

            # Clean up after previous plot

            global current_atoms_violin

            # catplot is a figure-level function, doesn't accept target axes, can't use it for this usage
            if args.s:
                violin_df = parse_violin(chunk_size * indexer, ((chunk_size * indexer) + chunk_size), vectors1=vectors1,
                                         traj1_name=traj1_name, vectors2=vectors2,
                                         traj2_name=traj2_name)  # start from atom 0, end with atom 3 (excluding)
                current_atoms_violin_arr = violin_df['atom_id'].unique()
                current_atoms_violin = ''
                for atom in current_atoms_violin_arr:
                    current_atoms_violin += atom
                current_atoms_violin.replace(' ', '_')

                # order = violin_df['atom_id'].unique()
                # hue_order = violin_df['traj_id'].unique()
                # sb.catplot(data=violin_df, x='datapoint_id', y='coords', hue='traj_id', col='atom_id', kind='violin', split='traj_id', sharex=False, inner='quartiles')
            else:
                violin_df = parse_violin(chunk_size * indexer, ((chunk_size * indexer) + chunk_size), vectors1=vectors1,
                                         traj1_name=traj1_name)
                # sb.catplot(data=violin_df, x='datapoint_id', y='coords', hue='axis_id', cut=0, inner='quartiles')

            for i, (n, grp) in enumerate(violin_df.groupby("atom_id")):
                ax1 = fig.add_subplot(2, 2, 1)
                print('adding subplot')

                # 1-row, 3-cols, index
                # sb.countplot(x="sex", hue="smoker", data=grp,
                #              order=order, hue_order=hue_order, ax=ax)
                if args.s:
                    sb.violinplot(x='datapoint_id', y='coords', hue='traj_id', data=grp, split='True', cut=0,
                                                inner='quartile', ax=ax1)
                else:
                    axs_list[i] = sb.violinplot(x='datapoint_id', y='coords', hue='axis_id', data=grp, cut=0,
                                                inner='quartile', ax=ax1)
                ax1.set_title(f"atom = {n}")
                ax1.get_legend().remove()
                # axs_list[i].get_xaxis().set_visible(False)
            ax1.legend()


            ax2 = fig.add_subplot(2, 2, 2)
            ax3 = fig.add_subplot(2, 2, 3)
            ax4 = fig.add_subplot(2, 2, 4, projection='3d')

            trajs_vol = VOLDF.iloc[indexer, :]
            #traj1_volume_name = VOLDF[traj1_name][indexer]
            #traj2_volume_name = VOLDF[traj2_name][indexer]
            trajs_grid_rms = GRIDF[f'{traj1_name}/{traj2_name}']

            print(trajs_vol)

            # Volume explored by atom plot
            # x x
            # o x
            x = [traj1_name, traj2_name]
            y = [trajs_vol[f'V({traj1_name})/step'], trajs_vol[f'V({traj2_name})/step']]
            sb.barplot(x=x, y=y, ax=ax3)
            #axs_list[2].set_title(f"atom {indexer+1} explored volume per traj. step")

            # RMS plot
            rms_x = [i for i in range(1, len(GRIDF.iloc[:, 0])+1)]
            rms_y = GRIDF.iloc[:, 0]
            ax2.plot(rms_x, rms_y, zorder=1)
            ax2.scatter(rms_x[indexer], rms_y[indexer], s=40, color='red', marker='x', zorder=2)
            ax2.xaxis.set_major_locator(ticker.MultipleLocator(15))
            ax2.set_xticklabels(ax2.get_xticks(), rotation=90)

            print(GRIDF)



            # 3D plot for the atom
            atomic_grid_keys = list(atomic_grid.keys())
            atomic_grid_keys_2 = list(atomic_grid_2.keys())
            # print(atomic_grid_uniq[f'atom {current_plot+1}'])

            """ Get (X,Y,Z) unique grid coordinates for TRAJ1"""
            grids1 = tuple(map(tuple, atomic_grid_uniq[atomic_grid_keys[indexer]]))
            """ Get (X,Y,Z) unique grid coordinates for TRAJ2"""
            grids2 = tuple(map(tuple, atomic_grid_2_uniq[atomic_grid_keys_2[indexer]]))
            """ Have to keep order, so can't use sets with intersections etc. """

            grids1_unq = [grid for grid in grids1 if grid not in grids2]  # Coords where only grids 1 reside
            grids2_unq = [grid for grid in grids2 if grid not in grids1]  # Coords where only grids 2 reside
            grids_intersect = [grid for grid in grids1 if grid in grids2]  # Intersection of the two

            if len(grids1_unq) > 0:
                plot_cubes(data=grids1_unq, color='blue',
                           label=f'{traj1_name} only', ax=ax4)  # Plot cubes which are only present in dataset 1
            if len(grids2_unq) > 0:
                plot_cubes(data=grids2_unq, color='purple',
                           label=f'{traj2_name} only', ax=ax4)  # Plot cubes which are only present in dataset 2
            if len(grids_intersect) > 0:
                plot_cubes(data=grids_intersect, color='orange',
                           label=f'intersection', ax=ax4)  # Intersections (i.e. both 1 and 2 have a point here)

            ax4.legend(loc='upper left')
            set_ax_lims_3d(grids1, grids2, ax=ax4)
            #axs_list[3].set_title(f'Atom {indexer + 1}\ngridsize {args.grid}, pop. thresh. {args.pop_threshold}')

            #fig.tight_layout()
            plt.draw()

        def save_all_violin(_):
            for i in range(0, math.ceil(len(vector1_keys) / chunk_size)):
                plot_violin(i)
                plt.savefig(f'{args.o}/violinplot_{current_atoms_violin}.png')

        # Setup the plots
        fig, axs = plt.subplots(figsize=(15, 12))
        plt.axis('off')
        sb.set_theme(style="whitegrid")
        sb.despine(offset=10)
        plt.subplots_adjust(bottom=0.25)
        axbox = fig.add_axes([0.3, 0.05, 0.6, 0.075])
        save_axbox = fig.add_axes([0.8, 0.075, 0.2, 0.05])
        save_all_axbox = fig.add_axes([0.8, 0.025, 0.2, 0.05])
        ###

        global chunk_size


        if str.lower(args.method) == 'all' and args.s:

            chunk_size = 1

            #print(VOLDF)
            #print(GRIDF)

            #slider = Slider(ax=axbox, label='Atom #', valmin=0, valmax=len(vector1_keys) - 1,
            #                valinit=0, valstep=1)
            #slider.on_changed(plot_violin_complex)

            text_box = TextBox(axbox, "Atom", textalignment="center")
            text_box.on_submit(plot_violin_complex)
            text_box.set_val(0)

            save_button = Button(ax=save_axbox, label='Save current')
            save_button.on_clicked(lambda x: plt.savefig(f'{args.o}/violinplot_{current_atoms_violin}.png'))

            save_all_button = Button(ax=save_all_axbox, label='Save all')
            save_all_button.on_clicked(save_all_violin)
            # FIX CLEANING IN PLOT_VIOLIN (NOT COMPLEX)
            global axs_list
            axs_list = []

            plot_violin_complex(0)  # Initial plot

            plt.show()

        else:

            chunk_size = 3

            #slider = Slider(ax=axbox, label='Atom/Residue #', valmin=0, valmax=(len(vector1_keys) / chunk_size) - 1,
            #                valinit=0, valstep=1)
            #slider.on_changed(plot_violin)
            save_button = Button(ax=save_axbox, label='Save current')
            save_button.on_clicked(lambda x: plt.savefig(f'{args.o}/violinplot_{current_atoms_violin}.png'))

            save_all_button = Button(ax=save_all_axbox, label='Save all')
            save_all_button.on_clicked(save_all_violin)


            text_box = TextBox(axbox, "Atom", textalignment="center")
            text_box.on_submit(plot_violin_complex)
            text_box.set_val(0)

            plot_violin(0)  # Initial plot

            plt.show()

    # Add possibility to multiprocess the scan
    if str.lower(args.method) == 'grid_scan':

        print('Will perform a scan for an optimal value of --grid parameter. This may take a while.')
        fig, axs = plt.subplots()
        for gridsize in [3, 2, 1.5, 1, 0.5, 0.25, 0.10]:
            args.grid = gridsize
            atomic_grid = grid(vectors1)
            atomic_grid_uniq, grid_count = return_unique(atomic_grid)
            atomic_grid_2 = grid(vectors2)
            atomic_grid_2_uniq, grid_count_2 = return_unique(atomic_grid_2)

            ### Check statistics ###
            vector1_keys = list(vectors1.keys())
            vectors_num = sum([len(vectors1[key]) for key in vector1_keys])
            grid1_keys = list(atomic_grid.keys())
            grid2_keys = list(atomic_grid.keys())
            grid1_size = sum([len(atomic_grid[i]) for i in grid1_keys])
            grid2_size = sum([len(atomic_grid[i]) for i in grid2_keys])
            grid1_unq_size = sum([len(atomic_grid_uniq[i]) for i in grid1_keys])
            grid2_unq_size = sum([len(atomic_grid_2_uniq[i]) for i in grid2_keys])

            print(f'gs={gridsize} >Total amount of vectors {vectors_num * 2}')
            print(f'gs={gridsize} >Total amount of gridpoints {grid1_size + grid2_size}')
            print(f'gs={gridsize} >Total amount of unique gridpoints {grid1_unq_size + grid2_unq_size}')

            rms = grid_rms(atomic_grid_uniq, atomic_grid_2_uniq, grid_count, grid_count_2)

            rms_out = pd.Series(rms).T





            rms_keys = list(rms.keys())

            x = [x for x in range(0, len(rms_keys))]
            y = []
            for key in rms_keys:
                y.append(rms[key])

            norm_y1 = [float(i) / max(y) for i in y]

            rms_out = pd.Series(norm_y1).T
            rms_out = pd.DataFrame(rms_out, columns=[f'{traj1_name}/{traj2_name}'])
            rms_out.to_csv(f'{args.o}/{args.o}_grid_g{args.grid}_p{args.pop_threshold}.csv')


            axs.set_xlabel('Atom #')
            axs.set_ylabel('Apparent RMS / grid units^2')
            axs.set_title(f'norm. appRMS with --grid {args.grid}')
            axs.plot(x, norm_y1, label=args.grid)

        axs.legend()
        plt.show()





if __name__ == "__main__":
    main()

