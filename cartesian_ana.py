"""
Combined cartesian_diff and cartesian_grid modules into cartesian_ana
"""

import sys, argparse, os
# ujson is recommended but not needed, speeds up json vectors loading
"""
try:
    import ujson as json
except ImportError:
    try:
        import simplejson as json
    except ImportError:
        import json
"""

# orjson could be used for speedup of json operations

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
from cartesian_diff import write_to_pdb_beta
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from cartesian import __prepare_matplotlib
from cartesian_diff import write_to_pdb_beta
from matplotlib.widgets import TextBox, Slider, Button
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

__prepare_matplotlib()


class NumpyEncoder(json.JSONEncoder): # Encodes numpy ndarrays into serializable lists for json
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def analyse_space(vector):
    """
    :param vector: Cartesian.py produced MD trajectory atomic vectors
    :return: Pandas DataFrame of volumes explored by atoms, approximated with Convex Hull algorithm in 3D-space
    """
    output_df = pd.DataFrame()
    atom_keys = list(vector.keys())
    with alive_bar(len(atom_keys)) as bar:
        print('Calculating vector explored volumes')
        for atom_key in atom_keys:
            # Histogramming for plots
            x = [i[0] for i in vector[atom_key]]
            y = [i[1] for i in vector[atom_key]]
            z = [i[2] for i in vector[atom_key]]

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

            out_series = pd.Series({'atom' : atom_key,
                                    'vol': outvol,
                                    'volstep': outvolstep
                                    })
            output_df = output_df.append(out_series, True)
            bar()
    return(output_df)

def assign_into_grid_member(coordinate):
    # print(coordinate)
    if (coordinate / args.grid).is_integer():
        # if the atom coordinate is a whole number, we want to make it a 50 % chance it will be
        # assigned to the lower grid or 50 % to upper grid (so grids are, with enough samples, of the same size)
        if coordinate % 2 == 0:  # had if atomX % 2 ==0 ; WEIRD? !!!!!
            grid_x = (
                             coordinate / args.grid) - 1  # If the atomX is an even number, assign to "the left" (into lower grid member)
        else:
            grid_x = coordinate / args.grid  # If odd, assign to higher grid member
    else:
        grid_x = math.floor(coordinate / args.grid)  # Assign into the grid member
    return (int(grid_x))

def atom_atomic_grids_asynch(key):
    new_grids = []
    atomXs = [vector[0] for vector in glob_vectors[key]]
    atomYs = [vector[1] for vector in glob_vectors[key]]
    atomZs = [vector[2] for vector in glob_vectors[key]]
    for atomX, atomY, atomZ in zip(atomXs, atomYs, atomZs):
        new_grid_tuple = [None, None, None]
        new_grid_tuple[0] = assign_into_grid_member(atomX)
        new_grid_tuple[1] = assign_into_grid_member(atomY)
        new_grid_tuple[2] = assign_into_grid_member(atomZ)
        new_grids.append(tuple(new_grid_tuple))
    return_dict = {}
    return_dict[key] = tuple(new_grids)

    return (return_dict)

def grid(vectors):
    vectors_keys = list(vectors.keys())

    def scan_through_grid(vector):
        vectors_grid_positions = {}

        """
        Scan for each of the atoms in vector dict
        If there are grid points (0.0, 0.5, 1.0, 1.5) and there's two coordinates (0.5, 1.0),
        the first coordinate=0.5 will be assigned to the first grid (0.0-0.5), the second coordinate=1.0
        will be assigned to 3rd grid (1.0-1.5). This is to solve uncertainty in assigning and with big enough
        dataset, statistically the probabilities should both be 50 % and error caused by this should "cancel out"
        somewhat.        

        As compared to the previous calculation method, where grids were pre-assigned and then several nested for loops were used,
        we're using simple math to assign coordinates into correct brackets. This is much faster and much more elegant.        
        """

        global glob_vectors
        glob_vectors = vectors
        glob_vectors_keys = tuple(glob_vectors.keys())

        """
        Improve with numpy or pandas, apply the function to the whole array,
        or use Python mapping, don't iterate       

        """

        atomic_grids = {}
        if args.mp:

            print(f'Using multiprocessing with {args.mp}')
            # Assign coordinates to a dictionary
            pool = mp.Pool(args.mp)

            atomic_grids = pool.map_async(atom_atomic_grids_asynch, glob_vectors_keys).get()
            # atomic_grids[atom] = result.get()
            pool.close()
            pool.join()
            atomic_grids = {list(grid_tuple.keys())[0]: tuple(grid_tuple.values())[0] for grid_tuple in atomic_grids}

        else:
            print('NOT using multiprocessing')
            for atom in vectors_keys:
                atomic_grids[atom] = []

                atomXs = [vector[0] for vector in vector[atom]]
                atomYs = [vector[1] for vector in vector[atom]]
                atomZs = [vector[2] for vector in vector[atom]]

                for atomX, atomY, atomZ in zip(atomXs, atomYs, atomZs):
                    new_grid_tuple = [None, None, None]
                    new_grid_tuple[0] = assign_into_grid_member(atomX)
                    new_grid_tuple[1] = assign_into_grid_member(atomY)
                    new_grid_tuple[2] = assign_into_grid_member(atomZ)
                    atomic_grids[atom].append(new_grid_tuple)

        return (atomic_grids)

    """
    Return unique returns only unique grid members
    and the "population" of them (how many times the vector
    visited these grid members during the trajectory)
    """

    # Scan through the grid space for vectors
    grid_pop_vectors = scan_through_grid(vectors)

    return (grid_pop_vectors)

def grid_rms(grid1, grid2, grid_count1, grid_count2):
    atom1_keys = list(grid1.keys())
    atom2_keys = list(grid2.keys())
    if len(atom1_keys) != len(atom2_keys):
        print(
            'WARNING grid-datasets used for RMS calculation are of different size, RMS values will suffer some uncertainty, atomic indices will be assigned from the second trajectory!')

    rmsd_lst = {}

    # two lens, ceil and /2 because datasets can theoretically have different size
    # this way the bar really calculates progress in such a case
    with alive_bar(int(math.ceil(len(atom1_keys) + len(atom2_keys)) / 2)) as bar:
        print(f'Calculating inter-grid RMS using gridsize {args.grid} nm')
        for key1, key2 in zip(atom1_keys, atom2_keys):
            rmsd_sum = 0
            atom1_xs = [x[0] for x in grid1[key1]]
            atom1_ys = [x[1] for x in grid1[key1]]
            atom1_zs = [x[2] for x in grid1[key1]]
            atom2_xs = [x[0] for x in grid2[key2]]
            atom2_ys = [x[1] for x in grid2[key2]]
            atom2_zs = [x[2] for x in grid2[key2]]

            current_atom_grid_count = grid_count1[key1]
            current_atom_grid_count2 = grid_count2[key2]
            # current_atom_internal_rms = int_rms_1[key1]
            # current_atom_2_internal_rms = int_rms_2[key2]

            # For i-th entry in atom1(K)
            i = 0
            j = 0
            while i < len(atom1_xs):

                if j == len(atom2_xs):
                    i += 1
                    j = 0
                    continue

                rmsd_sum += (current_atom_grid_count[i] * current_atom_grid_count2[j]) * (
                            atom1_xs[i] - atom2_xs[j]) ** 2 + (atom1_ys[i] - atom2_ys[j]) ** 2 + (
                                        atom1_zs[i] + atom2_zs[j]) ** 2
                j += 1

            rmsd = math.sqrt((1 / sum(current_atom_grid_count)) * rmsd_sum)
            # rmsd += -(current_atom_internal_rms + current_atom_2_internal_rms) Subtract internal RMS of the atoms
            # print(current_atom_internal_rms)
            rmsd_lst[f'{key1}-{key2}'] = rmsd
            bar()
    return (rmsd_lst)

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

def return_unique(grid_pop_vectors):
    grid_keys = list(grid_pop_vectors.keys())
    new_dict = {}
    count_dict = {}

    backend = 1  # Numpy backend doesn't seem to work later when asking for dataset differences and intersection

    if backend == 0:
        for grid_key in grid_keys:
            new_dict[grid_key] = []
            count_dict[grid_key] = []
            print(f'Returning unique {grid_key}')
            # For each grid visited by the atom
            for grid_coordinates in grid_pop_vectors[grid_key]:
                if not grid_coordinates in new_dict[grid_key]:

                    # Throw away bins that have populations below args.pop_threshold, append the rest
                    count = grid_pop_vectors[grid_key].count(grid_coordinates)
                    if count > args.pop_threshold:
                        new_dict[grid_key].append(grid_coordinates)
                        count_dict[grid_key].append(count)

    if backend == 1:
        # Improved and faster with numpy
        with alive_bar(len(grid_keys)) as bar:
            print('Finding unique bins')
            for grid_key in grid_keys:
                current_vectors = np.array(grid_pop_vectors[grid_key])
                uniqs, counts = np.unique(current_vectors, axis=0,
                                          return_counts=True)  # Returns unique tuples of X,Y,Z (axis=0) and their counts in the original array (return_counts=True)
                which_to_delete = []
                for i, count in enumerate(counts):
                    if count < args.pop_threshold:
                        which_to_delete.append(i)
                new_dict[grid_key] = np.delete(uniqs, which_to_delete, axis=0)
                count_dict[grid_key] = np.delete(counts, which_to_delete, axis=0)
                bar()
    """
    Returns something like this (new_dict)
    {'p1': [[1, 0, 0], [-1, 0, -1]], 'p2': [[2, 2, 5], [1, 2, 5], [2, 2, 6]], 'p3': [[-9, -3, 3], [-10, -4, 3], [-10, -3, 3]]}
    Read: particle 1 uniquely visited gridpoint (1,0,0) gridpoint (-1,0,-1).
          particle 2 uniquelyvisidet gridpoint (2,2,5), gridpoint (1,2,5), gridpoint (2,2,6)
    And this (count_dict)
    {'p1': [2, 1], 'p2': [1, 1, 1], 'p3': [1, 1, 1]}
    Read: particle 1 visited the first unique gridpoint in new_dict (1,0,0) 2x, the (-1,0,-1) gridpoint was visited once.

    """
    return (new_dict, count_dict)

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



def main(argv=sys.argv[1:]):

    ## General I/O
    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=str, help='Input vector.json for trajectory 1', required=True)
    parser.add_argument("--s", type=str, help='OPTIONAL Input vector.json for trajectory 2', default=False)
    parser.add_argument("--o", type=str, help='Output directory, defaults to names of trajectories used separated by _', required=False, default=False)
    parser.add_argument("--pdbs", type=str, help='OPTIONAL Input structure in .pdb format of the second trajectory', default=False)
    parser.add_argument("--resi", type=str, help='OPTIONAL atoms-to-residues assignment file.json. Will turn on residue mode', required=False, default=False)


    # cartesian grid unique arguments
    parser.add_argument('--method', type=str, help='Method - grid, grid_scan, violin, volume, all; defaults to grid', required=False, default='grid')

    ## Grid method specifics
    parser.add_argument('--grid', type=float, help='Grid size in nm, defaults to 0.5', required=False, default=0.5)
    parser.add_argument('--pop_threshold', type=int,
                        help='Bins with populations lower than this will be disregarded for plotting, defaults to 10 (good for throwing away SSAP caused artefacts)',
                        required=False, default=10)
    parser.add_argument('--mp', type=int, help='Nthreads to use for grid calculations, defaults to 1', required=False, default=False)

    ## Plotting
    parser.add_argument("--plot", type=str2bool, help='Plot spatial stuff, defaults to False', const=True, default=True, nargs='?')
    parser.add_argument('--plot_3d', type=str2bool, help='OPTIONAL plot results in a 3D-plot', const=True, required=False, default=True, nargs='?')
    parser.add_argument('--plot_positions', type=str2bool, help='OPTIONAL plot positional violin plots', const=True, required=False, default=True, nargs='?')
    parser.add_argument('--plot_diff', type=str2bool, help='OPTIONAL plot explored volume by atom', const=True, required=False, default=True, nargs='?')
    parser.add_argument('--plot_violin', type=str2bool, help='OPTIONAL plot violin plots of spatial positions', const=True, required=False, default=True, nargs='?')




    ##


    global args
    args = parser.parse_args(argv)

    # Set up default output directory, trajectory names
    traj1_name = args.f.replace('.json', '')
    if args.s:
        traj2_name = args.s.replace('.json', '')

    if not args.o:
        if args.s:
            args.o = f'{traj1_name}_{traj2_name}'
        else:
            args.o = traj1_name
    if not os.path.exists(args.o):
        os.makedirs(args.o)
    # Load vectors
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
                # Append to a file that is later used to make sense of output folders for cartesian_batch.py
                # Traj1, traj2, directory, dynamical data, conformational data

                """
                This part only works if you're using two trajectories, perhaps modify later, so it prints out an outputs.txt
                file for cartesian_batch.py even if we're only interested in obtaining distributions with a single batch of trajectories
                """

                with open('outputs.txt', 'at') as file:
                    if args.resi:
                        diffname = 'diff_resis.csv'
                    else:
                        diffname = 'diff_atom.csv'
                    confname = f'{args.o}_grid_g{args.grid}_p{args.pop_threshold}.csv'
                    grids1name = f'{traj1_name}_grid_g{args.grid}.json'
                    grids2name = f'{traj2_name}_grid_g{args.grid}.json'
                    grids1count = f'{traj1_name}_grid_g{args.grid}_count.json'
                    grids2count = f'{traj2_name}_grid_g{args.grid}_count.json'

                    file.write(f'{args.f},{args.s},{args.o},{diffname},{confname},{grids1name},{grids2name},{grids1count},{grids2count}\n')

        except:
            print(f'{args.s} not found, will proceed with only one vector file analysis, obtained from --f')
            args.s = False





    if str.lower(args.method) == 'volume' or str.lower(args.method) == 'all':
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
            output_df.loc[output_df.index[0], f'SUMV({traj1_name})/step'] = tot_explored_volume_1 / len(
                list(vectors1.keys()))
            output_df.loc[output_df.index[0], f'SUMV({traj2_name})'] = tot_explored_volume_2
            output_df.loc[output_df.index[0], f'SUMV({traj2_name})/step'] = tot_explored_volume_2 / len(
                list(vectors2.keys()))
        else:
            output_df = output_df1
            tot_explored_volume_1 = output_df[f'V({traj1_name})'].sum()
            output_df.loc[output_df.index[0], f'SUMV({traj1_name})'] = tot_explored_volume_1
            output_df.loc[output_df.index[0], f'SUMV({traj1_name})/step'] = tot_explored_volume_1 / len(
                list(vectors1.keys()))

        output_df.to_csv(f'{args.o}/diff_atom.csv', float_format='{:.12f}')

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
                min_atom = atoms_in_residue[0]
                max_atom = atoms_in_residue[-1]
                sum_vol1 = output_df.loc[f'atom {min_atom}':f'atom {max_atom}', f'V({traj1_name})'].sum()
                sum_vol1_step = output_df.loc[f'atom {min_atom}':f'atom {max_atom}', f'V({traj1_name})/step'].sum()

                if args.s:
                    sum_vol2 = output_df.loc[f'atom {min_atom}':f'atom {max_atom}', f'V({traj2_name})'].sum()
                    sum_vol2_step = output_df.loc[f'atom {min_atom}':f'atom {max_atom}', f'V({traj2_name})/step'].sum()
                    residue_df = pd.DataFrame({
                        f'residue': key,
                        f'V({traj1_name})': sum_vol1,
                        f'V({traj1_name})/step': sum_vol1_step,
                        f'V({traj2_name})': sum_vol2,
                        f'V({traj2_name})/step': sum_vol2_step
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

        if args.plot_diff:
            # Plot atom-wise diff
            fig, ax = plt.subplots(figsize=(15, 10))
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
                fig_resi, ax_resi = plt.subplots(figsize=(15, 10))
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
        global VOLDF
        VOLDF = output_df

    if str.lower(args.method) == 'grid' or str.lower(args.method) == 'all':

        atomic_grid = grid(vectors1)
        atomic_grid_uniq, grid_count = return_unique(atomic_grid)
        # Calculate original points amount
        vector1_keys = list(vectors1.keys())
        vectors_num = sum([len(vectors1[key]) for key in vector1_keys])
        # Calculate the total amount of gridpoints in the first datasets
        grid1_keys = list(atomic_grid.keys())
        grid1_size = sum([len(atomic_grid[i]) for i in grid1_keys])
        # Calculate total amount of unique gridpoints
        grid1_unq_size = sum([len(atomic_grid_uniq[i]) for i in grid1_keys])

        if args.s:
            atomic_grid_2 = grid(vectors2)
            atomic_grid_2_uniq, grid_count_2 = return_unique(atomic_grid_2)
            grid2_keys = list(atomic_grid.keys())
            grid2_size = sum([len(atomic_grid[i]) for i in grid2_keys])
            grid2_unq_size = sum([len(atomic_grid_2_uniq[i]) for i in grid2_keys])
            rms = grid_rms(atomic_grid_uniq, atomic_grid_2_uniq, grid_count, grid_count_2)
            rms_out = pd.Series(rms).T
            rms_out = pd.DataFrame(rms_out, columns=[f'{traj1_name}/{traj2_name}'])
            rms_out.to_csv(f'{args.o}/{args.o}_grid_g{args.grid}_p{args.pop_threshold}.csv')
        else:
            grid2_size = 0
            grid2_unq_size = 0

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

