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
import subprocess

__prepare_matplotlib()
caller_procs = 12 # Asynchronous processes for saving output files of 3D and violin plots


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


    divided = (vectors / args.grid)  # Divide all cells by args.grid #
    divided = divided.apply(np.floor).astype(np.int16) # Floor means that the grids start from 0,
    # not 0 (because then for example x=0.05 / grid=0.1 would end up in bin=0, which is the zeroth one)

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
    uniq_df.columns = [str(i) for i in uniq_df.columns]
    count_df.columns = [str(i) for i in count_df.columns]

    return(uniq_df, count_df)

def grid_rms(grid_unq1, grid_unq2, grid_count1, grid_count2, method='genius'):
    """
    Calculates the so-called "R-score", which is derived from classic molecular RMSD.
    We have two trajectories, where i-th and j-th indexes denote a unique gridpoint inside each trajectory
    (or a trajectory datapoint). We binned these (reducing accuracy) and identified unique triples (X,Y,Z).

    We have two datasets per trajectory now, one with unique triples and one with it's population per trajectory.
    The amount of samples is the sum of all counts (gets us the amount of points in original trajectories)
    Using matrix trickery, we calculate (computationally simplified) version of this equation, where 1, 2 denote
    trajectories.

    pairR_n-th_atom = sqrt[ (1/samples) * SUM_i SUM_j (X1i-X2j)^2 + (Y1i-Y2j)^2 + (Z1i-Z2j)^2  ]
    firstinternalR_n-th_atom = sqrt[ (1/samples) * SUM_i SUM_j (X1i-X1j)^2 + (Y1i-Y1j)^2 + (Z1i-Z1j)^2  ]
    secondinternalR_n-th_atom = sqrt[ (1/samples) * SUM_i SUM_j (X2i-X2j)^2 + (Y2i-Y2j)^2 + (Z2i-Z2j)^2  ]

    The R-score described such as this is unfortunately dependent on the distributions widths. If the two distributions
    (of the two trajectories) were identical, but both got wider, the R-score value would also rise, falsely indication
    an increasing change in the atom's positions. To solve this, we calculate "internal" R-scores for each of the distri-
    butions, which describe the widths of them.

    To obtain the final R-score for an atom, we calculate (2 x pairR_n-th_atom) - (firstinternalR_n-th_atom + secondinternalR_n-th_atom).

    The higher the R-score for the atom, the more it's positional distribution changed between the two trajectories.


    :param grid_unq1:
    :param grid_unq2:
    :param grid_count1:
    :param grid_count2:
    :param method:
    :return:
    """




    # Prepare RMSmetric list (one value for one compared atom)
    rms_lst = []

    if method == 'genius':

        # How many triples
        colsno_1 = len(grid_unq1.columns)
        colsno_2 = len(grid_unq2.columns)



        if colsno_1 != colsno_2:
            print('WARNING grid-datasets used for RMS calculation are of different size, RMS value may be faulty. Will iterate through'
                  'atoms until one of the datasets runs out.')

        # two lens, ceil and /2 because datasets can theoretically have different size
        # this way the bar really calculates progress in such a case


        # If we just straight up dropna, rows containing SOME None values would be dropped in entirety
        # This causes only some of the atom RMS to be calculated wrong, dependent on the atom with the least
        # unique values (the one where Nones start the earlies in the column)
        # This is solved by calculating column by column

        print(f'Calculating R-Score using gridsize {args.grid} nm')
        for col1, col2 in zip(grid_unq1.columns, grid_unq2.columns): # This works even if one of the datasets has less atoms
            rmsd_sum = 0

            triple1 = pd.DataFrame(grid_unq1[col1].dropna().tolist())
            triple2 = pd.DataFrame(grid_unq2[col2].dropna().tolist())
            triple1_counts = grid_count1[col1].dropna().tolist()
            triple2_counts = np.array(grid_count2[col2].dropna().tolist())  # has to be matrix multiplied by numpy, so has to be an array
            sum_counts1, sum_counts2 = sum(triple1_counts), sum(triple2_counts)
            triple1_counts = pd.DataFrame(triple1_counts)


            xs2, ys2, zs2 = tuple(triple2[0]), tuple(triple2[1]), tuple(triple2[2])

            samples = sum_counts1*sum_counts2

            def map_matrices(xyz):
                return(xyz[0] - xs2, xyz[1] - ys2, xyz[2] - zs2)

            def map_counts(xyz):
                return(np.array(xyz) * triple2_counts)

            # Prepare matrices
            matrix = pd.DataFrame(triple1.apply(map_matrices, axis=1))
            counts = triple1_counts.apply(map_counts, axis=1)

            # Calculate ij pairs (each row is a single i, each column inside the array is a single j) including their weights
            matrix = pd.DataFrame(matrix[0].to_list())
            matrix = matrix**2
            matrix = matrix.sum(axis=1)
            matrix = matrix * counts

            # Calculate sums and the final R-factor for current atom
            matrix_sum = matrix.sum().sum()
            rmsd = math.sqrt((1 / samples) * matrix_sum)

            # Append to list
            rms_lst.append(rmsd)

    elif method == 'smart':
        print(f'Calculating R-Score using gridsize {args.grid} nm')
        for col1, col2 in zip(grid_unq1.columns,
                              grid_unq2.columns):  # This works even if one of the datasets has less atoms
            rmsd_sum = 0

            triple1 = pd.DataFrame(grid_unq1[col1].dropna().tolist())
            triple2 = pd.DataFrame(grid_unq2[col2].dropna().tolist())
            triple1_counts = tuple(grid_count1[col1].dropna().tolist())
            triple2_counts = np.array(grid_count2[col2].dropna().tolist())  # has to be matrix multiplied by numpy, so has to be an array

            xs2, ys2, zs2 = tuple(triple2[0]), tuple(triple2[1]), tuple(triple2[2])
            sum_counts1, sum_counts2 = sum(triple1_counts), sum(triple2_counts)
            samples = sum_counts1 * sum_counts2

            # For i-th entry in first trajectory grids, for j-th entry in second trajectory grids
            i = 0

            matrix_x = []
            matrix_y = []
            matrix_z = []
            matrix_weight = []

            while i < len(triple1):  # This is faster than using list comprehension

                # X-term (single i with all j), Y-term, Z-term
                x_term_i = triple1.iloc[i, 0] - xs2  # Create a series for i=1 minus all possible j
                y_term_i = triple1.iloc[i, 1] - ys2
                z_term_i = triple1.iloc[i, 2] - zs2

                # Weights for all pairs
                i_weight = triple1_counts[i] * triple2_counts

                matrix_x.append(x_term_i)
                matrix_y.append(y_term_i)
                matrix_z.append(z_term_i)
                matrix_weight.append(i_weight)

                i += 1

            for ix, iy, iz, iweight in zip(matrix_x, matrix_y, matrix_z, matrix_weight):
                # print(ix)
                # print(iy)
                # print(iz)

                ix = ix ** 2
                iy = iy ** 2
                iz = iz ** 2

                ixyz = ix + iy + iz
                ixyz_w = ixyz * iweight

                rmsd_sum += ixyz_w.sum()

            rmsd = math.sqrt((1 / samples) * rmsd_sum)

            rms_lst.append(rmsd)

    elif method == 'brute-force':

        colsno_1 = len(grid_unq1.columns)
        colsno_2 = len(grid_unq2.columns)

        rms_lst_brute = []

        if colsno_1 != colsno_2:
            print(
                'WARNING grid-datasets used for RMS calculation are of different size, RMS value may be faulty. Will iterate through'
                'atoms until one of the datasets runs out.')
        for col1, col2 in zip(grid_unq1.columns, grid_unq2.columns): # This works even if one of the datasets has less atoms
            rmsd_sum = 0

            triple1 = pd.DataFrame(grid_unq1[col1].dropna().tolist())
            triple2 = pd.DataFrame(grid_unq2[col2].dropna().tolist())
            triple1_counts = tuple(grid_count1[col1].dropna().tolist())
            triple2_counts = tuple(grid_count2[col2].dropna().tolist())

            xs, ys, zs = tuple(triple1[0]), tuple(triple1[1]), tuple(triple1[2])
            xs2, ys2, zs2 = tuple(triple2[0]), tuple(triple2[1]), tuple(triple2[2])

            sum_counts1, sum_counts2 = sum(triple1_counts), sum(triple2_counts)

            samples = sum_counts1*sum_counts2

            i = 0
            j = 0

            rmsd_sum = 0
            while i < len(xs):
                j = 0
                while j < len(xs2):
                    rmsd_sum += ((xs[i]-xs2[j])**2 + (ys[i]-ys2[j])**2 + (zs[i]-zs2[j])**2) * (triple1_counts[i]*triple2_counts[j])
                    j += 1

                i += 1

            rmsd = math.sqrt((1/samples) * rmsd_sum)
            rms_lst_brute.append(rmsd)

    return(np.array(rms_lst))

def parse_violin(start, stop, vectors1, traj1_name, vectors2=pd.DataFrame(), traj2_name=False):

    def prep_half_violin(start, stop, vectors, traj_name):
        # Prepare base dataframes
        current_vectors = vectors[vectors.columns[start:stop]]

        melt_vect = current_vectors.melt()


        # Add axis_IDs for Seaborn
        len_curr_vect = len(current_vectors.iloc[:, 0])


        x_id = np.full(len_curr_vect, ['x'], dtype=str)
        y_id = np.full(len_curr_vect, ['y'], dtype=str)
        z_id = np.full(len_curr_vect, ['z'], dtype=str)
        axis_id = np.concatenate((x_id, y_id, z_id))

        how_many_atoms = int(len(melt_vect) / (3 * len_curr_vect))


        axis_ids = np.repeat(axis_id, how_many_atoms)
        axis_ids = pd.Series(axis_ids, name='axis_id')

        # Add atom_IDs
        atom_identifier = lambda atom: re.search("[a-zA-Z]+ [0-9]*", atom)[0].lower()
        atom_ids = melt_vect['variable'].apply(atom_identifier)

        # Add traj_id
        traj_ids = np.full(len(axis_ids), [traj_name])
        traj_ids = pd.Series(traj_ids, name='traj_id')

        # Put together to get a half of the violin df
        half_violin = pd.concat((melt_vect, axis_ids, traj_ids, atom_ids), axis=1)


        return(half_violin)

    violin_df = prep_half_violin(start, stop, vectors1, traj1_name)

    if not vectors2.empty:
        half_violin2 = prep_half_violin(start, stop, vectors2, traj2_name)
        violin_df = pd.concat((violin_df, half_violin2), axis=0)

    violin_df.columns = ('datapoint_id', 'coords', 'axis_id', 'traj_id', 'atom_id')
    violin_df['atom_id'] = violin_df['atom_id'].astype(str)
    violin_df['datapoint_id'] = violin_df['datapoint_id'].apply(str.lower) # Fixes a weird bug where vectors 2 atoms are lowercase, vectors 1 start with a capital A

    return (violin_df)

def _save_all_3d_plots(iter):
    global gfig, gax, gatomic_grid_uniq, gatomic_grid_uniq2, gtraj1_name, gtraj2_name

    _submit_3d_plot(gax, iter, gatomic_grid_uniq, gatomic_grid_uniq2, gtraj1_name, gtraj2_name)
    plt.savefig(f'{args.o}/gridplot_atom{iter}.png')

def plot_3d_handler(fig, ax, atomic_grid_uniq, atomic_grid_uniq2, traj1_name, traj2_name):

    global gfig, gax, gatomic_grid_uniq, gatomic_grid_uniq2, gtraj1_name, gtraj2_name
    gfig, gax, gatomic_grid_uniq, gatomic_grid_uniq2, gtraj1_name, gtraj2_name = fig, ax, atomic_grid_uniq, atomic_grid_uniq2, traj1_name, traj2_name


    def save_all_mp_caller(_):
        import multiprocessing as mp
        p = mp.Pool(caller_procs)
        r = p.map_async(_save_all_3d_plots, iterable=[ i for i in range(len(atomic_grid_uniq.columns)) ])

    global current_plot
    current_plot = 0
    _submit_3d_plot(ax, current_plot, atomic_grid_uniq, atomic_grid_uniq2, traj1_name, traj2_name)


    axbox = fig.add_axes([0.3, 0.05, 0.6, 0.075])
    textbox = fig.add_axes([0.15, 0.05, 0.10, 0.075])
    save_axbox = fig.add_axes([0.8, 0.075, 0.2, 0.05])
    save_all_axbox = fig.add_axes([0.8, 0.025, 0.2, 0.05])
    text_box = TextBox(textbox, "Atom #", textalignment="center")
    text_box.on_submit(lambda val: _submit_3d_plot(ax, val, atomic_grid_uniq, atomic_grid_uniq2, traj1_name, traj2_name))
    text_box.set_val("1")  # Trigger `submit` with the initial string.

    slider = Slider(ax=axbox, label='Atom #', valmin=1, valmax=len(atomic_grid_uniq.columns), valinit=1,
                    valstep=1)
    slider.on_changed(lambda val: _submit_3d_plot(ax, val, atomic_grid_uniq, atomic_grid_uniq2, traj1_name, traj2_name))

    save_button = Button(ax=save_axbox, label='Save current')
    save_button.on_clicked(lambda x: plt.savefig(f'{args.o}/gridplot_atom{current_plot}.png'))

    save_all_button = Button(ax=save_all_axbox, label='Save all')
    save_all_button.on_clicked(save_all_mp_caller)

    # Show
    plt.show()

def _submit_3d_plot(ax, expression, atomic_grid_uniq, atomic_grid_uniq2, traj1_name, traj2_name):
    global current_plot


    # Handle maximum and minimum possible plots and cases where indexer isn't int
    try:
        current_plot = int(expression)
    except ValueError:
        current_plot = 1

    if current_plot < 1:
        current_plot = 1
    elif current_plot > len(atomic_grid_uniq.columns):
        current_plot = len(atomic_grid_uniq.columns)
    current_plot -= 1

    """ Get (X,Y,Z) unique grid coordinates for TRAJ1"""

    tuplize = lambda xyz: tuple(xyz)

    current_atom_1 = atomic_grid_uniq.iloc[:, current_plot].dropna().apply(tuplize)
    current_atom_2 = atomic_grid_uniq2.iloc[:, current_plot].dropna().apply(tuplize)

    grids_intersect = current_atom_1[~current_atom_1.apply(tuple, 1).isin(current_atom_2.apply(tuple, 1))]
    grids1_unq = [grid for grid in current_atom_1 if grid not in current_atom_2]  # Coords where only grids 1 reside
    grids2_unq = [grid for grid in current_atom_2 if grid not in current_atom_1]  # Coords where only grids 2 reside
    # Intersection of the two

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
    _set_ax_lims_3d(current_atom_1, current_atom_2, ax=ax)
    ax.set_title(f'Atom {current_plot+1}\ngridsize {float(args.grid)}')

    plt.draw()
    return(ax)

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

def _save_all_violin(i):
    global gfig, gtraj1_name, gtraj2_name, gchunk_size

    _submit_violin(gfig, i, gchunk_size, gtraj1_name, gtraj2_name)
    plt.savefig(f'{args.o}/violinplot_{current_atoms_violin}.png')

def _submit_violin(fig, indexer, chunk_size, traj1_name, traj2_name, ax = False):
    global axs_list

    # Handle indexer being too low or too high or not an int
    try:
        indexer = int(indexer)
    except ValueError:
        indexer = 1

    indexer = int(indexer)
    if indexer < 1:
        indexer = 1
    elif indexer > len(vectors1.columns):
        indexer = len(vectors1.columns)
    indexer -= 1

    # Clean up after previous plot
    try:
        for i in axs_list:
            i.remove()
    except:
        axs_list = []

    global current_atoms_violin
    current_atoms_violin = ''

    # Get violin DFs
    # catplot is a figure-level function, doesn't accept target axes, can't use it for this usage
    if args.s:
        violin_df = parse_violin((chunk_size * indexer) * 3, ((chunk_size * indexer) + chunk_size) * 3,
                                 vectors1=vectors1,
                                 traj1_name=traj1_name, vectors2=vectors2,
                                 traj2_name=traj2_name)  # start from atom 0, end with atom 3 (excluding)

    else:
        violin_df = parse_violin((chunk_size * indexer) * 3, ((chunk_size * indexer) + chunk_size) * 3,
                                 vectors1=vectors1,
                                 traj1_name=traj1_name)


        # sb.catplot(data=violin_df, x='datapoint_id', y='coords', hue='axis_id', cut=0, inner='quartiles')
    #####

    # Get atom names for output files
    current_atoms_violin_arr = violin_df['atom_id'].unique()
    current_atoms_violin.replace(' ', '_')
    for atom in current_atoms_violin_arr:
        current_atoms_violin += atom
    #####


    # Iterate violin DFs atom by atom and plot
    for i, (n, grp) in enumerate(violin_df.groupby("atom_id")):
        _ = fig.add_subplot(1, chunk_size, i + 1)

        axs_list.append(_)
        if args.s:
            axs_list[i] = sb.violinplot(x='datapoint_id', y='coords', hue='traj_id', data=grp, split=True, cut=0,
                                        inner='quartile')
        else:
            axs_list[i] = sb.violinplot(x='datapoint_id', y='coords', hue='axis_id', data=grp, cut=0, inner='quartile')
        axs_list[i].set_title(f"atom = {n}")
        axs_list[i].get_legend().remove()
        # axs_list[i].get_xaxis().set_visible(False)
    axs_list[-1].legend()
    plt.draw()
    return(axs_list)
    #####



def pymol_dynamicity(file_name, traj1name, traj2name, type):

    # Coloring by B-factor in Pymol doesn't work when the factors are low. Maybe they only consider integers?
    # Solved by rescaling normalized datasets by a factor of 100x

    if type == 'diff':
        minimum_beta = -100
        session = 'visualized_diff.pse'
        spectrum = 'marine_gray70_raspberry'
        spectrum_ramp = '[marine, gray70, raspberry]'
        sphere_select = f'select more_in_{traj2name}, b > 0; select more_in_{traj1name}, b < 0; show_as sticks cartoon sphere,more_in_{traj1name};show_as sticks cartoon sphere,more_in_{traj2name};'

    else:
        minimum_beta = 0
        session = 'visualized_positions.pse'
        spectrum = 'gray70_pink_raspberry_purple'
        spectrum_ramp = '[gray70, pink, raspberry, purple]'
        sphere_select = ''

    pymol_sub2 = 'cmd.alter("*", "vdw=0.6")'
    pymol_sub3 = f'cmd.label("more_in_{traj1name}","str(ID)")'
    pymol_sub4 = f'cmd.label("more_in_{traj2name}","str(ID)")'
    pymol_command = f"set orthoscopic, on; bg_color white; ramp_new colorbar, none, [{minimum_beta}, 0, 100], {spectrum_ramp}; spectrum b, {spectrum}, minimum={minimum_beta}, maximum=100; {sphere_select} ;set seq_view; show lines; {pymol_sub2}; {pymol_sub3};{pymol_sub4}; set cartoon_discrete_colors, on; set valence, 1; set label_shadow_mode, 2; set label_size,-0.6; set label_font_id,7; set label_outline_color, black; set label_color, white; set label_position,(0,0,2); save {args.o}/{session}"

    p = subprocess.Popen(f"pymol {file_name} -d '{pymol_command}' &", stdout=subprocess.PIPE, shell=True, start_new_session=True)

def _set_ax_lims_3d(grids1, grids2, ax):


    triple1 = pd.DataFrame(grids1.tolist())
    triple2 = pd.DataFrame(grids2.tolist())
    triples = pd.concat((triple1, triple2), axis=0, ignore_index=True)


    minima = triples.min()
    widths = triples.max() - triples.min()

    width = widths.max() # find the "width" of the data in all axes and the widest one


    maxima = minima + width

    minima = minima*args.grid
    maxima = maxima*args.grid


    x_range, y_range, z_range = (minima[0], maxima[0]), (minima[1], maxima[1]), (minima[2], maxima[2])

    ax.set_xlabel('Bin(x) / nm')
    ax.set_ylabel('Bin(y) / nm')
    ax.set_zlabel('Bin(z) / nm')
    # ax.autoscale_view()
    ax.set(xlim=x_range, ylim=y_range, zlim=z_range)

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
    #parser.add_argument("--resi", type=str, help='OPTIONAL atoms-to-residues assignment file.json. Will turn on residue mode', required=False, default=False)


    # cartesian grid unique arguments
    parser.add_argument('--method', type=strlower, help='Method - grid, grid_scan, violin, volume, all; defaults to all', required=False, default='grid',
                        choices=('grid','grid_scan','violin','volume','all'))

    ## Grid method specifics
    parser.add_argument('--grid', type=float, help='Grid size in nm, defaults to 0.1 nm', required=False, default=0.1)
    #parser.add_argument('--pop_threshold', type=int,
                        #help='Bins with populations lower than this will be disregarded for plotting, defaults to 10 (good for throwing away SSAP caused artefacts)',
                        #required=False, default=10)
    #parser.add_argument('--mp', type=int, help='Nthreads to use for grid calculations, defaults to 1', required=False, default=False)
    parser.add_argument('--gridbackend', type=int, help='Backend to use for unique grid assignment', required=False, default=1, choices=(0,1))

    ## Plotting
    parser.add_argument("--p", type=str2bool, help='Disable all visualization and plotting, overrides all other plot parameters', default=False)
    parser.add_argument("--pdbshow", type=str2bool, help='Creates a PyMol visualization of dynamicity and position change, requires --pdbs. '
                                                    'Default=True', default=True)
    parser.add_argument("--plotrscore", type=str2bool, help='Plot R-scores from grid, defaults to True', const=True, default=True, nargs='?')


    parser.add_argument('--plot_3d', type=str2bool, help='OPTIONAL plot results in a 3D-plot', const=True, required=False, default=True, nargs='?')
    parser.add_argument('--plot_diff', type=str2bool, help='OPTIONAL plot explored volume by atom', const=True, required=False, default=True, nargs='?')
    parser.add_argument('--plot_violin', type=str2bool, help='OPTIONAL plot violin plots of spatial positions', const=True, required=False, default=True, nargs='?')
    ##

    global args
    args = parser.parse_args(argv)

    # Disable or enable all plotting with a single argument
    if args.p:
        args.pdbshow, args.plotrscore, args.plot_3d, args.plot_violin, args.plot_diff = False, False, False, False, False

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
                confname = f'{traj1_name}_{traj2_name}_grid_g{args.grid}.csv'
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
                delta = output_df[f'V({traj1_name})']*1000
                diff_pdb = 'REMARK Volume explored by atom in cartesian space PER STEP writen in Beta-factors in pm^3 (picometers)\n'
            else:
                print('(!!) Will print B-factors deltas between the two trajectories (s - f) in nm^3')
                delta = (output_df[f'V({traj2_name})'] - output_df[f'V({traj1_name})'])
                diff_pdb = 'REMARK Difference of volumes explored by atoms in two trajectories in cartesian space writen in Beta-factors in nm^3. Normalized -100 to 100 \n'


            # Normalize for better visualization between -1 and 1
            delta = (delta-delta.mean())/delta.std()
            delta = delta * 100 # Needed because of PyMol shenanigans


            diff_pdb += write_to_pdb_beta(pdb=args.pdbs, delta=delta)

            pdb_filename = f'{args.o}/{args.o}_diff.pdb'
            with open(pdb_filename, 'w') as file:
                file.write(diff_pdb)

            if args.pdbshow:
                pymol_dynamicity(pdb_filename, traj1_name, traj2_name, type='diff')


        if args.plot_diff:
            # Plot atom-wise diff
            fig, ax = plt.subplots(figsize=(15, 10))
            vol_atoms_1 = output_df[f'V({traj1_name})/step']
            name_atoms_1 = output_df[traj1_name]
            ax.plot(name_atoms_1, vol_atoms_1, color='blue', label=traj1_name)
            # Rework plotting with seaborn instead of matplotlib
            if args.s:
                vol_atoms_2 = output_df[f'V({traj2_name})/step']
                name_atoms_2 = output_df[traj2_name]
                ax.plot(name_atoms_2, vol_atoms_2, color='orange', label=traj2_name)
            ax.set_xlabel('Atom number')
            ax.set_ylabel('Explored volume per step / nm^3')
            ax.legend()
            plt.savefig(f'{args.o}/diff_plot_atom.png', dpi=300)

    if args.method == 'grid' or args.method == 'all':

        atomic_grid_uniq, atomic_grid_count = grid(vectors1)
        nongridded_vectors1_num = vectors1.count(axis=0).sum() # how many vectors before
        gridded_vectors1_num = atomic_grid_uniq.count(axis=0).sum() # how many vectors after binning
        nongridded_vectors2_num = 0
        gridded_vectors2_num = 0

        # Save grids for later use with cartesian_batch
        atomic_grid_uniq.to_parquet(f'{args.o}/{traj1_name}_grid_g{args.grid}.cart')
        atomic_grid_count.to_parquet(f'{args.o}/{traj1_name}_grid_g{args.grid}_count.cart')

        if args.s:
            atomic_grid_uniq2, atomic_grid_count2 = grid(vectors2)
            nongridded_vectors2_num = vectors2.count(axis=0).sum()
            gridded_vectors2_num = atomic_grid_uniq2.count(axis=0).sum()

            rms = grid_rms(atomic_grid_uniq, atomic_grid_uniq2, atomic_grid_count, atomic_grid_count2, method='genius')
            int_rms1 = grid_rms(atomic_grid_uniq, atomic_grid_uniq, atomic_grid_count, atomic_grid_count, method='genius')
            int_rms2 = grid_rms(atomic_grid_uniq2, atomic_grid_uniq2, atomic_grid_count2, atomic_grid_count2, method='genius')
            rms = (2 * rms) - (int_rms1+int_rms2) # Final R-score is a difference between twice the pairscore and a sum of internal R-scores of both distributions

            # Output results of grid method
            rms_out = pd.DataFrame(rms, columns=[f'{traj1_name}/{traj2_name}'])
            rms_out.to_csv(f'{args.o}/{traj1_name}_{traj2_name}_grid_g{args.grid}.csv')
            rms_out_norm = (rms_out-rms_out.min())/(rms_out.max()-rms_out.min())
            rms_out_norm.to_csv(f'{args.o}/{traj1_name}_{traj2_name}_grid_g{args.grid}_normalized.csv')

            atomic_grid_uniq2.to_parquet(f'{args.o}/{traj2_name}_grid_g{args.grid}.cart')
            atomic_grid_count2.to_parquet(f'{args.o}/{traj2_name}_grid_g{args.grid}_count.cart')

            if args.pdbs:
                delta = pd.Series(rms_out_norm.iloc[:, 0].values)
                delta = delta * 100 # Needed because of PyMol shenanigans
                print('(!!) Will print B-factors of conformational deltas between the two trajectories in arbitrary units. Normalized 0 to 100.')
                grid_pdb = 'REMARK Differences in exploring cartesian space (R-scores) per atom writen in Beta-factors in arbitrary units \n'
                grid_pdb += write_to_pdb_beta(pdb=args.pdbs, delta=delta)
                pdb_filename = f'{args.o}/{traj1_name}_{traj2_name}_grid{args.grid}.pdb'

                with open(pdb_filename, 'w') as file:
                    file.write(grid_pdb)

                if args.pdbshow:
                    pymol_dynamicity(pdb_filename, traj1_name, traj2_name, type='conf')

            if args.plotrscore:
                fig, ax = plt.subplots(figsize=(15, 12))
                x = rms_out.index + 1
                ax.plot(x.astype(int), rms_out.values.astype(float))
                ax.set_xlabel('Atom')
                ax.set_ylabel('R-score')
                ax.set_title(f'{traj1_name} vs {traj2_name} R-score with grid {args.grid}')
                ax.set(xlim=(1, len(x)))

        print(f'Total amount of coordinates: {nongridded_vectors1_num + nongridded_vectors2_num}, '
              f'after binning {gridded_vectors1_num + gridded_vectors2_num}')


        if args.plot_3d and args.s:

            # Adjust bottom to make room for Buttons
            fig, ax = plt.subplots(figsize=(15, 12))
            plt.axis('off')
            ax = fig.add_subplot(projection='3d')
            plt.subplots_adjust(bottom=0.25)

            plot_3d_handler(fig, ax, atomic_grid_uniq, atomic_grid_uniq2, traj1_name, traj2_name)

    if args.method == 'violin' and args.plot_violin and False:
        """
        Using non-unique pre-binning vector data, plot a violin plot to show x,y,z and density positions in a grid using violin plot       
        """

        def save_all_mp_caller(_):
            import multiprocessing as mp
            p = mp.Pool(caller_procs)
            r = p.map_async(_save_all_violin, iterable=[ i for i in range(1, math.ceil(len(vectors1.columns) / (3 * chunk_size))) ])

        # Control the amount of triples (xyz) plotted at once
        chunk_size = 1

        # Setup the plots
        fig, axs = plt.subplots(figsize=(15, 12))
        plt.axis('off')
        sb.set_theme(style="whitegrid")
        sb.despine(offset=10)
        plt.subplots_adjust(bottom=0.25)
        axbox = fig.add_axes([0.35, 0.05, 0.2, 0.075])
        save_axbox = fig.add_axes([0.6, 0.075, 0.2, 0.05])
        save_all_axbox = fig.add_axes([0.6, 0.025, 0.2, 0.05])
        ###

        global gfig, gtraj1_name, gtraj2_name, gchunk_size
        gfig, gtraj1_name, gtraj2_name, gchunk_size = fig, traj1_name, traj2_name, chunk_size

        # slider = Slider(ax=axbox, label='Atom/Residue #', valmin=0, valmax=(len(vector1_keys) / chunk_size) - 1,
        #                valinit=0, valstep=1)
        # slider.on_changed(plot_violin)
        save_button = Button(ax=save_axbox, label='Save current')
        save_button.on_clicked(lambda x: plt.savefig(f'{args.o}/violinplot_{current_atoms_violin}.png'))

        save_all_button = Button(ax=save_all_axbox, label='Save all')
        save_all_button.on_clicked(save_all_mp_caller)

        text_box = TextBox(axbox, "Atom #", textalignment="center")
        text_box.on_submit(lambda val: _submit_violin(fig, val, chunk_size, traj1_name, traj2_name))
        text_box.set_val('1')

        _submit_violin(fig, 1, chunk_size, traj1_name, traj2_name)  # Initial plot

        plt.show()

    if args.method == 'all' and args.s:
        def plot_violin_complex(indexer, fig):

            try:
                indexer = int(indexer)
            except:
                indexer = 1

            if indexer < 1:
                indexer = 1
            elif indexer > len(vectors1.columns):
                indexer = len(vectors1.columns)

            indexer -= 1

            global gax1
            global gax2
            global gax3
            global gax4
            try:
                gax1.remove()
                gax2.remove()
                gax3.remove()
                gax4.remove()
            except:
                pass


            gax1 = fig.add_subplot(2, 2, 1)
            gax2 = fig.add_subplot(2, 2, 2)
            gax3 = fig.add_subplot(2, 2, 3)
            gax4 = fig.add_subplot(2, 2, 4, projection='3d')

            # Clean up after previous plot
            global current_atoms_violin

            trajs_vol = output_df.iloc[indexer, :]  # Output_df is the output of differences from method=diff

            # Violin plot
            violin_df = parse_violin(indexer * 3, (indexer + 1) * 3,
                                     vectors1=vectors1,
                                     traj1_name=traj1_name, vectors2=vectors2,
                                     traj2_name=traj2_name)  # start from atom 0, end with atom 3 (excluding)
            current_atoms_violin = violin_df['atom_id'].unique()

            sb.violinplot(x='datapoint_id', y='coords', hue='traj_id', data=violin_df, split=True, cut=0,
                          inner='quartile', ax=gax1)



            # Volume explored by atom plot
            # x x
            # o x
            x = [traj1_name, traj2_name]
            y = [trajs_vol[f'V({traj1_name})/step'], trajs_vol[f'V({traj2_name})/step']]
            sb.barplot(x=x, y=y, ax=gax3)

            # RMS plot
            rms_x = rms_out.index
            rms_y = rms_out.values

            gax2.plot(rms_x, rms_y, zorder=1)
            gax2.scatter(rms_x[indexer], rms_y[indexer], s=40, color='red', marker='x', zorder=2)
            gax2.xaxis.set_major_locator(ticker.MultipleLocator(15))
            gax2.xaxis.set_ticks(gax2.get_xticks())
            gax2.set_xticklabels(gax2.get_xticks(), rotation=90)



            # 3D plot for the atom
            _submit_3d_plot(gax4, indexer, atomic_grid_uniq, atomic_grid_uniq2, traj1_name, traj2_name)
            #fig.tight_layout()

            #plt.draw()

            #slider = Slider(ax=axbox, label='Atom #', valmin=1, valmax=len(atomic_grid_uniq.columns),
                            #valinit=1, valstep=1)
            #slider.on_changed(lambda val: plot_violin_complex(val, fig))



        cfig, _ = plt.subplots(figsize=(15, 12))
        # Setup the plots
        plt.axis('off')
        sb.set_theme(style="whitegrid")
        sb.despine(offset=10)
        plt.subplots_adjust(bottom=0.25)

        axbox = cfig.add_axes([0.3, 0.05, 0.6, 0.075])
        save_axbox = cfig.add_axes([0.8, 0.075, 0.2, 0.05])
        save_all_axbox = cfig.add_axes([0.8, 0.025, 0.2, 0.05])
        save_button = Button(ax=save_axbox, label='Save current')
        save_button.on_clicked(lambda x: plt.savefig(f'{args.o}/allplot_{current_atoms_violin}.png'))

        save_all_button = Button(ax=save_all_axbox, label='Save all')
        save_all_button.on_clicked(lambda _: _)

        ###
        text_box = TextBox(axbox, "Atom #", textalignment="center", initial='1')
        text_box.on_submit(lambda val: plot_violin_complex(val, cfig))
        text_box.set_val(1)
        plot_violin_complex(1, cfig)  # Initial plot


        """
        Make an MP caller for saving en masse        
        
        """







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
    plt.show()  # use only one plt.show() at the end
