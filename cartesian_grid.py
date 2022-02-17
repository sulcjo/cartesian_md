import argparse
import sys, os
import json
import numpy as np
import math
import matplotlib
#matplotlib.use('GTK3Agg') #depends on pycairo, PyGObject packages
from alive_progress import alive_bar
from cartesian_diff import write_to_pdb_beta
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from cartesian import __prepare_matplotlib
from cartesian_diff import write_to_pdb_beta

# Prepare matplotlib parameters
__prepare_matplotlib()

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

        atomic_grids = {}

        def assign_into_grid_member(coordinate):
            #print(coordinate)
            if (coordinate / args.grid).is_integer():
                # if the atom coordinate is a whole number, we want to make it a 50 % chance it will be
                # assigned to the lower grid or 50 % to upper grid (so grids are, with enough samples, of the same size)
                if atomX % 2 == 0:
                    grid_x = (coordinate / args.grid) - 1  # If the atomX is an even number, assign to "the left" (into lower grid member)
                else:
                    grid_x = coordinate / args.grid  # If odd, assign to higher grid member
            else:
                grid_x = math.floor(coordinate / args.grid)  # Assign into the grid member
            return(int(grid_x))

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

        return(atomic_grids)

    """
    Return unique returns only unique grid members
    and the "population" of them (how many times the vector
    visited these grid members during the trajectory)
    """

    # Scan through the grid space for vectors
    grid_pop_vectors = scan_through_grid(vectors)

    return(grid_pop_vectors)

def return_unique(grid_pop_vectors):
    grid_keys = list(grid_pop_vectors.keys())
    new_dict = {}
    count_dict = {}
    # For each atom
    for grid_key in grid_keys:
        new_dict[grid_key] = []
        count_dict[grid_key] = []
        # For each grid visited by the atom
        for grid_coordinates in grid_pop_vectors[grid_key]:
            if not grid_coordinates in new_dict[grid_key]:
                new_dict[grid_key].append(grid_coordinates)
                count = grid_pop_vectors[grid_key].count(grid_coordinates)
                count_dict[grid_key].append(count)

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

            #print(atom_ys)
            #print(atom_zs)


            i = 0
            j = 0
            while i < len(atom_xs):
                #print(i,j)
                #print(rmsd)
                if j == len(atom_xs)-1:
                    i += 1
                    j = 0
                    continue
                if i != j:
                    rmsd_sum += (current_atom_grid_count[i]*current_atom_grid_count[j])*(atom_xs[i]-atom_xs[j])**2 + (atom_ys[i]-atom_ys[j])**2 + (atom_zs[i]+atom_zs[j])**2
                j += 1

            rmsd = math.sqrt((1/2)*(1/len(atom_xs))*rmsd_sum)
            rmsd_dict[atom]=rmsd
            bar()
    return(rmsd_dict)

def grid_rms(grid1, grid2, grid_count1, grid_count2):
    atom1_keys = list(grid1.keys())
    atom2_keys = list(grid2.keys())
    if len(atom1_keys) != len(atom2_keys):
        print('WARNING grid-datasets used for RMS calculation are of different size, RMS values will suffer some uncertainty, atomic indices will be assigned from the second trajectory!')

    rmsd_lst = {}

    # two lens, ceil and /2 because datasets can theoretically have different size
    # this way the bar really calculates progress in such a case
    with alive_bar(int(math.ceil(len(atom1_keys)+len(atom2_keys))/2)) as bar:
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
            #current_atom_internal_rms = int_rms_1[key1]
            #current_atom_2_internal_rms = int_rms_2[key2]

            # For i-th entry in atom1(K)
            i = 0
            j = 0
            while i < len(atom1_xs):

                if j == len(atom2_xs):
                    i += 1
                    j = 0
                    continue

                rmsd_sum += (current_atom_grid_count[i] * current_atom_grid_count2[j]) * (atom1_xs[i] - atom2_xs[j]) ** 2 + (atom1_ys[i] - atom2_ys[j]) ** 2 + (atom1_zs[i] + atom2_zs[j]) ** 2
                j += 1

            rmsd = math.sqrt((1 / sum(current_atom_grid_count)) * rmsd_sum)
            #rmsd += -(current_atom_internal_rms + current_atom_2_internal_rms) Subtract internal RMS of the atoms
            #print(current_atom_internal_rms)
            rmsd_lst[f'{key1}-{key2}'] = rmsd
            bar()
    return(rmsd_lst)

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', type=str, help='First set of vectors in .json', required=True, default=False)
    parser.add_argument('--s', type=str, help='Second set of vectors in .json', required=False, default=False)
    parser.add_argument('--method', type=str, help='Method - grid, grid_scan; defaults to grid', required=False, default='grid')
    parser.add_argument('--grid', type=float, help='Grid size in nm, defaults to 0.5', required=False, default=0.5)
    #parser.add_argument('--threshold', type=int, help='Threshold in gridpoints for perturbation detection, defaults to 1', required=False, default=1)
    parser.add_argument('--o', type=str, help='Output directory, defaults to names of trajectories used separated by _', required=False)
    parser.add_argument('--pdbs', type=str, help='OPTIONAl .pdb file to generate rms-based coloring', required=False, default=False)
    parser.add_argument('--plot', type=str, help='OPTIONAL plot results in a 3D-plot', required=False, default=False)
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

    # Load vectors generated by cartesian.py
    try:
        with open(args.f) as file:
            vectors1 = json.load(file)
        with open(args.s) as file:
            vectors2 = json.load(file)
    except FileNotFoundError:
        print(f'No {args.f} or {args.s} found, will quit')
        exit()
    ###

    if args.method == 'grid':

        """
        # TEST VECTORS
        vectors1 = {
            'p1': [(0,0,0), (0,0,0), (0,0,0)],
        }

        vectors2 = {
            'p1': [(1, 0, 0), (0, 1, 0), (1, 0, 0)],
            'p2': [(2.3, 2.4, 5.6), (1.9, 2.45, 5.93), (20.01, 2.47, 60.00)],
            'p3': [(-8.6, -2.5, 3.0), (-9.78, -3.15, 3.6), (-9.99, -2.5, 3.7)]
        }
         """
        atomic_grid = grid(vectors1)
        atomic_grid_uniq, grid_count = return_unique(atomic_grid)
        atomic_grid_2 = grid(vectors2)
        atomic_grid_2_uniq, grid_count_2 = return_unique(atomic_grid_2)

        # Calculate original points amount
        vector1_keys = list(vectors1.keys())
        vectors_num = sum([len(vectors1[key]) for key in vector1_keys])

        # Calculate the total amount of gridpoints in both datasets
        grid1_keys = list(atomic_grid.keys())
        grid2_keys = list(atomic_grid.keys())
        grid1_size = sum([len(atomic_grid[i]) for i in grid1_keys])
        grid2_size = sum([len(atomic_grid[i]) for i in grid2_keys])

        # Calculate total amount of unique gridpoints
        grid1_unq_size = sum([len(atomic_grid_uniq[i]) for i in grid1_keys])
        grid2_unq_size = sum([len(atomic_grid_2_uniq[i]) for i in grid2_keys])

        print(f'Total amount of vectors {vectors_num*2}')
        print(f'Total amount of gridpoints {grid1_size+grid2_size}')
        print(f'Total amount of unique gridpoints {grid1_unq_size+grid2_unq_size}')

        # Calculate internal RMSD of grid-set 1 and 2
        #int_rms_1 = internal_grid_rms(atomic_grid_uniq, grid_count)
        #int_rms_2 = internal_grid_rms(atomic_grid_2_uniq, grid_count_2)
        #print(int_rms_1)
        #print(int_rms_2)

        rms=grid_rms(atomic_grid_uniq, atomic_grid_2_uniq, grid_count, grid_count_2)
        rms_out=pd.Series(rms).T
        rms_out=pd.DataFrame(rms_out, columns=[f'{traj1_name}/{traj2_name}'])
        rms_out.to_csv(f'{args.o}/{args.o}_grid.csv')


    if args.pdbs:
        print('(!!) Will print B-factors of conformational deltas between the two trajectories in arbitrary units')
        delta = pd.Series(rms)

        grid_pdb = 'REMARK Differences in exploring cartesian space per atom writen in Beta-factors in arbitrary units \n'
        grid_pdb += write_to_pdb_beta(pdb=args.pdbs, delta=delta)

        with open(f'{args.o}/{args.o}_grid.pdb', 'w') as file:
            file.write(grid_pdb)
    if args.plot:

        atomic_grid_keys = list(atomic_grid.keys())
        atomic_grid_keys_2 = list(atomic_grid_2.keys())
        count_dict_keys = list(grid_count.keys())
        count_dict_keys_2 = list(grid_count_2.keys())
        fig = plt.figure()



        maxcols=8
        nrows = math.ceil(len(atomic_grid_keys) / maxcols)
        col = 0
        row = 0
        def coloring(pop_values):
            colors = []
            for value in pop_values:

                if 2 >= value > 0:
                    colors.append('cyan')
                elif 5 >= value > 2:
                    colors.append('blue')
                elif 10 >= value > 5:
                    colors.append('orange')
                else:
                    colors.append('red')
            return(colors)

        i = 0
        for key, key_2 in zip(atomic_grid_keys, atomic_grid_keys_2):

            xs = [x[0] for x in atomic_grid[key]]
            ys = [x[1] for x in atomic_grid[key]]
            zs = [x[2] for x in atomic_grid[key]]

            xs_2 = [x[0] for x in atomic_grid_2[key_2]]
            ys_2 = [x[1] for x in atomic_grid_2[key_2]]
            zs_2 = [x[2] for x in atomic_grid_2[key_2]]

            #colors = coloring(count_dict[key])
            ax = fig.add_subplot(nrows, maxcols, i+1, projection='3d')



            ax.scatter(xs=xs, ys=ys, zs=zs, s=10, c='blue', marker='s')
            ax.scatter(xs=xs_2, ys=ys_2, zs=zs_2, s=10, c='red', marker='s')
            ax.set_title(key+' '+key_2)
            ax.set(xlim=(-10,10), ylim=(-10,10), zlim=(-10,10))
            ax.set_xlabel('Grid X')
            ax.set_ylabel('Grid Y')
            ax.set_zlabel('Grid Z')

            col += 1
            i += 1
        fig.tight_layout()

        plt.show()
        """

    if args.method == 'grid_scan':
        ### TEST
        print('Will perform a scan for an optimal value of --grid parameter. This may take a while.')
        fig, axs = plt.subplots()
        for gridsize in [3, 2, 1.5, 1, 0.5, 0.25]:
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
            ### ###




            rms = grid_rms(atomic_grid_uniq, atomic_grid_2_uniq, grid_count, grid_count_2)

            rms_keys = list(rms.keys())

            x = [x for x in range(0, len(rms_keys))]
            y = []
            for key in rms_keys:
                y.append(rms[key])

            norm_y1 = [float(i) / max(y) for i in y]
            axs.set_xlabel('Atom #')
            axs.set_ylabel('Apparent RMS / grid units^2')
            axs.set_title(f'norm. appRMS with --grid {args.grid}')
            axs.plot(x, norm_y1, label=args.grid)

        axs.legend()
        plt.show()

        ### /TEST

if __name__ == '__main__':
    main()