import argparse
import sys
import json
import numpy as np
import math

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
        return(new_dict, count_dict)

    # Scan through the grid space for vectors
    grid_pop_vectors = scan_through_grid(vectors)
    grid_pop_vectors_uniq = return_unique(grid_pop_vectors)


    return(grid_pop_vectors_uniq)

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', type=str, help='First set of vectors in .json', required=True, default=False)
    parser.add_argument('--s', type=str, help='Second set of vectors in .json', required=False, default=False)
    parser.add_argument('--method', type=str, help='Method - grid, cog; defaults to grid', required=False, default='grid')
    parser.add_argument('--grid', type=float, help='Grid size in nm, defaults to 0.5', required=False, default=0.5)
    parser.add_argument('--threshold', type=int, help='Threshold in gridpoints for perturbation detection, defaults to 1', required=False, default=1)

    global args
    args = parser.parse_args(argv)


    try:
        with open(args.f) as file:
            vectors1 = json.load(file)
    except FileNotFoundError:
        print(f'No {args.f} found, will quit')
        exit()

    if args.s:
        try:
            with open(args.s) as file:
                vectors2 = json.load(file)
        except FileNotFoundError:
            print(f'No {args.s} found, will quit')
            exit()

    if args.method == 'grid':

        """
        # TEST VECTORS
        vectors1 = {
            'p1': [(1,0,0), (0,1,0), (1,0,0)],
            'p2': [(2.3,2.4,5.6), (1.9,2.45,5.93), (2.01,2.47,6.00)],
            'p3': [(-8.6,-2.5,3.0), (-9.78,-3.15,3.6), (-9.99,-2.5,3.7)]
        }

        vectors2 = {
            'p1': [(1, 1, 1), (1.1, 1.4, 1.5), (1.2, 1.6, 1.9)],
            'p2': [(2.3, 2.4, 5.6), (1.9, 2.45, 5.93), (2.01, 2.47, 6.00)],
            'p3': [(-8.6, -2.5, 3.0), (-9.78, -3.15, 3.6), (-11.6, -2.5, 3.7)]
        }
        """


        atomic_grid, count_dict = grid(vectors1)
        atomic_grid_2, count_dict_2 = grid(vectors2)


        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        atomic_grid_keys = list(atomic_grid.keys())
        atomic_grid_keys_2 = list(atomic_grid_2.keys())
        count_dict_keys = list(count_dict.keys())
        count_dict_keys_2 = list(count_dict.keys())
        fig = plt.figure()
        maxcols=8
        nrows = math.ceil(len(atomic_grid_keys[200:225]) / maxcols)
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
        for key, key_2 in zip(atomic_grid_keys[200:225], atomic_grid_keys_2[200:225]):

            xs = [x[0] for x in atomic_grid[key]]
            ys = [x[1] for x in atomic_grid[key]]
            zs = [x[2] for x in atomic_grid[key]]

            xs_2 = [x[0] for x in atomic_grid_2[key_2]]
            ys_2 = [x[1] for x in atomic_grid_2[key_2]]
            zs_2 = [x[2] for x in atomic_grid_2[key_2]]

            #colors = coloring(count_dict[key])
            ax = fig.add_subplot(nrows, maxcols, i+1, projection='3d')
            ax.scatter(xs=xs, ys=ys, zs=zs, label='key', s=10, c='blue', marker='s')
            ax.scatter(xs=xs_2, ys=ys_2, zs=zs_2, s=10, c='red', marker='s')
            ax.set_title(key+' '+key_2)
            ax.set(xlim=(-5,5), ylim=(-5,5), zlim=(-5,5))
            col += 1
            i += 1
        fig.tight_layout()
        plt.axis('off')
        plt.show()





if __name__ == '__main__':
    main()