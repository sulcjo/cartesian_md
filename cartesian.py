import sys, argparse
import re
import multiprocessing as mp
from alive_progress import alive_bar
import json
import sys


"""
python cartesian.py [--f file containing coordinates of atoms <cart.xvg>] [--s file containing coordinates of molecular COMs <mol.xvg>]

Cartesian requires an output file from gmx traj -f TRAJ.xtc -s TOP.tpr -ox CARTESIAN.xvg (vectors of atoms).
Optionally, it also requires an output from the same command, but with -com flag specified (vectors of COMs, beware - always 
compare comparable, that is COMs of the same part of the molecule for both systems!). Alternative to this
is just suplying a single coord.xvg for each of the trajectories, which were previously aligned using cartesian_prepare.sh.
In this case, vectors are compared in cartesian space after careful fitting of each frame using SSAP algorithm.

Please make sure that your trajectory is free of PBC and is fitted both for rotational as well as translational
movement of you molecule of interest. The molecule should also be (usually) exploring the equilibrium
part of the trajectory.

Consider omitting hydrogens out of the analysis, these are either constrained or too flexible and usually don't carry meaningful
information.
"""


def assign_to_dict(data, x_cart, y_cart, z_cart):
    times = []

    # Atoms don't start with a zero index, we'll start from the minimum value
    # They also don't neccessarily follow a linear progression of indexes
    first_xcart_key = list(x_cart.keys())[0]
    starting_atom = re.search('(?<=(atom)\s)[0-9]*', first_xcart_key)[0]
    starting_atom = int(starting_atom)
    x_keys = list(x_cart.keys())
    y_keys = list(y_cart.keys())
    z_keys = list(z_cart.keys())

    with alive_bar(len(data)) as bar:
        print('Assigning vectors to dictionaries')

        for line in data:
            ix, iy, iz = 0, 0, 0
            line_split = line.replace('\n', '').split('\t')
            times.append(line_split[0])

            curr = 'x'

            bar()
            for value in line_split[1:]:


                if curr == 'x':
                    x_cart[x_keys[ix]].append(value)
                    #print(f'Appending {value} to {x_keys[ix]} @ {line_split[0]}, ix={ix}')
                    ix += 1
                    curr = 'y'
                elif curr == 'y':
                    y_cart[y_keys[iy]].append(value)
                    #print(f'Appending {value} to {y_keys[iy]} @ {line_split[0]}, iy={iy}')
                    iy += 1
                    curr = 'z'
                elif curr == 'z':
                    z_cart[z_keys[iz]].append(value)
                    #print(f'Appending {value} to {z_keys[iz]} @ {line_split[0]}, iz={iz}')
                    iz += 1
                    curr = 'x'

    return(x_cart, y_cart, z_cart, times)

def parse_cartesian(path):
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
    comments = [line for line in lines if '@' in line]
    data = [line for line in lines if '@' not in line and '#' not in line]

    # Get atom numbers from the file using regex, create 3 dictionaries containing X Y Z values for atoms
    pattern = '[a-z]+\s+[0-9]+\s+[XYZ]*'

    """
    This doesn't work, generated dictionaries can't have one of the lists appended without modifying all the other lists,
    this is normal in Python    
    x_cart = dict.fromkeys([re.search(pattern, line)[0] for line in comments if 'atom' in line and 'X' in line], [])
    y_cart = dict.fromkeys([re.search(pattern, line)[0] for line in comments if 'atom' in line and 'Y' in line], [])
    z_cart = dict.fromkeys([re.search(pattern, line)[0] for line in comments if 'atom' in line and 'Z' in line], [])
    """

    x_keys = [re.search(pattern, line)[0] for line in comments if 'atom' in line and 'X' in line]
    y_keys = [re.search(pattern, line)[0] for line in comments if 'atom' in line and 'Y' in line]
    z_keys = [re.search(pattern, line)[0] for line in comments if 'atom' in line and 'Z' in line]
    x_cart, y_cart, z_cart = {}, {}, {}
    for x_key, y_key, z_key in zip(x_keys, y_keys, z_keys):
        x_cart[x_key] = []
        y_cart[y_key] = []
        z_cart[z_key] = []



    # Get indexes for parsing
    #pattern = '(?<=s)[0-9]*'
    #x_ind = [re.search(pattern, line)[0] for line in comments if 'atom' in line and 'X' in line]
    #y_ind = [re.search(pattern, line)[0] for line in comments if 'atom' in line and 'Y' in line]
    #z_ind = [re.search(pattern, line)[0] for line in comments if 'atom' in line and 'Z' in line]


    # Compare lengths of indexes and dictionary keys to be sure
    #if len(x_cart.keys())==len(x_ind) and len(y_cart.keys())==len(y_ind) and len(z_cart.keys())==len(z_ind):
    #    pass
    #else:
    #    print('Mismatch between indices and coordinates, will exit')
    #    exit()



    # Init multiprocessing.Pool()
    pool = mp.Pool(mp.cpu_count())
    # Assign coordinates to a dictionary
    x_cart, y_cart, z_cart, times = pool.apply(assign_to_dict, args=(data, x_cart, y_cart, z_cart))
    pool.close()

    return(x_cart, y_cart, z_cart, times)

def parse_com(path):
    """
    parse_com() parses a dataset of this type:
        # comments for .xvg reader
        @ notes about the plot for .xvg reader (names, legends etc.)
        25000	6.4015	6.40326	3.00702
        25010	6.40063	6.40136	3.00749
        25020	6.40089	6.39934	3.00821
        25030	6.40104	6.40203	3.01097
        25040	6.39872	6.40116	3.00877
        25050	6.39943	6.39993	3.00686

    The datalines (no @ or #) read as: at time 25000, the (X Y Z) of the center of mass of the molecule was (6.4015 6.40326 3.00702)


    :return: assigned dictionaries of COM coordinates x_cart, y_cart, z_cart
    """

    try:
        with open(path) as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f'{path} does not exist')
        exit()

    # Read only lines not containing '#'
    comments = [line for line in lines if '@' in line]
    data = [line for line in lines if '@' not in line and '#' not in line]
    cart_x = {'com X' : []}
    cart_y = {'com Y' : []}
    cart_z = {'com Z' : []}
    for dat in data:
        dat_split = dat.replace('\n','').split()
        cart_x['com X'].append(dat_split[1])
        cart_y['com Y'].append(dat_split[2])
        cart_z['com Z'].append(dat_split[3])

    return(cart_x, cart_y, cart_z)

def recalculate_vector(x, y, z, x_com=False, y_com=False, z_com=False):
    """
    Position vectors are by default in the Cartesian space, with (0 0 0) as the start. This means
    that two structures can't be directly compared. Imagine a multidomain system of a PDZ3-TrpCage and another
    with only PDZ3, thanks to centering of the structure in a different place in the cartesian space,
    the two vectors could not be compared. That's why we use the Center of Mass (but beware, always use the center of mass
    of the correct part of the molecule, if you're for example comparing domains!) as a reference point (0 0 0)/start
    of the atom position vector.

    :return: dictionary of lists of tuples of final vectors with starts in (COMx COMy COMz) and ends in (Atomx, atomy, atomz)
    """
    x_keys = list(x.keys())
    y_keys = list(y.keys())
    z_keys = list(z.keys())

    new_vectors = {}

    # Trigger if COM recalculation method was chosen
    if x_com and y_com and z_com:
        com_x = x_com[list(x_com.keys())[0]]
        com_y = y_com[list(y_com.keys())[0]]
        com_z = z_com[list(z_com.keys())[0]]

        with alive_bar(len(x_keys)) as bar:
            print('Recalculating vectors')
            for x_key, y_key, z_key in zip(x_keys, y_keys, z_keys):
                x_in_time = x[x_key]
                y_in_time = y[y_key]
                z_in_time = z[z_key]

                bar()
                new_vectors[x_key[:-2]] = []

                """
                Rewrite in NumPy, iterating is slow and stupid
                """
                for x_it, y_it, z_it, comx_it, comy_it, comz_it in zip(x_in_time, y_in_time, z_in_time, com_x, com_y, com_z):
                    #print('COM: ',comx_it, comy_it, comz_it, f' @ {time}', 'COORD', x_it, y_it, z_it)

                    new_x = float(x_it) - float(comx_it)
                    new_y = float(y_it) - float(comy_it)
                    new_z = float(z_it) - float(comz_it)
                    new_vectors[x_key[:-2]].append((new_x,new_y,new_z))
    else:
        with alive_bar(len(x_keys)) as bar:
            print('Assigning SSAP aligned vectors')
            for x_key, y_key, z_key in zip(x_keys, y_keys, z_keys):
                x_in_time = x[x_key]
                y_in_time = y[y_key]
                z_in_time = z[z_key]

                bar()
                new_vectors[x_key[:-2]] = []

                """
                Rewrite in NumPy, iterating is slow and stupid
                """
                for x_it, y_it, z_it in zip(x_in_time, y_in_time, z_in_time):
                    #print('COM: ',comx_it, comy_it, comz_it, f' @ {time}', 'COORD', x_it, y_it, z_it)

                    new_x = float(x_it)
                    new_y = float(y_it)
                    new_z = float(z_it)
                    new_vectors[x_key[:-2]].append((new_x,new_y,new_z))

    return(new_vectors)

def get_vectors(x_cart, y_cart, z_cart, x_com=False, y_com=False, z_com=False):
    """
    Multi-cpu caller for the recalculate_vector function.

    :param x_cart:
    :param y_cart:
    :param z_cart:
    :param x_com:
    :param y_com:
    :param z_com:
    :return: Vectors starting from COM and ending in the atomic position
    """
    pool = mp.Pool(mp.cpu_count())
    # Assign coordinates to a dictionary
    vectors = pool.apply(recalculate_vector, args=(x_cart, y_cart, z_cart, x_com, y_com, z_com))
    pool.close()
    return(vectors)

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=str, help='Input traj -ox filepath (aligned by cartesian_prepare or in conjuction with --s COM file)', required=True)
    parser.add_argument("--s", type=str, help='OPTIONAL Input traj -com -ox filepath', required=False)
    parser.add_argument("--o", type=str, help='Output.json (json)', required=True, default='cartesian_outfile')
    global args
    args = parser.parse_args(argv)
    x_cart, y_cart, z_cart, times = parse_cartesian(args.f)

    if args.s:
        print('######')
        print('(!!) Recalculating vectors with COM positions as (0,0,0)')
        print('######')
        x_com, y_com, z_com = parse_com(args.s)
        vectors = get_vectors(x_cart, y_cart, z_cart, x_com, y_com, z_com)
    else:
        print('######')
        print('(!!) Using SSAP aligned trajectories instead of COM fitting')
        print('######')
        vectors = get_vectors(x_cart, y_cart, z_cart)

    with open(args.o,'w') as fp:
        print(f'Outputting vector dictionary to {args.o}')
        json.dump(vectors, fp)



if __name__ == '__main__':
    main()






"""
exit()
# Plotting part, add to a function later
fig = plt.figure()
axs = []


max_cols = 7
#plots_total = len(vectors)
plots_total=10
ax_index = 1
axs = []



for ind, vector_key in enumerate(vectors):
    vectors_of_atom = vectors[vector_key]
    x = np.array([i[0] for i in vectors_of_atom])
    y = np.array([i[1] for i in vectors_of_atom])
    z = np.array([i[2] for i in vectors_of_atom])



    #ax = fig.add_subplot(math.ceil(plots_total/max_cols), max_cols, ax_index, projection='3d')
    #ax.set_title(vector_key)
    #ax.scatter(xs=x,ys=y,zs=z)
    #ax_index += 1
    #ax.set(xlim=(0, 10), ylim=(0, 10), zlim=(0, 10))
    #axs.append(ax)



    if ind > 10:
        break

#plt.show()

"""




