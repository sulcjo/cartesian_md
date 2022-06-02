import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import math
import re
import os
from operator import itemgetter
from scipy.optimize import curve_fit

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def get_simple_dataset(path):
    with open(path, "r") as file:
        read = file.readlines()
        array_lines = []
        info_lines = []
        axis_lines = []
        for line in read:
            if "#" not in line and '@' not in line:
                array_lines.append(line.replace("\n", ""))
            elif '#' in line:
                info_lines.append(line.replace("\n", ""))
            elif '@' in line:
                axis_lines.append(line.replace("\n", ""))
        x = []
        y = []

        for i in array_lines:
            holder = i.split()
            holder = [x for x in holder if x != '']
            try:
                x.append(float(holder[0]))
                y.append(float(holder[1]))
            except:
                continue

        return([x, y, info_lines, axis_lines])

def get_umbrella_profile(path):
    with open(path) as profile:
        profile_lines = profile.readlines()
        profile_x = []
        profile_y = []
        std_error = []
        for profile_line in profile_lines:
            if '@' not in profile_line and '#' not in profile_line:
                line_split = profile_line.split()
                profile_x.append(float(line_split[0]))
                profile_y.append(float(line_split[1]))
                if len(line_split) > 2:
                    std_error.append(float(line_split[2]))

    return([profile_x, profile_y, std_error])

def get_umbrella_histogram(path):

    with open(path) as histo_file:
            # The data form is kinda weird, X values (distances) don't correspond to runs
            # Every X value has N values in a line, where N corresponds to number of runs
            # And every run provided some counts into this X-value (weird histogram-style data)
            # So every line has N+1 values
        histo_x = []

            # Read histo.xvg line by line
        histo_lines = histo_file.readlines()
        run_y_values = None
        for i, line in enumerate(histo_lines):
                # In lines that are not comments
            if '@' not in line and '#' not in line:
                    # Split current line by whitespaces
                line_split = line.split()

                    # Now reconstruct data from specific runs (N runs, corresponding to N values
                    # in a line other than zeroth value - distance

                if not run_y_values:
                    number_runs = len(line_split[1:])
                    run_y_values = [[] for i in range(number_runs)]

                for j, value in enumerate(line_split[1:]):
                    if float(value) > 0:
                        run_y_values[j].append(float(value))
                    else:
                        run_y_values[j].append(None)
                    # Append first value (distance) to histo_x list (list of distances)
                histo_x.append(float(line_split[0]))


        return([histo_x, run_y_values])

def read_rama(path_to_file):
    with open(path_to_file,'r') as file:

        #split rama.xvg file into lines
        lines=file.readlines()
        data_lines = []

        #Iterate through lines, exclude comment lines, create array of data-carrying lines
        for line in lines:
            if '#' not in line and '@' not in line:
                data_lines.append(line.replace('\n',''))

    #Obtain names of residues into [residues], skip duplicities
        residues = {'Names' : [], 'Phis' : [], 'Psis' : [], 'Times' : []}

        for line in data_lines:
            split_line = line.split()
            residue = split_line[-1]

            already_in = 0
            for already_in_residue in residues['Names']:
                if residue == already_in_residue:
                    already_in = 1
            if already_in != 1:
                residues['Names'].append(residue)



    #Append an empty array for every residue in 'Names' into 'Phis' and 'Psis'
        for residue in residues['Names']:
            residues['Phis'].append([])
            residues['Psis'].append([])

        for line in data_lines:
            for i, residue in enumerate(residues['Names']):
                if residue == line.split()[-1]:

                    split_line = line.split()
                    residues['Phis'][i].append(float(split_line[0]))
                    residues['Psis'][i].append(float(split_line[1]))

    #Get times which are only implicitly available from rama.xvg, every 10th frame was saved, hence the 10x multiplier
        number_residues = len(residues['Names'])
        total_data_lines = len(data_lines)
        max_time = 10*(total_data_lines/number_residues)
        residues['Times'].append([int(time) for time in range(0, int(max_time), 10)])

        return(residues)

def anglea_coloring(psi_x, phi_y, phi_limits = (-90, -35), psi_limits = (-70, -15), colors = ('red', 'blue')):
    #If residue psi and phi lies in a certain Ramachandran plot area, it will be colored
    if psi_limits[0] <= psi_x <= psi_limits[1] and phi_limits[0] <= phi_y <= phi_limits[1]:
        return colors[0]
    else:
        return colors[1]

# Plot linker dihedrals
def plot_angles(dataset, type = 'Psi', ma = False, angle_limits = None, title = '', limit_resi = [0,100], limres = None):
    times = []
    for time in dataset['Times'][0]:
        times.append(time)

    if limres:
        new_dataset = {'Names': [], 'Phis': [], 'Psis': []}
        for i, dataset_name in enumerate(dataset['Names']):
            if dataset_name in limres:

                new_dataset['Names'].append(dataset['Names'][i])
                new_dataset['Phis'].append(dataset['Phis'][i])
                new_dataset['Psis'].append(dataset['Psis'][i])

        dataset = new_dataset





    fig, axs = plt.subplots(nrows=1, ncols=len(dataset['Names'][limit_resi[0] : limit_resi[1]]), figsize=(16,4))
    fig.suptitle(f'{type} angle for {title}', fontsize=24)
    plt.subplots_adjust(left=0.06,
                        bottom=0.15,
                        right=0.988,
                        top=0.80,
                        wspace=0.405,
                        hspace=0.285)
    #Limres = ['ARG-16','GLN-46', 'ARG-49', 'LYS-72']





    for i, residue in enumerate(dataset['Names'][limit_resi[0] : limit_resi[1]]):

        if type == 'Phi':
            y = dataset['Phis'][i]
        elif type == 'Psi':
            y = dataset['Psis'][i]

        print(f'{residue} len {len(y)}')

        axs[i].scatter(times, y, s = 1.4)

        if ma:
            df = pd.DataFrame(y)
            ma_dataset = df.rolling(window=100).mean()

            axs[i].plot(times, ma_dataset, color='black', linewidth = 1, alpha=0.8)



        if angle_limits:
            #Plot line for lower and upper angle limits
            axs[i].plot(times, [angle_limits[0] for time in times], color='red', linewidth=2, linestyle='--')
            axs[i].plot(times, [angle_limits[1] for time in times], color='red', linewidth=2, linestyle='--')



        axs[i].set_title(residue, fontsize=20)
        axs[i].tick_params(axis='y', labelsize=15)
        axs[i].tick_params(axis='x', labelsize=15)
        axs[i].set_xlabel('Time / ns', fontsize=16)
        axs[i].set_ylabel(f'{type}', fontsize=16)

def plot_3D_scatter(dataset, angle_coloring = False, limit_resi = [0,10], title='3D Rama'):

    #Convert ps times to ns times for improved graph clarity
    times = []
    for time in dataset['Times'][0]:
        times.append(time)

    fig = plt.figure(figsize=(8*(len(dataset['Names'][limit_resi[0] : limit_resi[1]])), 10))
    fig.suptitle(f'{title}', fontsize=24)
    for i, residue in enumerate(dataset['Names'][limit_resi[0] : limit_resi[1]]):



        ax = fig.add_subplot(1,len(dataset['Names'][limit_resi[0] : limit_resi[1]]), i+1, projection='3d')

        z = np.array(times)
        x = (np.array(dataset['Phis'][i]))
        y = (np.array(dataset['Psis'][i]))

        if angle_coloring:
            color_list = [anglea_coloring(phi_x, psi_y) for phi_x, psi_y in zip(x, y)]
            ax.scatter(x, y, z, color=color_list, alpha=0.5, marker='x')
        else:
            ax.scatter(x, y, z, alpha=0.5)

        ax.set_title(residue, fontsize=20)
        ax.tick_params(axis='y', labelsize=15)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='z', labelsize=12)
        ax.set_zlabel('Time / ns', fontsize=16, labelpad=15)
        ax.set_xlabel('Phi', fontsize=16, labelpad=15)
        ax.set_ylabel('Psi', fontsize=16, labelpad=15)
        ax.view_init(45, 270)
        ax.set_xlim((-180, 180))
        ax.set_ylim((-180, 180))
    plt.subplots_adjust(left=0.024,
                        bottom=0.050,
                        right=0.900,
                        top=0.90,
                        wspace=0.098,
                        hspace=0.285)

def plot_ramachandran(dataset, title = '', limit_resi = [0,100], limres = None, time_coloring = False):
    times = []
    for time in dataset['Times'][0]:
        times.append(time)

    if limres:
        new_dataset = {'Names': [], 'Phis': [], 'Psis': []}
        for i, dataset_name in enumerate(dataset['Names']):
            if dataset_name in limres:

                new_dataset['Names'].append(dataset['Names'][i])
                new_dataset['Phis'].append(dataset['Phis'][i])
                new_dataset['Psis'].append(dataset['Psis'][i])

        dataset = new_dataset

    import math

    how_many_plots_pre_row = 3
    how_many_resis = len(dataset['Names'][limit_resi[0] : limit_resi[1]])
    rows = int(math.ceil(how_many_resis)/how_many_plots_pre_row)
    cols = int(math.ceil(how_many_resis)/rows)
    print(rows)
    print(cols)
    fig, axs = plt.subplots(nrows=rows, ncols=cols+1, figsize=(16,4))
    fig.suptitle(f'Dihedrals for {title}', fontsize=24)
    plt.subplots_adjust(left=0.06,
                        bottom=0.15,
                        right=0.988,
                        top=0.80,
                        wspace=0.405,
                        hspace=0.285)
    #Limres = ['ARG-16','GLN-46', 'ARG-49', 'LYS-72']

    col = 0
    row = 0
    for i, residue in enumerate(dataset['Names'][limit_resi[0] : limit_resi[1]]):

        if col > cols:
            col = 0
            row += 1



        x = dataset['Phis'][i]
        y = dataset['Psis'][i]

        print(f'{residue} len {len(y)}')
        print('Current col: ' + str(col))
        print('Current rows: ' + str(row))
        data = pd.DataFrame([x, y, times])
        data = data.T
        data.columns = ('x', 'y', 'times')

        if time_coloring:


            norm = plt.Normalize(data['times'].min(), data['times'].max())
            sm = plt.cm.ScalarMappable(cmap="RdBu", norm=norm)
            sm.set_array([])

            if rows == 1:
                sns.scatterplot(data=data, x='x', y='y', size=3, hue='times', palette="RdBu", ax=axs[i], legend=False)
                fig.colorbar(sm, ax=axs[i], label='ns')
            else:
                sns.scatterplot(data=data, x='x', y='y', size=3, hue='times', palette="RdBu", ax=axs[row][col], legend=False)
                fig.colorbar(sm, ax=axs[row][col], label='ns')



        else:

            if rows == 1:
                sns.scatterplot(data=data, x='x', y='y', s=5, ax=axs[i])
            else:
                sns.scatterplot(data=data, x='x', y='y', s=5, ax=axs[row][col])

        if rows == 1:
            axs[i].set(xlim=(-180,180),ylim=(-180,180))
            axs[i].set_title(residue, fontsize=20)
            axs[i].tick_params(axis='y', labelsize=15)
            axs[i].tick_params(axis='x', labelsize=15)
            axs[i].set_xlabel('Phi', fontsize=16)
            axs[i].set_ylabel('Psi', fontsize=16)
        else:
            axs[row][col].set(xlim=(-180, 180), ylim=(-180, 180))
            axs[row][col].set_title(residue, fontsize=20)
            axs[row][col].tick_params(axis='y', labelsize=15)
            axs[row][col].tick_params(axis='x', labelsize=15)
            axs[row][col].set_xlabel('Phi', fontsize=16)
            axs[row][col].set_ylabel('Psi', fontsize=16)
        col += 1

def plot_angles_chi(dataset_array, title='', plot='scatter', corr=False, rama=False):

    """
    Included dataset array looks like this:
    [
        [ 'chi1ARG53' ],
        [ 'time1, time2...'],
        [ 'value1, value2 ...'],

        [ 'chi2ARG53' ],
        [ 'time1, time2...'],
        [ 'value1, value2 ...'],
        ...
    ]
    """

    # Array of datasets is provided, this checks how many individual datasets are present

    num_plots = len(dataset_array)


    print(num_plots)
    if num_plots > 1:
        # Maximum number of chi angles in any canonical aminoacid is 6, plus 2 for psi and chi
        nrows = math.ceil(num_plots / 8)
        fig, axs = plt.subplots(ncols=8, nrows=nrows)
        for ax in fig.axes:
            ax.set_axis_off()
        # Now iterate through the dataset array and see if the current iteration is the same aminoacid as before
        # In that case keep it in the current row, else (different AA) skip to next row
        # v this splits the aminoacid name chain after the first CAPITAL letter (file names such as phiARG44, chi1GLY65)
        last_aminoacid = re.findall('[A-Z]{3}[1-9]*', dataset_array[0][0])
        row = 0
        col = 0
        for dataset in dataset_array:
            if rama:
                dataset = [dataset[0], dataset[1], dataset[2]]


            current_aminoacid = re.findall('[A-Z]{3}[1-9]*', dataset[0])

            if  current_aminoacid != last_aminoacid and nrows > 1:
                row += 1
                col = 0
            last_aminoacid = re.findall('[A-Z]{3}[1-9]*', dataset[0])

            if nrows > 1:
                print(f'populating row:{row} col:{col} with {current_aminoacid}')
                axs[row][col].set_axis_on()
                if plot == 'scatter':

                    axs[row][col].scatter(dataset[1], dataset[2], s=1)
                elif plot == 'line':
                    axs[row][col].plot(dataset[1], dataset[2])
                axs[row][col].set_title(dataset[0])
                # A quick and dirty fix to see if it's plotting correlation or angle data (corr doesn't go over 1)
                if not corr:
                    axs[row][col].set(ylim=(-180, 180))
                    if col == 0:
                        axs[row][col].set_xlabel('Time / ns')
                    if row == nrows-1:
                        axs[row][col].set_ylabel('Angle / deg')
                    if rama:
                        axs[row][col].set(xlim=(-180, 180))
                        if col == 0:
                            axs[row][col].set_xlabel('Chi1')
                        if row == nrows - 1:
                            axs[row][col].set_ylabel('Chi2')
                else:
                    axs[row][col].set(ylim=(-1,1))
                    if col == 0:
                        axs[row][col].set_ylabel('C(t)')
                    if row == nrows-1:
                        axs[row][col].set_xlabel('Time / ns')
            elif nrows == 1:
                print(f'populating col:{col} with {current_aminoacid}')
                axs[col].set_axis_on()
                if plot == 'scatter':
                    axs[col].scatter(dataset[1], dataset[2], s=1)
                elif plot == 'line':
                    axs[col].plot(dataset[1], dataset[2])

                axs[col].set_title(dataset[0])
                if not corr:
                    axs[col].set(ylim=(-180, 180))
                    axs[col].set_xlabel('Time / ns')
                    axs[col].set_ylabel('Angle / deg')
                    if rama:
                        axs[col].set(xlim=(-180,180))
                        axs[col].set_xlabel('Chi1')
                        axs[col].set_ylabel('Chi2')
                else:
                    axs[col].set(ylim=(-1,1))
                    axs[col].set_ylabel('C(t)')
                    axs[col].set_xlabel('Time / ns')
            col += 1
        plt.suptitle(title)

def obtain_dataset_array(path):

    dataset_array = []
    files = []


    for file in os.listdir(path):
        if file.endswith('.xvg'):
            files.append(file)
        # This sorts the lists inside dataset_array according to the first value (angle type + residue name) in each dataset

    # This sorts all the lists, main criterion is the aminoacid identifier, secondary is angle-type
    def sort_list(files_list):

        # First sort the files list according to aminoacids, yields list like "[angleAA1, angleAA1, angleAA1, angleAA2, angleAA2, angleAA2]" where angles are not sorted
        def get_sorting_aminoacid(file_from_list):
            return(re.findall('[A-Z]{3}[1-9]*', file_from_list))

        sort_list = sorted(files_list, key=get_sorting_aminoacid)
        # Now split the list into lists made up of files for individual aminoacids (angles still not sorted)
        last_AA_identifier = re.findall('[A-Z]{3}[1-9]*', sort_list[0])

        list_of_sorted_lists = [[]]
        i = 0
        for sort_value in sort_list:
            current_AA_identifier = re.findall('[A-Z]{3}[1-9]*', sort_value)

            # Split the list into several lists, each containing file names for only one aminoacid
            if last_AA_identifier == current_AA_identifier:
                list_of_sorted_lists[i].append(sort_value)
                print(f'Appending {current_AA_identifier} to list index {i}')
            else:
                list_of_sorted_lists.append([])
                i += 1
                list_of_sorted_lists[i].append(sort_value)
            last_AA_identifier = current_AA_identifier

        master_list = []
        for AA_files_list in list_of_sorted_lists:
            master_list.append(sorted(AA_files_list))

        return(master_list)





    files = sort_list(files)

    # Flatten the list of lists into a list
    files = [item for sublist in files for item in sublist]

    for file in files:

        # Simply read the dataset and append stuff like "chi1ASP17" to the beggining of the dataset list
        dataset = get_simple_dataset(path+file)
        dataset.insert(0, file.replace('.xvg',''))

        dataset_array.append(dataset)





    return(dataset_array)

def get_rdf_time(path, filename_contain='cn_rdf.xvg', index = 2.5):
    """
    This function is used to plot RDF / specified time index to see the solvatation dynamics
    Index is customizable
    Sigmoid fitting also available (upper plateau - solvation max at time x)

    :param path:
    :param filename_contain:
    :param index:
    :return:
    """

    files = []
    for file in os.listdir(path):
        if filename_contain in file:
            files.append(file)
    dataset_array = []

    def get_sorting_numbers(file_from_list):
        list_re = re.findall('\d+', file_from_list)
        return(list(map(int, list_re)))



    files = sorted(files, key=get_sorting_numbers)


    for file in files:
        dataset_array.append(get_simple_dataset(path+file))


    fig, axs = plt.subplots()
    dataset_at_index = [[], []]
    for time, dataset in enumerate(dataset_array):
        for x, y in zip(dataset[0], dataset[1]):
            if x == index:
                dataset_at_index[0].append(time*10)
                dataset_at_index[1].append(y)
    return(dataset_at_index)

def get_gmxmmpbsa_dataset(path):

    with open(path) as datafile:
        lines = datafile.readlines()

        # Split into parts (PB + GB)
        reading_pb = False
        reading_gb = False
        reading_diff = False
        pb_total, pb_total_std, gb_total, gb_total_std = None,None,None,None

        for line in lines:
            if 'GENERALIZED BORN' in line:
                reading_gb = True
            elif 'POISSON BOLTZMANN' in line:
                reading_pb = True
            elif 'Delta' in line:
                reading_diff = True

            if reading_pb and reading_diff and 'TOTAL' in line:
                pb_total = line.split()[2]
                pb_total_std = line.split()[3]
                reading_pb = False
                reading_diff = False
            elif reading_gb and reading_diff and 'TOTAL' in line:
                gb_total = line.split()[1]
                gb_total_std = line.split()[2]
                reading_gb = False
                reading_diff = False


    return ({'binding': gb_total, 'bindinge': gb_total_std})

    #self.datasets[f'{dataset_name}_pb'] = {'binding': pb_total, 'bindinge': pb_total_std}
    #self.datasets[f'{dataset_name}_gb'] = {'binding': gb_total, 'bindinge': gb_total_std}
def plot_mmpbsa(dataset_array, dataset = 'gmxmmpbsa_gb', modelIndexes = None, title = ''):

    def add_value_labels(ax, spacing=5):
        """Add labels to the end of each bar in a bar chart.

        Arguments:
            ax (matplotlib.axes.Axes): The matplotlib object containing the axes
                of the plot to annotate.
            spacing (int): The distance between the labels and the bars.
        """

        # For each bar: Place a label
        for rect, variant, error in zip(ax.patches, variants_plotted, binding_energies_errors):
            # Get X and Y placement of label from rect.
            y_value = rect.get_height()
            x_value = rect.get_x() + rect.get_width() / 2

            # Number of points between bar and label. Change to your liking.
            space = 15
            # Vertical alignment for positive values
            va = 'bottom'

            # If value of bar is negative: Place label below bar
            if y_value < 0:
                # Invert space to place label below
                space *= -1
                # Vertically align label at top
                va = 'top'

            # Use Y value as label and format number with one decimal place
            label = "{:.1f}".format(y_value)
            error_label = "{:.1f}".format(error)

            # Create annotation
            ax.annotate(
                str(variant) + " : \n" + label + " +/- \n" + str(error_label),  # Use `label` as label
                (x_value, y_value),  # Place label at end of the bar
                fontsize=18,
                xytext=(0, space),  # Vertically shift label by `space`
                textcoords="offset points",  # Interpret `xytext` as offset in points
                ha='center',  # Horizontally center label
                va=va)  # Vertically align label differently for
            # positive and negative values.

        # Call the function above. All the magic happens there.

    binding_energies = []
    binding_energies_errors = []

    variants_plotted = []
    clrs = []

    for index in modelIndexes:

        # Types: gmmpbsa, gmxmmpbsa_pb, gmxmmpbsa_gb

        variant_dataset = dataset_array[index]
        binding_energy = float(variant_dataset['binding'])
        binding_error = float(variant_dataset['bindinge'])

        binding_energies.append(binding_energy)
        binding_energies_errors.append(binding_error)
        variants_plotted.append(index)

        if binding_energy + binding_error < 0 and binding_energy - binding_error < 0:
            clrs.append("g")
        elif binding_energy + binding_error > 0 and binding_energy - binding_error > 0:
            clrs.append('r')
        else:
            clrs.append('y')


    fig, ax = plt.subplots(figsize=(20,15))
    ax.bar(variants_plotted, binding_energies, color=clrs, yerr=binding_energies_errors, alpha=0.8, capsize=4,
            width=0.5, align='edge')
    plt.xticks(variants_plotted)

    ax.set_ylabel(r"$\Delta$G  / kcal/mol ", size=16)
    ax.set_title(f"{title}", size=18)
    ax.yaxis.grid(True)
    y_ticks = ax.get_yticks()
    ax.set_yticklabels(labels=y_ticks, fontsize=16)
    ax.set_xticklabels([])
    # plt.tight_layout()
    add_value_labels(ax)
    legend_entries = {'> 0': 'r', '?': 'y', '< 0': 'g'}
    labels = list(legend_entries.keys())
    handles = [plt.Rectangle((0, 0), 1, 1, color=legend_entries[label]) for label in labels]
    plt.legend(handles, labels, fontsize=16, loc='upper right')

"""
mmgbsa_arrays = [{'bla':'ble'}]
base_path = '/run/timeshift/backup/IOCB/docking/FDs/jama_12_md/fd3'
for i in range(1,11):
    mmgbsa_arrays.append(get_gmxmmpbsa_dataset(f'{base_path}/model_{i}/FINAL_RESULTS_MMPBSA.dat'))
print(mmgbsa_arrays)

plot_mmpbsa(mmgbsa_arrays, modelIndexes=[1,2,3,4,5,6,7,8,9,10], title='FD3A+JAMA-12\n1999 frames MM/GB(8)SA')
plt.show()

"""
#base_path = '/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_i/chi/pyplot/'
#base_path_corr = '/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_i/chi/corrs/pyplot/'
#base_path_ramas = '/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_i/chi/ramas/pyplot/'
#angles_dataset = obtain_dataset_array(base_path)
#correlation_dataset = obtain_dataset_array(base_path_corr)
#ramas_dataset = obtain_dataset_array(base_path_ramas)


#plot_angles_chi(angles_dataset, title='TrpCage-GGGGGG-PDZ3 closed start')
#plot_angles_chi(correlation_dataset, title='TrpCage-GGGGGG-PDZ3 closed start self-correlation', plot='line', corr=True)
#plot_angles_chi(ramas_dataset,title='TrpCage-GGGGGG-PDZ3 closed start RAMA', plot='scatter', rama=True)

"""
# Load RMSD
#rsmd = get_simple_dataset('/home/sulcjo/IOCB/md/trp_gggggg_pdz_closed_400ns/xmgrace/')
#rg = get_simple_dataset('/home/sulcjo/IOCB/md/trp_gggggg_pdz_closed_400ns/xmgrace')
#domain_distance = get_simple_dataset('/home/sulcjo/IOCB/md/trp_gggggg_pdz_closed_400ns/xmgrace')
"""
"""
# Load umbrella contacts, histograms and pmf curve
us_contacts = get_simple_dataset('/run/timeshift/backup/IOCB/md/FDs/umbrella/centroids_comparison/fd4a_2/wham_results/contacts_pulling.xvg')
pmf_curve = get_umbrella_profile('/run/timeshift/backup/IOCB/md/FDs/umbrella/centroids_comparison/fd4a_2/wham_results/profile_errors.xvg')
pmf_histograms = get_umbrella_histogram('/run/timeshift/backup/IOCB/md/FDs/umbrella/centroids_comparison/fd4a_2/wham_results/histo.xvg')
#solvation_curve_25 = get_rdf_time('/run/media/sulcjo/sulcjo-data/IOCB/md/umbrella/cluster5/cluster5/wham_results/', index=2.0)

## Plot US
convergence_lowlim=2.6
convergence_highlim=4.5
converged_vals = [val for i,val in enumerate(pmf_curve[1]) if convergence_highlim < float(pmf_curve[0][i]) > convergence_lowlim ]
print(np.mean(converged_vals))

fig_us, axs_us = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios' : [3,2,2]})
setLineColor = 'blue'
setFontSizeLarge = 18
setFontSizeMedium = 14

axs_us[0].errorbar(pmf_curve[0], pmf_curve[1], yerr=pmf_curve[2], capsize=5, ecolor='black',color=setLineColor)
axs_us[0].set_title('PMF Curve for FD4A representative cut II', fontsize=setFontSizeLarge)
axs_us[0].set_ylabel('PMF / kJ/mol', fontsize=setFontSizeMedium)


#axs_us[1].scatter(solvation_curve_25[0],solvation_curve_25[1])
#axs_us[1].set_title('SOL molecules around com of TrpCage vs time @ r=2.0 nm')
#axs_us[1].set_ylabel('SOL molecules')
#rdf_time = pd.DataFrame(solvation_curve_25[1]).rolling(30).mean()
#axs_us[1].plot(solvation_curve_25[0],rdf_time)

for run in pmf_histograms[1]:
    axs_us[1].plot(pmf_histograms[0], run, color=setLineColor)
axs_us[1].set_title('Histograms', fontsize=setFontSizeLarge)
axs_us[1].set_ylabel('Counts', fontsize=setFontSizeMedium)
axs_us[1].set_xlabel('COM-COM distance / nm', fontsize=setFontSizeMedium)

axs_us[2].plot(us_contacts[0], us_contacts[1], color=setLineColor)
axs_us[2].set_title('< 0.6 nm contacts', fontsize=setFontSizeLarge)
axs_us[2].set_xlabel('Time / ns', fontsize=setFontSizeMedium)
axs_us[2].set_ylabel('Count', fontsize=setFontSizeMedium)
"""
"""
"""

#plot_3D_scatter(rama, angle_coloring=True, title='Linker (TrpCage-GGGGGG-PDZ3)', limit_resi=(24,30))
#plot_angles(rama, 'Phi', ma=True, angle_limits = (-90, -35), title='US Linker (TrpCage-GGGGGG-PDZ3)', limit_resi=(24,30))
#plot_angles(rama, 'Psi', ma=True, angle_limits = (-70, -15), title='US Linker (TrpCage-GGGGGG-PDZ3)', limit_resi=(24,30))
"""
fig_md_closed, axs_md_closed = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios' : [3,1,1]})
rama_md = read_rama('/run/timeshift/backup/IOCB/md/trp_pdz_open_amber/rama.xvg')
plot_angles(rama_md, 'Phi', ma=True, angle_limits = None, title='MD Linker (TrpCage-GGGGGG-PDZ3 open start AMBER)', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'])
plot_angles(rama_md, 'Psi', ma=True, angle_limits = None, title='MD Linker (TrpCage-GGGGGG-PDZ3 open start AMBER)', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'])
plot_ramachandran(rama_md, title='MD Linker (TrpCage-GGGGGG-PDZ3 open start AMBER)', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'], time_coloring=True)
"""
#fig_md_closed, axs_md_closed = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios' : [3,1,1]})
#rama_md = read_rama('/run/media/sulcjo/sulcjo-data/IOCB/md/pdz_gggggg_trp_150ns/rama.xvg')
#plot_angles(rama_md, 'Phi', ma=True, angle_limits = None, title='MD Linker (PDZ3-GGGGGG-TrpCage)', limit_resi=(0,100), limres = ['GLY-83' ,'GLY-84', 'GLY-85', 'GLY-86', 'GLY-87', 'GLY-88'])
#plot_angles(rama_md, 'Psi', ma=True, angle_limits = None, title='MD Linker (PDZ3-GGGGGG-TrpCage)', limit_resi=(0,100), limres = ['GLY-83' ,'GLY-84', 'GLY-85', 'GLY-86', 'GLY-87', 'GLY-88'])
#plot_ramachandran(rama_md, title='MD Linker (PDZ3-GGGGGG-TrpCage)', limit_resi=(0,100), limres = ['GLY-83' ,'GLY-84', 'GLY-85', 'GLY-86', 'GLY-87', 'GLY-88'], time_coloring=True)

"""
fig_md_closed, axs_md_closed = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios' : [3,1,1]})
rama_md = read_rama('/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_ii/rama.xvg')
plot_angles(rama_md, 'Phi', ma=True, angle_limits = None, title='MD Linker (TrpCage-GGGGGG-PDZ3 closed start II) 400 ns', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'])
plot_angles(rama_md, 'Psi', ma=True, angle_limits = None, title='MD Linker (TrpCage-GGGGGG-PDZ3 closed start II) 400 ns', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'])
plot_ramachandran(rama_md, title='MD Linker (TrpCage-GGGGGG-PDZ3 closed start II) 400 ns', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'], time_coloring=True)
"""
"""
fig_us, axs_us = plt.subplots(nrows=3, ncols=1)

domain_distance = get_simple_dataset('/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_ii/xmgrace/distance_domains.xvg')
#domain_distance_open = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/xmgrace/distance_domains.xvg')
axs_us[0].plot(domain_distance[0], domain_distance[1])
#axs_us[0].plot(domain_distance_open[0], domain_distance_open[1], color='red',label='Open Start')
axs_us[0].set_title('TrpCage-GGGGGG-PDZ3 COM distances (closed start II)')
#axs_us[0].set_title('TrpCage-GGGGGG-PDZ3 COM distances')
axs_us[0].set_ylabel('d / nm')
axs_us[0].legend()

domain_contacts = get_simple_dataset('/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_ii/xmgrace/domain_contacts.xvg')
#domain_contacts_open = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/xmgrace/domain_contacts.xvg')
axs_us[1].plot(domain_contacts[0], domain_contacts[1])
#axs_us[1].plot(domain_contacts_open[0], domain_contacts_open[1],color='red')
axs_us[1].set_title('TrpCage-GGGGGG-PDZ3 COM contacts (closed start II)')
axs_us[1].set_ylabel('count')

domain_histograms = get_simple_dataset('/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_ii/xmgrace/domains_distance_histo.xvg')
#domain_histograms_open = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/xmgrace/domains_distance_histo.xvg')
axs_us[2].plot(domain_histograms[0], domain_histograms[1])
#axs_us[2].plot(domain_histograms_open[0], domain_histograms_open[1],color='red')
axs_us[2].set_title('TrpCage-GGGGGG-PDZ3 COM distances frequency (closed start II)')
axs_us[2].set_ylabel('Count')

#domain_distances_autocorr = get_simple_dataset('/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_i/xmgrace/autocorr_domain_distances.xvg')
#axs_us[3].plot(domain_distances_autocorr[0], domain_distances_autocorr[1])
#axs_us[3].set_title('PDZ3-TrpCage (Closed start) COM distances autocorrelation')
#axs_us[3].set_ylabel('G(t)')
#axs_us[3].set_xlabel('t / ps')
#axs_us[3].set(ylim=(-1,1))
"""
"""
rmsf_pdztrp = get_simple_dataset('/home/sulcjo/IOCB/md/pdz_gggggg_trp_150ns/xmgrace/rmsf_all_atom.xvg')
rmsf_pdz = get_simple_dataset('/home/sulcjo/IOCB/md/pdz/xmgrace/rmsf_all_atom.xvg')
def delta_rmsf(dataset1, dataset2):
    delta_rmsf = []
    for d1, d2 in zip(dataset1[1],dataset2[1]):

        try:
            delta_rmsf.append(abs(d1-d2))
        except:
            delta_rmsf.append(0)
    return(delta_rmsf)

delta_rmsf = delta_rmsf(rmsf_pdz, rmsf_pdztrp)
fig_rmsf, axs_rmsf = plt.subplots(nrows=2,ncols=1, sharex=True)
axs_rmsf[0].plot(rmsf_pdztrp[0], rmsf_pdztrp[1], label='PDZ3-GGGGGG-TrpCage RMSF')
axs_rmsf[0].set_ylabel('nm', fontsize=setFontSizeMedium)
axs_rmsf[0].set_xlabel('atom', fontsize=setFontSizeMedium)
axs_rmsf[0].plot(rmsf_pdz[0], rmsf_pdz[1], label='PDZ3 RMSF')
axs_rmsf[1].plot(rmsf_pdz[0], delta_rmsf, color='red', label='delta')
axs_rmsf[0].legend()
axs_rmsf[1].legend()

fig_us, axs_us = plt.subplots(nrows=3, ncols=1)
setLineColor = 'blue'
setFontSizeLarge = 18
setFontSizeMedium = 14
domain_distance = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_closed_amber/xmgrace/distance_domains.xvg')
axs_us[0].plot(domain_distance[0], domain_distance[1],color=setLineColor, label='Closed start')
axs_us[0].set_title('TrpCage...PDZ3 no linker COM distances (closed start)', fontsize=setFontSizeLarge)
axs_us[0].set_title('TrpCage-GGGGGG-PDZ3 closed start AMBER COM distances', fontsize=setFontSizeLarge)
axs_us[0].set_ylabel('d / nm', fontsize=setFontSizeMedium)
axs_us[0].legend()

domain_contacts = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_closed_amber/xmgrace/domain_contacts.xvg')
axs_us[1].plot(domain_contacts[0], domain_contacts[1],color=setLineColor)
axs_us[1].set_title('TrpCage-GGGGGG-PDZ3 closed start AMBER COM contacts', fontsize=setFontSizeLarge)
axs_us[1].set_ylabel('count', fontsize=setFontSizeMedium)

domain_histograms = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_closed_amber/xmgrace/domains_distance_histo.xvg')
axs_us[2].plot(domain_histograms[0], domain_histograms[1],color=setLineColor)
axs_us[2].set_title('TrpCage-GGGGGG-PDZ3 closed start AMBER COM distances frequency', fontsize=setFontSizeLarge)
axs_us[2].set_ylabel('Count', fontsize=setFontSizeMedium)

us_contacts = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/umbrella/linker_amber/wham_results/contacts_pulling.xvg')
pmf_curve = get_umbrella_profile('/run/media/sulcjo/sulcjo-data/IOCB/umbrella/linker_amber/wham_results/profile_errors.xvg')
pmf_histograms = get_umbrella_histogram('/run/media/sulcjo/sulcjo-data/IOCB/umbrella/linker_amber/wham_results/histo.xvg')
fig_us, axs_us = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios' : [3,1,1]})
setLineColor = 'blue'
setFontSizeLarge = 18
setFontSizeMedium = 14

axs_us[0].errorbar(pmf_curve[0], pmf_curve[1], yerr=pmf_curve[2], capsize=5, ecolor='black',color=setLineColor)
axs_us[0].set_title('PMF Curve for TrpCage-GGGGGG-PDZ3 AMBER pulling', fontsize=setFontSizeLarge)
axs_us[0].set_ylabel('PMF / kJ/mol', fontsize=setFontSizeMedium)

for run in pmf_histograms[1]:
    axs_us[1].plot(pmf_histograms[0], run, color=setLineColor)
axs_us[1].set_title('Histograms', fontsize=setFontSizeLarge)
axs_us[1].set_ylabel('Counts', fontsize=setFontSizeMedium)
axs_us[2].set_xlabel('COM-COM distance / nm', fontsize=setFontSizeMedium)

axs_us[2].plot(us_contacts[0], us_contacts[1], color=setLineColor)
axs_us[2].set_title('< 0.6 nm contacts', fontsize=setFontSizeLarge)
axs_us[2].set_xlabel('Time / ns', fontsize=setFontSizeMedium)
axs_us[2].set_ylabel('Count', fontsize=setFontSizeMedium)
"""

"""
fig, ax = plt.subplots(ncols=2)
pca_linker = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/umbrella/linker/pca_pdz/2d_proj.xvg')
pca_nolinker = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/umbrella/nolinker/pca_pdz3/2d_proj.xvg')
pca_linker_amber = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/umbrella/linker_amber/pca/2d_proj.xvg')
pca_closed_simulations = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_closed_400ns/pca_pdz3/2d_proj.xvg')
pca_closed_simulation_amber = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_closed_amber/pca_pdz3/2d_proj.xvg')

ax[0].scatter(pca_linker[0], pca_linker[1], label = 'Harm. pot. with linker CHARMM')
ax[0].scatter(pca_nolinker[0], pca_nolinker[1], label = 'Harm. pot. no linker')
ax[0].scatter(pca_linker_amber[0], pca_linker_amber[1], label = 'Harm. pot. with linker AMBER')
ax[1].scatter(pca_closed_simulations[0], pca_closed_simulations[1], label = 'STD MD with linker CHARMM')
ax[1].scatter(pca_closed_simulation_amber[0], pca_closed_simulation_amber[1], label = 'STD MD with linker AMBER')
ax[0].legend()
ax[1].legend()
"""
"""
d12 = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/md_d12.xvg')
d13 = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/md_d13.xvg')
d12_histo = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/md_d12_histo.xvg')
d13_histo = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/md_d13_histo.xvg')

d12_330 = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/annealing_d12.xvg')
d13_330 = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/annealing_d13.xvg')
d12_histo_330 = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/annealing_d12_histo.xvg')
d13_histo_330 = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/annealing_d13_histo.xvg')

fig_d12, axs_12 = plt.subplots(nrows=2, ncols=2, gridspec_kw={'width_ratios': [3, 1], 'height_ratios': [3, 1]}, sharex='col', sharey='row')
axs_12[0][0].scatter(d12[1], d13[1], s=100, marker='x', label='310 K')
axs_12[0][0].scatter(d12_330[1], d13_330[1], s=100, marker='x', label='330 K', alpha=0.5)
axs_12[0][0].set_xlabel('COM-COM dist. D1-2 / nm')
axs_12[0][0].set_ylabel('COM-COM dist. D1-3 / nm')
axs_12[0][0].set(xlim=(2.50, 3.50), ylim=(2.50, 3.50))
axs_12[0][0].legend()

axs_12[1][0].plot(d12_histo[0], d12_histo[1], label='310 K')
axs_12[1][0].plot(d12_histo_330[0], d12_histo_330[1], label='330 K')
axs_12[1][0].set_xlabel('probability dist. D1-2')
axs_12[1][0].set_ylabel('COM-COM dist. D1-2 / nm')
axs_12[1][0].legend()

axs_12[0][1].plot(d13_histo[1], d13_histo[0], label='310 K')
axs_12[0][1].plot(d13_histo_330[1], d13_histo_330[0], label='330 K')
axs_12[0][1].set_xlabel('probability dist. D1-3')
axs_12[0][1].set_xlabel('COM-COM dist. D1-3 / nm')
axs_12[0][1].legend()

"""
"""
base_path = '/run/timeshift/backup/IOCB/MSM/'

rmsf_dataset_pdz = []
rmsf_dataset_fd3a = []
rmsf_dataset_fd4a = []
rmsf_dataset_pdz_amber = []
rmsf_dataset_fd3a_amber = []
rmsf_dataset_fd4a_amber = []

for run in [f'run_{i}' for i in range(1,19)]:
    new_dataset_pdz = get_simple_dataset(f'{base_path}/pdz/charmm_chwater/analyses/complexrmsf_{run}.xvg')
    new_dataset_fd3a = get_simple_dataset(f'{base_path}/pdz_l_trp/analyses_pdz/pdzrmsf_{run}.xvg')
    new_dataset_fd4a = get_simple_dataset(f'{base_path}/trp_l_pdz_closed/analyses_pdz/pdzrmsf_{run}.xvg')

    rmsf_dataset_pdz.append(new_dataset_pdz[1])
    rmsf_dataset_fd3a.append(new_dataset_fd3a[1])
    rmsf_dataset_fd4a.append(new_dataset_fd4a[1])
rmsf_dataset_pdz_df = pd.DataFrame(rmsf_dataset_pdz)
rmsf_dataset_fd3a_df = pd.DataFrame(rmsf_dataset_fd3a)
rmsf_dataset_fd4a_df = pd.DataFrame(rmsf_dataset_fd4a)


for run in [f'run_{i}' for i in range(20,29)]:
    new_dataset_pdz_amber = get_simple_dataset(f'{base_path}/pdz_amber/analyses/complexrmsf_{run}.xvg')
    new_dataset_fd3a_amber = get_simple_dataset(f'{base_path}/pdz_l_trp_amber/analyses_pdz/pdzrmsf_{run}.xvg')
    new_dataset_fd4a_amber = get_simple_dataset(f'{base_path}/trp_l_pdz_closed_amber/analyses_pdz/pdzrmsf_{run}.xvg')

    rmsf_dataset_pdz_amber.append(new_dataset_pdz[1])
    rmsf_dataset_fd3a_amber.append(new_dataset_fd3a[1])
    rmsf_dataset_fd4a_amber.append(new_dataset_fd4a[1])
rmsf_dataset_pdz_df = pd.DataFrame(rmsf_dataset_pdz)
rmsf_dataset_fd3a_df = pd.DataFrame(rmsf_dataset_fd3a)
rmsf_dataset_fd4a_df = pd.DataFrame(rmsf_dataset_fd4a)

rmsf_dataset_pdz_amber_df = pd.DataFrame(rmsf_dataset_pdz_amber)
rmsf_dataset_fd3a_amber_df = pd.DataFrame(rmsf_dataset_fd3a_amber)
rmsf_dataset_fd4a_amber_df = pd.DataFrame(rmsf_dataset_fd4a_amber)


atoms = [i for i in range(1, len(rmsf_dataset_pdz[0]))]
fig, axs = plt.subplots(ncols=2, nrows=4, figsize=(20,15), sharex='all', gridspec_kw={'height_ratios': [4, 1, 1, 1 ]})
axs[0][0].set_title('CHARMM36m average RMSF from 20x 400 ns run')
axs[0][0].plot(rmsf_dataset_pdz_df.mean(), label='PDZ3')
axs[0][0].plot(rmsf_dataset_fd3a_df.mean(), label='FD3A')
axs[0][0].plot(rmsf_dataset_fd4a_df.mean(), label='FD4A')
axs[0][0].legend()


delta_fd3a_pdz = rmsf_dataset_fd3a_df.mean() - rmsf_dataset_pdz_df.mean()
axs[1][0].set_title('d(FD3A - PDZ3) CHARMM36m')
axs[1][0].plot(delta_fd3a_pdz)

delta_fd4a_pdz = rmsf_dataset_fd4a_df.mean() - rmsf_dataset_pdz_df.mean()
axs[2][0].set_title('d(FD4A - PDZ3) CHARMM36m')
axs[2][0].plot(delta_fd4a_pdz)


delta_fd43a_pdz = rmsf_dataset_fd4a_df.mean() - rmsf_dataset_fd3a_df.mean()
axs[3][0].set_title('d(FD4A - FD3A) CHARMM36m')
axs[3][0].plot(delta_fd43a_pdz)
axs[3][0].set_xlabel('atom')





axs[0][1].set_title('SB99-ILDN average RMSF from 20x 400 ns run')
axs[0][1].plot(rmsf_dataset_pdz_amber_df.mean(), label='PDZ3')
axs[0][1].plot(rmsf_dataset_fd3a_amber_df.mean(), label='FD3A')
axs[0][1].plot(rmsf_dataset_fd4a_amber_df.mean(), label='FD4A')
axs[0][1].legend()

delta_fd3a_pdz_amber = rmsf_dataset_fd3a_amber_df.mean() - rmsf_dataset_pdz_amber_df.mean()
axs[1][1].set_title('d(FD3A - PDZ3) SB99-ILDN')
axs[1][1].plot(delta_fd3a_pdz_amber)

delta_fd4a_pdz_amber = rmsf_dataset_fd4a_amber_df.mean() - rmsf_dataset_pdz_amber_df.mean()
axs[2][1].set_title('d(FD4A - PDZ3) SB99-ILDN')
axs[2][1].plot(delta_fd4a_pdz_amber)

delta_fd43a_pdz_amber = rmsf_dataset_fd4a_amber_df.mean() - rmsf_dataset_fd3a_amber_df.mean()
axs[3][1].set_title('d(FD4A - FD3A) SB99-ILDN')
axs[3][1].set_xlabel('atom')
axs[3][1].plot(delta_fd43a_pdz_amber)

plt.show()




runs = [f'run_{i}' for i in range(20,29)]
rg_dataset_pdz = []
sasa_dataset_pdz = []
rg_dataset_pdzltrp = []
sasa_dataset_pdzltrp = []
rg_dataset_trplpdz = []
sasa_dataset_trplpdz = []

for run in runs:
    rg_dataset_pdz.append(get_simple_dataset(f'{base_path}pdz_amber/analyses/complexgyrate_{run}.xvg')[1])
    sasa_dataset_pdz.append(get_simple_dataset(f'{base_path}pdz_amber/analyses/complexsasa_{run}.xvg')[1])

    rg_dataset_pdzltrp.append(get_simple_dataset(f'{base_path}pdz_l_trp_amber/analyses/complexgyrate_{run}.xvg')[1])
    sasa_dataset_pdzltrp.append(get_simple_dataset(f'{base_path}pdz_l_trp_amber/analyses/complexsasa_{run}.xvg')[1])

    rg_dataset_trplpdz.append(get_simple_dataset(f'{base_path}trp_l_pdz_closed_amber/analyses/complexgyrate_{run}.xvg')[1])
    sasa_dataset_trplpdz.append(get_simple_dataset(f'{base_path}trp_l_pdz_closed_amber/analyses/complexsasa_{run}.xvg')[1])

fig, axs = plt.subplots(ncols=3, figsize=(20,15))

areas_pdz = []
areas_pdzltrp = []
areas_trplpdz = []


for rg_pdz, sasa_pdz, run, rg_pdzltrp, sasa_pdzltrp, rg_trplpdz, sasa_trplpdz in zip(rg_dataset_pdz, sasa_dataset_pdz, runs,
                                                                                     rg_dataset_pdzltrp, sasa_dataset_pdzltrp, rg_dataset_trplpdz, sasa_dataset_trplpdz):
    axs[0].scatter(rg_pdz, sasa_pdz, label=run, s=5)
    axs[1].scatter(rg_pdzltrp, sasa_pdzltrp, label=run, s=5)
    axs[2].scatter(rg_trplpdz, sasa_trplpdz, label=run, s=5)

    hist_pdz, xedges_pdz, yedges_pdz = np.histogram2d(np.asarray(rg_pdz), np.asarray(sasa_pdz), bins=(20, 20))
    hist_pdzltrp, xedges_pdzltrp, yedges_pdzltrp = np.histogram2d(np.asarray(rg_pdzltrp), np.asarray(sasa_pdzltrp), bins=(20, 20))
    hist_trplpdz, xedges_trplpdz, yedges_trplpdz = np.histogram2d(np.asarray(rg_trplpdz), np.asarray(sasa_trplpdz), bins=(20, 20))

    over_threshold_pdz = hist_pdz > 1
    over_threshold_pdzltrp = hist_pdzltrp > 1
    over_threshold_trplpdz = hist_trplpdz > 1

    areas_pdz.append(over_threshold_pdz.sum() * (xedges_pdz[1]-xedges_pdz[0]) * (yedges_pdz[1] - yedges_pdz[0]))
    areas_pdzltrp.append(over_threshold_pdzltrp.sum() * (xedges_pdzltrp[1]-xedges_pdzltrp[0]) * (yedges_pdzltrp[1] - yedges_pdzltrp[0]))
    areas_trplpdz.append(over_threshold_trplpdz.sum() * (xedges_trplpdz[1] - xedges_trplpdz[0]) * (
                yedges_trplpdz[1] - yedges_trplpdz[0]))

areas_pdz = list(map(lambda x: np.round(x, decimals=2), areas_pdz))
areas_pdzltrp = list(map(lambda x: np.round(x, decimals=2), areas_pdzltrp))
areas_trplpdz = list(map(lambda x: np.round(x, decimals=2), areas_trplpdz))

print(f'Areas PDZ3 (SB99-ILDN)  {areas_pdz} , AVG {np.round(np.mean(areas_pdz), decimals=3)} +- STD{np.round(np.std(areas_pdz), decimals=3)} (nm^3)\n')
print(f'Areas PDZ3-l-TrpCage (SB99-ILDN) {areas_pdzltrp} , {np.round(np.mean(areas_pdzltrp), decimals=3)} +- STD{np.round(np.std(areas_pdzltrp), decimals=3)} (nm^3)\n')
print(f'Areas TrpCage-l-PDZ3 (SB99-ILDN) {areas_trplpdz} , {np.round(np.mean(areas_trplpdz), decimals=3)} +- STD{np.round(np.std(areas_trplpdz), decimals=3)} (nm^3)\n')



#axs[0].legend(loc='upper right')
#axs[1].legend(loc='upper right')
#axs[2].legend(loc='upper right')

#for ax in axs:
    #ax.set(xlim=(1.15, 1.30), ylim=(46, 60))

axs[0].set_title('PDZ3 (SB99-ILDN)')
axs[1].set_title('FD3A (SB99-ILDN)')
axs[2].set_title('FD4A (SB99-ILDN)')

axs[0].set_xlabel('R(g) / nm')
axs[0].set_ylabel('SASA / nm^2')
axs[1].set_xlabel('R(g) / nm')
axs[1].set_ylabel('SASA / nm^2')
axs[2].set_xlabel('R(g) / nm')
axs[2].set_ylabel('SASA / nm^2')

"""
"""
base_path='/run/timeshift/backup/IOCB/docking/FDs/jama_12_md/fd3'
rmsd_datasets = []
distance_histograms = []

for model in range(1,11):
    rmsd_datasets.append(get_simple_dataset(f'{base_path}/ligand_rms_model_{model}.xvg'))
    distance_histograms.append(get_simple_dataset(f'{base_path}/histogram_complex_distance_model_{model}.xvg'))
###
fig = plt.figure(figsize=(20,15))
cols = 3
rows = math.ceil(len(rmsd_datasets)/cols)
for i, dataset in enumerate(rmsd_datasets):
    ax = plt.subplot(rows, cols, i+1)
    ax.plot(dataset[0], dataset[1])
    ax.set_title(f'FD3A+JAMA-12\nJAMA-12 backbone, model_{i+1}')
    ax.set_xlabel('Time / ps')
    ax.set_ylabel('RMSD / nm')

fig.tight_layout()
plt.subplots_adjust(hspace=0.7)

###

fig = plt.figure(figsize=(20,15))

cols = 3
rows = math.ceil(len(distance_histograms)/cols)
for i, dataset in enumerate(distance_histograms):
    ax = plt.subplot(rows, cols, i+1)
    ax.plot(dataset[0], dataset[1])
    ax.set_title(f'FD3A+JAMA-12\nCOM-COM distance, model_{i+1}')
    ax.set_xlabel('d / nm')
    ax.set_ylabel('frequency')
    ax.set(xlim=(1,4.0))
fig.tight_layout()
plt.subplots_adjust(hspace=0.7)
plt.show()
"""

def read_dat(path):
    with open(path) as file:
        lines = file.readlines()
    return [line.replace('\n','') for line in lines]

pdz3_cmd = read_dat("/run/media/sulcjo/sulcjo-data/IOCB/docking/FDs/jama_12/jama_12_gamd/pdz_conventional/xtc_trajectory/dat/com_prot_lig.dat")
pdz3_gamd = read_dat("/run/media/sulcjo/sulcjo-data/IOCB/docking/FDs/jama_12/jama_12_gamd/pdz_bound/xtc_trajectory/dat/com_prot_lig.dat")
pdz3_pepgamd = read_dat("/run/media/sulcjo/sulcjo-data/IOCB/docking/FDs/jama_12/jama_12_gamd/pdz_pepgamd/xtc_trajectory/dat/com_prot_lig.dat")
pdz3_hrexgamd = read_dat("/run/media/raid/IOCB_archive/docking/FDs/jama_12/jama_12_gamd/hrexgamd_pdz_jama12_model10/xtc_trajectory/dat/com_prot_lig.dat")
fd3_gamd = read_dat("/run/media/sulcjo/sulcjo-data/IOCB/docking/FDs/jama_12/jama_12_gamd/fd3_bound/xtc_trajectory/dat/com_prot_lig.dat")
fd4_gamd = read_dat("/run/media/sulcjo/sulcjo-data/IOCB/docking/FDs/jama_12/jama_12_gamd/fd4_bound/xtc_trajectory/dat/com_prot_lig.dat")
labels=["PDZ3+JAMA12 MD", "PDZ3+JAMA12 GaMD", "PDZ3+JAMA12 pep-GaMD", "PDZ3+JAMA12 HREX-GaMD", "FD3A+JAMA12 GaMD", "FD4A+JAMA12 GaMD"]
datasets=[pdz3_cmd, pdz3_gamd, pdz3_pepgamd, pdz3_hrexgamd, fd3_gamd, fd4_gamd]
fig, axs = plt.subplots(nrows=2, ncols=3)
col=-1
row=0
for data, label in zip(datasets, labels):
    col += 1
    if col > 2:
        col = -1
        row += 1

    axs[row][col].plot([x for x in range(0, len(data))], [float(y) for y in data], label=label)
    axs[row][col].legend()



"""
base_path = '/run/timeshift/backup/IOCB/md/FDs/'
open_i_rg = get_simple_dataset(f'{base_path}trp_gggggg_pdz_open_i/xmgrace/gyrate.xvg')
open_i_sasa = get_simple_dataset(f'{base_path}trp_gggggg_pdz_open_i/xmgrace/sasa.xvg')
open_ii_rg = get_simple_dataset(f'{base_path}trp_gggggg_pdz_open_ii/trp_pdz_open_ii/xmgrace/gyrate.xvg')
open_ii_sasa = get_simple_dataset(f'{base_path}trp_gggggg_pdz_open_ii/trp_pdz_open_ii/xmgrace/sasa.xvg')
open_amber_rg = get_simple_dataset(f'{base_path}trp_pdz_open_amber/xmgrace/gyrate.xvg')
open_amber_sasa = get_simple_dataset(f'{base_path}trp_pdz_open_amber/xmgrace/sasa.xvg')

closed_i_rg = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_i/xmgrace/gyrate.xvg')
closed_i_sasa = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_i/xmgrace/sasa.xvg')
closed_ii_rg = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_ii/xmgrace/gyrate.xvg')
closed_ii_sasa = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_ii/xmgrace/sasa.xvg')
closed_amber_rg = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_amber/xmgrace/gyrate.xvg')
closed_amber_sasa = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_amber/xmgrace/sasa.xvg')

remd_rg = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/remd_aur/remd/cold_replica/gyrate.xvg')
remd_sasa = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/remd_aur/remd/cold_replica/sasa.xvg')

#rg_AMBER = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/amd/trp_l_pdz_closed_i/rg.rms')
#sasa_AMBER = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/amd/trp_l_pdz_closed_i/sasa.rms')
#sasa_AMBER[1] = list(map(lambda x: x/100, sasa_AMBER[1]))
#rg_AMBER[1] = list(map(lambda x: x/10, rg_AMBER[1]))



fig, axs = plt.subplots()
whole_dataset = [ [open_i_rg[1], open_i_sasa[1]], [open_ii_rg[1], open_ii_sasa[1]],
                 [open_amber_rg[1], open_amber_sasa[1]], [closed_i_rg[1], closed_i_sasa[1]],
                 [closed_ii_rg[1], closed_ii_sasa[1]], [closed_amber_rg[1], closed_amber_sasa[1]],
                  [remd_rg[1], remd_sasa[1]]]
dataset_labels = ['Open I', 'Open II', 'Open AMBER', 'Closed I',  'Closed II', 'Closed AMBER', '150 ns REMD']
dataset_colors = ['tab:blue', 'tab:blue', 'tab:blue', 'tab:blue', 'tab:blue', 'tab:blue', 'red']

for dataset, label, color in zip(whole_dataset, dataset_labels, dataset_colors):
    print(f'plotting {label} in {color}')
    if label == 'Amber-16 ff14sb aMD':
        axs.scatter(dataset[0], dataset[1], color=color, s=5, zorder=1, marker='o', alpha=0.5)
    else:
        axs.scatter(dataset[0], dataset[1], color=color, s=5, zorder=1, marker='x')
    #axs.scatter(0,0, color=color, s=60, label=f'{label} ({int(math.floor(len(dataset[0])/100))}) ns')
    axs.scatter(0, 0, color=color, s=60, label=f'{label}')
    axs.scatter(dataset[0][0], dataset[1][0], marker='>', s=500, color=color, edgecolors='black', zorder=2, linewidth=5)
    axs.scatter(dataset[0][-1], dataset[1][-1], marker='s', s=500, color=color, edgecolors='black', zorder=2, linewidth=5)

axs.set_xlabel('R(g) / nm')
axs.set_ylabel('SASA / nm^2')
axs.scatter(0,0, marker='>', color='white', edgecolors='black', linewidth=8 ,label='1st frame', s=100)
axs.scatter(0,0, marker='s', color='white', edgecolors='black', linewidth=5 ,label='Last frame', s=100)
axs.set(xlim=(1.2, 2.8), ylim=(63, 110))
axs.legend()
"""
"""
open_i_distance = get_simple_dataset(f'{base_path}trp_gggggg_pdz_open_i/sexmgrace/domains_distance_histo.xvg')
open_ii_distance = get_simple_dataset(f'{base_path}trp_gggggg_pdz_open_ii/trp_pdz_open_ii/xmgrace/domains_distance_histo.xvg')
open_amber_distance = get_simple_dataset(f'{base_path}trp_pdz_open_amber/xmgrace/domains_distance_histo.xvg')
closed_i_distance = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_i/xmgrace/domains_distance_histo.xvg')
closed_ii_distance = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_ii/xmgrace/domains_distance_histo.xvg')
closed_amber_distance = get_simple_dataset(f'{base_path}trp_gggggg_pdz_closed_amber/xmgrace/domains_distance_histo.xvg')

distance_dataset = [open_i_distance, open_ii_distance, open_amber_distance, closed_i_distance, closed_ii_distance, closed_amber_distance]

fig_dist, axs_dist = plt.subplots()
for dataset, label, color in zip(distance_dataset, dataset_labels, dataset_colors):
    axs_dist.plot(dataset[0], dataset[1], label=label, color=color)

axs_dist.set_xlabel('COM TrpCage - COM PDZ3 distance / nm')
axs_dist.set_ylabel('Probability')
axs_dist.set(xlim=(1.5,3.5))
axs_dist.legend()


hiv_rg = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/md_gyrate.xvg')
hiv_sasa = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/md_sasa.xvg')

hiv_rg_330 = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/annealing_gyrate.xvg')
hiv_sasa_330 = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/analyses/annealing_sasa.xvg')

fig_hiv, axs_hiv = plt.subplots(ncols=2)
axs_hiv[0].scatter(hiv_rg[1], hiv_sasa[1], s=100, marker='x')
axs_hiv[0].set_xlabel('R(g) / nm')
axs_hiv[0].set_ylabel('SASA / nm^2')
axs_hiv[0].set_title('HmII repr., 310 K')

axs_hiv[1].scatter(hiv_rg_330[1], hiv_sasa_330[1], s=100, marker='x')
axs_hiv[1].set_xlabel('R(g) / nm')
axs_hiv[1].set_ylabel('SASA / nm^2')
axs_hiv[1].set_title('HmII repr., 310->330 K')




ax3d = plt.axes(projection='3d')
for dataset, distance, label, color in zip(whole_dataset, distance_dataset, dataset_labels, dataset_colors):
    ax3d.scatter(dataset[0][::20], dataset[1][::20], distance[::20])
    ax3d.scatter(0,0,0, color=color, s=60, label=f'{label} ({int(math.floor(len(dataset[0])/100))}) ns')
    ax3d.scatter(dataset[0][0], dataset[1][0], distance[0], marker='>', s=500, color=color, edgecolors='black', zorder=2, linewidth=5)
    ax3d.scatter(dataset[0][-1], dataset[1][-1], distance[-1], marker='s', s=500, color=color, edgecolors='black', zorder=2, linewidth=5)
ax3d.set_xlabel('R(g) / nm')
ax3d.set_ylabel('SASA / nm^2')
ax3d.set_zlabel('No. interdomain contacts')
ax3d.scatter(0,0,0, marker='>', color='white', edgecolors='black', linewidth=8 ,label='1st frame', s=100)
ax3d.scatter(0,0,0, marker='s', color='white', edgecolors='black', linewidth=5 ,label='Last frame', s=100)
ax3d.set(xlim=(1.2, 2.8), ylim=(63, 93), zlim=(0, 5000))
ax3d.legend()



rama_md = read_rama('/run/timeshift/backup/IOCB/md/trp_gggggg_pdz_closed_amber/rama.xvg')
times = []
for time in rama_md['Times'][0]:
    times.append(time)
limres = ['GLY-24', 'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30']

new_dataset = {'Names': [], 'Phis': [], 'Psis': []}
for i, dataset_name in enumerate(rama_md['Names']):
    if dataset_name in limres:

        new_dataset['Names'].append(rama_md['Names'][i])
        new_dataset['Phis'].append(rama_md['Phis'][i])
        new_dataset['Psis'].append(rama_md['Psis'][i])

histograms = {'Names': [], 'Phis': [], 'Psis': []}
for name, phi, psi in zip(new_dataset['Names'], new_dataset['Phis'], new_dataset['Psis']):
    if name in limres:

        phistogram = np.histogram(phi, bins=360, density=True)
        psistogram = np.histogram(psi, bins=360, density=True)

        histograms['Names'].append(name)
        histograms['Phis'].append(phistogram)
        histograms['Psis'].append(psistogram)

fig_angle_hist, axs_angle_hist = plt.subplots(nrows=2,ncols=len(limres), subplot_kw={'projection': 'polar'})
fig_angle_hist_normal, axs_angle_hist_normal = plt.subplots(nrows=2,ncols=len(limres))
i = 0


for name, phi_hist, psi_hist in zip(histograms['Names'], histograms['Phis'], histograms['Psis']):

    theta_phi = [math.radians(x) for x in phi_hist[1][:-1]]
    r_phi = [y for y in phi_hist[0]]

    theta_psi = [math.radians(x) for x in psi_hist[1][:-1]]
    r_psi = [y for y in psi_hist[0]]

    axs_angle_hist[0][i].plot(theta_phi,r_phi, label=f'{name} Phi')
    axs_angle_hist[0][i].legend()
    #axs_angle_hist[0][i].set_xticklabels([i for i in range(-180,180, 1)])
    #axs_angle_hist[0][i].set_theta_offset(np.pi / 2)
    axs_angle_hist[1][i].plot(theta_psi,r_psi,  label=f'{name} Psi')
    #axs_angle_hist[1][i].set_xticklabels([i for i in range(-180, 180, 1)])
    #axs_angle_hist[0][i].set_theta_offset(np.pi / 2)
    axs_angle_hist[1][i].legend()

    axs_angle_hist_normal[0][i].plot(phi_hist[1][:-1], phi_hist[0])
    axs_angle_hist_normal[1][i].plot(psi_hist[1][:-1], psi_hist[0])

    i+=1
"""
#plt.tight_layout()
plt.show()
