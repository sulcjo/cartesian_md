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

        if time_coloring:
            data = pd.DataFrame([x,y, times])
            data = data.T
            data.columns = ('x','y','times')

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

def plot_angles_chi(dataset_array, title='', plot='scatter'):

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

        # Now iterate through the dataset array and see if the current iteration is the same aminoacid as before
        # In that case keep it in the current row, else (different AA) skip to next row
        # v this splits the aminoacid name chain after the first CAPITAL letter (file names such as phiARG44, chi1GLY65)
        last_aminoacid = re.findall('[A-Z]{3}[1-9]*', dataset_array[0][0])
        row = 0
        col = 0
        for dataset in dataset_array:
            current_aminoacid = re.findall('[A-Z]{3}[1-9]*', dataset[0])

            if  current_aminoacid != last_aminoacid:
                row += 1
                col = 0
            last_aminoacid = re.findall('[A-Z]{3}[1-9]*', dataset[0])

            if nrows > 1:
                print(f'populating row:{row} col:{col} with {current_aminoacid}')

                if plot == 'scatter':
                    axs[row][col].scatter(dataset[1], dataset[2], s=1)
                elif plot == 'line':
                    axs[row][col].plot(dataset[1], dataset[2])
                axs[row][col].set_title(dataset[0])
                # A quick and dirty fix to see if it's plotting correlation or angle data (corr doesn't go over 1)
                if max(dataset[2]) > 1:
                    axs[row][col].set(ylim=(-180, 180))
                else:
                    axs[row][col].set(ylim=(-1,1))
            else:
                if plot == 'scatter':
                    axs[col].scatter(dataset[1], dataset[2], s=1)
                elif plot == 'line':
                    axs[col].plot(dataset[1], dataset[2])

                axs[col].set_title(dataset[0])
                if max(dataset[2]) > 1:
                    axs[col].set(ylim=(-180, 180))
                else:
                    axs[col].set(ylim=(-1,1))
            col += 1
        plt.suptitle(title)

base_path = '/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/chi/pyplot/'
base_path_corr = '/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/chi/corrs/pyplot/'
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

angles_dataset = obtain_dataset_array(base_path)
correlation_dataset = obtain_dataset_array(base_path_corr)

plot_angles_chi(angles_dataset, title='TrpCage-GGGGGG-PDZ3 open start')
plot_angles_chi(correlation_dataset, title='TrpCage-GGGGGG-PDZ3 open start', plot='line')

"""
# Load RMSD
#rsmd = get_simple_dataset('/home/sulcjo/IOCB/md/trp_gggggg_pdz_closed_400ns/xmgrace/')
#rg = get_simple_dataset('/home/sulcjo/IOCB/md/trp_gggggg_pdz_closed_400ns/xmgrace')
#domain_distance = get_simple_dataset('/home/sulcjo/IOCB/md/trp_gggggg_pdz_closed_400ns/xmgrace')
"""

# Load umbrella contacts, histograms and pmf curve
#us_contacts = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/umbrella/nolinker/wham_results/contacts_pulling.xvg')
#pmf_curve = get_umbrella_profile('/run/media/sulcjo/sulcjo-data/IOCB/md/umbrella/nolinker/wham_results/profile_errors.xvg')
#pmf_histograms = get_umbrella_histogram('/run/media/sulcjo/sulcjo-data/IOCB/md/umbrella/nolinker/wham_results/histo.xvg')


## Plot US
"""
fig_us, axs_us = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios' : [3,1,1]})
setLineColor = 'blue'
setFontSizeLarge = 18
setFontSizeMedium = 14

axs_us[0].errorbar(pmf_curve[0], pmf_curve[1], yerr=pmf_curve[2], capsize=5, ecolor='black',color=setLineColor)
axs_us[0].set_title('PMF Curve for TrpCage...PDZ3 pulling', fontsize=setFontSizeLarge)
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

#plot_3D_scatter(rama, angle_coloring=True, title='Linker (TrpCage-GGGGGG-PDZ3)', limit_resi=(24,30))
#plot_angles(rama, 'Phi', ma=True, angle_limits = (-90, -35), title='US Linker (TrpCage-GGGGGG-PDZ3)', limit_resi=(24,30))
#plot_angles(rama, 'Psi', ma=True, angle_limits = (-70, -15), title='US Linker (TrpCage-GGGGGG-PDZ3)', limit_resi=(24,30))

#fig_md_closed, axs_md_closed = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios' : [3,1,1]})
#rama_md = read_rama('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_closed_400ns/rama.xvg')
#plot_angles(rama_md, 'Phi', ma=True, angle_limits = None, title='MD Linker (TrpCage-GGGGGG-PDZ3 closed start)', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'])
#plot_angles(rama_md, 'Psi', ma=True, angle_limits = None, title='MD Linker (TrpCage-GGGGGG-PDZ3 closed start)', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'])
#plot_ramachandran(rama_md, title='MD Linker (TrpCage-GGGGGG-PDZ3 open start)', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'])

fig_md_closed, axs_md_closed = plt.subplots(nrows=3, ncols=1, gridspec_kw={'height_ratios' : [3,1,1]})
#rama_md = read_rama('/run/media/sulcjo/sulcjo-data/IOCB/md/pdz_gggggg_trp_150ns/rama.xvg')
#plot_angles(rama_md, 'Phi', ma=True, angle_limits = None, title='MD Linker (PDZ3-GGGGGG-TrpCage)', limit_resi=(0,100), limres = ['GLY-83' ,'GLY-84', 'GLY-85', 'GLY-86', 'GLY-87', 'GLY-88'])
#plot_angles(rama_md, 'Psi', ma=True, angle_limits = None, title='MD Linker (PDZ3-GGGGGG-TrpCage)', limit_resi=(0,100), limres = ['GLY-83' ,'GLY-84', 'GLY-85', 'GLY-86', 'GLY-87', 'GLY-88'])
#plot_ramachandran(rama_md, title='MD Linker (PDZ3-GGGGGG-TrpCage)', limit_resi=(0,100), limres = ['GLY-83' ,'GLY-84', 'GLY-85', 'GLY-86', 'GLY-87', 'GLY-88'], time_coloring=True)

"""
plot_angles(rama_md, 'Phi', ma=True, angle_limits = None, title='MD Linker (TrpCage-GGGGGG-PDZ3 closed start)', limit_resi=(0,100), limres = ['ARG-97', 'GLU-98'])
plot_angles(rama_md, 'Psi', ma=True, angle_limits = None, title='MD Linker (TrpCage-GGGGGG-PDZ3 closed start)', limit_resi=(0,100), limres = ['ARG-97', 'GLU-98'])

rama_md = read_rama('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/rama.xvg')
plot_ramachandran(rama_md, title='MD Linker (TrpCage-GGGGGG-PDZ3 open start)', limit_resi=(0,100), limres = ['GLY-24' ,'GLY-25', 'GLY-26', 'GLY-27', 'GLY-28', 'GLY-29', 'GLY-30'], time_coloring=True)

fig_us, axs_us = plt.subplots(nrows=3, ncols=1)
setLineColor = 'blue'
setFontSizeLarge = 18
setFontSizeMedium = 14

domain_distance = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_pdz_closed_nolinker/xmgrace/distance_domains.xvg')
domain_distance_open = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/xmgrace/distance_domains.xvg')
axs_us[0].plot(domain_distance[0], domain_distance[1],color=setLineColor, label='Closed start')
axs_us[0].plot(domain_distance_open[0], domain_distance_open[1], color='red',label='Open Start')
axs_us[0].set_title('TrpCage...PDZ3 no linker COM distances (closed start)', fontsize=setFontSizeLarge)
axs_us[0].set_title('TrpCage-GGGGGG-PDZ3 COM distances', fontsize=setFontSizeLarge)
axs_us[0].set_ylabel('d / nm', fontsize=setFontSizeMedium)
axs_us[0].legend()

domain_contacts = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_pdz_closed_nolinker/xmgrace/domains_contacts.xvg')
domain_contacts_open = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/xmgrace/domain_contacts.xvg')
axs_us[1].plot(domain_contacts[0], domain_contacts[1],color=setLineColor)
axs_us[1].plot(domain_contacts_open[0], domain_contacts_open[1],color='red')
axs_us[1].set_title('TrpCage-GGGGGG-PDZ3 COM contacts', fontsize=setFontSizeLarge)
axs_us[1].set_ylabel('count', fontsize=setFontSizeMedium)

domain_histograms = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_pdz_closed_nolinker/xmgrace/domains_distance_histo.xvg')
domain_histograms_open = get_simple_dataset('/run/media/sulcjo/sulcjo-data/IOCB/md/trp_gggggg_pdz_open/xmgrace/domains_distance_histo.xvg')
axs_us[2].plot(domain_histograms[0], domain_histograms[1],color=setLineColor)
axs_us[2].plot(domain_histograms_open[0], domain_histograms_open[1],color='red')
axs_us[2].set_title('TrpCage-GGGGGG-PDZ3 COM distances frequency', fontsize=setFontSizeLarge)
axs_us[2].set_ylabel('Count', fontsize=setFontSizeMedium)

#domain_distances_autocorr = get_simple_dataset('/home/sulcjo/IOCB/md/trp_gggggg_pdz_closed_400ns/xmgrace/autocorr_domain_distances.xvg')
#axs_us[3].plot(domain_distances_autocorr[0], domain_distances_autocorr[1],color=setLineColor)
#axs_us[3].set_title('PDZ3-TrpCage (Closed start) COM distances autocorrelation', fontsize=setFontSizeLarge)
#axs_us[3].set_ylabel('G(t)', fontsize=setFontSizeMedium)
#axs_us[3].set_xlabel('t / ps', fontsize=setFontSizeMedium)
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
"""



#plt.tight_layout()
plt.show()
