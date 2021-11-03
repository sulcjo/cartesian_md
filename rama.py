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

        for line in data_lines[0:10]:
            split_line = line.split()
            residue = split_line[-1]
            if residue not in residues['Names']:
                residues['Names'].append(residue)

    #Append an empty array for every residue in 'Names' into 'Phis' and 'Psis'
        for residue in residues['Names']:
            residues['Phis'].append([])
            residues['Psis'].append([])

        for line in data_lines:
            for i, residue in enumerate(residues['Names']):
                if residue in line:
                    try:
                        split_line = line.split()
                        residues['Phis'][i].append(float(split_line[0]))
                        residues['Psis'][i].append(float(split_line[1]))
                    except:
                        continue
    #Get times which are only implicitly available from rama.xvg, every 10th frame was saved, hence the 10x multiplier
        number_residues = len(residues['Names'])
        total_data_lines = len(data_lines)
        max_time = 10*(total_data_lines/number_residues)
        residues['Times'].append([int(time) for time in range(0, int(max_time), 10)])

        return(residues)

#Add option for customizing output
def anglea_coloring(psi_x, phi_y, phi_limits = (-90, -35), psi_limits = (-70, -15), colors = ('red', 'blue')):
    #If residue psi and phi lies in a certain Ramachandran plot area, it will be colored
    if psi_limits[0] <= psi_x <= psi_limits[1] and phi_limits[0] <= phi_y <= phi_limits[1]:
        return colors[0]
    else:
        return colors[1]

def plot_3D_scatter(dataset, angle_coloring = False):

    #Convert ps times to ns times for improved graph clarity
    ns_times = []
    for ps_time in dataset['Times'][0]:
        ns_times.append(ps_time / 1000)

    fig = plt.figure(figsize=(8*(len(dataset['Names'])), 10))
    fig.suptitle(f'Ramachandran plot time dependence for 7Trp, {forcefield}', fontsize=24)
    for i, residue in enumerate(dataset['Names']):



        ax = fig.add_subplot(1,len(dataset['Names']), i+1, projection='3d')

        z = np.array(ns_times)
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

def plot_angles(dataset, type = 'Psi', ma = False, angle_limits = None):
    ns_times = []
    for ps_time in dataset['Times'][0]:
        ns_times.append(ps_time / 1000)

    fig, axs = plt.subplots(nrows=1, ncols=len(dataset['Names']), figsize=(16,4))
    fig.suptitle(f'{type} angle for 7Trp, {forcefield}', fontsize=24)
    plt.subplots_adjust(left=0.06,
                        bottom=0.15,
                        right=0.988,
                        top=0.80,
                        wspace=0.405,
                        hspace=0.285)



    x = dataset['Times'][0]
    for i, residue in enumerate(dataset['Names']):

        if type == 'Phi':
            y = dataset['Phis'][i]
        elif type == 'Psi':
            y = dataset['Psis'][i]

        axs[i].scatter(ns_times, y, s = 1.4)

        if ma:
            df = pd.DataFrame(y)
            ma_dataset = df.rolling(window=100).mean()
            axs[i].plot(ns_times, ma_dataset, color='black', linewidth = 1, alpha=0.8)

        if angle_limits:
            #Plot line for lower and upper angle limits
            axs[i].plot(ns_times, [angle_limits[0] for time in ns_times], color='red', linewidth=2, linestyle='--')
            axs[i].plot(ns_times, [angle_limits[1] for time in ns_times], color='red', linewidth=2, linestyle='--')

        axs[i].set_title(residue, fontsize=20)
        axs[i].tick_params(axis='y', labelsize=15)
        axs[i].tick_params(axis='x', labelsize=15)
        axs[i].set_xlabel('Time / ns', fontsize=16)
        axs[i].set_ylabel(f'{type}', fontsize=16)

forcefield = 'AMBER ff99sb-ILDN'
dataset = read_rama(path+'7Trp_amber/rama.xvg')
plot_3D_scatter(dataset, angle_coloring=True)
plot_angles(dataset, 'Phi', ma=True, angle_limits = (-90, -35))
plot_angles(dataset, 'Psi', ma=True, angle_limits = (-70, -15))

# psi_limits = (-90, -35), phi_limits = (-70, -15)

plt.show()