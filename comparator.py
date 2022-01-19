import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import math
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import operator
from matplotlib import gridspec
from numpy.linalg import eig
from sklearn.preprocessing import MinMaxScaler, StandardScaler


class Comparator:

    # Comparator takes proteinModel instances as an argument
    # It has access to all the datasets present in the models
    def __init__(self, proteinModels = None):
        plt.rcParams.update({'font.sans-serif': 'Verdana'})
        models = proteinModels
        if isinstance(proteinModels,list):
            self.proteinModels = models
        else:
            self.proteinModels = [models]

        self.number_models = len(self.proteinModels)

        self.setMaxColumns = 5
        self.setFontSizeLarge = 24
        self.setLineWidth = 0.8
        self.setLineColor = 'blue'
        self.setFigSize = (20, 15)
        self.setFontSizeMedium = 16
        self.setFontSizeSmall = 14


        print(f'COMPARATOR LOADED MODELS {[model.annotation for model in self.proteinModels]}')

        plt.rc('xtick', labelsize=self.setFontSizeMedium)
        plt.rc('ytick', labelsize=self.setFontSizeMedium)

    def load_model(self, proteinModel):
        self.proteinModels.append(proteinModel)
        print(f'MANUAL LOADING {[model.annotation for model in proteinModel]}')

    def __get_model_indexes(self):
        return([i for i in range(self.number_models)])

    def plot_simple_dataset(self, dataset = 'rmsd', modelIndexes = None, subtitle = ''):
        # modelIndexes are a tuple containing indexes of proteinModels to be plotted
        # If none are provided, all the loaded models will be used
        if not modelIndexes:
            modelIndexes = self.__get_model_indexes()

        # Split plots, resulting plots will be stacked vertically, but at max 5 per column
        ncols = math.ceil(len(modelIndexes) / self.setMaxColumns)
        nrows = math.ceil(len(modelIndexes) / ncols)

        # Prepare matplotlib objects
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=self.setFigSize)
        plt.suptitle(subtitle, size=self.setFontSizeLarge)
        fig.tight_layout()
        plt.subplots_adjust(wspace=0.53, hspace=0.357, top=0.940, right=0.920, bottom=0.050, left=0.050)
        #plt.subplots_adjust(top=0.940)
        # Fill rows and columns of subplots
        for i, model in enumerate(modelIndexes):
            if i == 0:
                row = 0
                col = 0

            if row == nrows:
                row = 0
                col += 1

            # Obtain datasets from model __dict__

            try:
                variant_dataset = self.proteinModels[model].datasets[dataset]
            except:
                print(f"ERROR PLOTTING {dataset} OF {self.proteinModels[model].annotation}")

                continue

            if ncols == 1:

                axs[row].plot(variant_dataset[0], variant_dataset[1], linewidth=self.setLineWidth, color=self.setLineColor)
                axs[row].set_title(variant_dataset[3][0].replace("@    title \"", f"{self.proteinModels[model].annotation} ").replace("\"", ""))
                axs[row].set_xlabel(variant_dataset[3][1].replace("@    xaxis  label \"", "").replace("\"", ""), size=self.setFontSizeMedium)
                axs[row].set_ylabel(variant_dataset[3][2].replace("@    yaxis  label \"", "").replace("\"", ""), size=self.setFontSizeMedium)
                axs[row].grid(True)

            else:

                axs[row][col].plot(variant_dataset[0], variant_dataset[1], linewidth=self.setLineWidth, color=self.setLineColor)
                axs[row][col].set_title(variant_dataset[3][0].replace("@    title \"", f"{self.proteinModels[model].annotation} ").replace("\"", ""))
                axs[row][col].set_xlabel(variant_dataset[3][1].replace("@    xaxis  label \"", "").replace("\"", ""), size=self.setFontSizeMedium)
                axs[row][col].set_ylabel(variant_dataset[3][2].replace("@    yaxis  label \"", "").replace("\"", ""), size=self.setFontSizeMedium)
                axs[row][col].grid(True)

            row += 1

    def convert_units(self, dataset, orig, new):
        if orig == 'kcal' and new == 'kJ':
            return([value*4.184 for value in dataset] )

    def plot_mmpbsa(self, dataset = 'gmmpbsa', modelIndexes = None, title = ''):
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
                    fontsize=self.setFontSizeMedium,
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

        if not modelIndexes:
            modelIndexes = self.__get_model_indexes()

        for index in modelIndexes:
            try:
                # Types: gmmpbsa, gmxmmpbsa_pb, gmxmmpbsa_gb

                variant_dataset = self.proteinModels[index].datasets[dataset]
                binding_energy = float(variant_dataset['binding'])
                binding_error = float(variant_dataset['bindinge'])

                binding_energies.append(binding_energy)
                binding_energies_errors.append(binding_error)
                variants_plotted.append(self.proteinModels[index].annotation)
            except:
                print(f"ERROR PLOTTING MM/PBSA FOR {self.proteinModels[index].annotation}")
                continue

            if binding_energy + binding_error < 0 and binding_energy - binding_error < 0:
                clrs.append("g")
            elif binding_energy + binding_error > 0 and binding_energy - binding_error > 0:
                clrs.append('r')
            else:
                clrs.append('y')


        fig, ax = plt.subplots(figsize=self.setFigSize)
        ax.bar(variants_plotted, binding_energies, color=clrs, yerr=binding_energies_errors, alpha=0.8, capsize=4,
               width=0.5, align='edge')
        plt.xticks(variants_plotted)
        if 'gmxmmpbsa' in dataset:
            binding_energies = self.convert_units(binding_energies, 'kcal', 'kJ')
            binding_energies_errors = self.convert_units(binding_energies_errors, 'kcal', 'kJ')


        ax.set_ylabel(r"$\Delta$G  / kJ/mol ", size=self.setFontSizeMedium)
        ax.set_title(f"{title}", size=self.setFontSizeLarge)
        ax.yaxis.grid(True)
        y_ticks = ax.get_yticks()
        ax.set_yticklabels(labels=y_ticks, fontsize=self.setFontSizeLarge)
        ax.set_xticklabels([])
        # plt.tight_layout()
        add_value_labels(ax)
        legend_entries = {'> 0': 'r', '?': 'y', '< 0': 'g'}
        labels = list(legend_entries.keys())
        handles = [plt.Rectangle((0, 0), 1, 1, color=legend_entries[label]) for label in labels]
        plt.legend(handles, labels, fontsize=self.setFontSizeMedium, loc='upper right')

    def plot_umbrella(self, modelIndexes = None, fit = False, stderror = False, check_sampling = False, check_sampling_limit = 100, print_values = True, max_traj = None):
        """
        :param modelIndexes: Which models (zero indexing) to plot from the models included in this Comparator class (modelIndexes)
        :param fit: Do an exponential fit of the curve and calculate max. PMF
        :param stderror: Add error-bars, requires a file containing errors calculated by GMX WHAM
        :param check_sampling: Calculates total sampling for all x-values and plots them as a curve. Warns if any X-value has < check_sampling_limit counts.
        :param check_sampling_limit: Warning limit of counts for sampling, 100 by default.
        :return: Doesn't return anything. Creates an ax instance of matplotlib.
        """

        if not modelIndexes:
            modelIndexes = self.__get_model_indexes()

        for index in modelIndexes:
            try:
                profile_x = self.proteinModels[index].datasets['umbrella_profile'][0]
                profile_y = self.proteinModels[index].datasets['umbrella_profile'][1]
                histo_x = self.proteinModels[index].datasets['umbrella_histogram'][0]
                histo_y = self.proteinModels[index].datasets['umbrella_histogram'][1]






                fig, axs = plt.subplots(nrows=2, ncols=1, figsize=self.setFigSize, sharex=True)
                for run in histo_y:
                    axs[1].plot(histo_x, run, color=self.setLineColor)
                    axs[1].set_title('Histograms', fontsize=self.setFontSizeLarge)
                    axs[1].set_xlabel('COM-COM distance', fontsize=self.setFontSizeMedium)
                    axs[1].set_ylabel('Counts', fontsize=self.setFontSizeMedium)




                if check_sampling:
                    if 'Cumbrella_histogram_total_sampling' in self.proteinModels[index].Cdatasets:
                        sampling = self.proteinModels[index].Cdatasets['Cumbrella_histogram_total_sampling']
                    else:
                        sampling = []
                        for i, x in enumerate(histo_x):
                            y_values_for_x = [value[i] for value in histo_y]
                            y_values_for_x = filter(None, y_values_for_x)
                            sampling.append(sum(y_values_for_x))
                        self.proteinModels[index].Cdatasets['Cumbrella_histogram_total_sampling'] = sampling

                        # y = filter(None, y)
                    axs[1].plot(histo_x, sampling,color='red')
                    inadequate = []
                    for x, sample in zip(histo_x, sampling):
                        if sample < 10:
                            inadequate.append(x)
                    print(f'{self.proteinModels[index].annotation} WARNING INADEQUATE SAMPLING IN {inadequate}')



                if stderror:
                    profile_stderror = self.proteinModels[index].datasets['umbrella_profile'][2]
                    axs[0].errorbar(profile_x, profile_y, yerr=profile_stderror, capsize=5, ecolor='black')
                axs[0].plot(profile_x, profile_y, color=self.setLineColor)
                axs[0].set_title(r'$\zeta$ trajectory energy profile', fontsize=self.setFontSizeLarge)
                axs[0].set_xlabel(r'$\zeta$ distance', fontsize=self.setFontSizeMedium)
                axs[0].set_ylabel('Energy kJ/mol', fontsize=self.setFontSizeMedium)
                plt.suptitle(self.proteinModels[index].annotation, fontsize=self.setFontSizeLarge)



                """
                Method for getting minimum of PMF curve (ignoring outliers caused by bad sampling) and maximum
                (through exponential fitting to defeat oscillations and bad convergence) to calculated DeltaG
                
                """




                # If limits are set, this ignores them, so it can find a minimum outside the limits of profile_x, fix it
                def get_minimum(profile_y, histo_y, histo_x):
                    sampling = []
                    for i, x in enumerate(histo_x):
                        y_values_for_x = [value[i] for value in histo_y]
                        y_values_for_x = filter(None, y_values_for_x)
                        sampling.append(sum(y_values_for_x))
                    self.proteinModels[index].Cdatasets['Cumbrella_histogram_total_sampling'] = sampling

                    new_sampling = []
                    new_profile_y = []

                    for i, value in enumerate(sampling):
                        if value > 50:
                            new_sampling.append(sampling[i])
                            new_profile_y.append(profile_y[i])

                    return (min(new_profile_y))
                """
                def monoExp(x, m, t, b):
                    return m * np.exp(-t * x) + b

                p0 = (0, 2, 20)
                if max_traj:
                    limited_profile_x = [x for x in profile_x if float(x) < max_traj]
                    limited_profile_y = profile_y[:len(limited_profile_x)]
                    params, cv = scipy.optimize.curve_fit(monoExp, limited_profile_x, limited_profile_y, p0, maxfev=10000)
                else:
                    params, cv = scipy.optimize.curve_fit(monoExp, profile_x, profile_y, p0, maxfev=10000)
                m, t, b = params

                exp_y = []
                for value in profile_x:
                    exp_y.append(monoExp(value, m, t, b))

                
                maximum = b
                """

                if max_traj:
                    limited_profile_x = [x for x in profile_x if float(x) < max_traj]
                    limited_profile_y = profile_y[:len(limited_profile_x)]
                    maximum = max(limited_profile_y)
                else:
                    maximum = max(profile_y)
                minimum = get_minimum(profile_y, histo_y, histo_x)

                delta = round(maximum - minimum, 2)

                print(f'{self.proteinModels[index].annotation} min {minimum} upper {maximum}, delta {delta}')
                axs[0].plot(profile_x, [minimum for i in profile_x], color='red', label=f'{delta} kJ/mol')
                axs[0].plot(profile_x, [maximum for i in profile_x], color='red')
                axs[0].legend()
                if max_traj:
                    axs[0].set_xlim(right=max_traj)
                    axs[1].set_xlim(right=max_traj)


                """
                END OF METHOD
                """

                #plt.savefig(f'/run/media/sulcjo/Data/in silico/transfer/OBRAZKY/sirah/atomistic_{self.proteinModels[index].annotation}.png')

                if fit:

                    axs[0].plot(profile_x, exp_y, '--', label=fr"Monoexp. fit $\Delta$PMF = {round(b, 2)} kJ/mol",
                                color='r')





                    # kd = calculateKd(round(b, 2), 300)
                    # axs[0].scatter(profile_x, exp_y, s=0, label=f'Kd = {format(kd,".1E")} M')
                    axs[0].legend(fontsize=self.setFontSizeMedium)
            except KeyError:
                print(f'ERROR PLOTTING UMBRELLA SAMPLING PLOTS FOR {self.proteinModels[index].annotation}')



    def compare_IEM(self, modelIndexes, dataset = 'total_IEM'):
        """
        Both matrices are normalized first, this is done by obtaining a matrix mean, subtracting it from all the values
        and the dividing all values by the maximum value. This creates a normalized matrix, with values ranging from -1 to 1.
        These matrices are plotted alongside the delta matrix.

        IEM matrices from GRINN have a baseline for some reason. This baseline is obtained from each matrix so:
        1) Calculate mode of matrix values by column (returns a list of values)
        2) Calculate mean of this list (should return a list of means of modes) - this is mostly for making Pandas work properly
        3) Calculate mean of this list, returns a scalar
        4) Do this for both matrices, obtain two means of modes. Subtract this from each considered value
        5) Add 0.01 to baselines (empirical)

        Baselines are problematic, because no value is really zero, but a small number. This means that when comparing pairs
        present in matrix 1 with residues not present in matrix 2, some number between 0-1 (energy from matrix 1) would get divided by a small baseline.
        This would give a nonsensical result and effectively "imprints" matrix 1 into the difference matrix (where there's any value in matrix 1,
        it will be in matrix 2 as well). Even when subtracting the baseline, the problem remains, so a value of 0.01 was chosen to add to the
        baselines, so these new added baselines are used to subtract from all compared values. This may result in missing some insignificant differences however.

        Calculate abs(value1/value2) differences, there's several possible cases:
        A) Value1 > Value2: append 1/delta (so you get similarity between 0.00 and 1.00)
        B) Value1 < Value2: append delta
        c) Value2 == 0 (division by zero): append 0 (zero similarity)

        Then only plot these similarity values larger than 0.01 (1 %) to further crush the baseline and improve clarity and
        informational value

        :param modelIndexes: proteinModel instances with initialized IEM datasets to compare
        :param dataset: type of dataset (from proteinModel instances), default='total_IEM'
        :return: a difference dataset and a plt object
        """

        dataframe1 = self.proteinModels[modelIndexes[0]].datasets[dataset]
        #.replace(0, np.nan, inplace=True)
        dataframe2 = self.proteinModels[modelIndexes[1]].datasets[dataset]
        #dataframe1 = pd.DataFrame(MinMaxScaler().fit_transform(dataframe1))
        #dataframe2 = pd.DataFrame(MinMaxScaler().fit_transform(dataframe2))

        # Dataset normalization by all values, not rows or columns
        dataframe1 = dataframe1 - dataframe1.mean().mean()
        dataframe1 = dataframe1 / dataframe1.max().max()
        dataframe2 = dataframe2 - dataframe2.mean().mean()
        dataframe2 = dataframe2 / dataframe2.max().max()
        sequence = self.proteinModels[modelIndexes[0]].seq_1

        delta_dataframe = []
        y = 0
        baseline1 = dataframe1.mode().mean().mean() + 0.01
        baseline2 = dataframe2.mode().mean().mean() + 0.01

        similarity_non_zero = []
        while y < len(dataframe1):
            delta_dataframe.append([])
            for x in range(len(dataframe1)):
                value1 = float(dataframe1.iloc[x,y])
                value2 = float(dataframe2.iloc[x,y])
                if value2 <= baseline2 or value1 <= baseline1:
                    delta_value = 0
                elif value1 > value2:
                    delta_value = 1/(abs(value1/value2))
                elif value1 < value2:
                    delta_value = abs(value1 / value2)
                elif value1 == value2:
                    delta_value = 1
                else:
                    delta_value = 0

                if delta_value > 0.01:
                    delta_dataframe[y].append(delta_value)
                else:
                    delta_dataframe[y].append(0)


            y += 1

        delta_dataframe = pd.DataFrame(delta_dataframe)
        delta_dataframe.replace(0, np.nan)
        fig1, ax1 = plt.subplots(figsize=self.setFigSize)
        ax1 = sns.heatmap(dataframe1, cmap='RdBu', center=0, vmin=-1, vmax=1, cbar=True, cbar_kws={'label' : 'Lower = more stabilizing'})
        fig1.suptitle(f'Normalized IEM heatmap of {self.proteinModels[modelIndexes[0]].annotation}', size=self.setFontSizeLarge)

        fig2, ax2 = plt.subplots(figsize=self.setFigSize)
        ax2 = sns.heatmap(dataframe2, cmap='RdBu', center=0, vmin=-1, vmax=1, cbar=True, cbar_kws={'label' : 'Lower = more stabilizing'},)
        fig2.suptitle(f'Normalized IEM heatmap of {self.proteinModels[modelIndexes[1]].annotation}', size=self.setFontSizeLarge)

        fig3, ax3 = plt.subplots(figsize=self.setFigSize)
        #ax3.imshow(delta_dataframe)
        ax3 = sns.heatmap(delta_dataframe, cmap='Greys', robust=True, vmin=0.01, vmax=1, cbar_kws={'label' : '|[x1,y1]/[x2,y2]|'}, cbar=True)
        ax3.figure.axes[-1].yaxis.label.set_size(self.setFontSizeMedium)

        fig3.suptitle(str(f"$\Delta$-matrix of standardized matrices \n {self.proteinModels[modelIndexes[0]].annotation} {self.proteinModels[modelIndexes[1]].annotation} \n 1.00 means absolute similarity"), size=self.setFontSizeLarge)
        if self.proteinModels[modelIndexes[0]].split:
            protein_range = self.proteinModels[modelIndexes[0]].split[0]
            ligand_range = self.proteinModels[modelIndexes[0]].split[1]
        else:
            protein_range = (0, len(sequence))
            ligand_range = protein_range


        ticks_step = 10
        for plot in [ax1,ax2,ax3]:
            plot.set_xticks(ticks=[i for i in range(0, len(sequence), ticks_step)])
            plot.set_yticks(ticks=[i for i in range(0, len(sequence), ticks_step)])
            plot.set_xticklabels([i for i in sequence[::ticks_step]])
            plot.set_yticklabels([i for i in sequence[::ticks_step]])
            #plot.set(xlim=(ligand_range[0]+30, ligand_range[1]-40), ylim=(protein_range[0]+80, protein_range[1]-40))
            plot.set(xlim=(ligand_range[0], ligand_range[1]),
                     ylim=(protein_range[0], protein_range[1]))


            # Bottom frame line
            plot.axhline(y=protein_range[0], color='k', linewidth=2)

            # Top frame line
            plot.axhline(y=protein_range[1], color='k', linewidth=2)

            # Left frame line
            plot.axvline(x=ligand_range[0], color='k', linewidth=2)

            # Right frame line
            plot.axvline(x=ligand_range[1], color='k', linewidth=2)


        return(delta_dataframe)
        #values, vectors = eig(dataframe)
        #print(sum(values))
        #print(sum(vectors))


    def plot_IEM(self, modelIndexes = None, dataset = 'total_IEM',
                 vh_lines = True, sec_str=True, mark_resis=None, write_best_pairs=True):









        def get_best_pairs(dataframe, pairs=20):
            statistics_dataframe = dataframe.describe().T
            best = statistics_dataframe.nsmallest(n=pairs, columns='mean')
            best = best.drop_duplicates()

            return(best)

        if not modelIndexes:
            modelIndexes = self.__get_model_indexes()

        def obtain_sum_interactions(dataframe):
            sum_interactions = []
            for row_name, row in dataframe.items():
                sum_interactions.append(sum(row))

            return (sum_interactions)

        def get_sec_str_colors(ss):
            colors = []
            color_assignments = {'-': 'grey', 'E': 'blue', 'T': 'purple', 'S': 'purple', 'H': 'orange'}
            for struc in ss:
                try:
                    colors.append(color_assignments[struc])
                except:
                    colors.append('white')
            return (colors)

        for index in modelIndexes:
            # Build datasets
            type = dataset.split('_')[0]
            dataframe = self.proteinModels[index].datasets[dataset]




            """
            Temporary to show only some protein parts (only experimentally identified key interaction residues of the antibody
            
            
            print(self.proteinModels[index].seq_1)
            antibody_key_residues = ['W33', 'R50', 'E53', 'Y99', 'F100A', 'G100D', 'P100F', 'P100G', 'E100I', 'E100J', ]
            antibody_key_residues_indexes = []
            for i, aminoacid in enumerate(self.proteinModels[index].seq_1):
                if aminoacid in antibody_key_residues:
                    antibody_key_residues_indexes.append(i)
            #Non-key columns have to be removed 2x (once as rows in transposed, once as columns in normal - matrix is diagonally symetrical)
            #dataframe = dataframe.T
            for col in dataframe.columns:
                if col not in antibody_key_residues_indexes and col < 239:
                    dataframe[col].values[:] = 0
            dataframe = dataframe.T
            for col in dataframe.columns:
                if col not in antibody_key_residues_indexes and col < 239:
                    dataframe[col].values[:] = 0
            dataframe = dataframe.T

            print(self.proteinModels[index].seq_1)

"""
            excel_dataframe = pd.DataFrame(dataframe)
            excel_dataframe.columns = self.proteinModels[index].seq_3_chains
            excel_dataframe = excel_dataframe.loc[:, (excel_dataframe != 0).any(axis=0)]
            excel_dataframe = excel_dataframe.T
            excel_dataframe.columns = self.proteinModels[index].seq_3_chains
            excel_dataframe = excel_dataframe.loc[:, (excel_dataframe != 0).any(axis=0)]
            excel_dataframe = excel_dataframe.T

            """
            if dataset == 'total_IEM':
                excel_dataframe.to_excel(f'/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/md/md_310k_grinn_docker_2/python/all_{self.proteinModels[index].annotation}_total.xlsx')
            elif dataset == 'elec_IEM':
                excel_dataframe.to_excel(f'/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/md/md_310k_grinn_docker_2/python/all_{self.proteinModels[index].annotation}_elec.xlsx')
            elif dataset == 'vdw_IEM':
                excel_dataframe.to_excel(f'/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/md/md_310k_grinn_docker_2/python/all_{self.proteinModels[index].annotation}_vdw.xlsx')
            """







            if write_best_pairs:
                best_pairs = get_best_pairs(dataframe)
            else:
                best_pairs = ''

            sum_interactions = obtain_sum_interactions(dataframe)
            print(sum_interactions)

            sequence = self.proteinModels[index].seq_1_chains
            pairs_df = self.proteinModels[index].datasets[f'{type}_pairwise_IEM']
            title = f'{self.proteinModels[index].annotation} {dataset}'
            dataframe.replace(0, np.nan, inplace=True)
            #dataframe = pd.DataFrame(StandardScaler().fit_transform(dataframe))














            # Handle matrix splitting, if proteinModel doesnt have a split attribute, it will plot the complete (diagonally symmetrical matrix)
            if self.proteinModels[index].split:
                protein_range = self.proteinModels[index].split[0]
                ligand_range = self.proteinModels[index].split[1]
            else:
                protein_range = (0, len(sequence))
                ligand_range = protein_range












            # Save calculated values to proteinModel
            self.proteinModels[index].Cdatasets[f'Cbest_pairs_{dataset}'] = best_pairs
            self.proteinModels[index].Cdatasets[f'Csum_interactions_{dataset}'] = sum_interactions

            # Change step of ticks in heatmap
            ticks_step = 5




            # Prepare figure as a grid of 4 subplots
            fig = plt.figure(figsize=(15, 15))
            # plt.tight_layout()
            spec = gridspec.GridSpec(ncols=3, nrows=2, width_ratios=(4, 40, 1), height_ratios=(20, 3))
            ax_cbar = fig.add_subplot(spec[2])
            colors = ['red' if x < 0 else 'grey' if x == 0 else 'blue' for x in sum_interactions]
            # Plot heatmap using Seaborn
            ax3 = fig.add_subplot(spec[1])

            ax3 = sns.heatmap(dataframe, cmap='RdBu', robust=True, center=0,
                              cbar_kws={'label': 'kcal/mol'}, cbar=True, cbar_ax=ax_cbar, linewidths=0.25)
            # center=0
            # ax3 = plt.imshow(df, cmap='RdBu')

            # Plot ax1, upper-right plot showing ligand residues energies
            ax1 = fig.add_subplot(spec[4], sharex=ax3)

            ax1.scatter(sequence, sum_interactions, color=colors, marker='x')

            # Because heatmap is an image, it doesn't fill the whole subplot frame, therefore ax1 wouldn't correspond to heatmap axis
            # Add a little bit to the xlim (beyond 110 residues of MYO) so it compacts the frame
            # Hide x-axis ticks (resis)
            plt.setp(ax1.get_xticklabels(), visible=False)
            # Hide right, top and bottom parts of the frame
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.set_ylabel('kcal/mol')
            # Similar process for the ax2, showing Ab residues left of the heatmap
            ax2 = fig.add_subplot(spec[0], sharey=ax3)
            ax2.scatter(sum_interactions, sequence, color=colors, marker='x')
            ax2.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            ax2.set_xlabel('kcal/mol')
            plt.setp(ax2.get_yticklabels(), visible=False)

            ax2.tick_params(axis='both', which='both', length=0)
            ax1.tick_params(axis='both', which='both', length=0)

            ax3.set_xticks(ticks=[i for i in range(0, len(sequence), ticks_step)])
            ax3.set_yticks(ticks=[i for i in range(0, len(sequence), ticks_step)])
            ax3.set(xlim=ligand_range, ylim=protein_range)

            # Additional stuff
            if vh_lines:
                ax1.vlines(sequence, ymin=[0 for i in sequence], ymax=sum_interactions, color=colors)
                ax2.hlines(sequence, xmin=[0 for i in sequence], xmax=sum_interactions, color=colors)

            if mark_resis:
                for mark_resi in mark_resis:
                    if int(mark_resi) < ligand_range[0]:
                        ax3.scatter(ligand_range[0] + 2, mark_resi, color='red', marker='s', s=150, linewidths=0.25,
                                    edgecolors='black')

                    elif int(mark_resi) >= ligand_range[0]:
                        ax3.scatter(mark_resi, protein_range[0] + 2, color='red', marker='s', s=150, linewidths=0.25,
                                    edgecolors='black')
                    else:
                        print('Mark residue values out of range')

            if sec_str:
                if 'dssp' not in self.proteinModels[index].datasets:
                    self.proteinModels[index].DSSP_assign()
                secondary_structure = self.proteinModels[index].datasets['dssp']

                # Plot ligand SS (x-axis)
                ax3.scatter(dataframe.index.values[ligand_range[0]:ligand_range[1]],
                            [-2 for i in dataframe.index.values[ligand_range[0]:ligand_range[1]]],
                            color=get_sec_str_colors(secondary_structure[ligand_range[0]:ligand_range[1]]), s=150,
                            marker='s', linewidths=0.25, edgecolors='black')

                # Plot protein  SS (y-axis)
                ax3.scatter([(ligand_range[0] - 0.5) for i in dataframe.index.values], dataframe.index.values,
                            color=get_sec_str_colors(secondary_structure), s=150, marker='s', linewidths=0.25,
                            edgecolors='black')

                """ PLOTTING LIMITS """


                ax3.set(ylim=(protein_range[0] - 2.5, protein_range[1]))



                ax3.set(xlim=(ligand_range[0] - 1, ligand_range[1]))

            if mark_resis or sec_str:
                ax_ss_leg = fig.add_subplot(spec[3])
            else:
                ax_ss_leg = None

            if sec_str:
                # Construct a legend - make a non-visible scatter plot to obtain labels, put it into lower left subplot and hide everything plot-related, then create a legend
                for ss_color, ss_type in zip(['grey', 'blue', 'purple', 'orange'],
                                             ['Coil', r'$\beta$-sheet', 'Turn/Bend', r'$\alpha$-Helix']):
                    ax_ss_leg.scatter(-1000, -1000, color=ss_color, marker='s', label=ss_type, s=150)

            if mark_resis:
                ax_ss_leg.scatter(-1000, 1000, color='red', marker='s', label='Marked residues', s=150)

            if ax_ss_leg:
                ax_ss_leg.legend(ncol=1, loc='upper center', fontsize=20, frameon=False)
                ax_ss_leg.set(xlim=(0, 1), ylim=(0.1))
                ax_ss_leg.set_xticks([])
                ax_ss_leg.set_yticks([])
                for spine_pos in ax_ss_leg.spines:
                    ax_ss_leg.spines[spine_pos].set_visible(False)

            # Bottom frame line
            ax3.axhline(y=protein_range[0], color='k', linewidth=2)

            # Top frame line
            ax3.axhline(y=protein_range[1], color='k', linewidth=2)

            # Left frame line
            ax3.axvline(x=ligand_range[0], color='k', linewidth=2)

            # Right frame line
            ax3.axvline(x=ligand_range[1], color='k', linewidth=2)

            if write_best_pairs:
                best_pairs = get_best_pairs(pairs_df, pairs=14).T
                best_pairs_text = ''
                for row in best_pairs:
                    best_pairs_text += f'{row} {best_pairs[row]["mean"].round(decimals=1)} kcal/mol \n'

                ax3.scatter(x=-100, y=-100, label=best_pairs_text, marker='x', s=0)
                ax3.legend(loc='upper right', fontsize=self.setFontSizeSmall, frameon=False)



            plt.suptitle(title, size=self.setFontSizeLarge)
            """
            if dataset == 'total_IEM':
                plt.savefig(f'/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/md/md_310k_grinn_docker_2/python/all_{self.proteinModels[index].annotation}_totalIEMsquare.png')
            elif dataset == 'elec_IEM':
                plt.savefig(f'/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/md/md_310k_grinn_docker_2/python/all_{self.proteinModels[index].annotation}_elecIEMsquare.xlsx')
            elif dataset == 'vdw_IEM':
                plt.savefig(f'/run/media/sulcjo/sulcjo-data/IOCB/md/HIV/HmII_annealing/md/md_310k_grinn_docker_2/python/all_{self.proteinModels[index].annotation}_vdwIEMsquare.xlsx')
            """


    def show(self):




        plt.show()



"""
TO DO:

Plot saving, both in an image form and as an object so it can be loaded in the future (saved in the same way as proteinModels are going to be saved)

Conversion of units into every method (kJ <-> kcal)


"""