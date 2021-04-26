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


class Comparator:

    # Comparator takes proteinModel instances as an argument
    # It has access to all the datasets present in the models
    def __init__(self, proteinModels = None):
        self.proteinModels = proteinModels
        self.number_models = len(self.proteinModels)
        self.setMaxColumns = 5
        self.setFontSizeLarge = 24
        self.setLineWidth = 0.8
        self.setLineColor = 'blue'
        self.setFigSize = (20, 15)
        self.setFontSizeMedium = 16
        print(f'LOADED MODELS {[model.annotation for model in proteinModels]}')
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
        plt.subplots_adjust(wspace=0.140, hspace=0.450, top=0.950, right=0.920, bottom=0.050, left=0.050)

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
                    "Var. " + str(variant) + " : \n" + label + " +/- \n" + str(error_label),  # Use `label` as label
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

    def plot_umbrella(self, modelIndexes = None, fit = False, stderror = False, check_sampling = False, check_sampling_limit = 100):
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


                if fit:

                    # Fitting exponential function to obtain deltaG
                    def monoExp(x, m, t, b):
                        return m * np.exp(-t * x) + b

                    p0 = (0, 2, 60)
                    params, cv = scipy.optimize.curve_fit(monoExp, profile_x, profile_y, p0)
                    m, t, b = params

                    exp_y = []
                    for value in profile_x:
                        exp_y.append(monoExp(value, m, t, b))

                    axs[0].plot(profile_x, exp_y, '--', label=fr"Monoexp. fit $\Delta$PMF = {round(b, 2)} kJ/mol",
                                color='r')
                    # kd = calculateKd(round(b, 2), 300)
                    # axs[0].scatter(profile_x, exp_y, s=0, label=f'Kd = {format(kd,".1E")} M')
                    axs[0].legend(fontsize=self.setFontSizeMedium)
            except:
                print(f'ERROR PLOTTING UMBRELLA SAMPLING PLOTS FOR {self.proteinModels[index].annotation}')

    def plot_IEM(self, modelIndexes = None, dataset = 'total_IEM'):

        def get_best_pairs(dataframe, pairs=20):
            statistics_dataframe = dataframe.describe().T
            best = statistics_dataframe.nsmallest(n=pairs, columns='mean')
            best = best.drop_duplicates()

            return(best)

        if not modelIndexes:
            modelIndexes = self.__get_model_indexes()



        for index in modelIndexes:
            # Build datasets

            dataframe = self.proteinModels[index].datasets[dataset]
            best_pairs = get_best_pairs(dataframe)

            # Save calculated best pairs to proteinModel
            self.proteinModels[index].Cdatasets[f'Cbest_pairs_{dataset}'] = best_pairs









    def show(self):
        plt.show()