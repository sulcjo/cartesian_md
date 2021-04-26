import pickle
import types
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import operator
from matplotlib import gridspec
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os
from functools import wraps


def describeLoad(loadFunction):
    """
    Wrapper function that describes the loading process. All the loading methods of proteinModels take two arguments:
    path and dataset_name.

    Specific load function is selected by proteinModel.get_dataset using a dictionary of key:function pairs.
    proteinModel.get_dataset is wrapped by this function
    :param loadFunction: A specific function used to load a dataset, provided by proteinModel class
    :return: itself
    """

    @wraps(loadFunction)
    def function_wrapper(self, path, dataset_name):
        try:
            loadFunction(self, path, dataset_name)
            print(f' -----> LOADING {path} WITH {dataset_name}')
        except:
            print(f' XXXXXX UNSUCCESSFUL LOAD {path} OF {dataset_name}')

    return (function_wrapper)

class proteinModel:

    def __init__(self, annotation = 'Protein model', split = None, path=''):

        # If a pickle file is provided via path, this instance will set it's attributes according to it
        # This effectively returns an instance saved previously to a pickle file
        try:
            self.__load(path)

        except:

            self.annotation = annotation
            self.datasets = {}
            self.split = split
            self.Cdatasets = {}


    @staticmethod
    def is_iterable(set):
        if isinstance(set, tuple) or isinstance(set, list) or isinstance(set, types.GeneratorType) or isinstance(
                set, range):
            return True
        else:
            return False

    @staticmethod
    def parse_wildcards(dictionary, handle_parsing=False):

        # Dissassemble dictionary into key : value pairs
        replaced_all = []
        for key in dictionary.keys():

            # Get values from dictionary
            wildcarded_string = key  # Original string to be replaced

            # Check for extra selection algebra ^xxx
            # For example, path/***/***^2 : generator1
            # This tells the parser that generator1 is used to replace 2 wildcards with same values
            if '^' in key:
                extra_arguments = key.split('^')[1]
                extra_arguments = [int(i) for i in extra_arguments.split(",")]

                wildcarded_string = key.split('^')[0]
            else:
                extra_arguments = None

            # If handles, no dataset argument is passed in the last place in dict value
            if handle_parsing:
                replacement_generators = [replacement for replacement in dictionary[
                    key]]  # Lists/tuples/generators of things to replace wildcards with
            else:
                replacement_generators = [replacement for replacement in dictionary[key][
                                                                         :-1]]  # Lists/tuples/generators of things to replace wildcards with
                replacement_identifiers = dictionary[key][
                    -1]  # This identifies what replacements are supposed to mean (i.e. dataset type)

            # Check if number of wildcard matches the number of replacement generators
            if wildcarded_string.count('***') != len(
                    replacement_generators) and not handle_parsing and not extra_arguments:
                exit('NUMBER OF WILDCARDS DOES NOT MATCH THE NUMBER OF VALUES')

            # Now replace wildcards one by one
            # If dictionary values consist of iterable lists/tuples/generators i.e. it's nested such as ((2,4,6), (xxx))
            if proteinModel.is_iterable(replacement_generators[0]):

                for iteration, replacement_generator in enumerate(replacement_generators):
                    if iteration == 0:
                        replaced = [wildcarded_string]

                    replaced = proteinModel.replace_wildcard(replaced, replacement_generator, extra_arguments)
                    if extra_arguments:
                        extra_arguments.pop(0)


            # If it's just one list/tuple of numbers/strings for example (2,6,8)
            else:
                replaced = proteinModel.replace_wildcard([wildcarded_string], replacement_generators, extra_arguments)
                if extra_arguments:
                    extra_arguments.pop(0)

            replaced_all.append(replaced)

            # replaced_all = flatten(replaced_all)

        return (replaced_all)

    @staticmethod
    def replace_wildcard(string_holder, replacement_generator, occurences):
        replacement_holder = []
        if occurences == None:
            occurences = [1]

        for string in string_holder:
            if proteinModel.is_iterable(replacement_generator):
                replacements = [string.replace('***', str(single_replacement), occurences[0]) for single_replacement in
                                replacement_generator]
            else:
                replacements = [string.replace('***', str(replacement_generator), occurences[0])]
            replacement_holder.append(replacements)

            # Check if list is nested, if so, flatten it
            replacement_holder = proteinModel.flatten(replacement_holder)
        return (replacement_holder)

    @staticmethod
    def flatten(list_of_lists):
        if len(list_of_lists) == 0:
            return list_of_lists
        if isinstance(list_of_lists[0], list):
            return proteinModel.flatten(list_of_lists[0]) + proteinModel.flatten(list_of_lists[1:])
        return list_of_lists[:1] + proteinModel.flatten(list_of_lists[1:])

    @staticmethod
    def quick_load(arguments, handles):

        # Get all the neccessary values
        names = proteinModel.flatten(proteinModel.parse_wildcards(handles, handle_parsing=True))
        paths = proteinModel.parse_wildcards(arguments)
        dataset_names = [arguments[key][-1] for key in arguments]


        models = []



        # Generate every model according to handles/names
        for i, name in enumerate(names):

            new_model = proteinModel(annotation = name)
            print(f'{name.upper()} \n |')

            # Zip through dataset names (i.e. (('rmsd', 'mpd', 'rg'), 'gmmmpbsa', 'gmxmmpbsa') and paths
            # This assigns dataset types to paths
            for dataset_name, path in zip(dataset_names, paths):


                # Some dataset names can be lists/tuples, so iterate through them
                if proteinModel.is_iterable(dataset_name):

                    # For each name in a list or tuple, load the path into the model and remove the path from paths
                    for sub_name in dataset_name:

                        try:

                            new_model.get_dataset(path[0], sub_name)
                            #print(f' -----> LOADED {sub_name} FROM {path[0]}')

                        except:
                            #print(f' xxxxxx FILE {sub_name} PATH {path[0]} NOT FOUND')
                            pass
                        path.pop(0)


                # Same, but for cases where name is a string or an int (no real reason for it to be a float)
                elif type(dataset_name) is str or type(dataset_name) is int:

                    try:
                        new_model.get_dataset(path[0], dataset_name)
                        #print(f' -----> LOADED {dataset_name} FROM {path[0]}')
                    except:
                        pass
                        #print(f' XXXXXX FILE {dataset_name} PATH {path[0]} NOT FOUND')
                    path.pop(0)

            models.append(new_model)

        return(models)

    def set_special_residues(self, marked_residues):
        self.spec_residues = marked_residues

    def save(self, path):
        with open(path,'wb') as fileObject:
            fileObject.write(pickle.dumps(self.__dict__))

    def __load(self, path):

        with open(path, 'rb') as fileObject:
            data_pickle = fileObject.read()
            self.__dict__ = pickle.loads(data_pickle)
            print(f'LOADING INSTANCE NAME {self.annotation} FROM {path}')

    def DSSP_assign(self):


        with open('_temp_pdb.pdb', 'w') as temp_pdb_file:
            temp_pdb_file.write(self.pdb_file)

            p = PDBParser()
            structure = p.get_structure('MOD', '_temp_pdb.pdb')
            model = structure[0]
            dssp = DSSP(model, '_temp_pdb.pdb')
            os.remove('_temp_pdb.pdb')

        ss = []
        for i, residue in enumerate(self.seq_3):
            try:
                a_key = list(dssp.keys())[i]
                ss.append(dssp[a_key][2])
            except:
                ss.append('x')

        self.datasets['dssp'] = ss

    @describeLoad
    def get_dataset(self, path, dataset_name):

        # This finds out which dataset loading function to use according to data type
        dataset_sorting_dict = {
            'rmsd'                  : self.__get_simple_dataset,
            'rmsf'                  : self.__get_simple_dataset,
            'mpd'                   : self.__get_simple_dataset,
            'rg'                    : self.__get_simple_dataset,
            'gmmpbsa'               : self.__get_gmmpbsa_dataset,
            'gmxmmpbsa'             : self.__get_gmxmmpbsa_dataset,
            'umbrella_histogram'    : self.__get_umbrella_histogram,
            'umbrella_profile'      : self.__get_umbrella_profile,
            'contacts'              : self.__get_simple_dataset,
            'sequence'              : self.__get_sequence_pdb,
            'total_IEM'             : self.__get_interaction_energy_matrix,
            'elec_IEM'              : self.__get_interaction_energy_matrix,
            'vdw_IEM'               : self.__get_interaction_energy_matrix,
            'total_pairwise_IEM'    : self.__get_pairwise_matrix,
            'elec_pairwise_IEM'     : self.__get_pairwise_matrix,
            'vdw_pairwise_IEM'      : self.__get_pairwise_matrix,
            'total_energies_list'   : self.__get_energy_per_residue_dataset
        }


        dataset_sorting_dict[dataset_name](path, dataset_name)

        #return(dataset_sorting_dict[dataset_name])




    def get_GRINN_datasets(self, path, dataset_names = []):
        # Add slash if not provided at the end
        if path[-1] != '/':
            path += '/'

        wildcarded_paths = {f'{path}***' :
                                (('energies_intEnMeanTotal.dat','energies_intEnTotal.csv','energies_intEnMeanTotalList.dat',
                                  'energies_intEnMeanElec.dat', 'energies_intEnElec.csv',
                                  'energies_intEnMeanVdW.dat', 'energies_intEnVdW.csv'),
                                 ('total_IEM','total_pairwise_IEM','total_energies_list',
                                  'elec_IEM', 'elec_pairwise_IEM',
                                  'vdw_IEM', 'vdw_pairwise_IEM'))}
        paths = proteinModel.flatten(proteinModel.parse_wildcards(wildcarded_paths))
        dataset_names = wildcarded_paths[f'{path}***'][-1]

        for path, dataset_name in zip(paths, dataset_names):
            self.get_dataset(path,dataset_name)

    def call_pymol(self):
        with open('_temp_pdb.pdb', 'w') as temp_pdb_file:
            temp_pdb_file.write(self.pdb_file)
            os.system('pymol _temp_pdb.pdb')
            os.remove('_temp_pdb.pdb')

    def __get_energy_per_residue_dataset(self, path, dataset_name):
        with open(path) as file:
            lines = file.readlines()
        # Obtain list of unique residues (first column from the file)
        unique_residues = []
        for line in lines:
            residue_candidate = line.split()[0]
            holder = residue_candidate if residue_candidate not in unique_residues else None
            if holder:
                unique_residues.append(holder)
        # Calculate total interaction energies for each unique residue
        residue_dict = {}
        for unique_residue in unique_residues:
            holder = []

            for line in lines:
                if unique_residue == line.split()[0]:
                    holder.append(float(line.split()[2]))
            residue_dict[unique_residue] = sum(holder)
        self.datasets[dataset_name] = residue_dict

    def __get_pairwise_matrix(self, path, dataset_name):
        dataframe = pd.read_csv(path)
        self.datasets[dataset_name] = dataframe

    def __get_interaction_energy_matrix(self, path, dataset_name):

        with open(path) as file:
            lines = file.readlines()
            lines_split = [line.split() for line in lines]
            df = pd.DataFrame(lines_split).astype(float)
        self.datasets[dataset_name] = df

    def __get_simple_dataset(self, path, dataset_name):
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

            self.datasets[dataset_name] = [x, y, info_lines, axis_lines]

    def __get_gmmpbsa_dataset(self, path, dataset_name):
        with open(path) as file:
            read = file.readlines()
            dict_lines = {}
            for line in read:
                if "van der Waal energy" in line:
                    indexer = "VDW"
                elif "Electrostattic energy" in line:
                    indexer = "ELE"
                elif "Polar Solvation energy" in line:
                    indexer = "POL"
                elif "SASA energy" in line:
                    indexer = "SASA"
                elif "Binding energy" in line:
                    indexer = "binding"
                else:
                    indexer = False
                if indexer:
                    dict_lines[indexer] = line.split("=")[1].split("+/-")[0].replace(" ", "")
                    dict_lines[indexer + "e"] = line.split("=")[1].split("+/-")[1].replace("kJ/mol\n", "").replace(" ", "")
        self.datasets[dataset_name] = dict_lines

    def __get_gmxmmpbsa_dataset(self, path, dataset_name):

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
                elif 'Differences' in line:
                    reading_diff = True

                if reading_pb and reading_diff and 'DELTA TOTAL' in line:
                    pb_total = line.split()[2]
                    pb_total_std = line.split()[3]
                    reading_pb = False
                    reading_diff = False
                elif reading_gb and reading_diff and 'DELTA TOTAL' in line:
                    gb_total = line.split()[2]
                    gb_total_std = line.split()[3]
                    reading_gb = False
                    reading_diff = False
        self.datasets[f'{dataset_name}_pb'] = {'binding': pb_total, 'bindinge': pb_total_std}
        self.datasets[f'{dataset_name}_gb'] = {'binding': gb_total, 'bindinge': gb_total_std}

    def __get_sequence_pdb(self, path, dataset_name):

        def three_to_one(sequence):
            assignment_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

            modified_sequence = []
            for residue in sequence:
                modified_sequence.append(residue[0] + assignment_dict[residue[1:4]] + residue[4:])
            return(modified_sequence)

        with open(path) as pdb_file:


            residues = []
            pdb_lines = pdb_file.readlines()
            remember_last_resindex = 0
            for line in pdb_lines:
                split_line = line.split()

                if 'ATOM' in split_line[0] and remember_last_resindex != split_line[5]:
                    residues.append(split_line[4] + split_line[3] + split_line[5])
                    remember_last_resindex = split_line[5]
            pdb_file.seek(0)
            self.pdb_file = pdb_file.read()


        self.seq_3_chains = residues
        self.seq_3 = [residue[1:] for residue in residues]

        modified_sequence = three_to_one(self.seq_3_chains)
        self.seq_1_chains = modified_sequence
        self.seq_1 = [residue[1:] for residue in modified_sequence]

    def __get_umbrella_histogram(self, path, dataset_name):

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


        self.datasets[dataset_name] = [histo_x, run_y_values]

    def __get_umbrella_profile(self, path, dataset_name):
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

        self.datasets[dataset_name] = [profile_x, profile_y, std_error]



"""
To do:
Manual for methods and classes so it can be used properly from terminal

Add possibility to get energy matrices (vdw, elec, total), sequences,

Maybe add a possibility of getting other stuff from gmx_MMPBSA files, not just the summary

Comparator class that can visualize dataset of input models or a single model
4) Grinn heatmap

Convert kJ to kcal and v/v in plots

Pass through of arguments from quick_loader to proteinModels (like splits etc.)

Heatmapper class to visualize IEMs

Better loading of GRINN files, just add the folder path and it will load everything

Saving models to a database/file, then read from it later

Ultimate goal: BASH or GUI access to saved models, quick visualization. Use BASH or GUI to load new datasets etc.
"""
