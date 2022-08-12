"""
Prepare and analyze with new cartesian

### start with .xvg files
for i in {1..29}; do
cartesian.py --f run_${i}.xvg --o pdz_${i}.cart # ~8 min to finish (old one took 2 m per vector...)
done

### analyse vectors one by one
for i in {1..19}; do mkdir pdz_${i}; cartesian_ana.py --f pdz_${i}.cart --method all --p False --o pdz_${i}; done
for i in {1..19}; do mkdir fd4_${i}; cartesian_ana.py --f fd4_${i}.cart --method all --p False --o fd4_${i}; done



"""

"""
EXAMPLE USAGE

cartesian_batch.py --f outputs.txt --id "pdz,fd4" --pval 0.01 --pdbs PDB.pdb --pdbshow 2 --o output_dir --stripplot False
"""

import sys, argparse, os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sb
from scipy.stats import mannwhitneyu as mwu
from cartesian import __prepare_matplotlib
from cartesian_ana import write_to_pdb_beta, str2bool
import subprocess
import re
__prepare_matplotlib()

def parse_paths(output_file):
    # file.write(f'{args.f},{args.s},{args.o},{diffname},{confname},{grids1name},{grids2name},{grids1count},{grids2count}\n')

    with open(output_file) as file:
        lines=file.readlines()
        lines=[line.strip('\n').strip('.json') for line in lines]
    keys1=[]
    keys2=[]
    diff_paths=[]
    rms_out_paths=[]

    grids1_paths=[]
    grids2_paths=[]
    grids1_count_paths=[]
    grids2_count_paths=[]

    # Case 1 is when cartesian_ana.py was used to do pairwise comparisons (i.e. using both --f and --s flags, then data for both trajectories are on a single line)
    for line in lines:
        try:
            split_line = line.split(',')

            keys1.append(split_line[0])
            keys2.append(split_line[1])
            dir = split_line[2]
            diff_paths.append(f'{dir}/{split_line[3]}')
            rms_out_paths.append(f'{dir}/{split_line[4]}')
            grids1_paths.append(f'{dir}/{split_line[5]}')
            grids2_paths.append(f'{dir}/{split_line[6]}')
            grids1_count_paths.append(f'{dir}/{split_line[7]}')
            grids2_count_paths.append(f'{dir}/{split_line[8]}')

            grids1_paths.append(grids1_paths)
            grids2_paths.append(grids2_paths)
            grids1_count_paths.append(grids1_count_paths)
            grids2_count_paths.append(grids2_count_paths)
        except:
            print(f'One of the lines in {output_file} could not be parsed')
    # Case 2 is when cartesian_ana.py was used to only prepare datasets for cartesian_batch.py using only a --f flag, in this case the datapaths are split among blocks
    # in outputs.txt and we need twice the amount of .csv files. This is handled ok by the original cartesian_batch parser



    return(keys1, keys2, diff_paths, rms_out_paths, grids1_paths, grids2_paths, grids1_count_paths, grids2_count_paths)

def aggregate_volumes(paths): # add option to read two sets from two .csv
    # Explored volumes are not pairwise, only grab unique values



    df = pd.read_csv(paths[0])

    for path in paths[1:]:
        try:
            new_df = pd.read_csv(path, float_precision='high')
            df = pd.concat([df, new_df], axis=1) #pd.join fails if there's two columns of the same name?
        except:
            pass

    df = df.loc[:,~df.columns.duplicated()] # remove duplicates, StackOverflow solution
    df = df.drop(columns=df.columns[0]) # drop first column
    df = df[df.columns.drop(list(df.filter(regex='SUMV')))] #drop SUMV column
    df = df.T.drop_duplicates(keep='first').T # Transpose, drop non-unique rows, transpose back


    # Rename first column >> change later, improve in cartesian_ana.py
    df = df.replace({'atom ': ''}, regex=True)
    names = df.columns.tolist()
    names[0] = 'atom'
    df.columns = names
    ##


    #volumes_traj1 = pd.DataFrame(df.iloc[:, 0]) # new dataframes, insert atom ID column
    volumes_traj1 = pd.DataFrame() #df.iloc[:, 0]
    volsteps_traj1 = pd.DataFrame()
    volumes_traj2 = pd.DataFrame()  # new dataframes, insert atom ID column
    volsteps_traj2 = pd.DataFrame()



    for column in df:
        if '/step' in column:
            if args.id[0] in column:
                volsteps_traj1 = pd.concat([volsteps_traj1, df[column]], axis=1)
            elif args.id[1] in column:
                volsteps_traj2 = pd.concat([volsteps_traj2, df[column]], axis=1)

        elif 'V' in column:
            if args.id[0] in column:
                volumes_traj1 = pd.concat([volumes_traj1, df[column]], axis=1)
            if args.id[1] in column:
                volumes_traj2 = pd.concat([volumes_traj2, df[column]], axis=1)

    volumes_traj1 = volumes_traj1.astype(np.float64)
    volumes_traj2 = volumes_traj2.astype(np.float64)
    volsteps_traj1 = volsteps_traj1.astype(np.float64)
    volsteps_traj2 = volsteps_traj2.astype(np.float64)




    return(volumes_traj1, volumes_traj2, volsteps_traj1, volsteps_traj2)

def build_total_df(grids, counts):
    """
    Takes aggregate dataframes of grids and counts and makes one total dataframe.
    Unique grids are identified from unique grids and their counts. The counts are then summed up
    in corresponding grids.

    :param grids:
    :param counts:
    :return:
    """



    total_grids = {}
    total_counts = {}

    triplize = lambda tr: f'{tr[0]},{tr[1]},{tr[2]}'
    detriplize = lambda tr: tuple(np.array(tr.split(',')).astype(np.int32))

    # Multiprocess this
    for column in grids.columns:
        grid_col = grids[column].dropna().apply(triplize)
        count_col = counts[column].dropna()
        grid_count_col = pd.concat([grid_col, count_col], axis=1) # make a two column DF, one column grids, one columns counts
        grid_count_col.columns = ('grid', 'count') # group-by columns, two identical grid rows will be below each other
        total_column = grid_count_col.groupby(['grid']).sum()  # sum them up, obtaining an aggregate dataset of unique columns with correct populations
        total_grids[column], total_counts[column] = pd.Series(total_column.index).apply(detriplize), total_column['count'].values.astype(np.int32)



    total_grids_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in total_grids.items() ]))
    total_counts_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in total_counts.items() ]))

    return(total_grids_df, total_counts_df)

def aggregate_conformation(paths):
    """
    Reads grid and count dataframes from parquet and aggregate them into concatenated dataframes.
    Handles two cases

    :param paths:
    :return:
    """

    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    expression = "[\w-]+?((?=\.)|(?=$))"

    for path in paths:

        try:
            filename = re.search(expression, str(path)).group(0)

            new_df = pd.read_parquet(path)
            if args.id[0] in filename:
                print(f'{path} ID={args.id[0]}')
                df1 = pd.concat([new_df, df1], axis=0, ignore_index=True)
            elif args.id[1] in filename:
                print(f'{path} ID={args.id[1]}')
                df2 = pd.concat([new_df, df2], axis=0, ignore_index=True)
            else:
                print(f"Can't identify a trajectory to which dataset {path} belongs")

        except:
            pass

    return(df1, df2)

def parse_identifiers(string):
    if ',' in string:
        identifiers=string.split(',')
    else:
        raise(argparse.ArgumentTypeError('Identifiers need to be of a form "ID1,ID2" separated by a colon'))
    return(identifiers)

# not used .v.
def volume_df_describe(df):

    describe = pd.DataFrame()

    for row in df.iterrows():
        clean_row = row[1][1:].astype(float)
        describe_row = clean_row.describe()
        describe = pd.concat([describe, describe_row], axis=1)


    atom_ids = [i for i in range(1, len(describe.T)+1)] # Improve later
    atom_ids_series = pd.DataFrame(atom_ids, columns=['atom'])
    describe = pd.concat([describe.T, atom_ids_series], axis=1) # atom indexing starts from 1

    return(describe)

def mann_whitney_u_volumes(df1, df2):
    """
    Performs Mann Whitney U testing on the two provided datasets. This is equivalent to a two sample t-test
    for datasets not following normal distribution.
    Assumptions: Dependent variable ordinal/continuous, independent variable should be two independent categorical groups, independent observations (no relationship between the two groups or within each group.
    Null-Hypothesis: The distributions of both samples are the same

    :param df1:
        atom            0   1
        V(run1_vect1)   a   b
        V(run2_vect1)   c   d

    :param df2:
        atom            0   1
        V(run1_vect2)   e   f
        V(run2_vect2)   g   h
    :return:
    """
    df1, df2 = df1.T, df2.T # Each atom has its own column, observations (simulations) are rows

    atom_index, statistics, p_vals = [], [], []
    for column1 in df1:
        statistic, p_val = mwu(df1[column1], df2[column1])
        atom_index.append(column1)
        statistics.append(statistic)
        p_vals.append(p_val)

    return(atom_index, statistics, p_vals)

def drop_above_pval(mwu_pvals, mwu_scores, volsteps_traj1, volsteps_traj2):
    """
    :param mwu_pvals:
    :param volsteps_traj1:
    :param volsteps_traj2:
    :return: Three dataframes. 1st+2nd are dataframes containing observations only for atoms where perturbation was detected (below pval),
    3rd dataframe is a binary perturbation-1/0 with all the atoms (detected or not) which is useful for coloring pdb structures
    """

    # Add p-vals to dataframes of volsteps
    mwu_pvals_df = pd.DataFrame(mwu_pvals, columns=['pvals'])
    volsteps_traj1_pvals = pd.concat([volsteps_traj1, mwu_pvals_df], axis=1)
    volsteps_traj2_pvals = pd.concat([volsteps_traj2, mwu_pvals_df], axis=1)

    # Create new DF of yes/no values (1,0) - was the atom detected as perturbed
    delta_dict = {}
    for atom, pval in zip(volsteps_traj1_pvals.index, volsteps_traj1_pvals['pvals']):
        new_val = 1 if pval<=args.pval else 0
        delta_dict[atom] = new_val
    delta_df = pd.Series(delta_dict)

    # Create new DFs containing only atoms with p_val below the limit
    volsteps_traj1_pvals.drop(volsteps_traj1_pvals.index[volsteps_traj1_pvals['pvals'] >= args.pval],
                              inplace=True, axis=0)  # Drop any row where the p_value is below the limit
    volsteps_traj2_pvals.drop(volsteps_traj2_pvals.index[volsteps_traj2_pvals['pvals'] >= args.pval],
                              inplace=True, axis=0)  # Drop any row where the p_value is below the limit

    return(volsteps_traj1_pvals, volsteps_traj2_pvals, delta_df)

def plot_dynamical_distributions(volsteps_traj1_pvals, volsteps_traj2_pvals, stacking=True):
    import time
    start = time.time()

    if stacking:

        # Add identifiers (which trajectory)
        high1 = ((len(volsteps_traj1_pvals.columns)-1)*len(volsteps_traj1_pvals.index))
        high2 = ((len(volsteps_traj2_pvals.columns)-1)*len(volsteps_traj2_pvals.index))
        identifier1 = [args.id[0] for i in range(0, high1)]
        identifier2 = [args.id[1] for i in range(0, high2)]
        identifier = identifier1 + identifier2
        index1 = list(volsteps_traj1_pvals.columns)
        index2 = list(volsteps_traj2_pvals.columns)
        identifier_df = pd.DataFrame(identifier, columns=['identifier'])
        #####

        # Add atom indices (which atom)
        atoms1 = list(volsteps_traj1_pvals.index)
        atoms2 = list(volsteps_traj1_pvals.index)
        all_atoms = []
        for i in range(len(volsteps_traj1_pvals.columns)-1):

            all_atoms.append(atoms1)
        for i in range(len(volsteps_traj2_pvals.columns)-1):
            all_atoms.append(atoms2)

        all_atoms = [item for sublist in all_atoms for item in sublist]
        all_atoms_df = pd.DataFrame(all_atoms, columns=['atom'])
        #####

        volsteps_total = pd.concat([volsteps_traj1_pvals.iloc[:,0:-1], volsteps_traj2_pvals.iloc[:, 0:-1]], axis=1)
        volsteps_total = pd.melt(volsteps_total)

        volsteps_total = pd.concat([volsteps_total, identifier_df, all_atoms_df], axis=1)

        ax1 = plt.subplot()
        sb.boxplot(data=volsteps_total, x='atom', y='value', ax=ax1, orient='v', hue='identifier', palette=['indianred','royalblue'])
        if args.stripplot:
            sb.stripplot(data=volsteps_total, x='atom', y='value', hue="identifier", edgecolor='gray', ax=ax1, palette=['indianred','royalblue'])
        handles, labels = ax1.get_legend_handles_labels()  # stripplot creates a lot of legend values, get all handles and labels of the ax
        ax1.legend(handles=[handles[0], handles[1]], labels=[labels[0], labels[1]])
        ax1.set_title(f'Explored volume per step distributions per-atom\nwhere p < {args.pval}')


    else:
        ax1 = plt.subplot()

        sb.boxplot(data=volsteps_traj1_pvals.T.iloc[0:-1, :], color='orange', ax=ax1, palette=['red','blue'])
        sb.stripplot(data=volsteps_traj1_pvals.T.iloc[0:-1, :], color="orange", edgecolor='gray', ax=ax1, label='Batch 1', palette=['indianred','royalblue'])

        sb.boxplot(data=volsteps_traj2_pvals.T.iloc[0:-1, :], color='blue', ax=ax1, palette=['red','blue'])
        sb.stripplot(data=volsteps_traj2_pvals.T.iloc[0:-1, :], color="blue", edgecolor='gray', ax=ax1, label='Batch 2', palette=['indianred','royalblue'])

        handles, labels = ax1.get_legend_handles_labels() # stripplot creates a lot of legend values, get all handles and labels of the ax
        # Only use the first and the last handles and labels
        ax1.legend(handles=[handles[0], handles[-1]], labels=[labels[0], labels[-1]])

        ax1.set_title(f'Explored volume per step distributions per-atom\nwhere p < {args.pval}')

def pymol_dynamicity():
        kind = args.pdbshow

        if kind == 1: # Shows which atoms have changed with probability (depends on p-values)
            pymol_sub1 = 'cmd.label("perturbed","str(ID)")'  # add "name" for specific atom name
            pymol_sub2 = 'cmd.alter("perturbed", "vdw=0.6")'
            pymol_command = f"set orthoscopic, on; bg_color white; spectrum b,gray70_gray70_raspberry,minimum=0,maximum=1;select perturbed,b>0;set seq_view; show lines; {pymol_sub2}; show_as sticks cartoon sphere,perturbed;{pymol_sub1}; set cartoon_discrete_colors, on; set valence, 1; set label_shadow_mode, 2; set label_size,-0.6; set label_font_id,7; set label_outline_color, black; set label_color, white; set label_position,(0,0,2)"
            # Add something like this to cartesian_ana in the future
            p = subprocess.Popen(f"pymol {args.o}/dynamically_perturbed_atoms.pdb -d '{pymol_command}'", stdout=subprocess.PIPE, shell=True)
        if kind == 2: # Shows which atoms explored more volume in which trajectory, but only for those where p-value is lower than args.pval. Others are represented by B-score=0
            pymol_sub1 = 'select zero, b "=" 0'
            pymol_sub2 = 'cmd.alter("*", "vdw=0.6")'
            pymol_sub3 = f'cmd.label("more_in_{args.id[0]}","str(ID)")'
            pymol_sub4 = f'cmd.label("more_in_{args.id[1]}","str(ID)")'
            pymol_command = f"set orthoscopic, on; bg_color white; spectrum b, marine_gray70_raspberry; {pymol_sub1}; select more_in_{args.id[0]}, b > 0; select more_in_{args.id[1]}, b < 0; color gray70, zero; set seq_view; show lines; {pymol_sub2}; show_as sticks cartoon sphere,more_in_{args.id[0]};show_as sticks cartoon sphere,more_in_{args.id[1]};{pymol_sub3};{pymol_sub4}; set cartoon_discrete_colors, on; set valence, 1; set label_shadow_mode, 2; set label_size,-0.6; set label_font_id,7; set label_outline_color, black; set label_color, white; set label_position,(0,0,2)"
            p = subprocess.Popen(f"pymol {args.o}/dynamically_perturbed_atoms_rankings.pdb -d '{pymol_command}'", stdout=subprocess.PIPE, shell=True)

def pymol_conformation():
    minimum_beta = 0
    global maximum_beta_conf
    session = 'visualized_positions.pse'
    spectrum = 'gray70_pink_raspberry_purple'
    spectrum_ramp = '[gray70, pink, raspberry, purple]'
    sphere_select = ''

    pymol_sub1 = 'cmd.label("changed","str(ID)")'
    pymol_sub2 = 'cmd.alter("*", "vdw=0.6")'
    pymol_command = f"set orthoscopic, on; bg_color white; ramp_new colorbar, none, [{minimum_beta}, 0, {maximum_beta_conf}], {spectrum_ramp}; " \
                    f"spectrum b, {spectrum}, minimum={minimum_beta}, maximum={maximum_beta_conf}; {sphere_select} ;set seq_view; show lines; {pymol_sub2};" \
                    f"set cartoon_discrete_colors, on; set valence, 1; set label_shadow_mode, 2; set label_size,-0.6; set label_font_id,7; set label_outline_color, black; " \
                    f"set label_color, white; set label_position,(0,0,2); select changed, b>0.5 ; show_as sticks cartoon sphere, changed; set cartoon_discrete_colors, on; set valence, 1; set label_shadow_mode, 2;" \
                    f" set label_size,-0.6; set label_font_id,7; set label_outline_color, black; set label_color, white; {pymol_sub1}; save {args.o}/{session}"

    p = subprocess.Popen(f"pymol {args.o}/conformation_perturbed_atoms.pdb -d '{pymol_command}' &", stdout=subprocess.PIPE, shell=True, start_new_session=True)

def mwu_score_ranking(volsteps_df, mwu_scores, m, n, mwu_pvals):
    baseline = m*n*0.5 # mwu score of two exactly equal distributions

    # Prepare DataFrame with all atoms, pvals and score columns
    mwu_pvals_df = pd.DataFrame(mwu_pvals, columns=['pvals'])
    mwu_scores_df = pd.DataFrame(mwu_scores, columns=['scores'])
    total_df = pd.concat([volsteps_df, mwu_pvals_df, mwu_scores_df], axis=1)
    volsteps_df = total_df.copy(deep=True) # otherwise it assigns a new name for the same object
    ###

    # Create new DF of MWU score values - positive means TRAJ1 more explored volume, negative TRAJ2 more
    volsteps_df['scores'] = volsteps_df['scores']-baseline # subtract baseline value
    where_pvals_above_lim = volsteps_df.index[volsteps_df['pvals'] >= args.pval] # Which rows/atoms have p_vals above limit
    volsteps_df.loc[where_pvals_above_lim, 'scores'] = 0 # set scores to zero where above p_value limit
    #print(volsteps_df.loc[where_pvals_above_lim, 'scores'])


    delta_df = pd.Series(volsteps_df['scores'], index=volsteps_df.index)
    return(delta_df, total_df)

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=str, help='outputs.txt file with folder assignments', required=True, default='outputs.txt')
    parser.add_argument("--id", type=parse_identifiers, help='Identifier of the two datasets, format "ID1,ID2"', required=True)
    parser.add_argument("--pval", type=float, help="P-value limit for Mann-Whitney U testing distributions, default=0.01", required=False, default=0.01)
    parser.add_argument("--pdbs", type=str, help='OPTIONAL Input structure in .pdb format of the molecule', default=False)
    parser.add_argument("--pdbshow", type=int, help='OPTIONAL Creates a PyMol visualization of dynamicity change, requires --pdbs. Default=1: 0 None, 1 all atoms below p-val, 2 all atoms below p-val, also less/more',
                        default=1, choices=(0,1,2))
    parser.add_argument("--o", type=str, help='Output directory, defaults to names of trajectories used separated by _', required=False, default='.')
    parser.add_argument("--stripplot", type=str2bool, help='OPTIONAL plot a stripplot alongside dynamicity distributions, takes some time', required=False, default=False)

    global args
    args = parser.parse_args(argv)

    #
    if not args.o:
        args.o = f'{args.id[0]}_{args.id[1]}'
    if not os.path.exists(args.o):
        os.makedirs(args.o)


    # Handle data loading
    keys1, keys2, diff_paths, rms_out_paths, grids1_paths, grids2_paths, gridscount1_paths, gridscount2_paths = parse_paths(args.f)


    """
    ########################################
    # Testing for perturbation of dynamics #
    ########################################
    """

    # Handle dynamicity aggregation and statistical testing
    volumes_traj1, volumes_traj2, volsteps_traj1, volsteps_traj2 = aggregate_volumes(diff_paths) # ==> volumes_traj1.T, volumes_traj2.T into boxplots or violins, compare the two. Use volsteps datasets, independent on traj length


    mwu_atoms, mwu_scores, mwu_pvals = mann_whitney_u_volumes(volsteps_traj1, volsteps_traj2) # ==> results of the Mann-Whitney U test. (atom_n, score, p-value)

    # Plot atoms where the distribution changed significantly (p_val < args.pval) as boxplots
    volsteps_traj1_pvals, volsteps_traj2_pvals, only_perturbed_atoms = drop_above_pval(mwu_pvals, mwu_scores, volsteps_traj1, volsteps_traj2) # ==> two datasets only containing atoms proven to be perturbed, a dataframe prepared for b-coloring with all atoms and binary yes/no for perturbed

    m = len(volumes_traj1.columns)  # amount of samples m
    n = len(volumes_traj2.columns)  # amount of samples n

    pdb_kind1_header = f'REMARK B-Factor 1 means that the atom was identified as perturbed by comparing explored\n' \
                       f'REMARK volume distributions using MWU-test with a p-value limit of {args.pval}\n' \
                       f'REMARK ANALYZED {m}x{args.id[0]},{n}x{args.id[1]} \n'

    pdb_kind2_header = f'REMARK B-Factor 0 means that the atom was NOT identified as perturbed by comparing explored\n' \
                       f'REMARK volume distributions using MWU-test with a p-value limit of {args.pval}. Positive B-Factor\n' \
                       f'REMARK means that explored volume was larger for the -{args.id[0]}- batch of trajectories and vice-versa.\n' \
                       f'REMARK B-factors are MWU scores with subtracted baselines (m*n*0.5) that represent the score in\n' \
                       f'REMARK case of both distributions being identical.\n' \
                       f'REMARK ANALYZED {m}x{args.id[0]},{n}x{args.id[1]} \n'

    if not volsteps_traj1_pvals.empty and not volsteps_traj2_pvals.empty:

        plot_dynamical_distributions(volsteps_traj1_pvals, volsteps_traj2_pvals, stacking=True) # Stacked or unstacked boxplot combined with stripplot

        # Rank change of dynamicity, which variant more/less explored volume

        mwu_score_delta_df, total_dynamicity_df = mwu_score_ranking(volumes_traj1, mwu_scores, m, n, mwu_pvals) # only one dataframe needed, they both contain scores and atom numbers

        # Save relevant dynamicity datasets
        total_dynamicity_df.to_csv(f'{args.o}/{m}x{args.id[0]},{n}x{args.id[1]}.csv')

        """
        Use PyMol API later    
        """
        if args.pdbs:

            # Prepare .pdb for plot of kind=1 (show perturbed atoms)
            pdb_kind1 = write_to_pdb_beta(args.pdbs, only_perturbed_atoms)

            pdb_kind2 = write_to_pdb_beta(args.pdbs, mwu_score_delta_df)

            with open(f'{args.o}/dynamically_perturbed_atoms.pdb', 'w') as file:
                file.write(pdb_kind1_header + pdb_kind1)
            with open(f'{args.o}/dynamically_perturbed_atoms_rankings.pdb', 'w') as file:
                file.write(pdb_kind2_header + pdb_kind2)

    else:
        print(f'For explored volume analysis, no atoms passed the MWU test with given p-value {args.pval}')
        if args.pdbs:

            # Prepare .pdb for plot of kind=1 (show perturbed atoms)


            zeros_series = pd.Series([0 for i in range(len(volumes_traj1))])
            pdb_kind1 = write_to_pdb_beta(args.pdbs, zeros_series)
            pdb_kind2 = write_to_pdb_beta(args.pdbs, zeros_series)

            with open(f'{args.o}/dynamically_perturbed_atoms.pdb', 'w') as file:
                file.write(pdb_kind1_header + pdb_kind1)
            with open(f'{args.o}/dynamically_perturbed_atoms_rankings.pdb', 'w') as file:
                file.write(pdb_kind2_header + pdb_kind2)

    if args.pdbs and args.pdbshow:
        pymol_dynamicity()




    """
    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    @ Testing for perturbation of conformation #
    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    """

    import time
    start = time.time()

    grids1a, grids2a = aggregate_conformation(grids1_paths)
    counts1a, counts2a = aggregate_conformation(gridscount1_paths)

    print(f'Aggregate 1 {time.time() - start}')

    grids1b, grids2b = aggregate_conformation(grids2_paths) # This is for handling a second case from cartesian_ana.py, where --f and --s flags were both used and
    counts1b, counts2b = aggregate_conformation(gridscount2_paths) # grids are now in two different columns
    print(f'Aggregate 2 {time.time() - start}')

    grids1 = pd.concat([grids1a, grids1b], axis=0, ignore_index=True)
    grids2 = pd.concat([grids2a, grids2b], axis=0, ignore_index=True)

    counts1 = pd.concat([counts1a, counts1b], axis=0, ignore_index=True)
    counts2 = pd.concat([counts2a, counts2b], axis=0, ignore_index=True)
    print(f'Join 1&2 {time.time() - start}')

    ###
    total_grids1, total_counts1 = build_total_df(grids1, counts1)
    total_grids2, total_counts2 = build_total_df(grids2, counts2)
    print(f'Build total {time.time() - start}')
    ###

    from cartesian_ana import grid_rms
    rms, rms_norm = grid_rms(total_grids1, total_grids2, total_counts1, total_counts2, method='genius')

    int_rms1, int_rms_norm1 = grid_rms(total_grids1, total_grids1, total_counts1, total_counts1, method='genius')
    int_rms2, int_rms_norm2 = grid_rms(total_grids2, total_grids2, total_counts2, total_counts2, method='genius')
    rms = (2 * rms) - (int_rms1 + int_rms2)
    print(f'Calculate R {time.time() - start}')

    rms_out = pd.DataFrame(rms, columns=[f'{args.id[0]}/{args.id[1]}'])
    rms_norm = pd.DataFrame(rms_norm, columns=[f'{args.id[0]}/{args.id[1]}'])

    rms_out.plot()

    if args.pdbs:
        # Prepare .pdb for plot of kind=1 (show perturbed atoms)
        delta = pd.Series(rms_out.iloc[:, 0].values)
        global maximum_beta_conf
        maximum_beta_conf = delta.max()


        pdb_conf = write_to_pdb_beta(args.pdbs, delta)

        pdb_conf_header = f'REMARK B-Factor describes the R-score of an identical atom between two sets of trajectories\n' \
                           f'REMARK larger score means more conformational change on the atom\n' \
                           f'REMARK ANALYZED {m}x{args.id[0]},{n}x{args.id[1]} \n'

        with open(f'{args.o}/conformation_perturbed_atoms.pdb', 'w') as file:
            file.write(pdb_conf_header + pdb_conf)

    if args.pdbs and args.pdbshow:
        pymol_conformation()

    plt.show()



    #print(grids1.iloc[:, 0].dropna().to_list())
    #print(counts1.iloc[:, 0].dropna().to_list())


    # Pass args.grid (add as parameter to batch) to R-score calculation call in cartesian_ana


    ##

    """
    Outputs into documents
    volumes_describe = volume_distribution(volumes_traj1) # volume distributions
    
    """



if __name__ == "__main__":
    main()
    plt.show()