"""
Prepare analyses of trajectories:

# Vectors should be named such as this with arbitrary order
<something><run_number_of_specific_variant><variant><something>.json
# You specify the <variant> identifiers so these need to be unique in the file name
###


mkdir pdz_fd4_batch
for i in {1..29}; do
    cp pdz/run_${i}_pdz_vectors.json fd4/run_${i}_fd4_vectors.json pdz_fd4_batch
done

cd pdz_fd4_batch

# Compare two trajectories replicas among each other (N^2 pairs)
for i in {1..29}; do
    for j in {1..29}; do
        echo $i $j
        cartesian_ana.py --f run_${i}_pdz_vectors.json --s run_${j}_fd4_vectors.json --plot False --plot_3d False --plot_positions False --plot_diff False --mp 5 --method all --grid 0.1 --plot_violin False --o out_${i}${j}
    done
done
"""

"""
EXAMPLE USAGE

cartesian_batch.py --f outputs.txt --id "pdz,fd4" --pval 0.01 --pdbs pdz_representative.pdb --pdbshow 2 --o output_dir

"""

import sys, argparse, os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sb
from scipy.stats import mannwhitneyu as mwu
from cartesian import __prepare_matplotlib
from cartesian_ana import write_to_pdb_beta
import subprocess
__prepare_matplotlib()

def parse_paths(output_file):
    # file.write(f'{args.f},{args.s},{args.o},{diffname},{confname},{grids1name},{grids2name},{grids1count},{grids2count}\n')

    with open(output_file) as file:
        lines=file.readlines()
        lines=[line.strip('\n').strip('.json') for line in lines]
    keys1=[]
    keys2=[]
    diff_paths=[]
    conf_paths=[]

    grids1_paths=[]
    grids2_paths=[]
    grids1_count_paths=[]
    grids2_count_paths=[]

    for line in lines:
        split_line = line.split(',')


        keys1.append(split_line[0])
        keys2.append(split_line[1])
        dir = split_line[2]
        diff_paths.append(f'{dir}/{split_line[3]}')
        conf_paths.append(f'{dir}/{split_line[4]}')
        grids1_paths.append(f'{dir}/{split_line[5]}')
        grids2_paths.append(f'{dir}/{split_line[6]}')
        grids1_count_paths.append(f'{dir}/{split_line[7]}')
        grids2_count_paths.append(f'{dir}/{split_line[8]}')

    return(keys1, keys2, diff_paths, conf_paths, grids1_paths, grids2_paths, grids1_count_paths, grids2_count_paths)

def aggregate_volumes(paths):
    # Explored volumes are not pairwise, only grab unique values
    df = pd.read_csv(paths[0])

    for path in paths[1:]:
        new_df = pd.read_csv(path, float_precision='high')
        df = pd.concat([df, new_df], axis=1) #pd.join fails if there's two columns of the same name?

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

def aggregate_conformation(paths):
    df = pd.read_csv(paths[0])
    for path in paths[1:]:
        new_df = pd.read_csv(path)
        df = pd.concat([df, new_df], axis=1) #pd.join fails if there's two columns of the same name?
    df = df.loc[:, ~df.columns.duplicated()]  # remove duplicates, StackOverflow solution
    df = df.drop(columns=df.columns[0])  # drop first column
    return(df)

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

        if stacking:

            # Add identifiers (which trajectory)
            high1 = ((len(volsteps_traj1_pvals.columns)-1)*len(volsteps_traj1_pvals.index))
            high2 = ((len(volsteps_traj2_pvals.columns)-1)*len(volsteps_traj2_pvals.index))
            identifier1 = [args.id[0] for i in range(0, high1)]
            identifier2 = [args.id[1] for i in range(0, high2)]
            identifier = identifier1 + identifier2

            index1 = list(volsteps_traj1_pvals.columns)
            index2 = list(volsteps_traj2_pvals.columns)
            index = index1[0:-1] + index2[0:-1]

            identifier_df = pd.DataFrame(identifier, columns=['identifier'])
            #####

            # Add atom indices (which atom)
            num_atoms = len(volsteps_traj1_pvals.index)
            atoms1 = list(volsteps_traj1_pvals.index)
            atoms2 = list(volsteps_traj1_pvals.index)
            all_atoms = []
            #print(volsteps_traj1_pvals)
            for i in range(len(volsteps_traj1_pvals.columns)-1):

                all_atoms.append(atoms1)
            for i in range(len(volsteps_traj2_pvals.columns)-1):
                all_atoms.append(atoms2)


            all_atoms = [item for sublist in all_atoms for item in sublist]
            all_atoms_df = pd.DataFrame(all_atoms, columns=['atom'])
            #####

            #volsteps_traj1_pvals = pd.concat([volsteps_traj1_pvals, identifier1_df], axis=0)
            #volsteps_traj2_pvals = pd.concat([volsteps_traj2_pvals, identifier2_df], axis=0)

            volsteps_total = pd.concat([volsteps_traj1_pvals.iloc[:,0:-1], volsteps_traj2_pvals.iloc[:, 0:-1]], axis=1)
            volsteps_total = pd.melt(volsteps_total)

            volsteps_total = pd.concat([volsteps_total, identifier_df, all_atoms_df], axis=1)
            #print(volsteps_total)

            #print(identifier_df)
            ax1 = plt.subplot()
            sb.boxplot(data=volsteps_total, x='atom', y='value', ax=ax1, orient='v', hue='identifier', palette=['indianred','royalblue'])
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
            p = subprocess.Popen(f"pymol dynamically_perturbed_atoms.pdb -d '{pymol_command}'", stdout=subprocess.PIPE, shell=True)
        if kind == 2: # Shows which atoms explored more volume in which trajectory, but only for those where p-value is lower than args.pval. Others are represented by B-score=0
            pymol_sub1 = 'select zero, b "=" 0'
            pymol_sub2 = 'cmd.alter("*", "vdw=0.6")'
            pymol_sub3 = f'cmd.label("more_in_{args.id[0]}","str(ID)")'
            pymol_sub4 = f'cmd.label("more_in_{args.id[1]}","str(ID)")'
            pymol_command = f"set orthoscopic, on; bg_color white; spectrum b, marine_gray70_raspberry; {pymol_sub1}; select more_in_{args.id[0]}, b > 0; select more_in_{args.id[1]}, b < 0; color gray70, zero; set seq_view; show lines; {pymol_sub2}; show_as sticks cartoon sphere,more_in_{args.id[0]};show_as sticks cartoon sphere,more_in_{args.id[1]};{pymol_sub3};{pymol_sub4}; set cartoon_discrete_colors, on; set valence, 1; set label_shadow_mode, 2; set label_size,-0.6; set label_font_id,7; set label_outline_color, black; set label_color, white; set label_position,(0,0,2)"
            p = subprocess.Popen(f"pymol dynamically_perturbed_atoms_rankings.pdb -d '{pymol_command}'", stdout=subprocess.PIPE, shell=True)

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
    parser.add_argument("--pdbshow", type=int, help='Creates a PyMol visualization of dynamicity change, requires --pdbs. Default=1: 0 None, 1 all atoms below p-val, 2 all atoms below p-val, also less/more',
                        default=1, choices=(0,1,2))
    parser.add_argument("--o", type=str, help='Output directory, defaults to names of trajectories used separated by _', required=False, default='.')

    global args
    args = parser.parse_args(argv)

    #
    if not args.o:
        args.o = f'{args.id[0]}_{args.id[1]}'
    if not os.path.exists(args.o):
        os.makedirs(args.o)



    # Handle data loading
    keys1, keys2, diff_paths, conf_paths, grids1, grids2, gridscount1, gridscount2 = parse_paths(args.f)


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
    plot_dynamical_distributions(volsteps_traj1_pvals, volsteps_traj2_pvals, stacking=True) # Stacked or unstacked boxplot combined with stripplot

    # Rank change of dynamicity, which variant more/less explored volume
    m = len(volumes_traj1.columns) # amount of samples m
    n = len(volumes_traj2.columns) # amount of samples n
    mwu_score_delta_df, total_dynamicity_df = mwu_score_ranking(volumes_traj1, mwu_scores, m, n, mwu_pvals) # only one dataframe needed, they both contain scores and atom numbers

    # Save relevant dynamicity datasets
    total_dynamicity_df.to_csv(f'{args.o}/{m}x{args.id[0]},{n}x{args.id[1]}.csv')

    """
    Use PyMol API later    
    """
    if args.pdbs:

        # Prepare .pdb for plot of kind=1 (show perturbed atoms)
        pdb_kind1 = write_to_pdb_beta(args.pdbs, only_perturbed_atoms)
        pdb_kind1_header = f'REMARK B-Factor 1 means that the atom was identified as perturbed by comparing explored\n' \
                           f'REMARK volume distributions using MWU-test with a p-value limit of {args.pval}\n' \
                           f'REMARK ANALYZED {m}x{args.id[0]},{n}x{args.id[1]} \n'
        pdb_kind2 = write_to_pdb_beta(args.pdbs, mwu_score_delta_df)
        pdb_kind2_header = f'REMARK B-Factor 0 means that the atom was NOT identified as perturbed by comparing explored\n' \
                           f'REMARK volume distributions using MWU-test with a p-value limit of {args.pval}. Positive B-Factor\n' \
                           f'REMARK means that explored volume was larger for the -{args.id[0]}- batch of trajectories and vice-versa.\n' \
                           f'REMARK B-factors are MWU scores with subtracted baselines (m*n*0.5) that represent the score in\n' \
                           f'REMARK case of both distributions being identical.\n' \
                           f'REMARK ANALYZED {m}x{args.id[0]},{n}x{args.id[1]} \n'
        with open(f'{args.o}/dynamically_perturbed_atoms.pdb', 'w') as file:
            file.write(pdb_kind1_header+pdb_kind1)
        with open(f'{args.o}/dynamically_perturbed_atoms_rankings.pdb', 'w') as file:
            file.write(pdb_kind2_header+pdb_kind2)

        if args.pdbshow:
            pymol_dynamicity()

    """
    #
    # Testing for perturbation of conformation
    #
    """




    conformations = aggregate_conformation(conf_paths)





    ##

    """
    Outputs into documents
    volumes_describe = volume_distribution(volumes_traj1) # volume distributions
    
    """



if __name__ == "__main__":
    main()
    plt.show()