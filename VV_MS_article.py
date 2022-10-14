import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import re

# Little helpers
def __prepare_matplotlib():
    import matplotlib.pyplot as plt
    # Prepare matplotlib parameters
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 16

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
__prepare_matplotlib()

def clean_sequence(line):
    if '(sub' in line: # Cases with substitutions, marked as XXXX(sub Y)XXXX(sub Z)XXX in the sequence
        list_line = list(line)
        indices = [i for i, x in enumerate(list_line) if x == "("] # Find indices of all ( brackets indicating substitutions or modifications

        for ind in indices:
            if list_line[ind+1] == 's' and list_line[ind+2] == 'u' and list_line[ind+3] == 'b': # Dirty, dirty, dirty
                list_line[ind-1] = list_line[ind+5] # Find the substitued AA before the bracket
                list_line[ind+1], list_line[ind+2], list_line[ind+3], list_line[ind+4], list_line[ind+5] = '0', '0', '0', '0', '0' # Ugly, but will get deleted by regex. Simple solutions work best :)
        line = ''.join(list_line)

    # Remove XXX(+42.36)XXX modifications from sequences
    pattern = re.compile(r"[a-zA-Z]+", re.IGNORECASE)
    matches = pattern.findall(line)

    return(''.join(matches))

def query_sequence(query, parent):

    output_dict = {'before': '', 'after': '', 'population': 0}

    # Query sequence can be present multiple times in the parental sequence
    start_ind_matches = [match.start() for match in re.finditer(query, parent)] # Starting indices for query matches
    parent_list = list(parent)
    output_dict['population'] = len(start_ind_matches)

    for i in start_ind_matches:
        try:
            start = f'{parent_list[i-1]}{i-1}_{parent_list[i]}{i}'
            output_dict['before'] = start
            if i-1 == -1:
                output_dict['idx_before'] = 0
            else:
                output_dict['idx_before'] = i # the space between aminoacids = hole
        except IndexError:
            output_dict['before'] = 'START'
        try:
            stop = f'{parent_list[i+len(query)]}{i+len(query)}_{parent_list[i+len(query)+1]}{i+len(query)+1}'
            output_dict['after'] = stop
            output_dict['idx_after'] = i+len(query)


        except IndexError:
            output_dict['after'] = 'END'
            output_dict['idx_after'] = len(parent)+1

    return(pd.Series(output_dict, name=query))

def number_sequence(sequence):
    new_sequence = ''
    sequence = list(sequence)
    for i, letter in enumerate(sequence):
        new_sequence += f'{letter}{i+1}_'
    # Remove last underscore _
    new_sequence = new_sequence[:-1]

    return(new_sequence)

def count_cleavage_sites(df):
    sites_count_before = pd.DataFrame(df['before'].value_counts())
    sites_count_after = pd.DataFrame(df['after'].value_counts())
    # Sum before and after
    sites_count = sites_count_before.join(sites_count_after, how='outer').sum(axis=1)
    return(pd.DataFrame(sites_count))

def prepare_dataset(path, negative_identifier = False, positive_identifier = True):
    df = pd.DataFrame(pd.read_csv(path)['Peptide'].apply(clean_sequence))
    df_id = pd.DataFrame(pd.read_csv(path)['Protein Accession'])
    df = pd.concat([df, df_id], axis=1)
    if negative_identifier:
        df = df[df['Protein Accession'].str.contains(negative_identifier) == False] # Needed because this dataset sucks, ISOI vs ISOII
    else:
        df = df[df['Protein Accession'].str.contains(positive_identifier) == True]
    df = df.drop('Protein Accession', axis=1)
    return(df)

def count_peptide_redundancy(df):
    df_a = df.groupby(['Peptide'], as_index=False).sum()
    df_b = df.groupby(['Peptide'], as_index=False).size()
    df = pd.concat([df_a, df_b], axis=1)
    df.columns = ['id', 'Peptide', 'pop']
    df = df.drop('id', axis=1)

    return(df)

def prepare_cleavage_df(df):
    cleavages_before = pd.concat([df['idx_before'], df['pop'], df['Peptide']], axis=1)
    cleavages_before.columns = ['site', 'pop', 'peptide']
    cleavages_after = pd.concat([df['idx_after'], df['pop'], df['Peptide']], axis=1)
    cleavages_after.columns = ['site', 'pop', 'peptide']

    cleavages = pd.concat([cleavages_before, cleavages_after], axis=0)
    cleavages = cleavages.sort_values('site')
    return(cleavages)


# Sequences
seq1 = 'GASVPFFPQQSGTPGMASLSLETMRQLGSLQRLNTLSQYSRYGFGKSFNSLWMHGLLPPHSSLPWMRPREHETQQYEYSLPVHPPPLPSQPSLKPQQPGLKPFLQSAAATTNQATALKEALQPPIHLGHLPLQEGELPLVQQQVAPSDKPPKPELPVDFADPQGPSLPGMDFPDPQGPSLPGLDFADPQGSTIFQIARLISHGPMPQNKQSPLYPGMLYVPFGANQLNAPVRLGIMSSEEVAGGREDPMAYGAMFPGFGGMRPGFEGMPHNPAMGGDFTLEFDSPVAATKGPENEEGGAQGSPMPEANPDNLENPAFLTELEPAPHAGLLALPKDDIPGLPRSPSGKMKGLPSVTPAAADPLMTPELADVYRTYDADMTTSVDFQEEATMDTTMAPNSLQTSMPGNKAQEPEMMHDAWHFQEPPREENLYFQ'
seq2 = 'GASVPFFPQQSGTPGMASLSLETMRQLGSLQRLNTLSQYSRYGFGKSFNSLWMHGLLPPHSSLPWMRPREHETQQPSLKPQQPGLKPFLQSAAATTNQATALKEALQPPIHLGHLPLQEGELPLVQQQVAPSDKPPKPELPVDFADPQGPSLPGMDFPDPQGPSLPGLDFADPQGSTIFQIARLISHGPMPQNKQSPLYPGMLYVPFGANQLNAPVRLGIMSSEEVAGGREDPMAYGAMFPGFGGMRPGFEGMPHNPAMGGDFTLEFDSPVAATKGPENEEGGAQGSPMPEANPDNLENPAFLTELEPAPHAGLLALPKDDIPGLPRSPSGKMKGLPSVTPAAADPLMTPELADVYRTYDADMTTSVDFQEEATMDTTMAPNSLQTSMPGNKAQEPEMMHDAWHFQEPPREENLYFQ'


# Load datasets from .csv and prepare dataframes with raw sub-sequences (dump non-sequence data), split into AMBN isoforms
csv1_path = '/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/protein-peptides ISO I+MMP20+KLK4.csv'
csv2_path = '/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/protein-peptides ISO II+MMP20+KLK4.csv'

csv1_df = prepare_dataset(csv1_path, negative_identifier='AMBN_IsoII')
csv2_df = prepare_dataset(csv2_path, positive_identifier='AMBN_IsoII')

# Return only unique peptides with populations
csv1_df = count_peptide_redundancy(csv1_df)
csv2_df = count_peptide_redundancy(csv2_df)

# Identify starting and ending aminoacids, see if peptide is in the parental sequence
csv1_matches = pd.DataFrame(csv1_df['Peptide'].apply(query_sequence, args=(seq1,)))
csv1 = pd.concat([csv1_df, csv1_matches], axis=1)
csv2_matches = pd.DataFrame(csv2_df['Peptide'].apply(query_sequence, args=(seq2,)))
csv2 = pd.concat([csv2_df, csv2_matches], axis=1)

# Count identified and missing queries
#missing1 = csv1[csv1['population'] < 1]
#missing2 = csv2[csv2['population'] < 1]
counts1 = (csv1['population'] < 1).value_counts()
counts2 = (csv2['population'] < 1).value_counts()
print(f'I) Identified: {counts1[False]}, II) Identified: {counts2[False]}')

# Drop peptides with 0 population
csv1 = csv1[csv1['population'] > 0] # Drop peptides with 0 population in sequence
csv2 = csv2[csv2['population'] > 0]
csv1 = csv1.drop('population', axis=1)
csv2 = csv2.drop('population', axis=1)

# Sort by N-terminal cleavage IDX
csv1 = csv1.sort_values('idx_before')
csv2 = csv2.sort_values('idx_before')

# ISOII has a missing sequence in the middle, realign with ISOI by adding +15 to every aminoacid after 75
realign_sequence = lambda x: x+15 if x > 75 else x
csv2['idx_before'] = csv2['idx_before'].apply(realign_sequence)
csv2['idx_after'] = csv2['idx_after'].apply(realign_sequence)

# Prepare cleavage datasets
cleavages1 = prepare_cleavage_df(csv1)
cleavages2 = prepare_cleavage_df(csv2)

# Plot alignment maps
fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True)
for j, (i, row) in enumerate(csv1.iterrows()):
    low, high = row[4], row[5]
    axs[0][0].barh(y=j, width=high-low, left=low, color='blue')
    axs[0][0].set_title('ISOI redundant alignments')

for j, (i, row) in enumerate(csv2.iterrows()):
    low, high = row[4], row[5]
    axs[0][1].barh(y=j, width=high-low, left=low, color='orange')
    axs[0][1].set_title('ISOII redundant alignments')



#print(sites1_plot.to_dict())

def calculate_populations(cleavages_df):
    sorted_df = cleavages_df.value_counts(sort=False)
    dictionary = sorted_df.to_dict()
    dictionary[(10000, 0, 'AAAA')] = 0 # Placeholder for last iteration

    outlist = []

    current_site = 0
    counter = 0
    peptides = []

    for line in dictionary:

        site = line[0]
        pop = line[1]
        peptide = line[2]

        if current_site != site and counter != 0: # Finished with site
            outlist.append([current_site, counter, peptides])
            current_site = site
            counter = 0
            peptides = []

            counter += pop
            peptides.append(peptide)

        else: # Currently iterating through site

            counter += pop
            peptides.append(peptide)


    df = pd.DataFrame(outlist, columns=['site', 'pop', 'peptides'])
    return(df)

sites1_plot = calculate_populations(cleavages1)
sites2_plot = calculate_populations(cleavages2)
print(sites1_plot)

axs[1][0].bar(sites1_plot['site'], sites1_plot['pop'], color='blue')
axs[1][1].bar(sites2_plot['site'], sites2_plot['pop'], color='orange')


plt.show()
exit()

"""
#
# Disregard populations, caused by MS methodology, nor in vivo reality
#


counts1 = []
for i in range(len(seq1)+1):
    i = float(i)
    _ = cleavages1.loc[cleavages1['site'] == i]
    _ = _['pop'].sum()
    counts1.append(_)

counts2 = []
for i in range(len(seq2)+1):
    i = float(i)
    _ = cleavages2.loc[cleavages2['site'] == i]
    _ = _['pop'].sum()
    counts2.append(_)

axs[1][0].bar([i+1 for i in range(len(counts1))], counts1, color='blue')
axs[1][1].bar([i+1 for i in range(len(counts2))], counts2, color='orange')
"""



#print(f'Seq1 cleavage sites: {list(sites1_plot.index)}')
#print(f'Seq2 cleavage sites: {list(sites2_plot.index)}')
cleavages1.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/cleavages_iso1.csv')
cleavages2.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/cleavages_iso2.csv')



fig, axs = plt.subplots()


###

probability1_temp_counts = np.zeros(len(seq1)+2, dtype=float) # 2 for either START or END
probability1_temp_df = pd.DataFrame([probability1_temp_counts])
probability2_temp_counts = np.zeros(len(seq2)+2, dtype=float)
probability2_temp_df = pd.DataFrame([probability2_temp_counts])


sites1 = probability1_temp_df.add(sites1_plot.T).T.fillna(0)
sites2 = probability2_temp_df.add(sites2_plot.T).T.fillna(0)



sites1_plot_probability = sites1.div(sites1.sum()) * 100
sites2_plot_probability = sites2.div(sites2.sum()) * 100

sites1_plot_probability.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/probabilities_iso1.csv')
sites1_plot_probability.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/probabilities_iso2.csv')

sites_delta = sites2_plot_probability - sites1_plot_probability

sites_delta.index = [int(x) for x in sites_delta.index]
sites_delta.plot(ax=axs, kind='bar')

plt.xticks([round(x) for x in sites_delta.index[::5]], rotation='vertical')
axs.set_xlabel('Residue after which cleavage occurs')
axs.set_ylabel(r'$\Delta$Probability / %')
axs.get_legend().remove()
#print(sites1_plot_probability.sum())
#print(sites2_plot_probability.sum())
#plt.show()

######################################################################################################
# Cleavage probability with considered redundant sequences
wredundant_sites1 = cleavages1.groupby(by='site').sum()
wredundant_sites2 = cleavages2.groupby(by='site').sum()
wredundant_sites1_plot_probability = wredundant_sites1.div(wredundant_sites1.sum()) * 100
wredundant_sites2_plot_probability = wredundant_sites2.div(wredundant_sites2.sum()) * 100

wredundant_sites1_plot_probability.columns = [0]
wredundant_sites2_plot_probability.columns = [0]

wredundant_sites1_plot_probability = wredundant_sites1_plot_probability.add(probability1_temp_df.T).fillna(0)
wredundant_sites2_plot_probability = wredundant_sites2_plot_probability.add(probability2_temp_df.T).fillna(0)

fig, axs = plt.subplots(ncols=2)
wredundant_sites1_plot_probability.plot(kind='bar', ax=axs[0])
wredundant_sites2_plot_probability.plot(kind='bar', ax=axs[1])
plt.show()
