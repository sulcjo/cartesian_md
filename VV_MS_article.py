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
        start = f'{parent_list[i-1]}{i-1}_{parent_list[i]}{i}'
        output_dict['before'] = start
        output_dict['idx_before'] = i # the space between aminoacids = hole
        try:
            stop = f'{parent_list[i+len(query)]}{i+len(query)}_{parent_list[i+len(query)+1]}{i+len(query)+1}'
            output_dict['after'] = stop

        except IndexError:
            output_dict['after'] = 'END'
        output_dict['idx_after'] = i+len(query)


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
pig  = 'SSAVPAFPRQPGTPGVASLSLETMRQLGSLQGLNMLSQYSRFGFGKSFNSLWMHGLLPPHSSFQWMRPREHETQQYEYSLPVHPPPLPSQPSLQPQQPGQKPFLQPTVVTSIQNPVQKGVPQPPIYQGHPPLQQVEGPMVQQQVAPSEKPPEAELPGLDFADPQDPSMFPIARLISQGPVPQDKPSPLYPGMFYMSYGANQLNSPARLGILSSEEMAGGRGGPLAYGAMFPGFGGMRPNLGGMPPNSAKGGDFTLEFDSPAAGTKGPEKGEGGAEGSPVAEANTADPESPALFSEVASGVLGGLLANPKGKIPNLARGPAGRSRGPPGVTPADADPLMTPGLADAYETYGADETTTLGLQEEMTMDSTATPYSEHTSMPGNKAQQPQIKRDAWRFQEP'
#pig_cut_nterm = 'MKDMVLILCLLKM'
#
#pig_cuts = [16, 18, 47, 146, 186, 238, 316, 335, 358] #Pig sequence has a deletion after 156th residues, add +10 after this resi
# also subtract 13 because of the cut N-terminus
#pig_cuts = [3, 5, 34, 133, 183, 235, 313, 332, 355] # MMP-20 aligned
#pig_cuts = [44, 59, 207, 233, 239, 240, 372, 417] # KLK4 unaligned
#pig_cuts = [31, 46, 204, 230, 236, 237, 369, 414] # KLK4 aligned
pig_cuts = [3, 5, 34, 133, 183, 235, 313, 332, 355, 31, 46, 204, 230, 236, 237, 369, 414] # Both aligned

# Load datasets from .csv and prepare dataframes with raw sub-sequences (dump non-sequence data), split into AMBN isoforms
#csv1_path = '/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/mmp20_klk4/protein-peptides ISO I+MMP20+KLK4.csv'
#csv2_path = '/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/mmp20_klk4/protein-peptides ISO II+MMP20+KLK4.csv'
csv1_path = '/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/mmp20/protein-peptides_ISO I.csv'
csv2_path = '/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/mmp20/protein-peptides_ISO II.csv'


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
    axs[0][0].set_ylabel('Peptide')

for j, (i, row) in enumerate(csv2.iterrows()):
    low, high = row[4], row[5]
    axs[0][1].barh(y=j, width=high-low, left=low, color='orange')
    axs[0][1].set_title('ISOII redundant alignments')
    axs[0][1].set_ylabel('Peptide')


#print(sites1_plot.to_dict())

def calculate_populations(cleavages_df, disregard_redundant = False):
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
            if disregard_redundant:
                counter += 1
            else:
                counter += pop
            peptides.append(peptide)


    df = pd.DataFrame(outlist, columns=['site', 'pop', 'peptides'])
    df['prob'] = (df['pop'] / df['pop'].sum()) * 100

    return(df)

sites1_plot = calculate_populations(cleavages1)
sites2_plot = calculate_populations(cleavages2)
#print(sites1_plot)

axs[1][0].bar(sites1_plot['site'], sites1_plot['prob'], color='blue')
axs[1][0].set_title(f'Redundant ISOI cleavage probabilities')
axs[1][0].set_xlabel(f'Aminoacid')
axs[1][0].set_ylabel('P / %')
axs[1][1].bar(sites2_plot['site'], sites2_plot['prob'], color='orange')
axs[1][1].set_title(f'Redundant ISOII cleavage probabilities')
axs[1][1].set_xlabel('Aminoacid after alignment to ISOI')
axs[1][1].set_ylabel('P / %')

# Save redundant table, but first return ISOII aminoacids to their original numbering
rerealign_sequence = lambda x: x-15 if x > 75 else x
sites2_plot['site'] = sites2_plot['site'].apply(rerealign_sequence)

sites1_plot.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/redundant_cleavages_prob_iso1.csv')
sites2_plot.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/redundant_cleavages__prob_iso2.csv')

sites2_plot['site'] = sites2_plot['site'].apply(realign_sequence)

#plt.show()

# Plot probabilities in differences between redundant ISOI and ISOII
# First create empty dataframes so all aminoacids have a row, even if its zero
def get_empty_df_probs(ln, left_df):
    occupied_sites = left_df['site'].to_list()

    sites = [i for i in range(ln) if i not in occupied_sites]
    pops = [0 for i in range(ln) if i not in occupied_sites]
    peptides = [[] for i in range(ln) if i not in occupied_sites]
    probs = [0 for i in range(ln) if i not in occupied_sites]
    df = pd.DataFrame([sites, pops, peptides, probs], index=['site', 'pop', 'peptides', 'prob']).T
    df.index = df['site']
    return(df)

# Use only the ISOI len, it's longer
sites1_plot = pd.concat([get_empty_df_probs(len(seq1)+2, sites1_plot), sites1_plot], axis=0).sort_values('site')
sites2_plot = pd.concat([get_empty_df_probs(len(seq1)+2, sites2_plot), sites2_plot], axis=0).sort_values('site')

delta_prob_np = sites2_plot['prob'].to_numpy() - sites1_plot['prob'].to_numpy() # Doesn't work with PD for some reason
fig, ax = plt.subplots()
ax.plot([i for i in range(len(delta_prob_np))], delta_prob_np)
ax.set_title('Redundant ISOI / ISOII cleavage differences')
ax.set_ylabel(r'$\Delta$P / %')
ax.set_xlabel('Aminoacid after alignment')
#plt.show()

# Redundant probabilities delta out
redundant_delta_prob_out = pd.DataFrame(delta_prob_np)
redundant_delta_prob_out = redundant_delta_prob_out.reset_index()
redundant_delta_prob_out.columns = ['site', 'p_delta']
redundant_delta_prob_out.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/redudant_cleavage_probability_deltas.csv')

# Joint plot for redundant datasets
fig, axj = plt.subplots(nrows=3, sharex=True)
for j, (i, row) in enumerate(csv1.iterrows()):
    low, high = row[4], row[5]
    axj[0].barh(y=j, width=high-low, left=low, color='blue', label='ISOI')
    axj[0].set_title('ISOI/ISOII redundant alignments')
    axj[0].set_ylabel('Peptide')

for j, (i, row) in enumerate(csv2.iterrows()):
    low, high = row[4], row[5]
    axj[0].barh(y=j, width=high-low, left=low, color='orange', label='ISOII')

csv1['id'] = ['ISOI' for i in csv1.iterrows()]
csv2['id'] = ['ISOII' for i in csv2.iterrows()]
csv_joint = pd.concat([csv1, csv2], axis=0)
import seaborn as sb
axj[1].bar(sites1_plot['site']-0.25, sites1_plot['prob'], color='blue', width=0.5, label='ISOI')
axj[1].bar(sites2_plot['site']+0.25, sites2_plot['prob'], color='orange', width=0.5, label='ISOII')
axj[2].plot([i for i in range(len(delta_prob_np))], delta_prob_np)

axj[1].legend()
axj[1].set_title('Probability of cleavage after aminoacid')
axj[1].set_ylabel('P / %')
axj[2].set_title('Differences in probabilities of cleavage')
axj[2].set_ylabel(r'$\Delta$P / %')
axj[2].set_xlabel('Aminoacid')
axj[0].set(xlim=(0,433))

# Are there any peaks in ISOI which are not present in ISOII?
x_new_peaks_1 = []
x_new_peaks_2 = []
for isoi, isoii in zip(sites1_plot.iterrows(), sites2_plot.iterrows()):

    if isoi[1]['prob'] > 0 and isoii[1]['prob'] == 0:
        x_new_peaks_1.append(isoi[1]['site'])
    if isoii[1]['prob'] > 0 and isoi[1]['prob'] == 0:
        x_new_peaks_2.append(isoi[1]['site'])

for new_peak_1, new_peak_2 in zip(x_new_peaks_1, x_new_peaks_2):
    axj[1].arrow(new_peak_1, -4, dx=0, dy=2, color='red', head_width=1)
    axj[1].arrow(new_peak_2, -4, dx=0, dy=2, color='green', head_width=1)

# Are there any peaks in ISOI which are not present in ISOII and pig AMBN at the same time?
x_new_peaks_1_depiged = []
x_new_peaks_2_depiged = []

for isoi, isoii in zip(sites1_plot.iterrows(), sites2_plot.iterrows()):

    if isoi[1]['prob'] > 0 and isoii[1]['prob'] == 0 and isoi[1]['site'] not in pig_cuts:
        x_new_peaks_1_depiged.append(isoi[1]['site'])
    if isoii[1]['prob'] > 0 and isoi[1]['prob'] == 0 and isoii[1]['site'] not in pig_cuts:
        x_new_peaks_2_depiged.append(isoi[1]['site'])

for new_peak_1, new_peak_2 in zip(x_new_peaks_1_depiged, x_new_peaks_2_depiged):
    axj[1].arrow(new_peak_1, -6, dx=0, dy=2, color='orange', head_width=1)
    axj[1].arrow(new_peak_2, -6, dx=0, dy=2, color='purple', head_width=1)

#print(sites1_plot['site'][sites1_plot['prob'] > 0])

print('ISOI:', x_new_peaks_1, len(x_new_peaks_1))
#print('ISOI-pig:', x_new_peaks_1_depiged, len(x_new_peaks_1_depiged))
print('ISOII:', x_new_peaks_2, len(x_new_peaks_2))
#print('ISOII-pig:', x_new_peaks_2_depiged, len(x_new_peaks_2_depiged))
#print('pig:', pig_cuts, len(pig_cuts))
#plt.show()





######
######
######
######
# Final table with thresholded sites
threshold = 1.5

print(f'ISOI nonuniq above {threshold}', len(sites1_plot['site'][sites1_plot['prob'] > threshold]))
print(f'ISOII nonuniq above {threshold}', len(sites2_plot['site'][sites2_plot['prob'] > threshold]))
sites1_thr = sites1_plot[sites1_plot['prob'] > threshold]
sites2_thr = sites2_plot[sites2_plot['prob'] > threshold]




iso1_thr_uniq = []
iso2_thr_uniq = []
for isoi, isoii in zip(sites1_plot.iterrows(), sites2_plot.iterrows()):

    if isoi[1]['prob'] > threshold and isoii[1]['prob'] == 0:
        iso1_thr_uniq.append(isoi[1])
    if isoii[1]['prob'] > threshold and isoi[1]['prob'] == 0:
        iso2_thr_uniq.append(isoii[1])

print('ISOI unq:', len(iso1_thr_uniq))
print('ISOII unq:', len(iso2_thr_uniq))

iso1_thr_uniq = pd.DataFrame(iso1_thr_uniq)
iso2_thr_uniq = pd.DataFrame(iso2_thr_uniq)

#print(sites1_thr)
#print(iso1_thr_uniq)


sites1_thr.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/mmp20_klk4/threshold/iso1_15_thr.csv')
sites2_thr.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/mmp20_klk4/threshold/iso2_15_thr.csv')
iso1_thr_uniq.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/mmp20_klk4/threshold/iso1_15_thr_uniq.csv')
iso2_thr_uniq.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/mmp20_klk4/threshold/iso2_15_thr_uniq.csv')
#
######
######
######
######
######
######





#
# Disregard populations, caused by MS methodology, not in vivo reality
#
sites1_plot_wr = calculate_populations(cleavages1, disregard_redundant=True)
sites2_plot_wr = calculate_populations(cleavages2, disregard_redundant=True)


fig, axs_wr = plt.subplots(ncols=2, sharex=True)
axs_wr[0].bar(sites1_plot_wr['site'], sites1_plot_wr['prob'], color='blue')
axs_wr[1].bar(sites2_plot_wr['site'], sites2_plot_wr['prob'], color='orange')

# Save redundant table, but first return ISOII aminoacids to their original numbering
sites2_plot_wr['site'] = sites2_plot_wr['site'].apply(rerealign_sequence)

sites1_plot_wr.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/nonredundant_cleavages_prob_iso1.csv')
sites2_plot_wr.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/nonredundant_cleavages__prob_iso2.csv')

sites2_plot_wr['site'] = sites2_plot_wr['site'].apply(realign_sequence)

# Use only the ISOI len, it's longer
sites1_plot_wr = pd.concat([get_empty_df_probs(len(seq1)+2, sites1_plot_wr), sites1_plot_wr], axis=0).sort_values('site')
sites2_plot_wr = pd.concat([get_empty_df_probs(len(seq1)+2, sites2_plot_wr), sites2_plot_wr], axis=0).sort_values('site')

delta_prob_np_wr = sites2_plot_wr['prob'].to_numpy() - sites1_plot_wr['prob'].to_numpy() # Doesn't work with PD for some reason
fig, ax_wr = plt.subplots()
ax_wr.plot([i for i in range(len(delta_prob_np_wr))], delta_prob_np_wr)
ax_wr.set_title('Non-redundant ISOI / ISOII cleavage differences')
ax_wr.set_ylabel(r'$\Delta$P / %')
ax_wr.set_xlabel('Aminoacid after alignment')
#plt.show()

# Non-edundant probabilities delta out
redundant_delta_prob_out_wr = pd.DataFrame(delta_prob_np_wr)
redundant_delta_prob_out_wr = redundant_delta_prob_out_wr.reset_index()
redundant_delta_prob_out_wr.columns = ['site', 'p_delta']
redundant_delta_prob_out_wr.to_csv('/run/media/sulcjo/sulcjo-data/IOCB/veronika/stepy_ms/align/nonredudant_cleavage_probability_deltas.csv')

plt.show()