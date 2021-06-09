import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


sec_files = pd.read_excel('/run/media/sulcjo/Data/BIOCEV/var92/SEC/300620 var92 SEC Akta/sec_data.xlsx')
empty_rows = 3
def clean_list(list):
    new_list = []

    for value in list:
        try:
            new_list.append(float(value))
        except:
            new_list.append(np.nan)
    return(new_list)

for column in sec_files.iteritems():

    column = list(column[1].values)
    if type(column) == list:
        if any('10_UV' in str(s) for s in column) and any('ml' in str(s) for s in column):
            list_vol = column
            list_vol = clean_list(list_vol)

        elif any('mAU' in str(s) for s in column):
            list_uv = column
            list_uv = clean_list(list_uv)

        elif any('_Fractions' in str(s) for s in column) and any('ml' in str(s) for s in column):
            list_fravol = column
            list_fravol = clean_list(list_fravol)
        elif any('(Fractions)' in str(s) for s in column):
            list_franum = column


plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
fig, axs = plt.subplots()
axs.plot(list_vol, list_uv, linewidth=2, color='black')
axs.set()
axs.set(xlim=(7,25))
axs.set_xlabel('V / ml', size=16)
axs.set_ylabel('A / mAU', size=16)
axs.set_title('mňaumesin', size=20)



vlines_x = [x for x in list_fravol if type(x) == float]
vlines_x = vlines_x[3:]
vlines_labels = [x for x in list_franum if type(x) != float and any(i.isdigit() for i in x)]
vlines_max = [x for x in list_uv if x is not np.nan]
vlines_max = max(vlines_max) + 1
axs.vlines(vlines_x, ymin = -0.1, ymax = vlines_max, color='red', linestyles='-', alpha=0.7)


even_odd = 1
for x, label in zip(vlines_x, vlines_labels):
    #label = "{:.2f}".format(label)
    #print(label)
    if even_odd %2 == 0:
        y_distance = -15
    else:
        y_distance = 0

    label = label.replace('"', '')
    plt.annotate(label, # this is the text
                 (x,vlines_max), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(12,y_distance), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 size=16)
    even_odd += 1
plt.show()
plt.savefig(f'{file_input}.png')




