import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#plt.rc('xtick', labelsize=18)
#plt.rc('ytick', labelsize=18)
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 18}
SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
#plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)



x_label = 'Temperature [°C]'
y_label = r'F$_{350nm}$ / F$_{330nm}$'
y_der_label = r'1st derivative  F$_{350nm}$ / F$_{330nm}$'

# MLD033, MLD068, MLD108
ratio_file = pd.read_excel('/run/media/sulcjo/Data/termostabilita_plot/wt_semssels.xlsx', sheet_name=1)
derivative_file = pd.read_excel('/run/media/sulcjo/Data/termostabilita_plot/wt_semssels.xlsx', sheet_name=2)

times = ratio_file.iloc[:,0]
temperatures = ratio_file.iloc[:,1]
cap_fluorescences = ratio_file.iloc[:,2:]

der_cap_fluorescences = derivative_file.iloc[:,2:]

fig, axs = plt.subplots(nrows=2, ncols=2,figsize=(20,15))
for i, column in enumerate(cap_fluorescences.iloc[0,:]):
        if i < 2:
            axs[0][0].plot(temperatures.iloc[2:], cap_fluorescences.iloc[2:, i],label=column)
        #axs[0][0].plot(0,0,label=column,linewidth=5)




axs[1][0].plot(temperatures.iloc[2:], der_cap_fluorescences.iloc[2:], linewidth=3)
axs[0][0].legend()

axs[0][0].set(xlabel=x_label, ylabel=y_label)
axs[1][0].set(xlabel=x_label, ylabel=y_der_label)


# MLB036, MLB041
ratio_file = pd.read_excel('/run/media/sulcjo/Data/termostabilita_plot/158_semssels.xlsx', sheet_name=1)
derivative_file = pd.read_excel('/run/media/sulcjo/Data/termostabilita_plot/158_semssels.xlsx', sheet_name=2)

times = ratio_file.iloc[:,0]
temperatures = ratio_file.iloc[:,1]
cap_fluorescences = ratio_file.iloc[:,2:]

der_cap_fluorescences = derivative_file.iloc[:,2:]

for i, column in enumerate(cap_fluorescences.iloc[0,:]):

    axs[0][1].plot(temperatures.iloc[2:], cap_fluorescences.iloc[2:, i],label=column)
    #axs[0][1].plot(0, 0, label=column, linewidth=5, color=colors[i])


import matplotlib.ticker as plticker

loc = plticker.MultipleLocator(base=5.0) # this locator puts ticks at regular intervals
minor_loc = plticker.MultipleLocator(base=2.5)
def fmt_xaxes(ax):
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_minor_locator(minor_loc)


fmt_xaxes(axs[0][0])
fmt_xaxes(axs[0][1])
fmt_xaxes(axs[1][0])
fmt_xaxes(axs[1][1])


axs[1][1].plot(temperatures.iloc[2:], der_cap_fluorescences.iloc[2:], linewidth=3)
axs[0][1].legend()

axs[0][1].set(xlabel=x_label, ylabel=y_label)
axs[1][1].set(xlabel=x_label, ylabel=y_der_label)



plt.suptitle('DSF measurements of protein thermal stability',size=24)
plt.subplots_adjust(hspace=0.2,wspace=0.25,top=0.930,left=0.08,bottom=0.07,right=0.950)
#plt.savefig('/run/media/sulcjo/Data/BIOCEV/prometheusgrafypropema/dsf.png')
plt.show()


