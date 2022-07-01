import json
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson, chi2, gamma, norm, lognorm, beta, rayleigh, weibull_min, skewnorm

def fit_function(k, a, b, c):
    '''poisson function, parameter lamb is the fit parameter'''
    #return beta.pdf(k, a, b)
    return(skewnorm.pdf(k, a, b, c))



with open('/run/timeshift/backup/IOCB/cartesian/rscore_testing/a/rscores_100000distributions_1000points_minus10toplus10all.json', 'r') as outfile:
    rscores = json.load(outfile)
#sns.histplot(rscores, stat='probability', bins=50)
with open('/run/timeshift/backup/IOCB/cartesian/rscore_testing/d/rscores.json', 'r') as outfile:
    rscores2 = json.load(outfile)
with open('/run/timeshift/backup/IOCB/cartesian/rscore_testing/e/rscores.json', 'r') as outfile:
    rscores3 = json.load(outfile)

# Calculate and plot the histogram
rscores = np.array(rscores)[0]
rscores2 = np.array(rscores2)
rscores3 = np.array(rscores3)

rscores = np.concatenate([rscores, rscores2, rscores3])




#weights = np.ones_like(rscores) / (len(rscores))
entries, bin_edges, patches = plt.hist(rscores, bins=100, label='Data', density=True)
bin_middles = 0.5 * (bin_edges[1:] + bin_edges[:-1])
#histo, bin_edges = np.histogram(rscores, bins=250, normed=False, density=True) # For chi-squared test


# calculate bin centres


# fit with curve_fit
parameters, cov_matrix = curve_fit(fit_function, bin_middles, entries)

# plot poisson-deviation with fitted parameter
x_plot = np.arange(0, 2, 0.01)

plt.plot(x_plot, fit_function(x_plot, *parameters), marker='', linestyle='-', label='Skew normal fit')
plt.legend()


plt.show()
f_obs = np.array(entries, dtype=np.float64)
f_exp = np.array([fit_function(i, a=parameters[0], b=parameters[1], c=parameters[2]) for i in bin_middles], dtype=np.float64) # Both need to have the same sum, needed for chisquare. It was off by a little bit
diff_sum = np.abs(f_obs.sum() - f_exp.sum()) / len(f_obs)
f_obs += diff_sum

import scipy.stats as stats


r, p = stats.pearsonr(f_obs, f_exp)

print(f'Parameters for Beta distribution are {parameters}, pearson R of fit/observed is {r}, p_val={p}')

plt.plot(f_obs)
plt.plot(f_exp)
plt.show()


probability_df = pd.DataFrame([x_plot, beta.cdf(x_plot, a=1.55239973, b=4.48648693)], index=['R-score', 'probability of null']).T
#probability_df.iloc[:, 1] = probability_df.iloc[:, 1] / probability_df.iloc[:, 1].sum() # Make it so probabilities integrate to 1

#probability_df.iloc[:, 1] = (probability_df.iloc[:, 1]-probability_df.iloc[:, 1].min())/(probability_df.iloc[:, 1].max()-probability_df.iloc[:, 1].min())

with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(probability_df)

probability_df.to_csv('/run/timeshift/backup/IOCB/cartesian/rscore_testing/rscore_table_p.csv')