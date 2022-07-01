import pandas as pd
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import time

start = time.time()

def calculateR_nonuniq(triple1, triple2):

    ######## R-score pair
    rmsd_sum = 0
    xs2, ys2, zs2 = tuple(triple2[0]), tuple(triple2[1]), tuple(triple2[2])
    samples = len(triple1) * len(triple2)

    # For i-th entry in first trajectory grids, for j-th entry in second trajectory grids
    i = 0
    matrix_x = []
    matrix_y = []
    matrix_z = []
    while i < len(triple1):  # This is faster than using list comprehension

        # X-term (single i with all j), Y-term, Z-term
        x_term_i = triple1.iloc[i, 0] - xs2  # Create a series for i=1 minus all possible j
        y_term_i = triple1.iloc[i, 1] - ys2
        z_term_i = triple1.iloc[i, 2] - zs2


        matrix_x.append(x_term_i)
        matrix_y.append(y_term_i)
        matrix_z.append(z_term_i)


        i += 1
    for ix, iy, iz in zip(matrix_x, matrix_y, matrix_z):
        # print(ix)
        # print(iy)
        # print(iz)

        ix = ix ** 2
        iy = iy ** 2
        iz = iz ** 2

        ixyz = ix + iy + iz

        rmsd_sum += ixyz.sum()
    rmsd = math.sqrt((1 / samples) * rmsd_sum)
    return(rmsd)
    ########

def getR(triple1, triple2):
    r = calculateR_nonuniq(triple1, triple2)
    rint1 = calculateR_nonuniq(triple1, triple1)
    rint2 = calculateR_nonuniq(triple2, triple2)
    return(( 2*r )-( rint1 + rint2 ))

"""
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$ Testing for R-score dependence on distribution widths
$ The corrected R-score is resistant to both distributions getting wider at the same rate
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
"""

# Generate randomly distributed datasets

grid = 0.1
first_low = -10
first_high = -9.8

second_low = -10
second_high = -9.99

########
x1 = pd.Series([random.uniform(first_low,first_high) for _ in range(100)]) / grid
x1 = x1.apply(np.ceil).astype(np.int16)
y1 = pd.Series([random.uniform(first_low,first_high) for _ in range(100)]) / grid
y1 = y1.apply(np.ceil).astype(np.int16)
z1 = pd.Series([random.uniform(first_low,first_high) for _ in range(100)]) / grid
z1 = z1.apply(np.ceil).astype(np.int16)
########
x2 = pd.Series([random.uniform(second_low,second_high) for _ in range(100)]) / grid
x2 = x2.apply(np.ceil).astype(np.int16)
y2 = pd.Series([random.uniform(second_low,second_high) for _ in range(100)]) / grid
y2 = y2.apply(np.ceil).astype(np.int16)
z2 = pd.Series([random.uniform(second_low,second_high) for _ in range(100)]) / grid
z2 = z2.apply(np.ceil).astype(np.int16)
########
triple1 = pd.concat([x1,y1,z1], axis=1)
triple2 = pd.concat([x2,y2,z2], axis=1)
######## !!!!!
triple1_counts = tuple([ 1 for i in range(len(triple1[0])) ])
triple2_counts = np.array([ 1 for i in range(len(triple1[0])) ])
######## !!!!!

#print(f'Total R: {getR(triple1, triple2)}')

###
from scipy.stats import norm
from scipy.stats import gamma as gm
from scipy.stats import beta as bt
from scipy.stats import skewnorm as sk

def generate_distribution(low, high, lowshift, highshift, datapoints=1000):

    def generate_gaussian(distribution):

        mean = random.randrange(low, high) # random mean
        sd = random.randint(1, 2) # random stdev originally 8 max, change back
        gaussians = norm.pdf(x, loc=mean, scale=sd)
        return(distribution + gaussians)

    def generate_gamma(distribution):
        a = random.randrange(2, 5)
        mean = random.randrange(low, high)
        gamma = gm.pdf(x, loc=mean, a=a)
        return(distribution + gamma)

    def generate_beta(distribution):
        a = random.randrange(1, 5)
        b = random.randrange(1, 5)
        mean = random.randrange(low, high)
        beta = bt.pdf(x, loc=mean, a=a, b=b)
        return(distribution + beta)

    def generate_skewed(distribution):
        mean = random.randrange(low, high)
        a = random.randrange(-10, 10) # Negative a skewed left, positive skewed right. 0 means normal.
        sd = random.randint(1, 2)  # random stdev, originally 6 max, change back
        skewed = sk.pdf(x, loc=mean, a=a, scale=sd)
        return (distribution + skewed)

    x = np.linspace(low-lowshift, high+highshift, datapoints)
    # Number of functions in the mixture is N
    seed = random.randint(0, 99)
    if seed >= 70:
        N = 1  # 30 % chance of generating a single function
    elif 60 <= seed < 70:
        N = 2
    elif 50 <= seed < 60:
        N = 3
    elif 40 <= seed < 50:
        N = 4 # 15 % four
    elif 30 <= seed < 40:
        N = 5
    elif 20 <= seed < 30:
        N = 6
    elif 10 <= seed < 20:
        N = 7
    else:
        N = 8

    distribution = np.zeros(datapoints)
    for function in range(0, N):
        # Decide on the distribution type
        seed = random.randint(0, 9)

        if seed >= 4:
            distribution = generate_gaussian(distribution) # 60 % chance of generating Gaussian (and Gaussian-mixed) distributions
        else:
            distribution = generate_skewed(distribution) # 40 % Skewed

    return(x, distribution)

####
####
####
####
Tot = 2500000 #
Datapoints = 250
####
####
####
####

Cols = 5
Rows = Tot // Cols
Rows += Tot % Cols
Position = range(1,Tot + 1)


triples1 = {}
triples2 = {}

for i in range(Tot):
    print(f'Generating {i}')
    lim1 = random.randrange(-10, 10)
    lim2 = random.randrange(-10, 10)
    low = min(lim1, lim2)
    high = max(lim1, lim2)
    ## Extends the distribution beyond low/high for gaussian and skewed mean calculation, so we can "see the whole curve" in the range
    lowshift = 5
    highshift = 5
    ##
    if int(low) == int(high):
        high = low + random.randrange(1, 5)

    #ax = fig.add_subplot(Rows, Cols, Position[i])
    dist1x = generate_distribution(low, high, lowshift, highshift, datapoints=Datapoints)[1]
    dist1y = generate_distribution(low, high, lowshift, highshift, datapoints=Datapoints)[1]
    dist1z = generate_distribution(low, high, lowshift, highshift, datapoints=Datapoints)[1]

    dist2x = generate_distribution(low, high, lowshift, highshift, datapoints=Datapoints)[1]
    dist2y = generate_distribution(low, high, lowshift, highshift, datapoints=Datapoints)[1]
    dist2z = generate_distribution(low, high, lowshift, highshift, datapoints=Datapoints)[1]

    triples1[f'x{i}'] = dist1x
    triples1[f'y{i}'] = dist1y
    triples1[f'z{i}'] = dist1z

    triples2[f'x{i}'] = dist2x
    triples2[f'y{i}'] = dist2y
    triples2[f'z{i}'] = dist2z

df1 = pd.DataFrame(triples1)
df2 = pd.DataFrame(triples2)
#df1.to_parquet('/run/timeshift/backup/IOCB/cartesian/rscore_testing/d/df1')
#df2.to_parquet('/run/timeshift/backup/IOCB/cartesian/rscore_testing/d/df2')

global rms_list
rms_list = []

def mp_caller_rscore(i):
    global rms_list
    print(i)

    low, high = i * 3, (i * 3) + 3
    triple1 = df1.iloc[:, low:high]
    triple2 = df2.iloc[:, low:high]
    triple1.columns = [0, 1, 2]
    triple2.columns = [0, 1, 2]



    return(getR(triple1, triple2))

def collector(result):
    global rms_list
    rms_list.append(result)

"""
for i in range(how_many_triples):
    print(f'Calculating R-score for {i}')
    low, high = i*3, (i*3)+3
    triple1 = df1.iloc[:, low:high]
    triple2 = df2.iloc[:, low:high]
    triple1.columns = [0,1,2]
    triple2.columns = [0,1,2]

    rms_list.append(getR(triple1, triple2))
"""

how_many_triples = int(len(df1.columns) / 3)
import multiprocessing as mp
p = mp.Pool(8)
r = p.map_async(mp_caller_rscore, iterable=[ i for i in range(how_many_triples) ], callback=collector)
r.wait()




print(f'{Tot} distributions with {Datapoints} datapoints, {time.time() - start} ')
#print(rms_list[0])
import seaborn as sns
import json
with open('/run/timeshift/backup/IOCB/cartesian/rscore_testing/e/rscores.json', 'w') as outfile:
    json.dump(rms_list[0], outfile)
sns.histplot(rms_list[0], stat='probability', bins=100)
plt.savefig('/run/timeshift/backup/IOCB/cartesian/rscore_testing/e/hist.png')

# 1000000 distributions with 500 datapoints, 22633.361141443253


plt.show()



