import pandas as pd
import numpy as np
import math
import random


grid = 0.1
first_low = -10
first_high = -9

second_low = 4
second_high = 400

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




######## R-score pair
rmsd_sum = 0
xs2, ys2, zs2 = tuple(triple2[0]), tuple(triple2[1]), tuple(triple2[2])
sum_counts1, sum_counts2 = sum(triple1_counts), sum(triple2_counts)
samples = sum_counts1 * sum_counts2
# For i-th entry in first trajectory grids, for j-th entry in second trajectory grids
i = 0
matrix_x = []
matrix_y = []
matrix_z = []
matrix_weight = []
while i < len(triple1):  # This is faster than using list comprehension

    # X-term (single i with all j), Y-term, Z-term
    x_term_i = triple1.iloc[i, 0] - xs2  # Create a series for i=1 minus all possible j
    y_term_i = triple1.iloc[i, 1] - ys2
    z_term_i = triple1.iloc[i, 2] - zs2

    # Weights for all pairs
    i_weight = triple1_counts[i] * triple2_counts

    matrix_x.append(x_term_i)
    matrix_y.append(y_term_i)
    matrix_z.append(z_term_i)
    matrix_weight.append(i_weight)

    i += 1
for ix, iy, iz, iweight in zip(matrix_x, matrix_y, matrix_z, matrix_weight):
    # print(ix)
    # print(iy)
    # print(iz)

    ix = ix ** 2
    iy = iy ** 2
    iz = iz ** 2

    ixyz = ix + iy + iz
    ixyz_w = ixyz * iweight

    rmsd_sum += ixyz_w.sum()
rmsd = math.sqrt((1 / samples) * rmsd_sum)
########

######## R-score internal 1
rmsd_sum = 0
xs2, ys2, zs2 = tuple(triple1[0]), tuple(triple1[1]), tuple(triple1[2])
sum_counts1, sum_counts2 = sum(triple1_counts), sum(triple1_counts)
samples = sum_counts1 * sum_counts2
# For i-th entry in first trajectory grids, for j-th entry in second trajectory grids
i = 0
matrix_x = []
matrix_y = []
matrix_z = []
matrix_weight = []
while i < len(triple1):  # This is faster than using list comprehension

    # X-term (single i with all j), Y-term, Z-term
    x_term_i = triple1.iloc[i, 0] - xs2  # Create a series for i=1 minus all possible j
    y_term_i = triple1.iloc[i, 1] - ys2
    z_term_i = triple1.iloc[i, 2] - zs2

    # Weights for all pairs
    i_weight = triple1_counts[i] * triple2_counts

    matrix_x.append(x_term_i)
    matrix_y.append(y_term_i)
    matrix_z.append(z_term_i)
    matrix_weight.append(i_weight)

    i += 1
for ix, iy, iz, iweight in zip(matrix_x, matrix_y, matrix_z, matrix_weight):
    # print(ix)
    # print(iy)
    # print(iz)

    ix = ix ** 2
    iy = iy ** 2
    iz = iz ** 2

    ixyz = ix + iy + iz
    ixyz_w = ixyz * iweight

    rmsd_sum += ixyz_w.sum()
rmsd_int1 = math.sqrt((1 / samples) * rmsd_sum)
########

######## R-score internal 2
rmsd_sum = 0
xs2, ys2, zs2 = tuple(triple2[0]), tuple(triple2[1]), tuple(triple2[2])
sum_counts1, sum_counts2 = sum(triple2_counts), sum(triple2_counts)
samples = sum_counts1 * sum_counts2
# For i-th entry in first trajectory grids, for j-th entry in second trajectory grids
i = 0
matrix_x = []
matrix_y = []
matrix_z = []
matrix_weight = []
while i < len(triple2):  # This is faster than using list comprehension

    # X-term (single i with all j), Y-term, Z-term
    x_term_i = triple2.iloc[i, 0] - xs2  # Create a series for i=1 minus all possible j
    y_term_i = triple2.iloc[i, 1] - ys2
    z_term_i = triple2.iloc[i, 2] - zs2

    # Weights for all pairs
    i_weight = triple2_counts[i] * triple2_counts

    matrix_x.append(x_term_i)
    matrix_y.append(y_term_i)
    matrix_z.append(z_term_i)
    matrix_weight.append(i_weight)

    i += 1
for ix, iy, iz, iweight in zip(matrix_x, matrix_y, matrix_z, matrix_weight):
    # print(ix)
    # print(iy)
    # print(iz)

    ix = ix ** 2
    iy = iy ** 2
    iz = iz ** 2

    ixyz = ix + iy + iz
    ixyz_w = ixyz * iweight

    rmsd_sum += ixyz_w.sum()
rmsd_int2 = math.sqrt((1 / samples) * rmsd_sum)
########





print(f'pairR : {rmsd}')
print(f'1st R: {rmsd_int1}')
print(f'2nd R: {rmsd_int2}')
print('#####')
print(f'Total R: {2*rmsd-(rmsd_int1+rmsd_int2)}')