# plot read depth ratio for all bins
import os
import sys
import numpy as np
import pandas as pd

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HG19_CHM_LEN = [-1, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392,
                90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]

file1 = sys.argv[1]
file2 = sys.argv[2]
merge_bins = 1
if len(sys.argv) > 3:
    merge_bins = int(sys.argv[3])

file1_name = os.path.basename(file1).split('.')[0]
file2_name = os.path.basename(file2).split('.')[0]

bins1 = pd.read_table(file1)
bins2 = pd.read_table(file2)

def process(bins):
    bins['Abs_Start'] = bins[['Chrm', 'Start']].apply(lambda (c, s): s + sum(HG19_CHM_LEN[1:c]), axis=1)
    bins = bins.groupby(bins.index // merge_bins).agg({'Abs_Start': 'mean', 'numTumor': 'sum', 'numNormal': 'sum'})

    tumor_total = bins['numTumor'].sum()
    normal_total = bins['numNormal'].sum()

    bins['normTC'] = bins['numTumor'] / tumor_total
    bins['normNC'] = bins['numNormal'] / normal_total
    bins['ratio'] = bins['normTC'].astype('float') / bins['normNC']
    return bins

bins1 = process(bins1)
bins2 = process(bins2)

fig = plt.figure(figsize=(10, 3))
ax = fig.add_subplot(111)

chm_len_cum = [x + sum(HG19_CHM_LEN[1:i]) for i, x in enumerate(HG19_CHM_LEN)]
ax.scatter(bins1['Abs_Start'], bins1['ratio'], s=.1, c='green', label=file1_name, clip_on=False)
ax.scatter(bins2['Abs_Start'], bins2['ratio'], s=.1, c='purple', label=file2_name, clip_on=False)
ax.legend(loc='upper left')

# all in one chromosome x-axis prep package
ax.set_xticks(chm_len_cum)
ax.set_xticks([x - HG19_CHM_LEN[i] // 2 for i, x in enumerate(chm_len_cum[1:], 1)], minor=True)
ax.set_xticklabels(map(str, range(1, 23)), minor=True)
ax.set_xticklabels([], minor=False)
ax.set_xlim(0, chm_len_cum[22])
ax.xaxis.grid(which='major')

ax.set_ylim(0.0, 5.0)
ax.set_xlabel('chm')
ax.set_ylabel('rdr')
fig.tight_layout()

fig.savefig("{}_{}_bin_rd_{}.png".format(file1_name,
                                         file2_name,
                                         merge_bins))
