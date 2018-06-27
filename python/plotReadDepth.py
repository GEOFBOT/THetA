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

file = sys.argv[1]
merge_bins = 1
if len(sys.argv) > 2:
    merge_bins = int(sys.argv[2])

bins = pd.read_table(file)
bins['Abs_Start'] = bins[['Chrm', 'Start']].apply(lambda (c,s): s + sum(HG19_CHM_LEN[1:c]), axis=1)
bins = bins.groupby(bins.index//merge_bins).agg({'Abs_Start': 'mean', 'numTumor': 'sum', 'numNormal': 'sum'})

tumor_total = bins['numTumor'].sum()
normal_total = bins['numNormal'].sum()

bins['normTC'] = bins['numTumor'] / sum(bins['numTumor'])
bins['normNC'] = bins['numNormal'] / sum(bins['numNormal'])
bins['ratio'] = bins['normTC'].astype('float') / bins['normNC']

fig = plt.figure(figsize=(10,3))
ax = fig.add_subplot(111)

chm_len_cum = [x + sum(HG19_CHM_LEN[1:i]) for i,x in enumerate(HG19_CHM_LEN)]
ax.scatter(bins['Abs_Start'], bins['ratio'], s=.1)
ax.set_xticks(chm_len_cum)
ax.set_xticklabels(range(23))
ax.set_ylim(0.0, 5,0)
ax.set_xlim(0, chm_len_cum[22])
ax.set_xlabel('chm')
ax.set_ylabel('rdr')
fig.tight_layout()

print(bins)

fig.savefig(os.path.basename(file).split('.')[0] + "_bin_rd_{}.png".format(merge_bins))
