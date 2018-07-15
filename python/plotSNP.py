import os
import sys
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from FileIO import read_snp_file, read_interval_file_BAF
from RunBAFModel import calculate_BAF

HG19_CHM_LEN = [-1, 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                               146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392,
                               90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566]

def plot_snp(tumor_snp, normal_snp, interval_file):
    tumorData = read_snp_file(tumor_snp)
    normalData = read_snp_file(normal_snp)
    chrmsToUse, intervalData = read_interval_file_BAF(interval_file)
    minSNP = 10
    gamma = 0.05
    print "Calculating BAFs"
    tumorBAF, normalBAF, tumorData, normalData = calculate_BAF(tumorData, normalData, chrmsToUse, minSNP, gamma,
                                                               16)

    snp_chm = [x[0] for x in tumorData]
    snp_pos = [x[1] for x in tumorData]

    snps = pd.DataFrame.from_dict({'chm': snp_chm,
                                   'pos': snp_pos,
                                   'baf': tumorBAF})

    snps['abs_pos'] = snps[['chm', 'pos']].apply(lambda (c,s): s + sum(HG19_CHM_LEN[1:c]), axis=1)

    fig = plt.figure(figsize=(10, 3))
    ax = fig.add_subplot(111)

    ax.scatter(snps['abs_pos'], snps['baf'], s=.5, clip_on=False)

    chm_len_cum = [x + sum(HG19_CHM_LEN[1:i]) for i, x in enumerate(HG19_CHM_LEN)]
    ax.set_xticks(chm_len_cum)
    ax.set_xticks([x  - HG19_CHM_LEN[i] // 2 for i, x in enumerate(chm_len_cum[1:], 1)], minor=True)
    ax.set_xticklabels(map(str, range(1, 23)), minor=True)
    ax.set_xticklabels([], minor=False)
    ax.set_xlim(0, chm_len_cum[22])
    ax.xaxis.grid(which='major')

    ax.set_ylim(0.0, 1.0)

    # from matplotlib.ticker import FormatStrFormatter
    # ax.xaxis.set_minor_formatter(FormatStrFormatter("%i"))

    ax.set_xlabel('chm')
    ax.set_ylabel('baf')
    fig.tight_layout()

    fig.savefig(os.path.basename(tumor_snp).split('.')[0] + "_snps.png")

if __name__ == '__main__':
    plot_snp(sys.argv[1], sys.argv[2], sys.argv[3])
