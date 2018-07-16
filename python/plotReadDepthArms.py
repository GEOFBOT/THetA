# plot read depth ratio, grouped by chromosome arms
import os
import sys
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HG19_CHM_LEN = [-1, 249250621, 243199373, 198022430, 191154276, 180915260,
                171115067, 159138663, 146364022, 141213431, 135534747,
                135006516, 133851895, 115169878, 107349540, 102531392,
                90354753, 81195210, 78077248, 59128983, 63025520,
                48129895, 51304566]

HG19_CENTROMERE = [-1, 121535434, 92326171, 90504854, 49660117, 46405641,
                   58830166, 58054331, 43838887, 47367679, 39254935, 51644205,
                   34856694, 16000000, 16000000, 17000000, 35335801, 22263006,
                   15460898, 24681782, 26369569, 11288129, 13000000]

HG19_CENTROMERE_END = [-1, 124535434, 95326171, 93504854, 52660117, 49405641,
                       61830166, 61054331, 46838887, 50367679, 42254935,
                       54644205, 37856694, 19000000, 19000000, 20000000,
                       38335801, 25263006, 18460898, 27681782, 29369569,
                       14288129, 16000000]

def get_abs_pos(chrm, rel_pos):
    """Get absolute coor. from relative coor."""
    return rel_pos + sum(HG19_CHM_LEN[1:chrm])

def get_read_depth_arms(file):
    bins = pd.read_table(file)
    bins['Arm'] = bins.apply(lambda r: 'p' if r['Start'] < HG19_CENTROMERE[r['Chrm']] else 'q', axis=1)

    bins['Abs_Start'] = bins[['Chrm', 'Start']].apply(lambda (c,s): s + sum(HG19_CHM_LEN[1:c]), axis=1)
    bins = bins.groupby(['Chrm', 'Arm']).agg({'Abs_Start': 'min', 'numTumor': 'sum', 'numNormal': 'sum'}).reset_index()

    tumor_total = bins['numTumor'].sum()
    normal_total = bins['numNormal'].sum()

    bins['normTC'] = bins['numTumor'] / sum(bins['numTumor'])
    bins['normNC'] = bins['numNormal'] / sum(bins['numNormal'])
    bins['ratio'] = bins['normTC'].astype('float') / bins['normNC']
    return bins

if __name__ == '__main__':
    file = sys.argv[1]

    fig = plt.figure(figsize=(10,3))
    ax = fig.add_subplot(111)
    bins = get_read_depth_arms(file)

    chm_len_cum = [get_abs_pos(i, x) for i,x in enumerate(HG19_CHM_LEN)]
    ctm_cum = [get_abs_pos(i, x) for i,x in enumerate(HG19_CENTROMERE)]
    ctm_cum_end = [get_abs_pos(i, x) for i,x in enumerate(HG19_CENTROMERE_END)]

    bins[['Abs_Start', 'ratio', 'Arm', 'Chrm']].apply(
        lambda (x_c,y,a,c): ax.plot(
            (chm_len_cum[c-1] if a == 'p' else ctm_cum_end[c], ctm_cum[c] if a == 'p' else chm_len_cum[c]),
            (y, y), linewidth=3, clip_on=False, c='b', solid_capstyle='butt'),
        axis=1)

    # all in one chromosome x-axis prep package
    ax.set_xticks(chm_len_cum)
    ax.set_xticks([x  - HG19_CHM_LEN[i] // 2 for i, x in enumerate(chm_len_cum[1:], 1)], minor=True)
    ax.set_xticklabels(map(str, range(1, 23)), minor=True)
    ax.set_xticklabels([], minor=False)
    ax.set_xlim(0, chm_len_cum[22])
    ax.xaxis.grid(which='major')

    ax.set_ylim(0.0, 5.0)
    ax.set_xlabel('chm')
    ax.set_ylabel('rdr')
    ax.set_title(os.path.basename(file).split('.')[0])
    fig.tight_layout()

    print(bins)

    fig.savefig(os.path.basename(file).split('.')[0] + "_bin_rd_arms.png")
