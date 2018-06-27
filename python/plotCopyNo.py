import os
import sys
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from FileIO import read_snp_file, read_interval_file_BAF
from RunBAFModel import calculate_BAF, generate_pi, calculate_interval

def plot_copy_no(tumor_snp, normal_snp, result_file, interval_file):
    tumorData = read_snp_file(tumor_snp)
    normalData = read_snp_file(normal_snp)
    chrmsToUse, intervalData = read_interval_file_BAF(interval_file)
    minSNP = 10
    gamma = 0.05
    print "Calculating BAFs"
    tumorBAF, normalBAF, tumorData, normalData = calculate_BAF(tumorData, normalData, chrmsToUse, minSNP, gamma,
                                                               16)

    m = len(intervalData)
    pi = generate_pi(intervalData)
    SNPToIntervalMap = [calculate_interval(pi, snp[0], snp[1]) for snp in tumorData]
    meanBAFs = [0 for i in range(m)]

    numSNPs = [0 for i in range(m)]
    for i in range(len(SNPToIntervalMap)):
        mapping = SNPToIntervalMap[i]
        if mapping is None: continue
        meanBAFs[mapping] += abs(tumorBAF[i] - 0.5)
        numSNPs[mapping] += 1.0
    meanBAFs = map(lambda (num, denom): num / denom if denom > 0 else -1, zip(meanBAFs, numSNPs))

    ifile = pd.read_table(interval_file)
    with open(result_file) as f:
        f.readline()
        C = map(int, f.readline().strip().split('\t')[2].split(':'))
    ifile['normTC'] = ifile['tumorCount'] / sum(ifile['tumorCount'])
    ifile['normNC'] = ifile['normalCount'] / sum(ifile['normalCount'])
    ifile['ratio'] = ifile['normTC'].astype('float') / ifile['normNC']
    ifile['C'] = C
    ifile['meanBAFs'] = meanBAFs
    ifile[['#ID', 'ratio', 'C']].sort_values(by='ratio')

    ifile = ifile.loc[ifile['meanBAFs'] != -1]
    print ifile

    print "Plotting clusters with copy numbers..."
    cmap = plt.get_cmap('gist_rainbow')
    maxCopy = np.max(ifile['C'])
    colors = [cmap(i) for i in  np.linspace(0, 1, maxCopy + 1)]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # copyNoAssignments = ifile['C']
    ifile['color'] = ifile['C'].map(int).map(lambda i: colors[i])
    # colorAssignment = map(lambda assignment: colors[int(assignment)], copyNoAssignments)

    ax.scatter(ifile['ratio'], ifile['meanBAFs'], c=ifile['color'])
    ifile[['ratio', 'meanBAFs', 'C']].apply(lambda (x,y,c): ax.annotate(str(int(c)), xy=(x,y)), axis=1)

    fig.savefig(os.path.basename(result_file).split('.')[0] + "_assignment_copyNo.png")

if __name__ == '__main__':
    plot_copy_no(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
