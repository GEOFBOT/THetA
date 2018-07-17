import sys

import pandas as pd
import seaborn as sns
sns.set_style('whitegrid')

runs_file = sys.argv[1]
chrm = sys.argv[2]
base_dir = sys.argv[3]

runs = pd.read_table(runs_file, header=None, names=['Pop','Run'])

def get_read_depth_arm(run):
    print('Processing {}'.format(run))
    import os.path
    from plotReadDepthArms import get_read_depth_arms
    file_name = os.path.join(base_dir, run) + '.intervals.no.m'

    arms = get_read_depth_arms(file_name)
    arms['ca'] = arms['Chrm'].map(str) + arms['Arm']

    rd = arms.loc[arms['ca'] == chrm]
    if rd.shape[0] == 0:
        return float('nan')
    else:
        return rd['ratio'].iloc[0]

runs['ratio'] = runs.apply(lambda (p,r): get_read_depth_arm(r), axis=1)
print(runs)
fig = sns.violinplot(x='Pop', y='ratio', data=runs, inner='stick').get_figure()
fig.savefig('{}.png'.format(chrm))
