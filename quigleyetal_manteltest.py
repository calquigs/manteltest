#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import math


#migration matrix
migmat = np.array(())

#pairwise Fst matrix
fst = np.array(())

#labels for sites
site_labs = np.array(())

#remove diagonals from both matrices
migmat_nodiag = migmat.copy()
for i in range(len(migmat_nodiag)):
	migmat_nodiag[i,i] = np.nan

fst_nodiag = fst.copy()
for i in range(len(fst_nodiag)):
	fst_nodiag[i,i] = np.nan


def mantel(xs, ys, n, ax):
	#xs = xs.flatten()
	#ys = ys.flatten()
	rs = np.zeros((n))
	for i in range(n):
		ys_p = np.random.permutation(ys)
		rs[i] = stats.pearsonr(xs, ys_p)[0]
	r = stats.pearsonr(xs, ys)[0]
	sns.kdeplot(rs, ax=ax, color = 'black')
	ax.axvline(x=r, color = 'red')
	ax.axvline(x=np.sort(rs)[math.floor(n*.025)], linestyle = ':', color = 'black')
	ax.axvline(x=np.sort(rs)[math.floor(n*.975)], linestyle = ':', color = 'black')
	ax.set_title(f'Distribution of correlation coefficients\nof permuted data (n = {n})')
	ax.set_xlabel('Pearson\'s R')



#create dataframes to plot matrices

fst_df = pd.DataFrame(data = fst_nodiag, index = site_labs, columns = site_labs)
mig_df = pd.DataFrame(data = migmat, index = site_labs, columns = site_labs)



#Plot matrices and results of mantel test

fig, axes = plt.subplots(2,2)
axes = axes.flatten()

x = sns.heatmap(mig_df, cbar_kws={'label': '% migration'}, ax = axes[0])
cbar = axes[0].collections[0].colorbar
#cbar.set_ticks([1,0,-1,-2,-3,-4])
#cbar.set_ticklabels([10,1,.1,.01,.001,.0001])
axes[0].invert_xaxis()
axes[0].set_xlabel("Destination Population")
axes[0].set_ylabel("Source Population")
axes[0].set_title("Oceanographic Connectivity")

cmap = sns.cm.rocket_r
x = sns.heatmap(fst_df, cbar_kws={'label': 'Fst'}, cmap=cmap, ax = axes[1])
cbar = axes[1].collections[0].colorbar
#cbar.set_ticks([1,0,-1,-2,-3])
#cbar.set_ticklabels([10,1,.1,.01,.001])
axes[1].invert_xaxis()
plt.xlabel("Destination Population")
plt.ylabel("Source Population")
axes[1].set_title("Genetic Connectivity")


m, b = np.polyfit(np.log10(migmat)).flatten(), fst.flatten(), 1)
r = stats.pearsonr(np.log10(migmat+.00001).flatten(), fst.flatten())

axes[2].scatter(np.log10(migmat+.00001).flatten(), fst.flatten(), color = 'black')
axes[2].plot(np.log10(migmat+.00001).flatten(), m*np.log10(migmat+.00001).flatten() + b, color = 'red')
axes[2].set_xlabel('Oceanographic connectivity (log10(% migrants))')
axes[2].set_ylabel('Genetic connectivity (Fst)')
axes[2].set_title(f'Pairwise ocean connectivitiy vs. Fst relationships\nR = {r[0]:.3f}, p = {r[1]:.3e}')
#axes[2].text(.9,.9, f'R = {r[0]:.3f}\np = {r[1]:.3e}', horizontalalignment='right')

mantel(np.log10(migmat+.00001), fst, 100000, axes[3])

plt.show()


