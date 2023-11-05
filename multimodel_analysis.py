# EBI-full cibersort results multimodel fitting
# in this script, we fit multi-model on cibersort results from EBI-full signature - TCGA. Multi-model fitting will be applied to per signature and per tumor type, we only keep the ones that have bi-modelity
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from smm.smm import SMM
import scipy as scp
import math
from scipy.stats import combine_pvalues
from scipy.stats import rankdata
#import torch
import copy
from scipy import optimize
import itertools
import warnings
from scipy.stats.stats import pearsonr
import scipy.stats
import glob
import os
random_seed = 1

#SMM model document: https://t-student-mixture-models.readthedocs.io/en/latest/_modules/smm/smm.html#SMM.bic

#run bimodality test on all 33 cancer cibersort results

import glob
import os

os.chdir("EBI_sig217_louvain/EBI_sig217_louvain_decon_results")
files = glob.glob('Cibersort_results_EBI-full_TCGA*.txt')

for file in files:
    print(file)
    cancer = file.split('_')[3].split('-')[1].split('.txt')[0]
    cibersort = pd.read_csv(file, sep='\t', index_col=0)
    sigs = cibersort.columns[0:-4]
    #fit in student-t mixture model
    fit_sig = []
    fit_bimodel_sig_boundary = dict()
    unimodel_sig = []
    warning = []

    sd_min = []
    mean_min = []
    sd_max = []
    mean_max = []
    warnings.filterwarnings("error")
    for sig in sigs:
        number_zeros = len(np.where(cibersort[sig].values == 0)[0])
        percent_zeros = len(np.where(cibersort[sig].values == 0)[0])/cibersort.shape[0]
        if percent_zeros <= 0.9:
            try:
                smm1 = SMM(n_components=1)
                smm1.fit(cibersort[sig].values[:,None])
                smm2 = SMM(n_components=2)
                smm2.fit(cibersort[sig].values[:,None])
                if smm1.bic(cibersort[sig].values[:,None]) > smm2.bic(cibersort[sig].values[:,None]):
                    #print("processing signature: ", sig)
                    #extract weights, degree, mean and covariance
                    i = np.argmin(smm2.means)
                    j = np.argmax(smm2.means)

                    wt1 = smm2.weights[i]
                    df1 = smm2.degrees[i]
                    mn1 = smm2.means[i,0]
                    sd1 = np.sqrt(smm2.covariances[i,0,0])

                    wt2 = smm2.weights[j]
                    df2 = smm2.degrees[j]
                    mn2 = smm2.means[j,0]
                    sd2 = np.sqrt(smm2.covariances[j,0,0])

                    # Find decision boundary
                    difBetDists = lambda x: wt1*scp.stats.t.pdf(x, df1, mn1, sd1) - wt2*scp.stats.t.pdf(x, df2, mn2, sd2)
                    x0 = np.median(cibersort[sig].values) # Initial Guess
                    x1 = np.mean(cibersort[sig].values)   # Second  Guess
                    decBoundary = optimize.root_scalar(f = difBetDists, x0 = x0, x1 = x1).root
                    if decBoundary >0 and decBoundary<1 and decBoundary > mn1 and decBoundary < mn2:
                        fit_bimodel_sig_boundary[sig] = decBoundary
                        fit_sig.append(sig)
                        mean_min.append(mn1)
                        mean_max.append(mn2)
                        sd_min.append(sd1)
                        sd_max.append(sd2)
                    else:
                        pass
                else:
                    unimodel_sig.append(sig)

            except RuntimeWarning:
                warning.append(sig)


    d = {'sigs': list(fit_bimodel_sig_boundary.keys()), 'boundary': list(fit_bimodel_sig_boundary.values()), 'mean1':mean_min, 'mean2':mean_max,'sd1':sd_min, 'sd2':sd_max}
    df = pd.DataFrame(data=d)
    df.to_csv('./EBI-full_'+cancer+'-bimodel-signatures_decision-boundary.tsv', sep='\t')
