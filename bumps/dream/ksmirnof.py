
"""
Convergence test statistic from Gelman and Rubin, 1992.
"""

from __future__ import division

from numpy import var, mean, ones, sqrt,sum,transpose,reshape,cov,corrcoef,array,floor
from scipy.stats import ks_2samp

def ks(seq,p=0.5):    
    #uses mean of parameters as per ROOT:
    #
    mean_seqs = mean(seq,axis=2)
    length,_ = mean_seqs.shape
    transposed = transpose(mean_seqs,(1,0))
    chlen = length*p
    return [ks_2samp(transposed[x,:chlen], transposed[x,-chlen:])[0] for x in range(len(transposed))]
    #print 'MEAN SHAPE',mean_seqs.shape,'SEQ1',seq1.shape,'SEQ2',seq2.shape,length,chains1.shape,chains2.shape


if __name__ == "__main__":
    test()
