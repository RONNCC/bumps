"""
Convergence test statistic from Gelman and Rubin, 1992.
"""

from __future__ import division

from numpy import var, mean, ones, sqrt,sum,transpose,reshape,array,log10,abs
from random import random
def geweke(sequences, portion=0.25):
    """
Calculates the Geweke convergence diagnostic

Refer to: 

    pymc-devs.github.com/pymc/modelchecking.html#informal-methods
    support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introbayes_sect008.html
    
"""

    # Find the size of the sample
    chain_len,Nchains,Nvar = sequences.shape
    # Only use the last portion of the sample
    Z_stat = 0
    if chain_len < 2:
        # Set the R-statistic to a large value
        Z_stat = -2 * ones(Nvar)
    else:
        new_len = int(chain_len*portion)
        #print "STARTING SHAPE",sequences.shape
        seq1 = sequences[:new_len,:,:]
        seq2 = sequences[-new_len:,:,:]
        #print "SEQ1",seq1.shape,'SEQ2',seq2.shape
        # Step 1: Determine the sequence means
        meanseq1 = mean(seq1, axis=0)
        meanseq2 = mean(seq2, axis=0)
        #print "SHAPEs",meanseq1.shape,meanseq2.shape
        var1 = var(seq1,axis=0)
        var2 = var(seq2, axis=0)
        Z_stat = (meanseq1 - meanseq2)/sqrt(var1 + var2)
        
        #Z_stat is now the Z score for every chain and parameter in that with shape (chains,vars)
        
        #To make it easier to look at, return the average for the vars.
        if 0:
            Avg_Z = mean(Z_stat,axis=0)
        Avg_Z = Z_stat
        #Print absolute value log so it looks cleaner
        LAvg_Z = log10(abs(Avg_Z))
    return LAvg_Z.tolist()
def test():
    raise NotImplementedError

if __name__ == "__main__":
    test()
