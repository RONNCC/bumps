"""
Convergence test statistic from Gelman and Rubin, 1992.
"""

from __future__ import division

from numpy import var, mean, ones, sqrt,sum,transpose,reshape,array

def geweke(sequences, portion=0.25):
    """
Calculates the Geweke convergence diagnostic

Refer to support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introbayes_sect008.htm
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
        var1 = var(meanseq1,axis=0)
        var2 = var(meanseq2, axis=0)
        Z_stat = (meanseq1 - meanseq2)/sqrt(var1 + var2)
        #print 'RETURNS',Z_stat.shape
    return Z_stat
def test():
    raise NotImplementedError

if __name__ == "__main__":
    test()
