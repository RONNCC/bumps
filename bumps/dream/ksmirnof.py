
"""
Convergence test statistic from Gelman and Rubin, 1992.
"""

from __future__ import division

from numpy import var, mean, ones, sqrt,sum,transpose,reshape,cov,corrcoef,array,floor,apply_along_axis
from random import random
from scipy.stats import ks_2samp

def ks(seq,p=0.5):    
    #uses mean of parameters as per ROOT:
    #
    chlen,nchains,nvars = seq.shape
    def ksm(chain):
        #only return the KS statistic value and not the 2 sided p tail.
        return ks_2samp(chain[:p*chlen],chain[-p*chlen:])[0]
    def ksmpt(chain):
    #only return the KS statistic value and not the 2 sided p tail.
        return ks_2samp(chain[:p*chlen],chain[-p*chlen:])[1]
    sampks = apply_along_axis(ksm,0,seq)
    samppt = apply_along_axis(ksmpt,0,seq)
    #print 'SAMP',samp.shape
    return (sampks.flatten().tolist(),samppt.flatten().tolist())
    #return [2,2*random()]
    #print 'MEAN SHAPE',mean_seqs.shape,'SEQ1',seq1.shape,'SEQ2',seq2.shape,length,chains1.shape,chains2.shape

def test():
    raise NotImplementedError

if __name__ == "__main__":
    test()
