"""
Tracks the median and confidence intervals
"""

from __future__ import division

from numpy import reshape,apply_along_axis,sort,average as avg,abs
from scipy.stats import ks_2samp        
                                                                                                                                                                                                                                                                                        
def med(seq, portion=.5 , level=.32):
    chlen,nchains,nvars = seq.shape
    front_len = int(portion*chlen)
    back_len = int(portion*chlen)
    fr,bk = level,1-level #front and back levels '
    seq1 = reshape(seq[:front_len,:,:],(front_len*nchains,nvars))
    seq2 = reshape(seq[-back_len:,:,:],(back_len*nchains,nvars))
    tmp = []
    for y in range(nvars):
        front = sort(seq1[:,y],axis=0)
        back = sort(seq2[:,y],axis=0)
        flen = len(front)
        blen = len(back)
        f_lo,f_mid,f_hi,b_lo,b_mid,b_hi = front[int(flen*fr)],front[flen//2], front[int(flen*bk)], back[int(blen*fr)],back[blen//2],back[int(blen*bk)]
        #tmp.append( [1.0*(f_lo - b_lo)/(b_mid - b_lo)] )
        #tmp.append([1.0*(f_hi-b_hi)/(b_hi-b_mid) ])
        #tmp.append( [1.0* avg( (f_hi - f_mid) - (b_hi- b_mid))/avg([f_hi - f_lo, b_hi-b_lo])  ] )
        #tmp.append( [1.0* avg( (f_mid - f_lo) - (b_mid- b_lo))/avg([f_hi - f_lo, b_hi-b_lo])  ] )
        tmp.append( [sum(abs([f_hi-b_hi, f_lo-b_lo, f_mid-b_mid]))/3.] )
        tmp.append([0])
    return tmp
