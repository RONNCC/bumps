from scipy.stats import *
import sys
from scipy.stats.mstats import find_repeats
from numpy import split,apply_along_axis,reshape,unique
ret  = {}
def foil(seq, **kargs):  
    # two sample cramer von mises test
    # assumes no repeats in the data - otherwise the midrank method can be used
    # right now it discards repeats
    chlen,nchains,nvars = seq.shape
    portion = kargs['portion'] if 'portion' in kargs else .5
    count = portion*chlen*nchains
    print 'STARTED'
    try:
        front_portion, back_portion = portion
    except TypeError:
        front_portion = back_portion = portion
    def cvm(chain):
        #cvm can use different sized samples - this currently uses same size. Change the front and back portion to change that.
        chain_len = len(chain)
        chain = unique(chain)
        front_len, back_len = int(chain_len*front_portion),int(chain_len*back_portion)
        seq1 = unique(chain[:front_len])
        seq2 = unique(chain[-back_len:])
        seq1r, seq2r = {},{}
        #print 'BLAH'
        #print find_repeats(chain)
        
        for i,r in enumerate(chain):
            #print seq1
            #print '--------------'
            #print seq2
            #print '>>>>>>>>>>>>>>>'
            #print 'R',r, r in seq1, r in seq2
            #assert( (r in seq1 and r in seq2)  == False), 'Elem in Both'
            if r in seq1:
                seq1r[r] = i
            elif r in seq2:
                seq2r[r] = i
            else:
                print 'DARNIT'
        #print 'DICT',seq1r
        N,M = len(seq1),len(seq2)
        U =         N*sum([(seq1r[r]-i)**2  for i,r in enumerate(seq1)]) + \
        M*sum([(seq2r[r]-i)**2  for i,r in enumerate(seq2)])
        T = U/(N*M*(N+M)) - (4*M*N-1)/(6*(M+N))
        return T
    ret['cramervonmises'] = apply_along_axis(cvm,0, reshape(seq, (chlen*nchains,nvars)))
    print 'RETURNED'
    return ret
    #ret['friedman'] = friedmanchisquare( )