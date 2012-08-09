from scipy.stats import *
import sys
from scipy.stats.mstats import find_repeats
from numpy import *
from numpy.random import *
ret  = {}
def foil(seq, **kargs):  
    # two sample cramer von mises test
    # assumes no repeats in the data - otherwise the midrank method can be used
    # right now it discards repeats
    chlen,nchains,nvars = seq.shape
    portion = kargs['portion'] if 'portion' in kargs else .5
    count = portion*chlen*nchains
    try:
        front_portion, back_portion = portion
    except TypeError:
        front_portion = back_portion = portion
    def cvm(chain):
        #cvm can use different sized samples - this currently uses same size. Change the front and back portion to change that.
        chain_len = len(chain)
        front_len, back_len = int(chain_len*front_portion),int(chain_len*back_portion)
        seq1 = sort(chain[:front_len])
        seq2 = sort(chain[-back_len:])
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
                print 'ERROR'
                #print 'DICT',seq1r
        N,M = len(seq1),len(seq2)
        U =         N*sum([(seq1r[r]-i)**2  for i,r in enumerate(seq1)]) + \
        M*sum([(seq2r[r]-i)**2  for i,r in enumerate(seq2)])
        T = U/(N*M*(N+M)) - (4*M*N-1)/(6*(M+N))
        return T
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # >  
    # >  Cramer Von Mises Two Sample needs to be fixed.
    # >
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #ret['cramervonmises'] = apply_along_axis(cvm,0, reshape(seq, (chlen*nchains,nvars)))
    
    def timeseries(chain):
        pass
    
    def cpoints(level=.32):
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
            tmp.append([f_lo,f_mid,f_hi,b_lo,b_mid,b_hi])
        ret['cpoints']= tmp
    def wilcox():
        #returns 0 until size is at least 20.
        #implementation uses value of sample size > 20 alghough usually n >= 10 is used.
        chlen,nchains,nvars = seq.shape
        def norm(W,p):
            #print 'W,p',W,p
            return (W/(1.0*count),p)
        def wxn(chain):
            front = len(chain)*portion
            back = len(chain)*portion
            #return ks_2samp(chain[:count], chain[-count:])
            #w = wilcoxon(chain[front:],chain[-back:])
            #print wilcoxon(chain[front:],chain[-back:])
            #print wilcoxon(rand(100),rand(100))
            #print 'WILCOX',w
            #return norm(wilcoxon(chain[front:],chain[-back:]))
            return wilcoxon(chain[front:],chain[-back:])
        W,p = apply_along_axis(wxn,0,reshape(seq, (chlen*nchains,nvars)))
        #print 'CHLEN',chlen
        ret['wilcoxon'] = [W,p]
        
    def ranktest():
        chlen,nchains,nvars = seq.shape
        #returns 0 until size is at least 20.
        #implementation uses value of sample size > 20 alghough usually n >= 10 is used.
        chlen,nchains,nvars = seq.shape
        def norm(W,p):
            #print 'W,p',W,p
            return (W/(1.0*count),p)
        def wxn(chain):
            N = len(chain)*portion
            back = sort(chain[-N:])
            front = sort(chain[:N])
            target = [int(N*p) for p in 0.32, 0.5, 0.68]
            actual = searchsorted(back, [front[ti] for ti in target])
            #return [exp(-0.5*(k-i)**2/k)/sqrt(2*pi*k) for k,i in zip(target,actual)]
            le=log10(e)
            return [(-0.5*(k-i)**2/k - 0.5*log(2*pi*k))/le for k,i in zip(target,actual)]
            #return ks_2samp(chain[:count], chain[-count:])
            #w = wilcoxon(chain[front:],chain[-back:])
            #print wilcoxon(chain[front:],chain[-back:])
            #print wilcoxon(rand(100),rand(100))
            #print 'WILCOX',w
            #return norm(wilcoxon(chain[front:],chain[-back:]))
            #return wilcoxon(chain[front:],chain[-back:])
        vars = apply_along_axis(wxn,0,reshape(seq, (chlen*nchains,nvars)))
        ret['rankstat'] = vars
        
        
        
    ranktest()
    wilcox()
    cpoints()
    return ret
    #ret['friedman'] = friedmanchisquare( )