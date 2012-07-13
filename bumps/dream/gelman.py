
"""
Convergence test statistic from Gelman and Rubin, 1992.
"""

from __future__ import division

from numpy import var, mean, ones, sqrt,sum,transpose,reshape,cov,corrcoef

def gelman(sequences, portion=0.5):
    """
Calculates the R-statistic convergence diagnostic

For more information please refer to: Gelman, A. and D.R. Rubin, 1992.
Inference from Iterative Simulation Using Multiple Sequences,
Statistical Science, Volume 7, Issue 4, 457-472.
doi:10.1214/ss/1177011136
"""

    # Find the size of the sample
    chain_len,Nchains,Nvar = sequences.shape

    # Only use the last portion of the sample
    chain_len = int(chain_len*portion)
    sequences = sequences[-chain_len:]

    if chain_len < 2:
        # Set the R-statistic to a large value
        R_stat = -2 * ones(Nvar)
    else:
        # Step 1: Determine the sequence means
        meanSeq = mean(sequences, axis=0)

        # Step 1: Determine the variance between the sequence means
        B = chain_len * var(meanSeq, axis=0, ddof=1)

        # Step 2: Compute the variance of the various sequences
        varSeq = var(sequences, axis=0, ddof=1)

        # Step 2: Calculate the average of the within sequence variances
        W = mean(varSeq,axis=0)

        # Step 3: Estimate the target mean
        #mu = mean(meanSeq)

        # Step 4: Estimate the target variance (Eq. 3)
        sigma2 = ((chain_len - 1)/chain_len) * W + (1/chain_len) * B

        # Step 5: Compute the R-statistic
        R_stat = sqrt((Nchains + 1)/Nchains * sigma2 / W - (chain_len-1)/Nchains/chain_len);

    return R_stat

def gelmanP(sequences, portion=0.5):
    """
Calculates the PSRF Refined Version convergence diagnostic

For more information please refer to: 
Brooks, S. P. and Gelman, A. (1997), 
"General Methods for Monitoring Convergence  of Iterative Simulations," 
Journal of Computational and Graphical Statistics, 7, 434-455. 
"""

    # Find the size of the sample
    chain_len,Nchains,Nvar = sequences.shape
    #useful for its relation to equation but only copies of the above variables
    N, M = chain_len, Nchains
    # Only use the last portion of the sample
    chain_len = int(chain_len*portion)
    sequences = sequences[-chain_len:]

    if chain_len < 2:
        # Set the R-statistic to a large value
        R_stat = -2 * ones(Nvar)
    else:
        meanSeq = mean(sequences, axis=0)
        mean_meanSeq = mean(meanSeq,axis=0)
        varSeq = var(sequences, axis=0, ddof=1)
        
        B = chain_len * var(meanSeq, axis=0, ddof=1)
        W = mean(varSeq,axis=0)
        
        sigma2 = ((chain_len - 1)/chain_len) * W + (1/chain_len) * B
        
        #Posterior Variance Estimate, V hat
        V = sigma2 + B/(M*N)
        
        if 0:
            # Simple PSRF
            PSRF = sqrt(V/W)
        else:
            pass
            # Refined version (Brooks and Gelman, 1998)
        # d degrees of freedom
        #d = 
        #Variance Hat ( V hat )
        #VarV = ( ((N-1)/N)**2) *(1/M)*varSeq + \
        #(((M+1)/(N*M))**2) * (2/M-1)*(B)**2  + \
        #2*( (M+1)*(N-1)/((N**2) *M)*(N/M))#(cov(varSeq,meanSeq**2) - 2*meanSeq*cov(varSeq**2, meanSeq ))
        #print VarV.shape
        #d hat
        #d = 2*V**2/VarV
        # Step 5: Compute the R-statistic
        #PSRF_stat = sqrt( ((d+3) / (d+1)) * VarV/W)
    return 1

def test():
        from numpy import reshape, arange, transpose
        from numpy.linalg import norm
        # Targe values computed from octave:
        #    format long
        #    S = reshape([1:15*6*7],[15,6,7]);
        #    R = gelman(S,struct('n',6,'seq',7))
        S = reshape(arange(1.,15*6*7+1)**-2, (15, 6, 7), order='F')
        S = transpose(S, (0,2,1))
        target = [1.06169861367116,   2.75325774624905,   4.46256647696399,
                  6.12792266170178,   7.74538715553575,   9.31276519155232]
        R = gelman(S, portion=1)
        #print R
        #print "target", array(target), "\nactual", R
        assert norm(R-target) < 1e-14
        R = gelman(S, portion=.1)
        assert norm(R - [-2, -2, -2, -2, -2, -2]) == 0
        original = reshape(arange(1,20,.5), (19,2,1))
        gelman(transpose( ),.5)

if __name__ == "__main__":
    test()
