"""
Sampling history for MCMC

MCMC keeps track of a number of things during sampling.

The results may be queried as follows::

    draws, generation, thinning
    sample(condition) returns draws, points, logp
    logp()            returns draws, logp
    acceptance_rate() returns draws, AR
    chains()          returns draws, chains, logp
    R_stat()          returns draws, R
    Z_stat()          returns draw, Z
    CR_weight()       returns draws, CR_weight
    best()            returns best_x, best_logp
    outliers()        returns outliers
    show()/save(file)/load(file)

Data is stored in circular arrays, which keeps the last N generations and
throws the rest away.

draws is the total number of draws from the sampler.

generation is the total number of generations.

thinning is the number of generations per stored sample.

draws[i] is the number of draws including those required to produce the
information in the corresponding return vector.  Note that draw numbers
need not be linearly spaced, since techniques like delayed rejection
will result in a varying number of samples per generation.

logp[i] is the set of log likelihoods, one for each member of the population.
The logp() method returns the complete set, and the sample() method returns
a thinned set, with on element of logp[i] for each vector point[i,:].

AR[i] is the acceptance rate at generation i, showing the proportion of
proposed points which are accepted into the population.

chains[i,:,:] is the set of points in the differential evolution population
at thinned generation i.  Ideally, the thinning rate of the MCMC process
is chosen so that thinned generations i and i+1 are independent samples
from the posterior distribution, though there is a chance that this may
not be the case, and indeed, some points in generation i+1 may be identical
to those in generation i.  Actual generation number is i*thinning.

points[i,:] is the ith point in a returned sample.  The i is just a place
holder; there is no inherent ordering to the sample once they have been
extracted from the chains.  Note that the sample may be from a marginal
distribution.

R[i] is the Gelman R statistic measuring convergence of the Markov chain.

CR_weight[i] is the set of weights used for selecting between the crossover
ratios available to the candidate generation process of differential
evolution.  These will be fixed early in the sampling, even when adaptive
differential evolution is selected.

outliers[i] is a vector containing the thinned generation number at which
an outlier chain was removed, the id of the chain that was removed and
the id of the chain that replaced it.  We leave it to the reader to decide
if the cloned samples, point[:generation, :, removed_id], should be included
in further analysis.

best_logp is the highest log likelihood observed during the analysis and
best_x is the corresponding point at which it was observed.

generation is the last generation number
"""
#TODO: state should be collected in files as we go
from __future__ import division

__all__ = ['MCMCDraw','load_state','save_state']

import re
import gzip

import numpy
from numpy import empty, sum, asarray, inf, argmax, hstack, dstack,array,append,vstack,concatenate,dtype
from numpy import savetxt,loadtxt, reshape
from .outliers import identify_outliers
from .util import draw, RNG

#EXT = ".mc.gz"
#CREATE = gzip.open
EXT = ".mc"
CREATE = open

def save_state(state, filename):
    # Build 2-D data structures
    draws, logp = state.logp(full=True)
    _, AR = state.acceptance_rate()
    chain = hstack((draws[:,None], AR[:,None], logp))

    _, point, logp = state.chains()
    Nthin,Npop,Nvar = point.shape
    point = dstack((logp[:,:,None], point))
    point = reshape(point, (point.shape[0]*point.shape[1],point.shape[2]))

    draws, R_stat = state.R_stat()
    _, CR_weight = state.CR_weight()
    _, Ncr = CR_weight.shape
    stats = hstack((draws[:,None], R_stat, CR_weight))

    #TODO: missing _outliers from save_state

    # Write convergence info
    file = CREATE(filename+'-chain'+EXT,'w')
    file.write('# draws acceptance_rate %d*logp\n'%Npop)
    savetxt(file,chain)
    file.close()

    # Write point info
    file=CREATE(filename+'-point'+EXT,'w')
    file.write('# logp point (Nthin x Npop x Nvar = [%d,%d,%d])\n'%(Nthin,Npop,Nvar))
    savetxt(file,point)
    file.close()

    # Write stats
    file=CREATE(filename+'-stats'+EXT,'w')
    file.write('# draws %d*R-stat %d*CR_weight\n'%(Nvar,Ncr))
    savetxt(file,stats)
    file.close()

IND_PAT = re.compile('-1#IND')
INF_PAT = re.compile('1#INF')

def loadtxt(file, report=0):
    """
    Like numpy.loadtxt, but adapted for windows non-finite numbers.
    """
    if not hasattr(file,'readline'):
        if file.endswith('.gz'):
            #print "opening with gzip"
            fh = gzip.open(file, 'r')
        else:
            fh = open(file, 'r')
    else:
        fh = file
    res = []
    section = 0
    lineno = 0
    for line in fh:
        lineno += 1
        if report and lineno%report==0:
            print "read",section*report
            section += 1
        IND_PAT.sub('nan', line)
        INF_PAT.sub('inf', line)
        line = line.split('#')[0].strip()
        values = line.split()
        if len(values) > 0:
            try:
                res.append([float(v) for v in values])
            except:
                print "Parse error:",values
    if fh != file:
        fh.close()
    return asarray(res)

def load_state(filename, skip=0, report=0):
    # Read chain file
    chain = loadtxt(filename+'-chain'+EXT)

    # Read point file
    file = open(filename+'-point'+EXT,'r')
    line = file.readline()
    point_dims = line[line.find('[')+1:line.find(']')]
    Nthin,Npop,Nvar = eval(point_dims)
    for _ in range(skip*Npop): file.readline()
    point = loadtxt(file,report=report*Npop)
    file.close()

    # Read stats file
    stats = loadtxt(filename+'-stats'+EXT)

    # Guess dimensions
    Ngen = chain.shape[0]
    thinning = 1
    Nthin -= skip
    Nupdate = stats.shape[0]
    #Ncr = stats.shape[1] - Nvar - 1

    # Create empty draw and fill it with loaded data
    state = MCMCDraw(0,0,0,0,0,0,thinning)
    #print "gen,var,pop",Ngen,Nvar,Npop
    state.draws = Ngen * Npop
    state.generation = Ngen
    state._gen_index = 0
    state._gen_draws = chain[:,0]
    state._gen_acceptance_rate = chain[:,1]
    state._gen_logp = chain[:,2:]
    state.thinning = thinning
    state._thin_count = Ngen//thinning
    state._thin_index = 0
    state._thin_draws = state._gen_draws[(skip+1)*thinning-1::thinning]
    state._thin_logp = point[:,0].reshape( (Nthin,Npop) )
    state._thin_point = reshape(point[:,1:], (Nthin,Npop,Nvar) )
    state._update_count = Nupdate
    state._update_index = 0
    state._update_draws = stats[:,0]
    state._update_R_stat = stats[:,1:Nvar+1]
    state._update_CR_weight = stats[:,Nvar+1:]
    state._outliers = []
    state._dof = None

    bestidx = numpy.argmax(point[:,0])
    state._best_logp = point[bestidx,0]
    state._best_x = point[bestidx,1:]

    return state

class MCMCDraw(object):
    """
    """
    _labels = None
    title = None
    @property
    def Nvar(self):
        """Number of parameters in the fit"""
        return self._thin_point.shape[2]

    def __init__(self, Ngen, Nthin, Nupdate, Nvar, Npop, Ncr, thinning):
        # Total number of draws so far
        self.draws = 0

        # Maximum observed likelihood
        self._best_x = None
        self._best_logp = -inf
        self._dof = None
        # Per generation iteration
        self.generation = 0
        self._gen_index = 0
        self._gen_draws = empty(Ngen, 'i')
        self._gen_logp = empty( (Ngen,Npop) )
        self._gen_acceptance_rate = empty(Ngen)

        # If we are thinning, we need to keep the current generation
        # separately. [Note: don't remember why we need both the _gen_*
        # and _thin_*]  [Note: the caller x vector is assigned to
        # _gen_current; this may lead to unexpected behaviour if x is
        # changed by the caller.
        self._gen_current = None

        # Per thinned generation iteration
        self.thinning = thinning
        self._thin_index = 0
        self._thin_count = 0
        self._thin_timer = 0
        self._thin_draws = empty(Nthin, 'i')
        self._thin_point = empty( (Nthin, Npop, Nvar) )
        self._thin_logp = empty( (Nthin, Npop) )

        # Per update iteration
        self._update_index = 0
        self._update_count = 0
        self._update_draws = empty(Nupdate, 'i')
        self._update_R_stat = empty( (Nupdate, Nvar) )
        self._update_CR_weight = empty( (Nupdate, Ncr) )

        #PR Stat
        self._updateP_index = 0
        self._updateP_count = 0
        self._updateP_draws = empty(Nupdate, 'i')
        self._updateP_stat = empty( (Nupdate, Nvar) )
        self._updateP_draws = empty(Nupdate, 'i')
        
        # Z stat
        self._updateZ_index = 0
        self._updateZ_count = 0
        self._updateZ_draws = empty(Nupdate, 'i')
        self._updateZ_Z_stat = []
        self._updateZ_CR_weight = empty( (Nupdate, Ncr) )
        
        #Ks Stat
        self._updateKs_index = 0
        self._updateKs_count = 0
        self._updateKs_draws = empty(Nupdate, 'i')
        self._updateKs_Ks_stat = []
        self._updateKs_Ks_ptail = []
        self._outliers = []
        
        ##Median Statistics
        
        self._updateMed_index = 0
        self._updateMed_count = 0
        self._updateMed_draws = empty(Nupdate, 'i')
        self._updateMed_Med = []
        
        self._updateOt_stat = {}
        

        # Query functions will not return outlier chains; initially, all
        # chains are marked as good.  Call mark_outliers to remove
        # outlier chains from the set.
        self._good_chains = slice(None,None)

    def resize(self, Ngen, Nthin, Nupdate, Nvar, Npop, Ncr, thinning):
        pNgen = self._gen_draws.shape[0]
        pNthin = self._thin_draws.shape[0]
        pNupdate = self._update_draws.shape[0]
        pNvar = self._thin_point.shape[2]
        pNpop = self._gen_logp.shape[1]
        pNcr = self._update_CR_weight.shape[1]

        if pNvar != Nvar or pNpop != Npop or pNcr != Ncr:
            raise ValueError("Cannot change Nvar, Npop or Ncr on resize")

        # For now, only handle the case where the we have one complete frame of data, such
        # as on reloading the state vector
        assert self._gen_index == 0 and self._update_index == 0 and self._thin_index == 0
        assert self.generation == pNgen and self._update_count == pNupdate and self._thin_count == pNthin

        self.thinning = thinning

        if Ngen > pNgen:
            self._gen_draws = numpy.resize(self._gen_draws, Ngen)
            self._gen_logp = numpy.resize(self._gen_logp,  (Ngen,Npop) )
            self._gen_acceptance_rate = numpy.resize(self._gen_acceptance_rate, Ngen)
            self._gen_index = pNgen
        elif Ngen < pNgen:
            self._gen_draws = self._gen_draws[-Ngen:].copy()
            self._gen_logp = self._gen_logp[-Ngen:,:].copy()
            self._gen_acceptance_rate = self._gen_acceptance_rate[-Ngen:].copy()

        if Nthin > pNthin:
            self._thin_draws = numpy.resize(self._thin_draws, Nthin)
            self._thin_point = numpy.resize(self._thin_point,  (Nthin, Npop, Nvar) )
            self._thin_logp = numpy.resize(self._thin_logp,  (Nthin, Npop) )
            self._thin_index = pNthin
        elif Nthin < pNthin:
            self._thin_draws = self._thin_draws[-Nthin:].copy()
            self._thin_point = self._thin_point[-Nthin:,:,:].copy()
            self._thin_logp = self._thin_logp[-Nthin:,:].copy()

        if Nupdate > pNupdate:
            self._update_draws = numpy.resize(self._update_draws, Nupdate)
            self._update_R_stat = numpy.resize(self._update_R_stat,  (Nupdate, Nvar) )
            self._update_CR_weight = numpy.resize(self._update_CR_weight,  (Nupdate, Ncr) )
            self._update_count = pNupdate
        elif Nupdate < pNupdate:
            self._update_draws = self._update_draws[-Nupdate:].copy()
            self._update_R_stat = self._update_R_stat[-Nupdate:,:].copy()
            self._update_CR_weight = self._update_CR_weight[-Nupdate:,:].copy()

    def save(self, filename):
        save_state(self,filename)

    def show(self, portion=1.0, figfile=None):
        from views import plot_all
        plot_all(self, portion=portion, figfile=figfile)

    def _last_gen(self):
        """
        Returns x, logp for most recent generation to dream.py.
        """
        # Note: if generation number has wrapped and _gen_index is 0
        # (the usual case when this function is called to resume an
        # existing chain), then this returns the last row in the array.
        return (self._thin_point[self._thin_index-1],
                self._thin_logp[self._thin_index-1])
    
    def _updateDof(self,pdof):#passed dof
        self._dof = pdof

    def _generation(self, new_draws, x, logp, accept, force_keep=False):
        """
        Called from dream.py after each generation is completed with
        a set of accepted points and their values.
        """
        # Keep track of the total number of draws
        # Note: this is first so that we tag the record with the number of
        # draws taken so far, including the current draw.
        self.draws += new_draws
        self.generation += 1

        # Record if this is the best so far
        maxid = argmax(logp)
        if logp[maxid] > self._best_logp:
            self._best_logp = logp[maxid]
            self._best_x = x[maxid,:]+0 # Force a copy

        # Record acceptance rate and cost
        i = self._gen_index 
        #print "generation",i,self.draws,"\n x",x,"\n logp",logp,"\n accept",accept
        self._gen_draws[i] = self.draws
        self._gen_acceptance_rate[i] = 100*sum(accept)/new_draws
        self._gen_logp[i] = logp
        i = i+1
        if i == len(self._gen_draws):
            i = 0
        self._gen_index = i

        # Keep every nth iteration
        self._thin_timer += 1
        if self._thin_timer == self.thinning or force_keep:
            self._thin_timer = 0
            self._thin_count += 1
            i = self._thin_index
            self._thin_draws[i] = self.draws
            self._thin_point[i] = x
            self._thin_logp[i] = logp
            i = i+1
            if i == len(self._thin_draws): i = 0
            self._thin_index = i
            self._gen_current = None
        else:
            self._gen_current = x+0 # force a copy


    def _update(self, R_stat, CR_weight):
        """
        Updates Gelman R Statistic
        Called from dream.py when a series of DE steps is completed and
        summary statistics/adaptations are ready to be stored.
        """
        self._update_count += 1
        i = self._update_index
        #print "update",i,self.draws,"\n Rstat",R_stat,"\n CR weight",CR_weight
        self._update_draws[i] = self.draws
        self._update_R_stat[i] = R_stat
        self._update_CR_weight[i] = CR_weight
        i = i+1
        if i == len(self._update_draws): i = 0
        self._update_index = i
        
 
    def _updateP(self, PR_stat):
        """
        Updates Gelman R Statistic
        Called from dream.py when a series of DE steps is completed and
        summary statistics/adaptations are ready to be stored.
        """
        self._updateP_count += 1
        i = self._updateP_index
        #print "update",i,self.draws,"\n Rstat",R_stat,"\n CR weight",CR_weight
        self._updateP_draws[i] = self.draws
        self._updateP_stat[i] = PR_stat
        i = i+1
        if i == len(self._updateP_draws): i = 0
        self._updateP_index = i
        
    def _updateZ(self,Z_stat):
        """
        Updates Geweke Z Statistic
        Called from dream.py when a series of DE steps is completed and
        summary statistics/adaptations are ready to be stored.
        """
        self._updateZ_count += 1
        i = self._updateZ_index
        #print "update",i,self.draws,"\n Rstat",R_stat,"\n CR weight",CR_weight
        self._updateZ_draws[i] =  self.draws
        #print 'CURRENTLY',self._updateZ_Z_stat.shape, Z_stat.shape
        self._updateZ_Z_stat.append(Z_stat)
        #print 'OVERALLZ',self._updateZ_Z_stat
        i = i+1
        if i == len(self._updateZ_draws): i = 0
        self._updateZ_index = i
        
    def _updateKs(self,Ks_stat,ptail):
        self._updateKs_count += 1
        i = self._updateKs_index
        self._updateKs_draws[i] = self.draws
        self._updateKs_Ks_stat.append(Ks_stat)
        self._updateKs_Ks_ptail.append(ptail)
        i = i+1
        if i == len(self._updateKs_draws): i = 0
        self._updateKs_index = i
        
    def _updateMed(self,Med):
        self._updateMed_count += 1
        i = self._updateMed_index
        self._updateMed_draws[i] = self.draws
        self._updateMed_Med.append(Med)    
        i = i+1
        if i == len(self._updateMed_draws): i = 0
        self._update_Med_index = i
        
    def _updateOt(self,kargs):
        for k,v in kargs.iteritems():
            if k not in self._updateOt_stat:
                self._updateOt_stat[k] = []
            if k in self._updateOt_stat:
                current = self._updateOt_stat[k]
                current.append(v)
                self._updateOt_stat[k] = current
        
    def _replace_outlier(self, old, new):
        """
        Called from outliers.py when a chain is replaced by the
        clone of another.
        """
        self._outliers.append((self._thin_index,old,new))

        self._gen_logp[:,old] = self._gen_logp[:,new]
        self._thin_logp[:,old] = self._thin_logp[:,new]
        self._thin_point[:,old,:] = self._thin_point[:,new,:]
        # PAK: shouldn't we reduce the total number of draws since we
        # are throwing way an entire chain?

    def _get_labels(self):
        if self._labels == None:
            return ["P%d"%i for i in range(self._thin_point.shape[2])]
        else:
            return self._labels
    def _set_labels(self, v):
        self._labels = v
    labels = property(fget=_get_labels,fset=_set_labels)

    def _draw_pop(self, Npop):
        """
        Generate a population from current generation and all history.
        """
        _, chains, _ = self.chains()
        Ngen,Nchain,Nvar = chains.shape
        points = reshape(chains,(Ngen*Nchain,Nvar))

        # There are two complications with the history buffer:
        # (1) due to thinning, not every generation is stored
        # (2) because it is circular, the cursor may be in the middle
        # If the current generation isn't in the buffer (but is instead
        # stored separately as _gen_current), then the entire buffer
        # becomes the history pool.
        # otherwise we need to exclude the current generation from
        # the pool.  If (2) happens, we need to increment everything
        # above the cursor by the number of chains.
        if self._gen_current != None:
            pool_size = Ngen*Nchain
            cursor = pool_size  # infinite
        else:
            pool_size = (Ngen-1)*Nchain
            k = len(self._thin_draws)
            cursor = Nchain*((k+self._thin_index-1)%k)

        # Make a return population and fill it with the current generation
        pop = empty((Npop,Nvar),'d')
        if self._gen_current != None:
            pop[:Nchain] = self._gen_current
        else:
            #print pop.shape, points.shape, chains.shape
            pop[:Nchain] = points[cursor:cursor+Nchain]

        if Npop > Nchain:
            # Find the remainder with unique ancestors.
            # Again, because this is a circular buffer, their may be random
            # numbers generated at or above the cursor.  All of these must
            # be shifted by Nchains to avoid the cursor.
            perm = draw(Npop-Nchain,pool_size)
            perm[perm>=cursor] += Nchain
            #print "perm",perm; raw_input('wait')
            pop[Nchain:] = points[perm]

        return pop

    def _unroll(self):
        """
        Unroll the circular queue so that data access can be done inplace.
        Call this when done stepping, and before plotting. Calls to
        logp, sample, etc. assume the data is already unrolled.
        """
        if self.generation > self._gen_index > 0:
            self._gen_draws[:] = numpy.roll(self._gen_draws,
                                            -self._gen_index, axis=0)
            self._gen_logp[:] = numpy.roll(self._gen_logp,
                                           -self._gen_index, axis=0)
            self._gen_acceptance_rate[:] = numpy.roll(self._gen_acceptance_rate,
                                                      -self._gen_index, axis=0)
            self._gen_index = 0

        if self._thin_count > self._thin_index > 0:
            self._thin_draws[:] = numpy.roll(self._thin_draws,
                                             -self._thin_index, axis=0)
            self._thin_point[:] = numpy.roll(self._thin_point,
                                             -self._thin_index, axis=0)
            self._thin_logp[:] = numpy.roll(self._thin_logp,
                                            -self._thin_index, axis=0)
            self._thin_index = 0

        if self._update_count > self._update_index > 0:
            self._update_draws[:] = numpy.roll(self._update_draws,
                                               -self._update_index, axis=0)
            self._update_R_stat[:] = numpy.roll(self._update_R_stat,
                                                -self._update_index, axis=0)
            self._update_CR_weight[:] = numpy.roll(self._update_CR_weight,
                                                   -self._update_index, axis=0)
            self._update_index = 0
            
    def remove_outliers(self, x, logp, test='IQR', portion=0.5):
        """
        Replace outlier chains with clones of good ones.  This should happen
        early in the sampling processes so the clones have an opportunity
        to evolve their own identity.

        *state* contains the chains, with log likelihood for each point
        *x*, *logp* is the current population and the corresponding log likelihoods
        *test* is the name of the test to use (one of IQR, Grubbs, Mahal or none).
        *portion* in (0,1] is the amount of the chain to use

        Updates *state*, *x* and *logp* to reflect the changes.

        See :module:`outliers` for details.
        """
        # Grab the last part of the chain histories
        _, chains = self.logp()
        chain_len, Nchains = chains.shape
        outliers = identify_outliers(test, chains[-chain_len:], x)

        # Loop over each outlier chain, replacing each with another
        for old in outliers:
            # Draw another chain at random, with replacement
            while True:
                new = RNG.randint(Nchains)
                if new not in outliers: break
            # Update the saved state and current population
            self._replace_outlier(old=old,new=new)
            x[old,:] = x[new,:]
            logp[old] = logp[new]

    def mark_outliers(self, test='IQR', portion=1.0):
        """
        Mark some chains as outliers but don't remove them.  This can happen
        after drawing is complete, so that chains that did not converge are
        not included in the statistics.

        *test* is 'IQR', 'Mahol' or 'none'.

        *portion* indicates what portion of the samples should be included
        in the outlier test.  The default is to include all of them.
        """
        _, chains, logp = self.chains()
        if test=='none':
            self._good_chains = slice(None,None)
        else:
            Ngen = chains.shape[0]
            start = int(Ngen*(1-portion)) if portion else 0
            outliers = identify_outliers(test, logp[start:], chains[-1])
            #print "outliers",outliers
            #print logp.shape, chains.shape
            self._good_chains = numpy.array([i for i in range(logp.shape[1])
                                             if i not in outliers])
            #print self._good_chains


    def logp(self, full=False):
        """
        Return the iteration number and the log likelihood for each point in
        the individual sequences in that iteration.

        For example, to plot the convergence of each sequence::

            draw, logp = state.logp()
            plot(draw, logp)

        Note that draw[i] represents the total number of samples taken,
        including those for the samples in logp[i].

        If full is True, then return all chains, not just good chains.
        """
        self._unroll()
        retval = self._gen_draws, self._gen_logp
        if self.generation == self._gen_index:
            retval = [v[:self.generation] for v in retval]
        draws,logp = retval
        return draws,(logp if full else logp[:,self._good_chains])

    def acceptance_rate(self):
        """
        Return the iteration number and the acceptance rate for that iteration.

        For example, to plot the acceptance rate over time::

            draw, AR = state.acceptance_rate()
            plot(draw, AR)

        """
        self._unroll()
        retval = self._gen_draws, self._gen_acceptance_rate
        if self.generation == self._gen_index:
            retval = [v[:self.generation] for v in retval]
        return retval

    def chains(self):
        """
        Returns the observed Markov chains and the corresponding likelihoods.

        The return value is a tuple (*draws*,*chains*,*logp*).

        *draws* is the number of samples taken up to and including the samples
        for the current generation.

        *chains* is a three dimensional array of generations X chains X vars
        giving the set of points observed for each chain in every generation.
        Only the thinned samples are returned.

        *logp* is a two dimensional array of generation X population giving
        the log likelihood of observing the set of variable values given in
        chains.
        """
        self._unroll()
        retval = self._thin_draws, self._thin_point, self._thin_logp
        if self._thin_count == self._thin_index:
            retval = [v[:self._thin_count] for v in retval]
        return retval

    def R_stat(self):
        #from .gelman import gelman
        """
        Return the R-statistics convergence statistic for each variable.

        For example, to plot the convergence of all variables over time::

            draw, R = state.R_stat()
            plot(draw, R)

        See :module:`dream.gelman` and references detailed therein.
        """
        self._unroll()
        #print 'CHAIN',self.chains()[1], self.chains()[1].shape
        retval = self._update_draws, self._update_R_stat
        if self._update_count == self._update_index:
            retval = [v[:self._update_count] for v in retval]
        #print retval, retval.shape
        return retval
    
    def R_stat2(self):
        #from .gelman import gelman
        """
        Return the R-statistics convergence statistic for each variable.
        
        For example, to plot the convergence of all variables over time::
        
        draw, R = state.R_stat()
        plot(draw, R)
        
        See :module:`dream.gelman` and references detailed therein.
        """
        self._unroll()
        #print 'CHAIN',self.chains()[1], self.chains()[1].shape
        retval = self._updateP_draws, self._updateP_stat
        if self._updateP_count == self._updateP_index:
            retval = [v[:self._updateP_count] for v in retval]
        #print retval, retval.shape
        return retval
    
        
    def Z_stat(self):
        self._unroll()
        #print 'CHAIN',self.chains()[1], self.chains()[1].shape
        retval = self._updateZ_draws, self._updateZ_Z_stat
        if self._updateZ_count == self._updateZ_index:
            retval = [v[:self._updateZ_count] for v in retval]
        #print retval, retval.shape
        return retval

    def Ks_stat(self):
        self._unroll()
        #print 'CHAIN',self.chains()[1], self.chains()[1].shape
        retval = self._updateKs_draws, self._updateKs_Ks_stat,self._updateKs_Ks_ptail
        if self._updateKs_count == self._updateKs_index:
            retval = [v[:self._updateKs_count] for v in retval]
        #print retval, retval.shape
        return retval
    
    def Ad_stat(self):
        self._unroll()
        #print 'CHAIN',self.chains()[1], self.chains()[1].shape
        retval = self._updateAd_draws, self._updateAd_stat
        #print retval, retval.shape
        return retval
    
    def Med(self):
        self._unroll()
        retval = self._updateMed_draws, self._updateMed_Med
        if self._updateMed_count == self._updateMed_index:
            retval = [v[:self._updateMed_count] for v in retval]
        return retval   

    def OtStat(self):
        return self._updateOt_stat

    def CR_weight(self):
        """
        Return the crossover ratio weights to be used in the next generation.

        For example, to see if the adaptive CR is stable use::

            draw, weight = state.CR_weight()
            plot(draw, weight)

        See :module:`dream.crossover` for details.
        """
        self._unroll()
        retval = self._update_draws, self._update_CR_weight
        if self._update_count == self._update_index:
            retval = [v[:self._update_count] for v in retval]
        return retval

    def outliers(self):
        """
        Return a list of outlier removal operations.

        Each outlier operation is a tuple giving the thinned generation
        in which it occurred, the old chain id and the new chain id.

        The chains themselves have already been updated to reflect the
        removal.

        Curiously, it is possible for the maximum likelihood seen so far
        to be removed by this operation.
        """
        return asarray(self._outliers, 'i')

    def best(self):
        """
        Return the best point seen and its log likelihood.
        """
        return self._best_x, self._best_logp

    def keep_best(self):
        """
        Place the best point at the end of the chain final good chain.

        Good chains are defined by mark_outliers.

        Because the Markov chain is designed to wander the parameter
        space, the best individual seen during the random walk may have
        been observed during the burn-in period, and may no longer be
        present in the chain.  If this is the case, replace the final
        point with the best, otherwise swap the positions of the final
        and the best.
        """

        # Get state as a 1D array
        _, chains, logp = self.chains()
        Ngen,Npop,Nvar = chains.shape
        points = reshape(chains,(Ngen*Npop,Nvar))
        logp = reshape(logp,Ngen*Npop)

        # Set the final position to the end of the last good chain.  If
        # mark_outliers has not been called, then _good_chains will
        # just be slice(None,None)
        if isinstance(self._good_chains, slice):
            final = -1
        else:
            final = self._good_chains[-1] - Npop

        # Find the location of the best point if it exists and swap with
        # the final position
        idx = numpy.where(logp==self._best_logp)[0]
        if len(idx) == 0:
            logp[final] = self._best_logp
            points[final,:] = self._best_x
        else:
            idx = idx[0]
            logp[final],logp[idx] = logp[idx], logp[final]
            points[final,:],points[idx,:] = points[idx,:],points[final,:]

    def sample(self, portion=1, vars=None, selection=None):
        """
        Return a sample from the posterior distribution.

        *portion* is the portion of each chain to use
        *vars* is a list of variables to return for each point
        *selection* sets the range for the returned marginal distribution

        *selection* is a dictionary of {variable: (low,high)} to set the
        range on each variable.  Missing variables default to the full
        range.

        To plot the distribution for parameter p1::

            points, logp = state.sample()
            hist(points[:,0])

        To plot the interdependence of p1 and p2::

            draw, points, logp = state.sample()
            plot(points[:,0],points[:,1],'.')

        """
        return _sample(self, portion, vars, selection)

    def derive_vars(self, fn, labels=None):
        """
        Generate derived variables from the current sample, adding columns
        for the derived variables to each sample of every chain.

        The new columns are treated as part of the sample.

        *fn* is a function taking points p[:,k] for k in 0 ... samples and
        returning a set of derived variables pj[k] for each sample k.  The
        variables can be returned as any kind of sequence including an
        array or a tuple with one entry per variable.  The caller uses
        asarray to convert the returned variables into a vars X samples array.
        For convenience, a single variable can be returned by itself.

        *labels* are the labels to use for the derived variables.

        The following example adds the new variable x+y = P[0] + P[1]::

            state.derive_vars(lambda p: p[0]+p[1], labels=["x+y"])
        """
        # Grab all samples as a set of points
        _, chains, _ = self.chains()
        Ngen,Npop,Nvar = chains.shape
        points = reshape(chains,(Ngen*Npop,Nvar))

        # Compute new variables from the points
        newvars = asarray(fn(points.T)).T
        Nnew = newvars.shape[1] if len(newvars.shape) == 2 else 1
        newvars.reshape((Ngen, Npop, Nnew))

        # Extend new variables to be the same length as the stored selection
        Nthin = self._thin_point.shape[0]
        newvars = numpy.resize(newvars, (Nthin, Npop, Nnew))

        # Add new variables to the points
        self._thin_point = dstack( (self._thin_point, newvars) )

        # Add labels for the new variables, if available.
        if labels != None:
            self.labels = self.labels + labels
        elif self._labels != None:
            labels = ["P%d"%i for i in range(Nvar,Nvar+Nnew)]
            self.labels = self.labels + labels
        else: # no labels specified, old or new
            pass


def _sample(state, portion, vars, selection):
    """
    Return a sample from a set of chains.
    """
    draw, chains, logp = state.chains()
    if portion is None:
        portion = 1
    start = int((1-portion)*len(draw))


    # Collect the subset we are interested in
    chains = chains[start:,state._good_chains,:]
    logp = logp[start:,state._good_chains]

    Ngen,Npop,Nvar = chains.shape
    points = reshape(chains,(Ngen*Npop,Nvar))
    logp = reshape(logp,(Ngen*Npop))
    if selection not in [None, {}]:
        idx = True
        for v,r in selection.items():
            if v == 'logp':
                idx = idx & (logp>=r[0]) & (logp<=r[1])
            else:
                idx = idx & (points[:,v]>=r[0]) & (points[:,v]<=r[1])
        points = points[idx,:]
        logp = logp[idx]
    if vars != None:
        points = points[:,vars]
    return points, logp


def test():
    from numpy.linalg import norm
    from numpy.random import rand
    from numpy import arange

    # Make some fake data
    Nupdate,Nstep = 3,5
    Ngen = Nupdate*Nstep
    Nvar,Npop,Ncr = 3,6,2
    xin = rand(Ngen,Npop,Nvar)
    pin = rand(Ngen,Npop)
    accept = rand(Ngen,Npop) < 0.8
    CRin = rand(Nupdate,Ncr)
    Rin = rand(Nupdate,1)
    #thinning = 2
    #Nthin = int(Ngen/thinning)

    # Put it into a state
    thinning = 2
    Nthin = int(Ngen/thinning)
    state = MCMCDraw(Ngen=Ngen, Nthin=Nthin, Nupdate=Nupdate,
                     Nvar=Nvar, Npop=Npop, Ncr=Ncr, thinning=thinning)
    for i in range(Nupdate):
        state._update(R_stat=Rin[i],CR_weight=CRin[i])
        for j in range(Nstep):
            gen = i*Nstep+j
            state._generation(new_draws=Npop, x=xin[gen],
                              logp=pin[gen], accept=accept[gen])

    # Check that it got there
    draws,logp = state.logp()
    assert norm(draws - Npop*arange(1,Ngen+1)) == 0
    assert norm(logp - pin) == 0
    draws,AR = state.acceptance_rate()
    assert norm(draws - Npop*arange(1,Ngen+1)) == 0
    assert norm(AR - 100*sum(accept,axis=1)/Npop) == 0
    draws,sample,logp = state.sample()
    assert norm(draws - thinning*Npop*arange(1,Nthin+1)) == 0
    assert norm(sample - xin[thinning-1::thinning]) == 0
    assert norm(logp - pin[thinning-1::thinning]) == 0
    draws,R = state.R_stat()
    assert norm(draws - Npop*Nstep*arange(Nupdate)) == 0
    assert norm(R-Rin) == 0
    draws,CR = state.CR_weight()
    assert norm(draws - Npop*Nstep*arange(Nupdate)) == 0
    assert norm(CR - CRin) == 0
    x,p = state.best()
    bestid = argmax(pin)
    i,j = bestid//Npop, bestid%Npop
    assert pin[i,j] == p
    assert norm(xin[i,j,:]-x) == 0

    # Check that outlier updates properly
    state._replace_outlier(1,2)
    outliers = state.outliers()
    draws,sample,logp = state.sample()
    assert norm(outliers -  asarray([[state._thin_index,1,2]])) == 0
    assert norm(sample[:,1,:] - xin[thinning-1::thinning,2,:]) == 0
    assert norm(sample[:,2,:] - xin[thinning-1::thinning,2,:]) == 0
    assert norm(logp[:,1] - pin[thinning-1::thinning,2]) == 0
    assert norm(logp[:,2] - pin[thinning-1::thinning,2]) == 0

if __name__ == "__main__":
    test()
