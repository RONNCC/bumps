from __future__ import with_statement
import os.path
import inspect
import sys

sys.path.append(os.path.abspath('..'))

MODULE_TEMPLATE=""".. Autogenerated by genmods.py -- DO NOT EDIT --

******************************************************************************
%(prefix)s%(module)s - %(title)s
******************************************************************************

.. currentmodule:: %(package)s.%(module)s

.. autosummary::
   :nosignatures:

   %(members)s

.. automodule:: %(package)s.%(module)s
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

"""

INDEX_TEMPLATE=""".. Autogenerated by genmods.py -- DO NOT EDIT --

.. _api-index:

##############################################################################
Reference
##############################################################################

.. only:: html

   :Release: |version|
   :Date: |today|

.. toctree::
   :hidden:

   %(rsts)s
**Modules defined within Bumps**

.. currentmodule:: %(package)s

.. autosummary::

   %(mods)s

"""

def getmembers(package, module):
    name = package+"."+module
    __import__(name)
    M = sys.modules[name]
    try:
        L = M.__all__
    except:
        L = [s for s in dir(M)
             if inspect.getmodule(getattr(M,s)) == M and not s.startswith('_')]
    return L

def genfiles(package, modules, dir='api', absolute=True):

    prefix = package+"." if absolute else ""
    if not os.path.exists(dir):
        os.makedirs(dir)
    for (module, title) in modules:
        members = "\n    ".join(getmembers(package, module))
        with open(os.path.join(dir,module+'.rst'), 'w') as f:
            f.write(MODULE_TEMPLATE%locals())
    rsts = "\n   ".join(module+'.rst' for module,_ in modules)
    mods = "\n   ".join(prefix+module for module,_ in modules)

    with open(os.path.join(dir,'index.rst'),'w') as f:
        f.write(INDEX_TEMPLATE%locals())


modules=[
    #('__init__', 'Top level namespace'),
    ('bspline', 'B-Spline interpolation library'),
    ('cheby', 'Freeform - Chebyshev'),
    ('cli', 'Command line interface'),
    #('corrtest', 'Test for residual structure'),
    ('errors','Plot sample profile uncertainty'),
    ('fitproblem', 'Interface between models and fitters'),
    ('fitservice', 'Remote job plugin for fit jobs'),
    ('fitters', 'Wrappers for various optimization algorithms'),
    ('initpop', 'Population initialization strategies'),
    #('interface', 'Models of interfacial roughness'),
    ('mapper', 'Parallel processing implementations'),
    ('mono', 'Freeform - Monotonic Spline'),
    ('numpyerrors', 'Decorator for function level error behaviour'),
    ('partemp', 'Parallel tempering optimizer'),
    #('plottable', 'Style-based plot definitions'),
    ('pytwalk', 'MCMC error analysis using T-Walk steps'),
    ('quasinewton', 'BFGS quasi-newton optimizer'),
    ('random_lines', 'random lines and particle swarm optimizers'),
    ('rebin','1D and 2D rebinning'),
    ('bumpsmodule','Low level calculations'),
    ('simplex','Nelder-Mead simplex optimizer (amoeba)'),
    ('stitch', 'Overlapping curve stitching'),
    ('support', 'Environment support'),
    ('util','Miscellaneous functions'),
    ('wsolve','Weighted linear and polynomial solver with uncertainty'),
    ('parameter', 'Parameters'),
    ('bounds', 'Bounds'),
    ('formatnum', 'Format numbers'),
]
package='bumps'
genfiles(package, modules)