"""
Exported names

In model definition scripts, rather than importing symbols one by one, you
can simply perform:

    from bumps.names import *

This is bad style for library and applications but convenient for
model scripts.

Numpy is available as np, and sys is available for sys.argv
"""
import sys
import numpy as np

from . import pmath
from .parameter import Parameter
from .modelfn import ModelFunction
from .fitproblem import preview, fit, mesh, FitProblem, MultiFitProblem

