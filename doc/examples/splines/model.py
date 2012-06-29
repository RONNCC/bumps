import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from bumps.mono import monospline
from bumps.names import FitProblem
from numpy import inf, ones_like

import splines

xvals = [1,2,3,4,5,6,7]
data = [2,5,10,17,26,37,50]
uncertainty = [1]*len(xvals)
def constant():
  M = splines.ConstantModel(C=5,x=range(len(data)), y=data, dy=uncertainty)
  M.C.range(0,10)
  return M
def linear():
    M = splines.LinearModel(A=0,B=5,x=range(len(data)), y=data, dy=ones_like(data))
    M.A.range(-10, 10) 
    M.B.range(-inf, inf)
    return M
def polynomial():
    M = splines.Polynomial([1,0,0],xvals,data,uncertainty)
    for p in M.parameters(): p.range(-20,20)
    return M
def monosp():
    pass
#M = polynomial()
#problem = FitProblem(M)
print monospline([1,2,3,4,5,6,7],[2,4,6,8,10,12,14],[8,9,10,11,12])
