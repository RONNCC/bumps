import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from bumps.mono import monospline
from bumps.names import FitProblem
from numpy import inf, ones_like,sort,random
import numpy
from itertools import izip
import splines
#sysarg 1 -> # of control points
#sysarg 2 -> seed
#sysarg 3 -> Model Type
#sysarg 4 -> number of sample points
num_points = int(sys.argv[1]) if len(sys.argv)>1 else 15
print num_points,sys.argv
if len(sys.argv)>2:
    seed = int(sys.argv[2])
else:
    seed = numpy.random.randint(10000)
model = sys.argv[3] if len(sys.argv)>3 else 'intmonosp'
print "Seed=",seed
random.seed(seed)

sample_points = sys.argv[4] if len(sys.argv)>4 else 10

xvals = [1,2,3,4,5,6,7]
data = [2,5,10,17,26,37,50]
uncertainty = [1]*len(xvals)
uncert = .01

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
    range = [0,10]
    x,y,Cx,Cy = splines.randomspline(num_points,sample_points, range[0],range[1])
    #print 'Values from randomly generated spline',rspline[0],rspline[1]
    #print 'Actual Control Points',sorted(list(izip(rspline[2],rspline[3])), cmp = lambda x,y: cmp(x[0],y[0]) )
    M = splines.Monospline(num_points,x,y,uncert*numpy.ones_like(x))
    for p in M.parameters():
        p[1].range(min(y),max(y))
    return M
def intmonosp():#interval monospline
    range = [0,10]
    x,y,Cx,Cy = splines.randomspline(num_points,sample_points, range[0],range[1])
    #x,y = splines.invx(1000)
    #print 'Values from randomly generated spline',rspline[0],rspline[1]
    #print 'Actual Control Points',sorted(list(izip(rspline[2],rspline[3])), cmp = lambda x,y: cmp(x[0],y[0]) )
    M = splines.MonosplineInterval(num_points,x,y,uncert*numpy.ones_like(x),targetpoints=(Cx,Cy))
    #for p in M.parameters():
    #    p[1].range(-2,2)
    return M
#M = polynomial()
M = eval(model)()
#M = monosp()
#M = intmonosp()
#M = polynomial()    
problem = FitProblem(M)

#print monospline([1,2,3,4,5,6,7],[2,4,6,8,10,12,14],[8,9,10,11,12])
