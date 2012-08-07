import pylab
from numpy import sin, cos, linspace, meshgrid, e, pi, sqrt, array, exp
from bumps.names import *

def prod(L):
    return reduce(lambda x,y: x*y, L, 1) 
def plot2d(fn,args,range=(-10,10)):
    def plotter(view=None, **kw):
        x,y = kw[args[0]],kw[args[1]]
        r = linspace(range[0],range[1],200)
        X,Y = meshgrid(x+r,y+r)
        kw['x'],kw['y'] = X,Y
        pylab.clf()
        pylab.pcolormesh(x+r,y+r,nllf(**kw))
        pylab.plot(x,y,'o',hold=True, markersize=6,markerfacecolor='red',markeredgecolor='black',markeredgewidth=1, alpha=0.7)
        pylab.colorbar()
    return plotter 
def fxy(fn):
    def cost(x,y): return fn((x,y))
    return cost
def fxyz(fn):
    def cost(x,y,z): return fn((x,y,z))
    return cost

def sin_plus_quadratic(x=0,y=0): 
    fx,fy = 2,3      # x,y frequency and between bowl barer height
    barrier = 2      # barrier height
    cx,cy = 3,1      # x,y center
    mx,my = 2,5      # x,y curvature
    width = 6        # size of the acceptance region
    return (barrier*(sin(fx*x) + sin(fy*y)+2) + ((x-cx)/mx)**2 + ((y-cy)/my)**2)/width

def ackley(x):
    n = len(x)
    return -20*exp(-0.2*sqrt(sum(xi**2 for xi in x)/n))-exp(sum(cos(2*pi*xi) for xi in x)/n) + 20 + e

def griewank(x):
    return 1 + sum(xi**2 for xi in x)**2/4000 - prod(cos(xi/sqrt(i+1)) for i,xi in enumerate(x))

def rastrigin(x):
    A = 10
    n = len(x)
    return A*n + sum(xi**2 - A*cos(2*pi*xi) for xi in x)

def uncoupled_rosenbrock(x):
    n = len(x)
    return sum(100*(x[2*i-1]**2 - x[2*i])**2 + (x[2*i-1] - 1)**2 for i in range(n/2))

def rosenbrock(x):
    n = len(x)
    return sum((1-x[i])**2 + 100*(x[i+1]-x[i]**2)**2 for i in range(n-1))

#nllf = sin_plus_quadratic
nllf = fxy(ackley)
#nllf = fxy(griewank)
#nllf = fxy(rastrigin)
#nllf = fxy(rosenbrock)
#nllf = fxyz(rosenbrock)
plot=plot2d(nllf,('x','y'),range=(-1,1))
M = ModelFunction(nllf,plot=plot)
for p in M.parameters().values():
    p.value = 200
    p.range(-200,200)
problem = FitProblem(M)
