import numpy
import bumps
from bumps import pmath
from bumps.parameter import Parameter,flatten
from bumps.mono import monospline
from random import randint,randrange, sample,uniform
#from logilab.common.compat import izip
from bumps.bounds import Bounds


class ConstantModel(object):
    def __init__(self, C, x, y, dy):
        self.C = Parameter(C,name="constant")
        self.x,self.y,self.dy = [numpy.array(v) for v in x,y,dy]

    def numpoints(self):
        return len(self.x)

    def parameters(self):
        return self.C

    def theory(self):
        #return self.parts[0](self.X,self.Y)
        #parts = [M(self.X,self.Y) for M in self.parts]
        #for i,p in enumerate(parts):
        #    if np.any(np.isnan(p)): print "NaN in part",i
        return numpy.ones_like(self.x)*self.C.value
    
    def residuals(self):
        #if np.any(self.err ==0): print "zeros in err"
        return (self.theory()-self.y)/(self.dy)

    def nllf(self):
        R = self.residuals()
        #if np.any(np.isnan(R)): print "NaN in residuals"
        return 0.5*numpy.sum(R**2)

    def plot(self, view='linear'):
        import pylab
        pylab.errorbar(self.x, self.y, yerr=self.dy, fmt='x')
        pylab.plot(self.x, self.theory(), '-')

    def simulate_data(self, noise):
        self.y = self.theory() + numpy.random.randn(*self.x.shape)*noise

    def save(self, basename):
        pass

    def update(self):
        pass

class LinearModel(object):
    def __init__(self, A, B, x, y, dy):
        self.A = Parameter(A,name="slope")
        self.B = Parameter(B,name="intercept")
        self.x,self.y,self.dy = [numpy.array(v) for v in x,y,dy]

    def numpoints(self):
        return len(self.x)

    def parameters(self):
        return dict(A=self.A,B=self.B)

    def theory(self):
        #return self.parts[0](self.X,self.Y)
        #parts = [M(self.X,self.Y) for M in self.parts]
        #for i,p in enumerate(parts):
        #    if np.any(np.isnan(p)): print "NaN in part",i
        return self.A.value*self.x + self.B.value
    
    def simulate_data(self, noise):
        self.y = self.theory() + numpy.random.randn(*self.x.shape)*noise

    def residuals(self):
        #if np.any(self.err ==0): print "zeros in err"
        return (self.theory()-self.y)/(self.dy)

    def nllf(self):
        R = self.residuals()
        #if np.any(np.isnan(R)): print "NaN in residuals"
        return 0.5*numpy.sum(R**2)

    def plot(self, view='linear'):
        import pylab
        pylab.errorbar(self.x, self.y, yerr=self.dy, fmt='x')
        pylab.plot(self.x, self.theory(), '-', hold=True)

        
    def save(self, basename):
        pass

    def update(self):
        pass


class Polynomial(object):
    r""" 
    the list of coefficients, where $A = [a_{i}, a_{i-1}, ... a_{0}]$, specifies the degree
    a polynomial of degree len(a)-1 is created with the initialization coefficients in A
    such that $a_i$ where $i = 0, ... , len(a)-1$ created a polynomial of form
    $a_{i}x^{len(A)-1}, a_{i}x^{len(A)-2}, ... , a_{0}x^{0}$
    
    Ex. given list [1,2,3,4,0,1] it would make the polynomial
    $1x^5 2x^4 3x^3 4x^2+ 0x+ 1$
    """
    def __init__(self, l, x, y, dy):
        self.li = l
        self.degree = len(l)-1
        self.x,self.y,self.dy = [numpy.array(v) for v in x,y,dy]
        self.polys = [Parameter(v, name='x{}'.format(self.degree-k)) for k,v in enumerate(self.li)]
    def numpoints(self):
        return len(self.x)

    def parameters(self):
        import pdb
        #pdb.set_trace()
        #print 'P',p
        return self.polys

    def theory(self):
        #return self.parts[0](self.X,self.Y)
        #parts = [M(self.X,self.Y) for M in self.parts]
        #for i,p in enumerate(parts):
        #    if np.any(np.isnan(p)): print "NaN in part",i
        return numpy.polyval([p.value for p in self.polys],self.x)
    
    def residuals(self):
        #if np.any(self.err ==0): print "zeros in err"
        return (self.theory()-self.y)/(self.dy)

    def nllf(self):
        R = self.residuals()    
        #if np.any(np.isnan(R)): print "NaN in residuals"
        return 0.5*numpy.sum(R**2)
    def __call__(self):
        raise NotImplementedError

    def plot(self, view='linear'):
        import pylab      
        pylab.errorbar(self.x, self.y, yerr=self.dy, fmt='x')
        pylab.plot(self.x, self.theory(), '-',hold=True)

    def simulate_data(self, noise):
        pass

    def save(self, basename):
        pass

    def update(self):
        pass
class Point(object):
    def __init__(self,x,y,name = "",xrange = [0,10],yrange=[0,10]):
        self.x = Parameter(x,name=name+'x', bounds = xrange)
        self.y = Parameter(y,name=name+'y',bounds=yrange)
        self.value = [self.x,self.y]
    def getParameters(self):
        return [self.x, self.y]
    def __str__(self):
        return str("Point({}, {})".format(self.x,self.y))
    
def randomrange(a,b,k):#range on [a,b) and number to pull from the range
    return sorted(sample(xrange(a,b)))
    
def invx(data_points):
    x = numpy.linspace(-1,1,data_points)
    y = 1/numpy.cos(x)
    return x,y
    
    
def randomspline(n,p,a,b,nsteps=100,distribution=numpy.random.random):    
    """ 
        increasing monotonically only
        n control points, x random points on the range [a,b) where c is the step
    """
    #x = numpy.arange(a,b,step)
    #xt = numpy.random.uniform(min(x),max(x),p)
    Cx = numpy.cumsum(numpy.random.exponential(size=n))
    Cx = Cx/Cx[-1]*(b-a)+a
    Cy = distribution(Cx.shape)
    xt = numpy.linspace(min(Cx),max(Cx),p)
    #print (xt,monospline(Cx,Cy,xt),Cx,Cy)
    return (xt,monospline(Cx,Cy,xt),Cx,Cy)

    #uniform(xr[z-1],xr[z]), uniform(miny,maxy)

class Monospline(object):
    """
        Monospline Interpolation for Control Points given a spline.
    """
    def __init__(self, n,xt,y,dy,name=''):
        """
        n = number of control points to use
        xt,y = input and output data for a spline p, in which p($x_i$) = $y_i$
        """
        self.xt = xt
        self.y = y
        self.dy = dy
        minx,maxx = min(xt),max(xt)
        miny,maxy = min(y),max(y)
        xr = numpy.linspace(minx,maxx,n+1)
        #print 'xr',xr
        self.cpoints = [Point(uniform(xr[z-1],xr[z]), uniform(miny,maxy),"C{}".format(z),
                              xrange=xr[z-1:z+1]) for z in range(1,n+1)]
        
        #print 'CPOINTS',[ (z.x.value, z.y.value) for z in self.cpoints]
        #self.cpoints = sorted(self.cpoints, cmp = lambda a,b: cmp(a.x.value, b.x.value))
        #print 'CPOINTS',[ (z.x.value, z.y.value) for z in self.cpoints]
    def parameters(self):
        return [z.getParameters() for z in self.cpoints]
    def flattened(self):
        return flatten(self.parameters())
    def numpoints(self):
        return len(self.y)
    def theory(self):
        return monospline([q.x.value for q in self.cpoints], [q.y.value for q in self.cpoints], self.xt)
    def residuals(self):
        return (self.theory()-self.y)/(self.dy)
    def nllf(self):
        R = self.residuals()
        #if np.any(np.isnan(R)): print "NaN in residuals"
        return 0.5*numpy.sum(((self.theory()-self.y)/(self.dy))**2)
#    def __call__(self):
#        raise NotImplementedError
    def plot(self, view='linear'):
        import pylab      
        xpts = [w.x.value for w in self.cpoints]
        ypts = [w.y.value for w in self.cpoints]
        #pylab.errorbar(xpts, ypts, yerr=self.dy, fmt='x')
        for i,(x,y) in enumerate(zip(xpts,ypts)):
            pylab.text(x,y,str(i+1))
        pylab.plot(self.xt,self.y,hold=True)
        pylab.plot(self.xt,self.theory(),hold=True)


class MonosplineInterval(object):
    """
        Monospline Interpolation for Control Points given a spline in which the control points vary as an interval
    """
    def __init__(self, n,xt,y,dy,name='',targetpoints=None):
        """
        n = number of control points to use
        xt,y = input and output data for a spline p, in which p($x_i$) = $y_i$
        """
        self.rmax = 100000
        self.rmin = 0
        offset=.01
        self.control_points = targetpoints
        self.n = n
        idx = numpy.argsort(xt)
        self.xt = numpy.asarray(xt)[idx]
        self.y = numpy.asarray(y)[idx]
        self.dy = numpy.asarray(dy)[idx]
        #print 'N',n
        #bound lower and bound higher
        self.bl =self.xt[0]
        self.bh = self.xt[-1]
        self.ymin = min(y)
        self.ymax = max(y)
        rs = numpy.random.uniform(0,1000,n-2)
        #print 'GENERATED RS',rs,sum(rs)
        #sumrs = 1+sum(rs)
        #rs /= sumrs
        #rs *= (self.bh - self.bl)
        #print 'RS',rs,sum(rs)
        #rs.sort()
        self.Cratios = [Parameter(rsi,name = 'ratio{}'.format(i+1), bounds = [self.rmin,self.rmax]) for i,rsi in enumerate(rs)]
        self.Cys = [Parameter(0,name='y{}'.format(w+1), bounds = [self.ymin,self.ymax]) for w in range(len(rs))]
#        self.cpoints = []
#        for r in range(len(self.ldx)):
#            #xval = sum([z.value for z in self.ldx[:r]])
#            yval = 0
#            self.cpoints.append(IPoint(yval,name="{}".format(r), vallist=self.ldx[:r]))
        
    def parameters(self):
        return zip(self.Cratios,self.Cys)
    def flattened(self):
        return flatten(self.parameters())
    def numpoints(self):
        return len(self.y)
    @property
    def Cx(self):
        xfirst, xlast = self.xt[0], self.xt[-1]
        ratios = numpy.hstack((1, [p.value for p in self.Cratios]))
        dx = ratios / sum(ratios) * (xlast-xfirst)
        return numpy.cumsum(numpy.hstack((xfirst, dx)))
        #print 'Q',[q.x.value.value for q in self.cpoints]
        #print 'LDX',[z.value for z in self.ldx],'\n','Qx',[q.x for q in self.cpoints],'\n','Qy',[q.y.value for q in self.cpoints],'\n'
    @property
    def Cy(self):
        return numpy.hstack((self.y[0],[y.value for y in self.Cys] ,self.y[-1]))

    def theory(self):
        return monospline(self.Cx, self.Cy, self.xt)
    def residuals(self):
        return (self.theory()-self.y)/(self.dy)
    def nllf(self):
        R = self.residuals()
        #if np.any(np.isnan(R)): print "NaN in residuals"
        return 0.5*numpy.sum(((self.theory()-self.y)/(self.dy))**2)
#    def __call__(self):
#        raise NotImplementedError
    def plot(self, view='linear'):
        import pylab      
        #pylab.errorbar(xpts, ypts, yerr=self.dy, fmt='x')
        #print 'PTS',xpts,ypts
        for i,(x,y) in enumerate(zip(self.Cx,self.Cy)):
            pylab.text(x,y,str(i))
        #print 'SELFXT',self.xt,'----------------','SELFTHEORY',self.theory()
        #print '\n'
        #print 'REAL CONTROLS AT',zip(self.control_points[0],self.control_points[1])
        #print 'CURRENT CONTROLS AT',zip(self.Cx, self.Cy)
        if self.control_points:
            pylab.plot(self.control_points[0],self.control_points[1],'bs',hold=True)
        pylab.plot(self.xt,self.y,hold=True)
        #pylab.plot(self.xt,self.theory(),hold=True)
        pylab.plot(self.xt,self.theory(),hold=True)
#

