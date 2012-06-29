import numpy
import bumps
from bumps.parameter import Parameter
from bumps.mono import monospline
from random import randint

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
    def __init__(self,x,y):
        if (isinstance(x, int) and isinstance(y,int)) == False:
            raise Not
        self.x = Parameter(x,name=name+'x')
        self.y = Parameter(y,name=name+'y')
        self.value = [self.x,self.y]
    def getParameters(self):
        return [self.x, self.y]
    
class Monospline(object):
    """
        Monospline Interpolation for Control Points given a spline.
    """
    def __init__(self, n,x,y,name=''):
        """
        n = number of control points to use
        x,y = input and output data for a spline p, in which p($x_i$) = $y_i$
        x and y are respectively lists of (x,y) points to pass to the monospline function. xt is a list
        of the control point x values to ask for back        """
        if all(isinstance(z,list) for z in (x,y,xt)):
            self.Points  = [Parameter(0,name=name+"x Values")]
        else:
            raise NotImplementedError('Only lists for the values for the x and y control points can be sent at the moment')
        
        def parameters(self):
            return [self.x,self.y]
        

#    def numpoints(self):
#        return len(self.x)
#
#    def parameters(self):
#        import pdb
#        #pdb.set_trace()
#        #print 'P',p
#        return self.polys
#
#    def theory(self):
#        #return self.parts[0](self.X,self.Y)
#        #parts = [M(self.X,self.Y) for M in self.parts]
#        #for i,p in enumerate(parts):
#        #    if np.any(np.isnan(p)): print "NaN in part",i
#        return numpy.polyval([p.value for p in self.polys],self.x)
#    
#    def residuals(self):
#        #if np.any(self.err ==0): print "zeros in err"
#        return (self.theory()-self.y)/(self.dy)
#
#    def nllf(self):
#        R = self.residuals()    
#        #if np.any(np.isnan(R)): print "NaN in residuals"
#        return 0.5*numpy.sum(R**2)
#    def __call__(self):
#        raise NotImplementedError
#
#    def plot(self, view='linear'):
#        import pylab      
#        pylab.errorbar(self.x, self.y, yerr=self.dy, fmt='x')
#        pylab.plot(self.x, self.theory(), '-',hold=True)
#
#    def simulate_data(self, noise):
#        pass
#
#    def save(self, basename):
#        pass
#
#    def update(self):
#        pass



