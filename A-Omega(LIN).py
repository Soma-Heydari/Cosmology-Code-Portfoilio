import scipy as sp
import numpy as np
import simpy as sm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.integrate import dblquad,nquad
from numpy import heaviside, log, sqrt, pi, inf, exp, cos, sin#, long,double

xa, ya = np.loadtxt('D:/Sci-F/my papers/paper6/revise-APJ/codes of referee case/natural-final codes-new reh&GWs&referee case/PkF.dat').T
##xb, yb =np.loadtxt('C:/Users/milad/Desktop/GI+natural/NEW_Natural/k_caseB.dat').T
##xc, yc =np.loadtxt('C:/Users/milad/Desktop/GI+natural/NEW_Natural/k_caseC558new.dat').T


pra = interp1d(xa,ya)
##prb = interp1d(xb,yb)
##prc = interp1d(xc,yc)
##prd = interp1d(xd,yd)


##plt.plot(xa,pra(xa),'-b' ,xb, prb(xb),'-k', xc,prc(xc),'-r')
##plt.xscale("log")
##plt.yscale("log")
##plt.legend(['case A', 'case B', 'case C'], loc='best')
##plt.title('Power spectrum in G-inflation model')
##plt.show()

def Si(x):
    from scipy.special import sici
    return sici(x)[0]

def Ci(x):
    from scipy.special import sici
    return sici(x)[1]

def H(x):
    return np.piecewise(x,[x<0,x>=0],[0,1])


##(-27*(-3 + u**2 + v**2))/(u**2*v**2) + (27*(-3 + u**2 + v**2)**2*(Ci((1 - (u - v)/sqrt(3))*x) + Ci((1 + (u - v)/sqrt(3))*x) - Ci((1 + (u + v)/sqrt(3))*x) - Ci(x*abs(1 - (u + v)/sqrt(3))) + log(abs((3 - (u + v)**2)/(3 - (u - v)**2)))))/(4.*u**3*v**3) + (27*(48*u*v*x**2*cos((u*x)/sqrt(3))*cos((v*x)/sqrt(3))*(-3*cos(x) + x*sin(x)) + 24*x*(6 - (3 - u**2 - v**2)*x**2)*sin(x)*sin((u*x)/sqrt(3))*sin((v*x)/sqrt(3)) - 48*sqrt(3)*x**2*sin(x)*(v*cos((v*x)/sqrt(3))*sin((u*x)/sqrt(3)) + u*cos((u*x)/sqrt(3))*sin((v*x)/sqrt(3))) + 8*sqrt(3)*x*cos(x)*(v*(18 - (3 + u**2 - v**2)*x**2)*cos((v*x)/sqrt(3))*sin((u*x)/sqrt(3)) + u*(18 - (3 - u**2 + v**2)*x**2)*cos((u*x)/sqrt(3))*sin((v*x)/sqrt(3))) + 24*(-18 + (3 + u**2 + v**2)*x**2*cos(x)*sin((u*x)/sqrt(3))*sin((v*x)/sqrt(3)))))/(8.*u**3*v**3*x**4)
def Ts(u,v,x):
    return (-27*(-3 + u**2 + v**2))/(u**2*v**2) + (27*(-3 + u**2 + v**2)**2*(Ci((1 - (u - v)/sqrt(3))*x) + Ci((1 + (u - v)/sqrt(3))*x) - Ci((1 + (u + v)/sqrt(3))*x) - Ci(x*abs(1 - (u + v)/sqrt(3))) + log(abs((3 - (u + v)**2)/(3 - (u - v)**2)))))/(4.*u**3*v**3) + (27*(48*u*v*x**2*cos((u*x)/sqrt(3))*cos((v*x)/sqrt(3))*(-3*cos(x) + x*sin(x)) + 24*(-18 + (3 + u**2 + v**2)*x**2)*cos(x)*sin((u*x)/sqrt(3))*sin((v*x)/sqrt(3)) + 24*x*(6 - (3 - u**2 - v**2)*x**2)*sin(x)*sin((u*x)/sqrt(3))*sin((v*x)/sqrt(3)) - 48*sqrt(3)*x**2*sin(x)*(v*cos((v*x)/sqrt(3))*sin((u*x)/sqrt(3)) + u*cos((u*x)/sqrt(3))*sin((v*x)/sqrt(3))) + 8*sqrt(3)*x*cos(x)*(v*(18 - (3 + u**2 - v**2)*x**2)*cos((v*x)/sqrt(3))*sin((u*x)/sqrt(3)) + u*(18 - (3 - u**2 + v**2)*x**2)*cos((u*x)/sqrt(3))*sin((v*x)/sqrt(3)))))/(8.*u**3*v**3*x**4)


##(-27*(-48*u*v*x**2*cos((u*x)/sqrt(3))*cos((v*x)/sqrt(3))*(x*cos(x) + 3*sin(x)) + 24*x*(-6 + (3 - u**2 - v**2)*x**2)*cos(x)*sin((u*x)/sqrt(3))*sin((v*x)/sqrt(3)) + 48*sqrt(3)*x**2*cos(x)*(v*cos((v*x)/sqrt(3))*sin((u*x)/sqrt(3)) + u*cos((u*x)/sqrt(3))*sin((v*x)/sqrt(3))) + 8*sqrt(3)*x*sin(x)*(v*(18 - (3 + u**2 - v**2)*x**2)*cos((v*x)/sqrt(3))*sin((u*x)/sqrt(3)) + u*(18 - (3 - u**2 + v**2)*x**2)*cos((u*x)/sqrt(3))*sin((v*x)/sqrt(3))) + 24*(-18 + (3 + u**2 + v**2)*x**2*sin(x)*sin((u*x)/sqrt(3))*sin((v*x)/sqrt(3)))))/(8.*u**3*v**3*x**4) - (27*(-3 + u**2 + v**2)**2*(Si((1 - (u - v)/sqrt(3))*x) + Si((1 + (u - v)/sqrt(3))*x) - Si((1 - (u + v)/sqrt(3))*x) - Si((1 + (u + v)/sqrt(3))*x)))/(4.*u**3*v**3)
def Tc(u,v,x):
    return (-27*(-48*u*v*x**2*cos((u*x)/sqrt(3))*cos((v*x)/sqrt(3))*(x*cos(x) + 3*sin(x)) + 24*x*(-6 + (3 - u**2 - v**2)*x**2)*cos(x)*sin((u*x)/sqrt(3))*sin((v*x)/sqrt(3)) + 24*(-18 + (3 + u**2 + v**2)*x**2)*sin(x)*sin((u*x)/sqrt(3))*sin((v*x)/sqrt(3)) + 48*sqrt(3)*x**2*cos(x)*(v*cos((v*x)/sqrt(3))*sin((u*x)/sqrt(3)) + u*cos((u*x)/sqrt(3))*sin((v*x)/sqrt(3))) + 8*sqrt(3)*x*sin(x)*(v*(18 - (3 + u**2 - v**2)*x**2)*cos((v*x)/sqrt(3))*sin((u*x)/sqrt(3)) + u*(18 - (3 - u**2 + v**2)*x**2)*cos((u*x)/sqrt(3))*sin((v*x)/sqrt(3)))))/(8.*u**3*v**3*x**4) - (27*(-3 + u**2 + v**2)**2*(Si((1 - (u - v)/sqrt(3))*x) + Si((1 + (u - v)/sqrt(3))*x) - Si((1 - (u + v)/sqrt(3))*x) - Si((1 + (u + v)/sqrt(3))*x)))/(4.*u**3*v**3)


##(((3*pi*(-3 + u**2 + v**2)**2*H(-sqrt(3) + u + v))/(4.*u**3*v**3) + Tc(u,v,1)/9.)**2 + ((27*(-3 + u**2 + v**2))/(u**2*v**2) - (27*(-3 + u**2 + v**2)**2*log(abs((3 - (u + v)**2)/(3 - (u - v)**2))))/(4.*u**3*v**3) + (Ts(u,v,1))**2)/81.)/2.
def I2RD(u,v):
    return 0.5*(((3*pi*(-3 + u**2 + v**2)**2*H(-sqrt(3) + u + v))/(4.*u**3*v**3) + Tc(u,v,1)/9.)**2 + ((27*(-3 + u**2 + v**2))/(u**2*v**2) - (27*(-3 + u**2 + v**2)**2*log(abs((3 - (u + v)**2)/(3 - (u - v)**2))))/(4.*u**3*v**3) + Ts(u,v,1))**2/81.)
##
def omega_A(u,v):
    return (zarib*(4*v**2 - (1 - u**2 + v**2)**2)**2*I2RD(u,v)*pra(xa[i]*u)*pra(xa[i]*v))/(96.*u**2*v**2)
##def omega_B(u,v):
##    return (zarib*(4*v**2 - (1 - u**2 + v**2)**2)**2*I2RD(u,v)*prb(xb[i]*u)*prb(xb[i]*v))/(96.*u**2*v**2)
##def omega_C(u,v):
##    return (zarib*(4*v**2 - (1 - u**2 + v**2)**2)**2*I2RD(u,v)*prc(xc[i]*u)*prc(xc[i]*v))/(96.*u**2*v**2)
##def omega_D(u,v):
##    return (zarib*(4*v**2 - (1 - u**2 + v**2)**2)**2*I2RD(u,v)*prd(xd[i]*u)*prd(xd[i]*v))/(96.*u**2*v**2)


kmax=xa[len(xa)-1]
kmin=xa[0]

def bound_v():
    return [0,np.inf]

def bound_u(v):
    return [np.abs(1-v),np.abs(1+v)]


beta=1.546*(10**(-15))
zarib=A=0.83*(9.93**(-1/3))*4.2*10**(-5)

b1=[]
b2=[]
b3=[]
b4=[]

c1=[]
c2=[]
c3=[]
c4=[]


##print(pra(xa[7]))
##i=1;
##print(omega_A(1,1))
##print(I2RD(1,1))


##print(len(xb))
##print(xc[0])
##i=35
##a1,e1=integrate.nquad(omega_A,[bound_u,bound_v])
##b1.append(a1)
##s1=beta*xa[i]
##c1.append(s1)
##print(s1,a1)


###PTA:50****LIGO:45****ogle,total:60-70
print("*************set A**********")
for i in range(91,190):
     # a1,e1=integrate.nquad(omega_A,[bound_u,bound_v])
    a1,e1=integrate.nquad(omega_A,[bound_u,[kmin/xa[i],xa[176]/xa[i]]])
    b1.append(a1)
    s1=beta*xa[i]
    c1.append(s1)
    print(s1,a1)
print("*********10****set B********90********")
##for i in range(10,12):
####    a2,e2=integrate.nquad(omega_B,[bound_u,bound_v])
##    b2.append(a2)
##    s2=beta*xb[i]
##    c2.append(s2)
##    print(s2,a2)
##print("*****15********set C******40**********")
##for i in range(20,80):
##    a3,e3=integrate.nquad(omega_C,[bound_u,bound_v])
##    b3.append(a3)
##    s3=beta*xc[i]
##    c3.append(s3)
##    print(s3,a3)


omegaf_A = interp1d(c1,b1)
##omegaf_B = interp1d(c2,b2)
##omegaf_C = interp1d(c3,b3)
##omegaf_D = interp1d(c4,b4)

##plt.plot(c1,omegaf_A(c1), '-b',c2,omegaf_B(c2), '-k',c3,omegaf_C(c3), '-r')
plt.plot(c1,omegaf_A(c1), '-b')
plt.xscale("log")
plt.yscale("log")
plt.axis([10**(-11), 10000, 10**(-19), 10**(-6)])
plt.legend(['case A', 'case B', 'case C'], loc='best')
plt.title('Energy spectra of secondary GWs')
plt.show()

