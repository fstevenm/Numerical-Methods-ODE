import numpy as np
import matplotlib.pyplot as plt

def f(t,y):
    return (t-y)/2

def RK4(f,a,b,h,y0):
    t = a + np.arange(4)*h
    y = np.zeros(4)
    k1 = np.zeros(4)
    k2 = np.zeros(4)
    k3 = np.zeros(4)
    k4 = np.zeros(4)
    y[0] = y0
    
    for j in range(1,4):
        k1[j] = f(t[j-1],y[j-1])
        k2[j] = f(t[j-1]+0.5*h, y[j-1]+0.5*k1[j]*h)
        k3[j] = f(t[j-1]+0.5*h, y[j-1]+0.5*k2[j]*h)
        k4[j] = f(t[j-1]+h, y[j-1]+k3[j]*h)
        y[j] = y[j-1] + (k1[j]+2*k2[j]+2*k3[j]+k4[j])*(h/6)
        
    return y

def abm(f,a,b,h,y0):
    n = (b-a)/h + 1
    
    T = a+np.arange(n)*h
    Y = np.zeros(len(T))
    P = np.zeros(len(T))
    
#    a, b, h, y0 = 0, 3, 0.125, 1
    
    y = RK4(f,a,b,h,y0)
    Y[0]=y[0]
    
    for j in range(1,4):
        Y[j] = y[j]
        
    for i in range(1, len(T)-3):
        j=i+3
        
        P[j] = Y[j-4] + (4*h/3)*(2*f(T[j-3],Y[j-3]) - f(T[j-2],Y[j-2]) + 2*f(T[j-1],Y[j-1]))
        Y[j] = Y[j-2] + (h/3)*(f(T[j-2],Y[j-2]) + 4*f(T[j-1],Y[j-1]) + f(T[j],P[j]))
    
    return T, Y, P

a, b, h, y0 = 0, 3, 0.125,1
T, Y, P = abm(f,a,b,h,y0)

#a,b,h,y0 = 0, 3, 0.125, 1
#t2, y2, p2 = abm(f,a,b,h,y0)

#Penyelesaian eksak
ye = 3*np.exp(-T/2)-2+T
#galat relatif
e = (Y-ye)/ye

plt.close('all')
plt.plot(T,Y,'-ob',label='ABM h=0.125')
plt.plot(T,ye,'r',label='eksak')
plt.xlabel('t')
plt.ylabel('y')
plt.axis([0,3,0.8,1.6])
plt.legend(loc='upper left')
plt.show()

print('Penyelesaian Metode Milne Simpson')
print(' n\t  t\t  pre\t y_num\ty_eksak\t e_rel')
print('----------------------------------------------')

for j,T in(enumerate(T)):
    print ("%2.0f"% (j),"\t%5.3f"% (T),"\t%5.4f"% (P[j]),"\t%5.4f"% (Y[j]),"\t%5.4f"% (ye[j]),"\t%5.4f"% (e[j])) 
print('----------------------------------------------')
