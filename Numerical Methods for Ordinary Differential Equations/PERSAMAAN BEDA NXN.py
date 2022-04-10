#import numpy as np
#from numpy.linalg import *
#
#A = [[1.9375,-0.75,0],
#     [-1.25, 1.9375, -0.75],
#     [0, -1.25, 1.9375]]
#A = np.array(A)
#
#B = [[6.3085],
#     [0.0469],
#     [7.5274]]
#
#print((inv(A)))

import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt
plt.style.use('seaborn-poster')

def p(t):
    return 2*t

def q(t):
    return -1

def r(t):
    return t**2 - 1

# Get A
def mdh(a,b,h,At,Bt):
    n = int((b-a)/h) + 1

    t = np.linspace(a,b,n)
    x = np.zeros(n)
    y = np.zeros(n)
    
    A = np.zeros((n-2, n-2))
    N = len(A)
    
    for i in range(0, N):
        if i == 0 :
            A[0, 0] = 2+(h**2)*q(t[i+1])
            A[0, 1] = (h/2)*p(t[i+1])-1
        elif i == N-1:
            A[N-1, N-2] = (-h/2)*p(t[i+1])-1
            A[N-1, N-1] = 2+(h**2)*q(t[i+1])
        else:
            A[i, i-1] = (-h/2)*p(t[i+1])-1
            A[i, i] = 2+(h**2)*q(t[i+1])
            A[i, i+1] = (h/2)*p(t[i+1])-1
        
    B = np.zeros((n-2,1))
    for i in range(0,N):
        if i == 0 :
            B[0, 0] = -(h**2)*r(t[i+1])-(-(h/2)*p(t[i+1])-1)*At
        elif i == N-1:
            B[N-1, 0] = -(h**2)*r(t[i+1])-((h/2)*p(t[i+1])-1)*Bt
        else:
            B[i, 0] = -(h**2)*r(t[i+1])
    
    print(A)
    print(B)
    
    X = np.dot(inv(A),B)
    return X
   

print(mdh(0,1,0.25,5,10))    
    
#plt.figure(figsize=(10,8))
#plt.plot(t, y)
#plt.plot(5, 50, 'ro')
#plt.xlabel('time (s)')
#plt.ylabel('altitude (m)')
#plt.show()

t = [0,0.25,0.5,0.75,1]
y = [5, 6.25, 7.25, 8.625,10]
plt.plot(t,y,'ro-')

plt.grid()
plt.show()