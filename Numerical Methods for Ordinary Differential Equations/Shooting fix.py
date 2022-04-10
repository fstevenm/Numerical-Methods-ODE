import sympy as sym
import numpy as np

def f(t,u,z):
    return ((2*t)/(1+t**2))*z - ((2)/(1+t**2))*u + 1
    
def g(t,v,z):
    return ((2*t)/(1+t**2))*z - ((2)/(1+t**2))*v
    
def kutta4_u(u_0,z_0,a,b,h):
    M = int((b-a)/h)
    t = a
    u = u_0
    z = z_0
    u_k4 = np.array([u])
    z_k4 = np.array([z])
    while t < b:
        k1 = z
        l1 = f(t,u,z)
        
        k2 = (z + ((1/2)*l1*h))
        l2 = f((t+(1/2)*h),(u+((1/2)*k1*h)),(z+((1/2)*l1*h)))
        
        k3 = (z + ((1/2)*l2*h))
        l3 = f((t+(1/2)*h),(u+((1/2)*k2*h)),(z+((1/2)*l2*h)))
        
        k4 = (z+l3*h)
        l4 = f((t+h),(u+k3*h),(z+l3*h))
        
        u += (1/6)*h*(k1 + 2*k2 + 2*k3 + k4)
        z += (1/6)*h*(l1 + 2*l2 + 2*l3 + l4)
#        print(z)
        u_k4 = np.append(u_k4,u)
        z_k4 = np.append(z_k4,z)
        
        t += h
        
    return u_k4
        
def kutta4_v(v_0,z_0,a,b,h):
    M = int((b-a)/h)
    t = a
    v = v_0
    z = z_0
    v_k4 = np.array([v])
    z_k4 = np.array([z])
    while t < b:
        k1 = z
        l1 = g(t,v,z)
        
        k2 = (z + ((1/2)*l1*h))
        l2 = g((t+(1/2)*h),(v+((1/2)*k1*h)),(z+((1/2)*l1*h)))
        
        k3 = (z + ((1/2)*l2*h))
        l3 = h*g((t+(1/2)*h),(v+((1/2)*k2*h)),(z+((1/2)*l2*h)))
        
        k4 = (z+l3*h)
        l4 = h*g((t+h),(v+k3*h),(z+l3*h))
        
        v += (1/6)*h*(k1 + 2*k2 + 2*k3 + k4)
        z += (1/6)*h*(l1 + 2*l2 + 2*l3 + l4)
        
        v_k4 = np.append(v_k4,v)
        z_k4 = np.append(z_k4,z)
        
        t += h
    return v_k4
        
def shot(a,b,h,A,B):
    u_k4 = kutta4_u(A,0,a,b,h)
    v_k4 = kutta4_v(0,1,a,b,h)
    t_s = np.arange(a,b+h,h)
    x_s = np.array([])
    w_k = np.array([])
    for i in range(len(t_s)):
        w = ((B-u_k4[-1])/v_k4[-1])*v_k4[i]
        w_k = np.append(w_k,w)
        x = u_k4[i] + w
        x_s = np.append(x_s,x)
    
    print('\n-------------------------SOLUSI-------------------------')
    print('--------------Metode Tembakan Linier-------------')
    print('t_k\tu_k\tv_k\tw_k\tx_k')
    print('--------------------------------------------------------------')
    for i in range(len(t_s)):
        print('%.2f\t%.3f\t%.3f\t%.3f\t%.3f'% (t_s[i],u_k4[i],v_k4[i],w_k[i],x_s[i]))
        
        
        
shot(0,4,0.2,1.25,-0.95)