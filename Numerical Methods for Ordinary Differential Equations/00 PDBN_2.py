import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

t,y = sym.symbols('t y')

def f(t,y):
    return (t-y)/2
    
def euler(y_0,a,b,h):
    M = int((b-a)/h)
    t = a
    y = y_0
    t_e = np.array([t])
    y_e = np.array([y])
    while t < b:
        y += h*f(t,y)
        t += h
        t_e = np.append(t_e,t)
        y_e = np.append(y_e,y)
        
    print('\n-----------SOLUSI-----------')
    print('-------Metode Euler---------')    
    print('k\tt_k\ty_k')
    print('----------------------------')
    for i in range(M+1):
        print('%.f\t%.1f\t%.7f'% (i,t_e[i],y_e[i]))
    
    plt.plot(t_e,y_e,"--o",label="Metode Euler h={}".format(h))
        
def heun(y_0,a,b,h):
    M = int((b-a)/h)
    t = a
    y = y_0
    t_h = np.array([t])
    y_h = np.array([y])
    while t < b:
        t += h
        y += (h/2)*(f(t-h,y)+f(t,(y+h*f(t-h,y))))
        t_h = np.append(t_h,t)
        y_h= np.append(y_h,y)

    print('\n-----------SOLUSI-----------')
    print('--------Metode Heun---------')    
    print('k\tt_k\ty_k')
    print('----------------------------')
    for i in range(M+1):
        print('%.f\t%.1f\t%.7f'% (i,t_h[i],y_h[i]))
        
    plt.plot(t_h,y_h,"--o",label="Metode Heun h={}".format(h))
        
def taylor(y_0,a,b,h,n):
    t_ = sym.symbols('t_')
    y_ = sym.Function('y_')(t_)
    
    M = int((b-a)/h)
    t = a
    t_t = np.array([t])
    y_t = np.array([y_0])
    
    for i in range(M):
        D = sym.symbols('D')
        d_t = np.array([y_t[i]])
        A = f(t_,y_t[i])
        
        for j in range(n):
            if j==0:
                d_t = np.append(d_t,float((A.subs(t_,t_t[i]))))
                A = f(t_,y_)
                
            else:
                A_ = sym.diff(A,t_)
                A__ = sym.lambdify((D,y_),A_.subs([(sym.diff(y_,t_),D)]))
                d_t = np.append(d_t,A__(d_t[1],y_))
                A = A__(f(t_,y_),y_)
        p=0
        for l in range(len(d_t)):
            p += (d_t[l]*h**l)/sym.factorial(l)
        y_t = np.append(y_t,p)
        t += h
        t_t = np.append(t_t,t)
        
    print('\n-------------SOLUSI-------------')
    print('------Metode Deret Taylor-------')
    print('k\tt_k\ty_k')
    print('--------------------------------')
    for i in range(len(t_t)):
        print('%.f\t%.1f\t%.7f'% (i,t_t[i],y_t[i]))

    plt.plot(t_t,y_t,"--o",label="Metode Taylor h={} orde {}".format(h,n))
    
def kutta2mid(y_0,a,b,h):
    M = int((b-a)/h)
    t = a
    y = y_0
    t_k2 = np.array([t])
    y_k2 = np.array([y])
    while t < b:
        k1 = f(t,y)
        k2 = f((t+(1/2)*h),(y+((1/2)*k1*h)))
        y += k2*h
        y_k2 = np.append(y_k2,y)
        t += h
        t_k2 = np.append(t_k2,t)
        
    print('\n-------------SOLUSI-------------')
    print('------Metode Rung-Kutta Orde 2-------\n --------------Titik Tengah--------------')
    print('k\tt_k\ty_k')
    print('--------------------------------')
    for i in range(len(t_k2)):
        print('%.f\t%.1f\t%.7f'% (i,t_k2[i],y_k2[i]))

    plt.plot(t_k2,y_k2,"--o",label="Metode Rung-Kutta Orde 2 Titik Tengah dengan h={}".format(h))
    
def anal(y_0,a,b):
    y = sym.Function('y')
    PDB = sym.Eq(sym.Derivative(y(t),t),f(t,y(t)))
    sol = sym.dsolve(PDB,y(t),ics={y(a):y_0})
    a_x = sol.rhs
    t_a = np.linspace(a,b,100)
    y_a = np.array([])
    for i in t_a:
        y_a = np.append(y_a,a_x.subs(t,i))
        
    plt.plot(t_a,y_a,label="Solusi Analitik",linewidth=3)
        
print(euler(1,0,3,0.5))
print(heun(1,0,3,0.5))
print(taylor(1,0,3,0.5,4))
print(kutta2mid(1,0,3,0.5))
print(anal(1,0,3))

plt.title("Grafik fungsi y' = {} \n Pada [{},{}] dengan y({}) = {}".format(f(t,y),0,3,0,1))
plt.grid()
plt.legend()
plt.show()