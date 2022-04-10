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
    t_k2t = np.array([t])
    y_k2t = np.array([y])
    while t < b:
        k1 = f(t,y)
        k2 = f((t+(1/2)*h),(y+((1/2)*k1*h)))
        y += k2*h
        y_k2t = np.append(y_k2t,y)
        t += h
        t_k2t = np.append(t_k2t,t)
        
    print('\n-------------SOLUSI-------------')
    print('------Metode Rung-Kutta Orde 2-------\n --------------Titik Tengah--------------')
    print('k\tt_k\ty_k')
    print('--------------------------------')
    for i in range(len(t_k2t)):
        print('%.f\t%.1f\t%.7f'% (i,t_k2t[i],y_k2t[i]))

    plt.plot(t_k2t,y_k2t,"--o",label="Metode Rung-Kutta Orde 2 Titik Tengah h={}".format(h))
    
def kutta2heun(y_0,a,b,h):
    M = int((b-a)/h)
    t = a
    y = y_0
    t_k2h = np.array([t])
    y_k2h = np.array([y])
    while t < b:
        k1 = f(t,y)
        k2 = f((t+h),(y+(k1*h)))
        y += (1/2)*(k1+k2)*h
        y_k2h = np.append(y_k2h,y)
        t += h
        t_k2h = np.append(t_k2h,t)
        
    print('\n-------------SOLUSI-------------')
    print('------Metode Rung-Kutta Orde 2-------\n --------------Heun--------------')
    print('k\tt_k\ty_k')
    print('--------------------------------')
    for i in range(len(t_k2h)):
        print('%.f\t%.1f\t%.7f'% (i,t_k2h[i],y_k2h[i]))

    plt.plot(t_k2h,y_k2h,"--o",label="Metode Rung-Kutta Orde 2 Heun h={}".format(h))
    
def kutta2ral(y_0,a,b,h):
    M = int((b-a)/h)
    t = a
    y = y_0
    t_k2r = np.array([t])
    y_k2r = np.array([y])
    while t < b:
        k1 = f(t,y)
        k2 = f((t+(3/4)*h),(y+((3/4)*k1*h)))
        y += ((1/3)*k1 + (2/3)*k2)*h
        y_k2r = np.append(y_k2r,y)
        t += h
        t_k2r = np.append(t_k2r,t)
        
    print('\n-------------SOLUSI-------------')
    print('------Metode Rung-Kutta Orde 2-------\n --------------Ralston--------------')
    print('k\tt_k\ty_k')
    print('--------------------------------')
    for i in range(len(t_k2r)):
        print('%.f\t%.1f\t%.7f'% (i,t_k2r[i],y_k2r[i]))

    plt.plot(t_k2r,y_k2r,"--o",label="Metode Rung-Kutta Orde 2 Ralston h={}".format(h))
    
def kutta3(y_0,a,b,h):
    M = int((b-a)/h)
    t = a
    y = y_0
    t_k3 = np.array([t])
    y_k3 = np.array([y])
    while t < b:
        k1 = f(t,y)
        k2 = f((t+(1/2)*h),(y+((1/2)*k1*h)))
        k3 = f((t+h),(y-k1*h+2*k2*h))
        y += (1/6)*(k1 + 4*k2 + k3)*h
        y_k3 = np.append(y_k3,y)
        t += h
        t_k3 = np.append(t_k3,t)
        
    print('\n-------------SOLUSI-------------')
    print('------Metode Rung-Kutta Orde 3-----')
    print('k\tt_k\ty_k')
    print('--------------------------------')
    for i in range(len(t_k3)):
        print('%.f\t%.1f\t%.7f'% (i,t_k3[i],y_k3[i]))

    plt.plot(t_k3,y_k3,"--o",label="Metode Rung-Kutta Orde 3 h={}".format(h))   

def kutta4(y_0,a,b,h):
    M = int((b-a)/h)
    t = a
    y = y_0
    t_k4 = np.array([t])
    y_k4 = np.array([y])
    while t < b:
        k1 = f(t,y)
        k2 = f((t+(1/2)*h),(y+((1/2)*k1*h)))
        k3 = f((t+(1/2)*h),(y+((1/2)*k2*h)))
        k4 = f((t+h),(y+k3*h))
        y += (1/6)*(k1 + 2*k2 + 2*k3 + k4)*h
        y_k4 = np.append(y_k4,y)
        t += h
        t_k4 = np.append(t_k4,t)
        
    print('\n-------------SOLUSI-------------')
    print('------Metode Rung-Kutta Orde 4-----')
    print('k\tt_k\ty_k')
    print('--------------------------------')
    for i in range(len(t_k4)):
        print('%.f\t%.1f\t%.7f'% (i,t_k4[i],y_k4[i]))

    plt.plot(t_k4,y_k4,"--o",label="Metode Rung-Kutta Orde 4 h={}".format(h))
    
def kuttafeh(y_0,a,b,h):
    M = int((b-a)/h)
    t = a
    y = y_0
    t_kfeh = np.array([t])
    y_kfeh = np.array([y])
    z_kfeh = np.array([y])
    while t < b:
        k1 = h*f(t,y)
        k2 = h*f((t+(1/4)*h),(y+((1/4)*k1)))
        k3 = h*f((t+(3/8)*h),(y+((3/32)*k1)+((9/32)*k2)))
        k4 = h*f((t+(12/13)*h),(y+((1932/2197)*k1)-((7200/2197)*k2)+((7296/2197)*k3)))
        k5 = h*f((t+h),(y+((439/216)*k1)-(8*k2)+((3680/513)*k3)-((845/4104)*k4)))
        k6 = h*f((t+(1/2)*h),(y-((8/27)*k1)+(2*k2)-((3544/2565)*k3)+((1859/4104)*k4)-((11/40)*k5)))
        z = y + ((16/135)*k1) + ((6656/12825)*k3) + ((28561/56430)*k4) - ((9/50)*k5) + ((2/55)*k6)
        z_kfeh = np.append(z_kfeh,z)
        y = y + ((25/216)*k1)+((1408/2565)*k3)+((2197/4101)*k4)-((1/5)*k5)
        y_kfeh = np.append(y_kfeh,y)
        t += h
        t_kfeh = np.append(t_kfeh,t)
        
    print('\n-------------SOLUSI-------------')
    print('------Metode Rung-Kutta Fehlberg-----')
    print('k\tt_k\ty_k\tz_k')
    print('--------------------------------')
    for i in range(len(t_kfeh)):
        print('%.f\t%.1f\t%.7f\t%.7f'% (i,t_kfeh[i],y_kfeh[i],z_kfeh[i]))

    plt.plot(t_kfeh,z_kfeh,"--o",label="Metode Rung-Kutta Fehlberg h={}".format(h))
    
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
        
#print(euler(1,0,3,0.5))
#print(heun(1,0,3,0.5))
#print(taylor(1,0,3,0.5,4))
#print(kutta2mid(1,0,3,0.5))
#print(kutta2heun(1,0,3,0.5))
#print(kutta2ral(1,0,3,0.5))
#print(kutta3(1,0,3,0.5))
#print(kutta4(1,0,3,0.5))
print(kuttafeh(1,0,3,0.5))
print(anal(1,0,3))

plt.title("Grafik fungsi penyelesaian y' = {} \n Pada [{},{}] dengan y({}) = {}".format(f(t,y),0,3,0,1))
plt.grid()
plt.legend()
plt.show()