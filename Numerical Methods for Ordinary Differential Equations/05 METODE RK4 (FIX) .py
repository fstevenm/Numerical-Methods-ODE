import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return (x-y)/2

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(0,10,0.01)
    y_eksak = lambda t:  t - 2 + (3/(np.e**(0.5*t)))
    y1 = y_eksak(t1)

def rk3(f,a,b,h,y_init,warna):
    n = int( (b-a)/h + 1 )
    print("Jumlah interval =",n,"dengan h =",h)
    
    xs = a + np.arange(n)*h
    ys = np.zeros(n)
    y = y_init
    
    print("---------------------------------------------------------------")
    print("| iterasi\t|  x\t |y\t\t|y eksak\t\t|Galat Relatif\t\t|")
    print("---------------------------------------------------------------")
    
    for j,x in enumerate(xs):
        ys[j] = y
        
        k1 = f(x,y)
        k2 = f(x+(1/2)*h, y+(1/2)*k1*h)
        k3 = f(x+(1/2)*h, y+(1/2)*k2*h)
        k4 = f(x+h, y+k3*h)
        
        error_relatif = ((ys[j]-y_eksak(x))/y_eksak(x))
        
        '''print hasil'''
        print("|  %.0f\t|"%(j)," %.4f\t|"%(xs[j]),"%.6f\t\t|"%(ys[j]),"%.6f\t\t|"%(y_eksak(x)),
              "%.6f\t\t|"%(error_relatif))
        
        y += (1/6)*(k1+2*k2+2*k3+k4)*h

    print("------------------------------------------------------------------------------")
    
    plt.plot(xs,ys, warna, label = 'Metode Runge-Kutta Orde 4 h={}'.format(h))
    
    return xs[-1],ys[-1]


eksak()
rk3(f,0,10,0.125,1,'bo--')
#rk3(f,0,10,0.5,1,'go--')
#rk3(f,0,10,0.1,1,'yo--')

plt.plot(t1,y1,'r',label = 'Penyelesaian eksak $y = e^{-x}$')
plt.title('Perbandingan Metode Runge-Kutta Orde 4\ndan Penyelesaian Eksaknya')
plt.xlabel('x')
plt.ylabel('Nilai y')
plt.axis([0,4,0,2])
plt.legend(loc = 'best')
plt.grid()
plt.show()

