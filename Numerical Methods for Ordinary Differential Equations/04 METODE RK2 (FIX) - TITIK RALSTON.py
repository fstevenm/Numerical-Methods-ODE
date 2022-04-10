import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return 3*y-4*(np.e**(-x))

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(0,10,0.01)
    y_eksak = lambda t: (np.e**(-t))
    y1 = y_eksak(t1)


def rk2_ralston(f,a,b,h,y_init,warna):
    n = int( (b-a)/h + 1 )
    print("Jumlah interval =",n,"dengan h =",h)
    
    xs = a + np.arange(n)*h
    ys = np.zeros(n)
    y = y_init
    
    print("------------------------------------------------------------------------------")
    print("| iterasi\t |    x\t   |     y hampiran\t\t\t     |   y eksak\t  | galat relatif |")
    print("-----------------------------------------------------------------------------")
    
    for j,x in enumerate(xs):
        y = float('%.5f'%(y))
        ys[j] = y
        
        if j == 0:
            k1 = 0
            k2 = 0
            
        error_relatif = ((ys[j]-y_eksak(x))/y_eksak(x))
        
        '''print hasil'''
        print("|    %.0f\t  |"%(j)," %.4f\t|"%(xs[j]),
              "  %.6f\t\t\t\t|"%(ys[j]),"  %.6f\t  |"%(y_eksak(x)),
              "      %.6f\t\t     |"%(error_relatif))
        k1 = f(x,y)
        k2 = f(x+(3/4)*h, y+(3/4)*k1*h)
        y += ((1/3)*k1+(2/3)*k2)*h

    print("------------------------------------------------------------------------------------------------------")
    
    plt.plot(xs,ys, warna, label = 'Metode Ralston - Runge-Kutta Orde 2 h={}'.format(h))
    
    return xs[-1],ys[-1]


eksak()
rk2_ralston(f,0,10,1,1,'bo--')
rk2_ralston(f,0,10,0.5,1,'go--')
rk2_ralston(f,0,10,0.1,1,'yo--')

plt.plot(t1,y1,'r',label = 'Penyelesaian eksak $y = e^{-x}$')
plt.title('Perbandingan Metode Ralston - Runge-Kutta Orde 2\ndan Penyelesaian Eksaknya')
plt.xlabel('x')
plt.ylabel('Nilai y')
plt.axis([0,10,-5,2])
plt.legend(loc = 'best')
plt.grid()
plt.show()

