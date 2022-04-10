import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return (x-y)/2

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(0,3,0.01)
    y_eksak = lambda t: t - 2 + (3/(np.e**(0.5*t)))
    y1 = y_eksak(t1)

def rkf45(f,a,b,h,y_init,warna):
    n = int( (b-a)/h + 1 )
    print("Jumlah interval =",n,"dengan h =",h)
    
    xs = a + np.arange(n)*h
    ys = np.zeros(n)
    zs = np.zeros(n)
    y = y_init
    z = y_init
    
    print("--------------------------------------------------------------------------------------------------------------------------------")
    print(" iterasi\t|  x\t  |k1\t    |k2\t        |k3\t        |k4\t         |k5\t        |k6\t        |y\t      |z\t        |y eksak\t  |  ERM  |")
    print("--------------------------------------------------------------------------------------------------------------------------------")
    
    for j,x in enumerate(xs):
#        y = float('%.5f'%(y))
        ys[j] = y
        zs[j] = z
#        print(x,y)
        k1 = h*f(x,y)
        k2 = h*f(x+(1/4)*h, y+(1/4)*k1)
        k3 = h*f(x+(3/8)*h, y+(3/32)*k1+(9/32)*k2)
        k4 = h*f(x+(12/13)*h, y+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3)
        k5 = h*f(x+h, y+(439/216)*k1-8*k2+(3680/513)*k3-(845/4104)*k4)
        k6 = h*f(x+(1/2)*h, y-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-(11/40)*k5)
        
#        error_relatif = abs((float('%.5f'%(ys[j])) - float('%.5f'%(y_eksak(x))))/float('%.5f'%(y_eksak(x))) )
        error_relatif = abs((zs[j]-y_eksak(x))/y_eksak(x))
        
        '''print hasil'''
        print("  %.0f\t|"%(j)," %.3f\t|"%(xs[j]),"%.5f\t|"%(k1)," %.5f\t|"%(k2)," %.5f\t|"%(k3),
              " %.5f\t|"%(k4)," %.5f\t|"%(k5)," %.5f\t|"%(k6),"%.5f\t|"%(ys[j]),
              "%.5f\t\t|"%(zs[j]),"%.5f\t\t|"%(y_eksak(x)),
              "%.5f\t"%(error_relatif))
        
        z = y + (16/135)*k1+(6656/12825)*k3+(28561/56430)*k4-(9/50)*k5+(2/55)*k6
        y = y + (25/216)*k1+(1408/2565)*k3+(2197/4101)*k4-(1/5)*k5

    print("--------------------------------------------------------------------------------------------------------------------------------")
    
    plt.plot(xs,zs, warna, label = 'Metode RKF45 h={}'.format(h))
    
    return xs[-1],zs[-1]


eksak()
rkf45(f,0,3,0.25,1,'bo--')

plt.plot(t1,y1,'r',label = 'Penyelesaian Eksak')
plt.title('Perbandingan Metode RKF45\ndan Penyelesaian Eksaknya')
plt.xlabel('t')
plt.ylabel('Nilai y')
plt.axis([0,3,0,2])
plt.legend(loc = 'best')
plt.grid()
plt.show()



#'''Euler definisi mudah'''
#def euler(y_0,a,b,h):
#    M = int((b-a)/h)
#    t = a
#    y = y_0
#    t_e = np.array([t])
#    y_e = np.array([y])
#    while t < b:
#        y += h*f(t,y)
#        t += h
#        t_e = np.append(t_e,t)
#        y_e = np.append(y_e,y)
#        
#    print('\n-----------SOLUSI-----------')
#    print('-------Metode Euler---------')    
#    print('k\tt_k\ty_k')
#    print('----------------------------')
#    for i in range(M+1):
#        print('%.f\t%.1f\t%.7f'% (i,t_e[i],y_e[i]))
#    
#    plt.plot(t_e,y_e,"--o",label="Metode Euler h={}".format(h))
#euler(1,0,3,0.5)



#b=3
#a=0
#h=0.6
#n = int((b-a)/h + 1)
#print(n)