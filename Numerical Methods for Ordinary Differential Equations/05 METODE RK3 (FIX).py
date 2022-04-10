import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return (x-y)/2

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(0,3,0.01)
    y_eksak = lambda t: t - 2 + (3/(np.e**(0.5*t)))
    y1 = y_eksak(t1)

def rk3(f,a,b,h,y_init,warna):
    n = int( (b-a)/h + 1 )
    print("Jumlah interval =",n,"dengan h =",h)
    
    xs = a + np.arange(n)*h
    ys = np.zeros(n)
    y = y_init
    
    print("----------------------------------------------------------------------------------------------------------------")
    print("| iterasi\t|  x\t  |k1\t\t|k2\t\t|k3\t\t|y\t\t|y eksak\t\t|Eror Relatif M\t\t|")
    print("----------------------------------------------------------------------------------------------------------------")
    
    for j,x in enumerate(xs):
#        y = float('%.5f'%(y))
        ys[j] = y
        
        k1 = f(x,y)
        k2 = f(x+(1/2)*h, y+(1/2)*k1*h)
        k3 = f(x+h, y-k1*h+2*k2*h)
        
#        error_relatif = abs((float('%.5f'%(ys[j])) - float('%.5f'%(y_eksak(x))))/float('%.5f'%(y_eksak(x))) )
        error_relatif = abs((ys[j]-y_eksak(x))/y_eksak(x))
        
        '''print hasil'''
        print("|  %.0f\t|"%(j)," %.4f\t|"%(xs[j]),"%.5f\t\t|"%(k1)," %.5f\t\t|"%(k2)," %.5f\t\t|"%(k3),
              "%.5f\t\t|"%(ys[j]),"%.5f\t\t|"%(y_eksak(x)),
              "%.5f\t\t|"%(error_relatif))
        
        y += (1/6)*(k1+4*k2+k3)*h

    print("----------------------------------------------------------------------------------------------------------------")
    
    plt.plot(xs,ys, warna, label = 'Metode Runge-Kutta Orde 3 h={}'.format(h))
    
    return xs[-1],ys[-1]


eksak()
rk3(f,0,3,0.25,1,'bo--')

plt.plot(t1,y1,'r',label = 'Penyelesaian Eksak')
plt.title('Perbandingan Metode Runge-Kutta Orde 3\ndan Penyelesaian Eksaknya')
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