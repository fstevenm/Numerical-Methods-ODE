import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return (x-y)/2

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(0,3,0.01)
    y_eksak = lambda t: t - 2 + (3/(np.e**(0.5*t)))
    y1 = y_eksak(t1)


def rk2_titik_tengah(f,a,b,h,y_init,warna):
    n = int( (b-a)/h + 1 )
    print("Jumlah interval =",n,"dengan h =",h)
    
    xs = a + np.arange(n)*h
    ys = np.zeros(n)
    y = y_init
    
    print("--------------------------------------------------------------------------------------------")
    print("| iterasi\t |    x\t   |    k1\t   |    k2\t   |     y\t\t   |   y eksak\t  | error relatif mutlak |")
    print("--------------------------------------------------------------------------------------------")
    
    for j,x in enumerate(xs):
#        y = float('%.5f'%(y))
        ys[j] = y
        
        k1 = f(x,y)
        k2 = f(x+(1/2)*h, y+(1/2)*k1*h)
        
#        error_relatif = abs((float('%.5f'%(ys[j])) - float('%.5f'%(y_eksak(x))))/float('%.5f'%(y_eksak(x))) )
        error_relatif = abs((ys[j]-y_eksak(x))/y_eksak(x))
        
        '''print hasil'''
        print("|    %.0f\t    |"%(j)," %.4f\t|"%(xs[j])," %.5f\t|"%(k1)," %.5f\t|"%(k2),
              "  %.5f\t\t|"%(ys[j]),"  %.5f\t  |"%(y_eksak(x)),
              "      %.5f\t\t     |"%(error_relatif))
        
        y += k2*h

    print("--------------------------------------------------------------------------------------------")
    
    plt.plot(xs,ys, warna, label = 'Metode Titik Tengah - Runge-Kutta Orde 2 h={}'.format(h))
    
    return xs[-1],ys[-1]


eksak()
rk2_titik_tengah(f,0,3,0.25,1,'bo--')

plt.plot(t1,y1,'r',label = 'Penyelesaian Eksak')
plt.title('Perbandingan Metode Titik Tengah - Runge-Kutta Orde 2\ndan Penyelesaian Eksaknya')
plt.xlabel('t')
plt.ylabel('Nilai y')
plt.axis([0,3,0,2])
plt.legend(loc = 'best')
plt.grid()
plt.show()



'''Titik tengah definisi mudah'''
#def kutta2mid(y_0,a,b,h):
#    M = int((b-a)/h)
#    t = a
#    y = y_0
#    t_k2 = np.array([t])
#    y_k2 = np.array([y])
#    while t < b:
#        k1 = f(t,y)
#        k2 = f((t+(1/2)*h),(y+((1/2)*k1*h)))
#        y += k2*h
#        y_k2 = np.append(y_k2,y)
#        t += h
#        t_k2 = np.append(t_k2,t)
#        
#    print('\n-------------SOLUSI-------------')
#    print('------Metode Rung-Kutta Orde 2-------\n --------------Titik Tengah--------------')
#    print('k\tt_k\ty_k')
#    print('--------------------------------')
#    for i in range(len(t_k2)):
#        print('%.f\t%.1f\t%.7f'% (i,t_k2[i],y_k2[i]))
#
#    plt.plot(t_k2,y_k2,"--o",label="Metode Rung-Kutta Orde 2 Titik Tengah dengan h={}".format(h))
#kutta2mid(1,0,3,0.5)



#b=3
#a=0
#h=0.6
#n = int((b-a)/h + 1)
#print(n)