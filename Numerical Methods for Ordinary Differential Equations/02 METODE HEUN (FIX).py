import numpy as np
import matplotlib.pyplot as plt

#def f(x,y):
#    return -x*y
#
#def eksak():
#    global t1,y1,y_eksak
#    t1 = np.arange(0,3,0.01)
#    y_eksak = lambda t: 1/(np.e**(0.5*t))
#    y1 = y_eksak(t1)

def f(x,y):
    return (x-y)/2

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(0,3,0.01)
    y_eksak = lambda t: t - 2 + (3/(np.e**(0.5*t)))
    y1 = y_eksak(t1)


def heun(f,a,b,h,y_init,warna):
    n = int( (b-a)/h + 1 )
    print("Jumlah interval =",n,"dengan h =",h)
    
    xs = a + np.arange(n+1)*h
    ps = np.zeros(n+1)
    ys = np.zeros(n+1)
    y = y_init
    p = 0

    print("----------------------------------------------------------------------------------")
    print("| iterasi\t |    x\t   |    p\t   |     y\t\t   |   y eksak\t  | error relatif mutlak |")
    print("----------------------------------------------------------------------------------")
    
    for j,x in enumerate(xs):
        if j < n:
#        y = (float('%.5f'%(y)))
            ys[j] = y
            ps[j] = p
            
            p = y + h*f(x,y)
            y += (h/2)*(f(x,y)+f(xs[j+1],p))
            
#        error_relatif = abs((float('%.5f'%(ys[j])) - float('%.5f'%(y_eksak(x))))/float('%.5f'%(y_eksak(x))) )
            error_relatif = abs((ys[j]-y_eksak(x))/y_eksak(x))
            
            '''print hasil'''
            print("|    %.0f\t    |"%(j)," %.4f\t|"%(xs[j])," %.5f\t|"%(ps[j]),
                  "  %.5f\t\t|"%(ys[j]),"  %.5f\t  |"%(y_eksak(x)),
                  "      %.5f\t\t     |"%(error_relatif))
#        else:
#            ys[j] = y
#            ps[j] = p
        
    print("----------------------------------------------------------------------------------")
    
    plt.plot(xs,ys, warna, label = 'Metode Heun h={}'.format(h))
    
    return xs[-1],ys[-1]


eksak()
heun(f,0,3,0.25,1,'bo--')

plt.plot(t1,y1,'r',label = 'Penyelesaian Eksak')
plt.title('Perbandingan Metode Heun\ndan Penyelesaian Eksaknya')
plt.xlabel('t')
plt.ylabel('Nilai y')
plt.axis([0,3,0,2])
plt.legend(loc = 'best')
plt.grid()
plt.show()



'''Heun definisi mudah'''
#def heun(y_0,a,b,h):
#    M = int((b-a)/h)
#    t = a
#    y = y_0
#    t_h = np.array([t])
#    y_h = np.array([y])
#    while t < b:
#        t += h
#        y += (h/2)*(f(t-h,y)+f(t,(y+h*f(t-h,y))))
#        t_h = np.append(t_h,t)
#        y_h= np.append(y_h,y)
#
#    print('\n-----------SOLUSI-----------')
#    print('--------Metode Heun---------')    
#    print('k\tt_k\ty_k')
#    print('----------------------------')
#    for i in range(M+1):
#        print('%.f\t%.1f\t%.7f'% (i,t_h[i],y_h[i]))
#        
#    plt.plot(t_h,y_h,"--o",label="Metode Heun h={}".format(h))
#heun(1,0,3,0.5)




#b=3
#a=0
#h=0.6
#n = int((b-a)/h + 1)
#print(n)