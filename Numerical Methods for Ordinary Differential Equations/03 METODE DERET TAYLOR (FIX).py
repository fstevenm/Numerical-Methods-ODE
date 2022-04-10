import numpy as np
import matplotlib.pyplot as plt

#def f(x,y):
#    return (x-y)/2
#
#def eksak():
#    global t1,y1,y_eksak
#    t1 = np.arange(0,3,0.01)
#    y_eksak = lambda t: t - 2 + (3/(np.e**(0.5*t)))
#    y1 = y_eksak(t1)

''' Metode Taylor dengan menggunakan definisi'''
#def d1(x,y):
#    return -x*y
#
#def d2(x,y):
#    return -y+(x**2)*y
#
#def d3(x,y):
#    return 3*x*y-(x**3)*y
#
#def eksak():
#    global t1,y1,y_eksak
#    t1 = np.arange(-3,3,0.01)
#    y_eksak = lambda t: 1/(np.e**(0.5*t**2))
#    y1 = y_eksak(t1)
def d1(x,y):
    return (x-y)/2

def d2(x,y):
    return (2-x+y)/4

def d3(x,y):
    return (-2+x-y)/8

def d4(x,y):
    return (2-x+y)/16

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(-3,3,0.01)
    y_eksak = lambda t: t - 2 + (3/(np.e**(0.5*t)))
    y1 = y_eksak(t1)

def taylor(d1,d2,d3,a,b,h,y_init,warna):
    n = int( (b-a)/h + 1 )
    print("Jumlah interval =",n,"dengan h =",h)
    
    xs = a + np.arange(n)*h
    ys = np.zeros(n)
    y = y_init
    
    print("--------------------------------------------------------------------------------------------------------")
    print("| iterasi\t |    x\t   |    d1\t   |    d2\t   |    d3\t   |    d4\t   |     y\t\t   |   y eksak\t  | error relatif mutlak |")
    print("--------------------------------------------------------------------------------------------------------")
    
    for j,x in enumerate(xs):
#        y = float('%.5f'%(y))
        ys[j] = y
        
#        error_relatif = abs((float('%.5f'%(ys[j])) - float('%.5f'%(y_eksak(x))))/float('%.5f'%(y_eksak(x))) )
        error_relatif = abs((ys[j]-y_eksak(x))/y_eksak(x))
        
        '''print hasil'''
        print("|    %.0f\t    |"%(j)," %.4f\t|"%(xs[j])," %.4f\t|"%(d1(x,y)),
              " %.4f\t|"%(d2(x,y))," %.4f\t|"%(d3(x,y))," %.4f\t|"%(d4(x,y)),
              "  %.5f\t\t|"%(ys[j]),"  %.5f\t  |"%(y_eksak(x)),
              "      %.5f\t\t     |"%(error_relatif))
        
        y += d1(x,y)*h + (d2(x,y)*(h**2)/2) + (d3(x,y)*(h**3)/6) + (d4(x,y)*(h**4)/24)

    print("--------------------------------------------------------------------------------------------------------")
    
    plt.plot(xs,ys, warna, label = 'Metode Deret Taylor h={}'.format(h))
    
    return xs[-1],ys[-1]


eksak()
taylor(d1,d2,d3,0,3,0.25,1,'bo--')

plt.plot(t1,y1,'r',label = 'Penyelesaian Eksak')
plt.title('Perbandingan Metode Deret Taylor\ndan Penyelesaian Eksaknya')
plt.xlabel('t')
plt.ylabel('Nilai y')
plt.axis([0,3,0,1.1])
plt.legend(loc = 'best')
plt.grid()
plt.show()



#b=3
#a=0
#h=0.6
#n = int((b-a)/h + 1)
#print(n)