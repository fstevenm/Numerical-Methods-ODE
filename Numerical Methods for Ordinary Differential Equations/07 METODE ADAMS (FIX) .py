import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return (x-y)/2

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(0,3,0.01)
    y_eksak = lambda t: (3*np.e**(-t/2))-2+t
    y1 = y_eksak(t1)

def ABM(f,a,b,h,y_init,warna):
    global xs, y_adams
    n = int( (b-a)/h)
    print("Jumlah interval =",n,"dengan h =",h)
    
    xs = a + np.arange(n+10)*h
    ys = np.zeros(n+10)
    ps = np.zeros(n+10)
    y = y_init

    y_adams = []
    p_adams = []
#     print(xs)
    xs_2 = xs[1:]
    xs_3 = xs[2:]
    xs_4 = xs[3:]
    xs_5 = xs[4:]
#     print(xs_2)
    print("----------------------------------------------------------------------------------------------------------------")
    print("|  iterasi  |     t     |      p      |    y adams    |  y eksak  |  Eror Relatif  |")
    print("----------------------------------------------------------------------------------------------------------------")
    
    for i,x in enumerate (xs):
        ys[i] = y
        
        k1 = f(x,y)
        k2 = f(x+(1/2)*h, y+(1/2)*k1*h)
        k3 = f(x+(1/2)*h, y+(1/2)*k2*h)
        k4 = f(x+h, y+k3*h)
        
        y += (1/6)*(k1+2*k2+2*k3+k4)*h
    

    for i in range (0,4):
        y_adams.append(ys[i])
        p_adams.append(0)
#     print(p_adams)
#     print(y_adams)
#     print(xs)
#     print(xs_5)

    for j in range (0,n+1):
        f0 = f(xs[j],y_adams[j])
        f1 = f(xs_2[j],y_adams[j+1])
        f2 = f(xs_3[j],y_adams[j+2])
        f3 = f(xs_4[j],y_adams[j+3])
        error_relatif = abs((y_adams[j]-y_eksak(xs[j]))/y_eksak(xs[j]))

        p_adams.append( y_adams[j+3] + (h/24)*(-9*f0+37*f1-59*f2+55*f3))
        y_adams.append( y_adams[j+3] + (h/24)*(f1-5*f2+19*f3+9*f(xs_5[j],p_adams[-1])))
        
        print("|    %.0f\t    |"%(j)," %.4f   |" %(xs[j])," %5.5f    |"  %(p_adams[j]),
              "   %.5f     |"%(y_adams[j])," %.4f |"%(y_eksak(xs[j])),"    %.4f    |"%(error_relatif))
    
    judul = 'h = '+str(h)
#    print(len(xs),len(y_adams))
    plot_gambar(warna,judul)
    
    return xs[-1],ys[-1]

def plot_gambar(warna,judul):
    '''------------ Plot grafik ------------------'''
    a=0
    b=3
    h = 1/8
    n = int( (b-a)/h)
    xs = a + np.arange(n+5)*h
    
    plt.plot(xs,y_adams, warna, label = 'Metode Adams Bashfroth Moulton dengan '+judul) 
    
eksak()
ABM(f,0,3,(1/8),1,'ro--')

plt.plot(t1,y1,'blue', label = 'Solusi Eksak',linewidth = 2)
plt.title('Perbandingan Hasil Metode Adams Bashfroth Moulton \ndan Penyelesaian y`=(t-y)/2')
plt.xlabel('t')
plt.ylabel('Nilai y')
plt.axis([-0,3,0,2.1])
plt.legend()

plt.grid()
plt.show()