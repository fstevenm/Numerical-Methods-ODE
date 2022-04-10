import numpy as np
import matplotlib.pyplot as plt

def dxdt(t,x,y):
    return y

def dydt(t,x,y):
    return -4*y-5*x

def eksak():
    global t1,y1,y_eksak
    t1 = np.arange(0,5,0.1)
    y_eksak = lambda t: (3*(np.e**(-2*t))*(np.cos(t)))+((np.e**(-2*t))*(np.sin(t)))
    y1 = y_eksak(t1)

def rk4(dxdt,dydt,a,b,h,x_init,y_init):
    
    n = int((b-a)/h) +1
    ts = a + np.arange(n)*h
    xs = np.zeros(n)
    ys = np.zeros(n)
       
    x = x_init
    y = y_init
    
    for j,t in enumerate(ts):
        xs[j] = x
        ys[j] = y

        
        k1x = dxdt(t,x,y)
        k1y = dydt(t,x,y)

        k2x = dxdt(t+(1/2)*h, x+(1/2)*k1x*h, y+(1/2)*k1y*h)
        k2y = dydt(t+(1/2)*h, x+(1/2)*k1x*h, y+(1/2)*k1y*h)

        k3x = dxdt(t+(1/2)*h, x+(1/2)*k2x*h, y+(1/2)*k2y*h)
        k3y = dydt(t+(1/2)*h, x+(1/2)*k2x*h, y+(1/2)*k2y*h)

        k4x = dxdt(t+h, x+k3x*h, y+k3y*h)
        k4y = dydt(t+h, x+k3x*h, y+k3y*h)
        
        x += (1/6)*(k1x+2*k2x+2*k3x+k4x)*h
        y += (1/6)*(k1y+2*k2y+2*k3y+k4y)*h


        '''print hasil'''
   
    print('\n---------------------SOLUSI----------------------------')
    print('-------------------------------------------------------')
    print('n\tt_k\t  x_k\t    eksak\t\terror relatif')
    print('-------------------------------------------------------')
    for j in range(n):
        error_relatif = ((y_eksak(ts[j])-xs[j])/y_eksak(ts[j]))
        print("|  %.0f\t|"%(j),"%.4f\t|"%(ts[j])," %.4f\t|"%(xs[j]),"%.5f\t\t|"%(y_eksak(ts[j])),
                "%.5f\t\t|"%(error_relatif))

    
    plt.plot(ts,xs, 'b')
    
    
    return ts[-1],xs[-1],ys[-1]

eksak()
rk4(dxdt,dydt,0,5,0.1,3,-5) #jurnal
# rk6(dsdt,dedt,dqdt,di1dt,di2dt,drdt,0,1000,1,0.6,0.2,0.2,0,0,0) #jika belum berkembang

plt.plot(t1,y1,'r',label = 'Penyelesaian Eksak')
plt.legend(["x(t)",eksak],loc = 'best')
plt.title('Grafik Penyelesaian Numeris Model Matematika Penyebaran COVID-19')
plt.xlabel('t')
plt.ylabel('kelompok(t)')
plt.axis([0,5,-0.1,3])
plt.grid()
plt.show()

