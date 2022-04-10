import numpy as np
import matplotlib.pyplot as plt
print ("")
print ("=================h=0.5================")
print ("")

#batas dan besar langkah dan yang diketahui
h1=0.5
y0=1
a=0
b=3

#untuk menghitung nilai eksak
t=np.arange(0,3,0.01)
y0=1
y=[]
f = lambda t: (3*np.e**(-t/2))-2+t
y=f(t)

#metode euler (nilai t, prediktor, korektor, eksak)
pk1=[]
tk1=[]
yk1=[]
yeksak=[]

t1=np.arange(a,b,h1)
print (" t      prediktor   korektor    eksak")
print ("----------------------------------------")
for i in range (len(t1)+1):
    if i==0:
        t0=a+(i*h1)
        tk1.append(t0)
        pk1.append(0)
        yk1.append(y0)
    
    else:
        tk1.append(tk1[i-1]+h1)
        pk1.append(yk1[i-1]+((h1/2)*(tk1[i-1]-yk1[i-1])))
        yk1.append(yk1[i-1]+((h1/2)*(((1/2)*(tk1[i-1]-yk1[i-1]))+((1/2)*(tk1[i]-pk1[i])))))
        
    print ("{:.3f}".format(tk1[i]),"  ", "{:.5f}".format(pk1[i]),"  ","{:.6f}".format(yk1[i]),"  ","{:.6f}".format(f(i-i/2)))

plt.plot(t,y,"black") #grafik y(t)
plt.plot(tk1,yk1,'bo--')

plt.xlabel("t")
plt.ylabel("y")
plt.title("$ y'=(t-y)/2 $")
plt.legend(["Solusi eksak","Heun dengan h=0.5"])

plt.axis([-0.1,3.1,0,2])
plt.grid()
plt.show()

