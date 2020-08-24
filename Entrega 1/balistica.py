import scipy as sp
from scipy.integrate import odeint
import matplotlib.pylab as plt

#Parámetros

#Unidades base
cm = 0.01
inch = 2.54*cm
g = 9.81

#Coeficiente de arrastre
ro = 1.225
cd = 0.47
D = 8.5*inch
r = D/2
A = sp.pi * r**2
CD = 0.5*ro*cd*A

#Masa
m = 15.       #Kg

#Viento
V = 0.       #m/s

# Función a integrar

#z es vector etado, z = [x, y, vx, vy]
#dx1/dt=z2

#      [z1      ]
#dz/dt=[        ]
#      [FD/m - g]

#Vector Estado
# z[0] = x
# z[1] = y
# z[2] = vx
# z[3] = vy




def bala(z, t,):
    zp= sp.zeros(4)
    
    zp[0] = z[2]
    zp[1] = z[3]
    
    v = z[2:4]   #ultimos dos componentes
    v[0] = v[0] - V
    v2 = sp.dot(v,v)

    vnorm = sp.sqrt(v2)
    
    FD = -CD * v2 * (v/vnorm)
    
    zp[2] = FD[0]/m
    zp[3] = FD[1]/m - g
    
    return zp

t = sp.linspace(0, 30, 1001)

#Parte en el origen con vx=vy=vi m/s

vi = 100*1000./3600.

z0 = sp.array([0, 0, vi, vi])


sol_0 = odeint(bala, z0, t)
    
x_0 = sol_0[:,0]    
y_0 = sol_0[:,1]

V = 10

sol_10 = odeint(bala, z0, t)

x_10 = sol_10[:,0]    
y_10 = sol_10[:,1]

V = 20

sol_20 = odeint(bala, z0, t)

x_20 = sol_20[:,0]    
y_20 = sol_20[:,1]
    
plt.figure(1)
plt.ylim(0, 50)
plt.xlim(0, 150)
plt.title("Trayectoria para distintos vientos")
plt.ylabel("Y (m)")
plt.xlabel("X (m)")
plt.grid(axis = 'both')
plt.plot(x_0,y_0, label='V = 0 m/s') 
plt.plot(x_10,y_10, label='V = 10.0 m/s') 
plt.plot(x_20,y_20, label='V = 20.0 m/s') 
plt.legend()
plt.tight_layout()
plt.savefig('Entrega_1.png', dpi=300)
plt.show()   
    
    
    