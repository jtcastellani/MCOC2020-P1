import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pylab as plt


SaH=1/(60*60) #Segundos a Hora

Vang = (7.27*10**-5)*SaH #rad/h
Mtierra = 5.972 * 10**24 #kg
Rtierra = 6371*1000 #m
G = 6.674*10**-11 #N*m^2/kg^2
Atmosfera = Rtierra + 80*1000 #metros


def satelite (z, t):
    zp= sp.zeros(6)
    teta = Vang*t
    
    #Se hacen los cambios de coordenadas 
    R = sp.array([[sp.cos(teta), -sp.sin(teta), 0],
                  [sp.sin(teta), sp.cos(teta), 0],
                  [0, 0, 1]])
    
    R_p = Vang * sp.array([[-sp.sin(teta), -sp.cos(teta), 0],
                          [sp.cos(teta), -sp.sin(teta), 0],
                          [0, 0, 0]])
                                    
    R_p_p = (Vang ** 2) * sp.array([[-sp.cos(teta), sp.sin(teta), 0],
                                    [-sp.sin(teta), -sp.cos(teta), 0],
                                    [0, 0, 0]])
    


    zp[0] = z[3]
    zp[1] = z[4]
    zp[2] = z[5]
      
    zp[3:6] = (-G*Mtierra)*z[0:3]/Rtierra**3 - R@(R_p_p@z[0:3] + 2*R_p@z[3:6])
    
    return zp

# Posicion inicial
x = Rtierra + 700*1000
y = 0
z = 0


# Velocidades iniciales m/s
vz = 0       
vx = 0
vy = 8020    

#Vector estado
z0 =  sp.array([x,y,z,vx,vy,vz])


t = sp.linspace(0,10000,1001)



sol = odeint(satelite,z0,t)  

x=sol[:,0]
y=sol[:,1]
z=sol[:,2]



plt.figure(1)

plt.plot(t, x, label="X(t)")
plt.plot(t, y, label="Y(t)")
plt.plot(t, z, label="Z(t)")

plt.title("Distancia satelite para dos orbitas")


plt.ylabel("Posicion (miles de Km)")
plt.xlabel("tiempo (s)")
plt.legend()  
plt.savefig("Grafico1.png", dpi=300)


plt.show()



circunferencia = np.linspace(0,2*3.14,360)

Xt = Rtierra*np.cos(circunferencia)
Yt = Rtierra*np.sin(circunferencia)

Xat = Atmosfera*np.cos(circunferencia)
Yat = Atmosfera*np.sin(circunferencia)


plt.figure(2)

d = sp.sqrt(x**2+y**2+z**2)
plt.plot(t, d, label = 'Satelite', c="y")
plt.xlabel('Tiempo (s)')
plt.ylabel('r(t) (miles de Km)')
plt.axhline(y=Atmosfera, c="c", label = "Atmosfera")
plt.axhline(y=Rtierra, c="green", label ="Tierra")

plt.grid()
plt.legend(loc = 1)
plt.tight_layout()
plt.title("Distancia al centro de la tierra")
plt.savefig("Grafico2.png", dpi=300)







plt.figure(3)

plt.plot(Xt, Yt, c="green", label ="Tierra")
plt.plot(Xat, Yat, c="c", label = "Atmosfera")

plt.title("Trayectoria satelite")

plt.plot(x,y, label = "Trayectoria a una velocidad de 8.020 m/s", c='y')
plt.ylabel("Y (miles de Km)")
plt.xlabel("X (miles de Km)")
plt.axis("Equal")
plt.legend()  
plt.savefig("Grafico3.png", dpi=300)


plt.show()
