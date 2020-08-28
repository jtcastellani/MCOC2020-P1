import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt

hr = 3600. #s
km = 1000. #m

Omg = 7.2921*10**-5 #rad/s
Mtierra = 5.972 * 10**24 #kg
Rtierra = 6371*km #m
G = 6.674*10**-11 #N*m^2/kg^2
Atmosfera = Rtierra + 80.*km #metros
hi = 700.* km # altura inicial

FgMax = G * Mtierra / Rtierra**2 

def satelite (z, t):
    zp= np.zeros(6)
    teta = Omg*t
    
    #Se hacen los cambios de coordenadas 
    R = np.array([[np.cos(teta), -np.sin(teta), 0],
                  [np.sin(teta), np.cos(teta), 0],
                  [0, 0, 1]])
    
    R_p = Omg * np.array([[-np.sin(teta), -np.cos(teta), 0],
                          [np.cos(teta), -np.sin(teta), 0],
                          [0, 0, 0]])
                                    
    R_p_p = (Omg ** 2) * np.array([[-np.cos(teta), np.sin(teta), 0],
                                   [-np.sin(teta), -np.cos(teta), 0],
                                   [0, 0, 0]])
    
    z1 = z[0:3]
    z2 = z[3:6]
    
    r2 = np.dot(z1, z1)
    r = np.sqrt(r2)
    
    Fg = (-G * Mtierra/r**2) * (R@(z1/r))
    
    zp[0:3] = z2
    zp[3:6] = R.T@(Fg - (2 * (R_p@z2) + (R_p_p@z1)))
   
    return zp

import datetime as dt

utc_EOF_format =            "%Y-%m-%dT%H:%M:%S.%f"
t1 = dt.datetime.strptime("2020-08-03T22:59:42.000000", utc_EOF_format)
t2 = dt.datetime.strptime("2020-08-05T00:59:42.000000", utc_EOF_format)

intervalo = t2 - t1
interv_en_s = intervalo.total_seconds()

#print(f"intervalo = {intervalo} s")
#print(f"Intervalo en Segundos = {interv_en_s} s")



# Estado inicial
x_i = 2054785.733796 #m
y_i = 6464390.044070 #m
z_i = 2008332.506688 #m

vx_i = 891.199859 #m/s
vy_i = 2505.711040 #m/s
vz_i = 7117.225886 #m/s



# Estado final
x_f = -950273.200762 #m
y_f = 340731.275287 #m
z_f = 6992937.840398 #m

vx_f = 1956.207151 #m/s
vy_f = 7325.242260 #m/s
vz_f = -90.860827 #m/s


#Vector tiempo
t = np.linspace(0, interv_en_s, 9360)


#Vector estado
z0 =  np.array([x_i,y_i,z_i,vx_i,vy_i,vz_i])


sol = odeint(satelite,z0,t) 

x = sol[:, :] 

pos_f = np.array([x_f,y_f,z_f,vx_f,vy_f,vz_f]) - sol[-1]


deriva = np.sqrt((x_f - sol[-1,0])**2 + (y_f - sol[-1, 1])**2 + (z_f - sol[-1, 2])**2) 
print(deriva)


   
Alt =np.sqrt(x[:, 0]**2 + x[:, 1]**2 + x[:, 2]**2) -Rtierra


def graficos ():
    plt.figure()

    for i in range(3):
        plt.subplot(3,3,1+i)
        plt.grid(True)
        plt.plot(t/hr, x[:,i])
    
    plt.figure()
    plt.grid(True)
    plt.plot(t/hr, Alt/km)
    plt.axhline(80.,  c="c")
    plt.axhline(0., c="green")

    plt.figure()
    plt.grid(True)
    plt.plot(x[:,0], x[:,1])

    th = np.linspace(0, 2 * np.pi, 500)
    plt.plot(Rtierra *np.cos(th), Rtierra * np.sin(th))
    plt.axis("equal")


    fig = plt.figure()
    ax= plt.axes(projection ='3d')
    ax.plot3D(x[:, 0], x[:, 1], x[:, 2])
    plt.show()
