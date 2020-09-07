import numpy as np
from scipy.integrate import odeint
import matplotlib.pylab as plt
from leer_eof import leer_eof
from time import perf_counter as cron

tti = cron()

archivo = "S1B_OPER_AUX_POEORB_OPOD_20200824T111152_V20200803T225942_20200805T005942.EOF"
at, ax, ay, az, avx, avy, avz = leer_eof(archivo)


hr = 3600. #s
km = 1000. #m

Omg = 7.2921*10**-5 #rad/s
Mtierra = 5.972 * 10**24 #kg
Rtierra = 6371*km #m
G = 6.674*10**-11 #N*m^2/kg^2
Atmosfera = Rtierra + 80.*km #metros
hi = 700.* km # altura inicial
FgMax = G * Mtierra / Rtierra**2 

J2 = (1.75553**10)*(km**5) #km5*s-2
J3 = (-2.61913**11)*(km**6) #km6*s-2



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
    
    x = z[0:3]
    xp = z[3:6]
    
    r = np.sqrt(np.dot(x, x))
    
    xstill = R@x
    rnorm = xstill/r
    Fg = -G * Mtierra/r**2 * rnorm
    
    z2 = xstill[2]**2
    rflat = xstill[0]**2 + xstill[1]**2
    FJ2 = J2*xstill/r**7
    FJ2[0] = FJ2[0] * (6*z2 - 1.5*rflat)
    FJ2[1] = FJ2[1] * (6*z2 - 1.5*rflat)
    FJ2[2] = FJ2[2] * (3*z2 - 4.5*rflat)
    
    FJ3 = np.zeros(3)
    FJ3[0] = J3 * xstill[0]*xstill[2]/r**9 * (10*z2 - 7.5*rflat)
    FJ3[1] = J3 * xstill[1]*xstill[2]/r**9 * (10*z2 - 7.5*rflat)
    FJ3[2] = J3                      /r**9 * ( 4*z2 *(z2 - 3*rflat) + 1.5*rflat**2) 
    
    
    zp[0:3] = xp
    
   
   
    #zp[3:6] = R.T@(Fg - (2 * (R_p@xp) + (R_p_p@x)))    se utilizó para las preguntas 1, 2 y 3

    zp[3:6] = R.T@(Fg + FJ2 + FJ3 - (2 * (R_p@xp) + (R_p_p@x)))   # se utiliza para la pregunta 4
   
    return zp

def eulerint(zp, z0, t, Nsubdivisiones):
	Nt = len(t)
	Ndim = len(z0)

	z = np.zeros((Nt, Ndim))
	z[0, :] = z0

	for i in range(1, Nt):
		t_anterior = t[i-1]
		dt = (t[i]- t[i-1])/Nsubdivisiones
		z_temp = z[i-1, :].copy()
		for k in range(Nsubdivisiones):
			z_temp += dt*zp(z_temp, t_anterior + k*dt)
		z[i, :] = z_temp
	return z


t = at
x0 = Rtierra + hi
vt = 6820.*3.6 #m/s

z0 = np.array([ax[0], ay[0], az[0], avx[0], avy[0], avz[0]])


''' Pregunta 1
sol = odeint(satelite, z0, t)

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]


plt.figure()
plt.subplot(3,1,1)
plt.suptitle("Posición")
plt.plot(at/hr, ax/km, label = "Real")
plt.plot(t/hr, x/km, label = "Analítica")
plt.ylabel("$X$ (Km)")
plt.subplot(3,1,2)
plt.plot(at/hr, ay/km, label = "Real")
plt.plot(t/hr, y/km, label = "Analítica")
plt.ylabel("$Y$ (Km)")
plt.subplot(3,1,3)
plt.plot(at/hr, az/km, label = "Real")
plt.plot(t/hr, z/km, label ="Analítica")
plt.ylabel("$Z$ (Km)")
plt.xlabel("Tiempo, $t$ (horas)")
plt.legend()
plt.savefig("GraficoP1.png", dpi=300)
plt.show()
'''


''' Pregunta 2
delta = at[-1]
t = np.linspace(0, delta, len(at))
t1=cron()

sol_odeint= odeint(satelite, z0, t)
t2=cron()

sol_eulerint = eulerint(satelite, z0, t, 1)
t3=cron()

tiempo_odeint = t2-t1
tiempo_euler = t3-t2

x_odeint = sol_odeint[:,0]
y_odeint = sol_odeint[:,1]
z_odeint = sol_odeint[:,2]
vx_odeint = sol_odeint[:,3]
vy_odeint = sol_odeint[:,4]
vz_odeint = sol_odeint[:,5]

x_euler=sol_eulerint[:,0]
y_euler=sol_eulerint[:,1]
z_euler=sol_eulerint[:,2]


print (f"Tiempo para Odeint = {tiempo_odeint} s")
print (f"Tiempo para Eulerint = {tiempo_euler} s")


grad_x = np.gradient(vx_odeint,at)
grad_y = np.gradient(vy_odeint,at)
grad_z = np.gradient(vz_odeint,at)

grad_ax = np.gradient(avx,at)
grad_ay = np.gradient(avy,at)
grad_az = np.gradient(avz,at)

delta_odeint = np.sqrt((x_odeint - ax)**2 + (y_odeint - ay)**2 + (z_odeint - az)**2)
delta_euler = np.sqrt((x_euler - ax)**2 + (y_euler - ay)**2 + (z_euler - az)**2)
print (f"Diferencia de deriva Eulerint y Odeint: {delta_euler[-1]/hr-delta_odeint[-1]/hr} km")


plt.figure(1)
plt.plot(t/hr ,delta_odeint/hr, label="Odeint")
plt.title(f"Diferencia entre posición real y analítica máxima = {delta_odeint[-1]/1000:.3f} (km)")
plt.ylabel(" $\\delta$ [KM]")
plt.xlabel("Tiempo[t/hr]")
plt.tight_layout()
plt.legend()
plt.savefig("Grafico2Odeint.png", dpi = 300)
plt.show()


plt.figure(2)
plt.plot(t/hr,delta_euler/hr, label="Eulerint")
plt.title(f"Diferencia entre posición real y analítica máxima = {delta_euler[-1]/hr:.3f} (km)")
plt.ylabel("Diferencia (km)")
plt.xlabel("Tiempo (hr)]")
plt.tight_layout()
plt.legend()
plt.savefig("Grafico2Eulerint.png", dpi = 300)
plt.show()

plt.figure(3)
plt.plot(t/hr,delta_odeint/hr, label="Odeint")
plt.plot(t/hr,delta_euler/hr, label="Eulerint")
plt.title(f" Diferencia para: Odeint = {delta_odeint[-1]/hr:.3f} (km), Eulerint = {delta_euler[-1]/hr:.3f} (km)")
plt.ylabel(" Diferencia (km)")
plt.xlabel("Tiempo (hr)")
plt.tight_layout()
plt.legend()
plt.savefig("GraficoP2DifEuleyOdeint.png", dpi = 300)
plt.show()

'''


''' Pregunta 3
delta = at[-1]
t = np.linspace(0, delta, len(at))
t1=cron()

sol_odeint= odeint(satelite, z0, t)
t2=cron()

sol_eulerint = eulerint(satelite, z0, t, 1500)
t3=cron()

tiempo_odeint = t2-t1
tiempo_euler = t3-t2

x_odeint = sol_odeint[:,0]
y_odeint = sol_odeint[:,1]
z_odeint = sol_odeint[:,2]
vx_odeint = sol_odeint[:,3]
vy_odeint = sol_odeint[:,4]
vz_odeint = sol_odeint[:,5]

x_euler=sol_eulerint[:,0]
y_euler=sol_eulerint[:,1]
z_euler=sol_eulerint[:,2]


print (f"Tiempo para Odeint = {tiempo_odeint} s")
print (f"Tiempo para Eulerint = {tiempo_euler} s")


grad_x = np.gradient(vx_odeint,at)
grad_y = np.gradient(vy_odeint,at)
grad_z = np.gradient(vz_odeint,at)

grad_ax = np.gradient(avx,at)
grad_ay = np.gradient(avy,at)
grad_az = np.gradient(avz,at)

delta_odeint = np.sqrt((x_odeint - ax)**2 + (y_odeint - ay)**2 + (z_odeint - az)**2)
delta_euler = np.sqrt((x_euler - ax)**2 + (y_euler - ay)**2 + (z_euler - az)**2)
print (f"Diferencia de deriva Eulerint y Odeint: {delta_euler[-1]/hr-delta_odeint[-1]/hr} km")


resta= np.zeros(len(at))
for i in range(len(at)):
    resta[i] = np.sqrt(np.dot((sol_odeint[i,:3] - sol_eulerint[i,:3]), (sol_odeint[i,:3] - sol_eulerint[i,:3])))
error = np.round_(resta[-1]/np.sqrt(np.dot(sol_odeint[-1,:3],sol_odeint[-1,:3])),1)
print (f"Error = {error*100} %")

plt.figure(1)
plt.plot(t/hr,delta_euler/hr, label="Eulerint")
plt.title(f"Diferencia entre posición real y analítica máxima = {delta_euler[-1]/hr:.3f} (km)")
plt.ylabel("Diferencia (km)")
plt.xlabel("Tiempo (hr)]")
plt.tight_layout()
plt.legend()
plt.savefig("Grafico3Eulerint.png", dpi = 300)
plt.show()
'''
sol = odeint(satelite, z0, t)

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]


plt.figure(1)
plt.subplot(3,1,1)
plt.suptitle("Posición")
plt.plot(at/hr, ax/km, label = "Real")
plt.plot(t/hr, x/km, label = "Analítica")
plt.ylabel("$X$ (Km)")
plt.subplot(3,1,2)
plt.plot(at/hr, ay/km, label = "Real")
plt.plot(t/hr, y/km, label = "Analítica")
plt.ylabel("$Y$ (Km)")
plt.subplot(3,1,3)
plt.plot(at/hr, az/km, label = "Real")
plt.plot(t/hr, z/km, label ="Analítica")
plt.ylabel("$Z$ (Km)")
plt.xlabel("Tiempo, $t$ (horas)")
plt.legend()
plt.savefig("GraficoP4J2J3.png", dpi=300)
plt.show()

delta = at[-1]
t = np.linspace(0, delta, len(at))
t1=cron()

sol_odeint= odeint(satelite, z0, t)
t2=cron()

sol_eulerint = eulerint(satelite, z0, t, 1)
t3=cron()

tiempo_odeint = t2-t1
tiempo_euler = t3-t2

x_odeint = sol_odeint[:,0]
y_odeint = sol_odeint[:,1]
z_odeint = sol_odeint[:,2]
vx_odeint = sol_odeint[:,3]
vy_odeint = sol_odeint[:,4]
vz_odeint = sol_odeint[:,5]

x_euler=sol_eulerint[:,0]
y_euler=sol_eulerint[:,1]
z_euler=sol_eulerint[:,2]


print (f"Tiempo para Odeint = {tiempo_odeint} s")
print (f"Tiempo para Eulerint = {tiempo_euler} s")


grad_x = np.gradient(vx_odeint,at)
grad_y = np.gradient(vy_odeint,at)
grad_z = np.gradient(vz_odeint,at)

grad_ax = np.gradient(avx,at)
grad_ay = np.gradient(avy,at)
grad_az = np.gradient(avz,at)

delta_odeint = np.sqrt((x_odeint - ax)**2 + (y_odeint - ay)**2 + (z_odeint - az)**2)
delta_euler = np.sqrt((x_euler - ax)**2 + (y_euler - ay)**2 + (z_euler - az)**2)
print (f"Diferencia de deriva Eulerint y Odeint: {delta_euler[-1]/hr-delta_odeint[-1]/hr} km")


plt.figure(2)
plt.plot(t/hr ,delta_odeint/hr, label="Odeint")
plt.title(f"Diferencia entre posición real y analítica máxima = {delta_odeint[-1]/1000:.3f} (km)")
plt.ylabel(" $\\delta$ [KM]")
plt.xlabel("Tiempo[t/hr]")
plt.tight_layout()
plt.legend()
plt.savefig("Grafico4OdeintJ2J3.png", dpi = 300)
plt.show()


plt.figure(3)
plt.plot(t/hr,delta_euler/hr, label="Eulerint")
plt.title(f"Diferencia entre posición real y analítica máxima = {delta_euler[-1]/hr:.3f} (km)")
plt.ylabel("Diferencia (km)")
plt.xlabel("Tiempo (hr)]")
plt.tight_layout()
plt.legend()
plt.savefig("Grafico4EulerintJ2J3.png", dpi = 300)
plt.show()

plt.figure(4)
plt.plot(t/hr,delta_odeint/hr, label="Odeint")
plt.plot(t/hr,delta_euler/hr, label="Eulerint")
plt.title(f" Diferencia para: Odeint = {delta_odeint[-1]/hr:.3f} (km), Eulerint = {delta_euler[-1]/hr:.3f} (km)")
plt.ylabel(" Diferencia (km)")
plt.xlabel("Tiempo (hr)")
plt.tight_layout()
plt.legend()
plt.savefig("GraficoP4DifEuleyOdeintJ2J3.png", dpi = 300)
plt.show()

ttf = cron()
ttotal = ttf-tti

print(f"Tiempo total de ejecución {ttotal} s")