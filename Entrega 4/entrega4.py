# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 22:48:40 2020

@author: Joseto
"""
  
from matplotlib.pylab import *
from scipy.integrate import odeint
import numpy as np

a = -1.

m = 1 #kg
f = 1.0 #Hz
e = 0.2
w = 2 * np.pi * f
k = m * w**2
c = 2 * e * w * m

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


def zp(z, t):
    zp = np.zeros(2)
    zp[0] = z[1]
    zp[1] = -(c/m) * z[1] - (k/m) * z[0]
    return zp

z0 = np.array([0, 1])
y = c/(2*m)
A = np.sqrt(z0[0]**2 + ((z0[1] + y*z0[0])/w)**2)
t = np.linspace(0, 4., 100)

z_odeint = odeint(zp, z0, t)
z_real = A * np.exp(-y*t) * np.sin(w*t)

z_euler1 = eulerint(zp, z0, t, 1)
z_euler10 = eulerint(zp, z0, t, 10)
z_euler100 = eulerint(zp, z0, t, 100)

figure(1)
plot(t, z_real,color = "black", linewidth = 2,label="real")
plot(t, z_odeint[:,0],"b",label="odeint")
plot(t,z_euler1[:,0], "g:",label="eulerint N=1")
plot(t,z_euler10[:,0], ":",color = "red", label="eulerint N=10")
plot(t,z_euler100[:,0], ":", color = "darkorange", label="eulerint N=100")
xlabel('t (s)')
ylabel('x(t) (m)')
legend()
grid()
savefig("Grafico4.png", dpi=300)

show()