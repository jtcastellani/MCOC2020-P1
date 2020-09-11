import numpy as np
from scipy.integrate import odeint
from time import perf_counter as cron
import xml
import xml.etree.ElementTree as ET
from numpy import zeros
import datetime as dt
from sys import argv



def utc2time(utc, ut1, EOF_datetime_format = "%Y-%m-%dT%H:%M:%S.%f"):
	t1 = dt.datetime.strptime(ut1,EOF_datetime_format)
	t2 = dt.datetime.strptime(utc,EOF_datetime_format)
	return (t2 - t1).total_seconds()


def leer_eof(fname):
	tree = ET.parse(fname)
	root = tree.getroot()

	Data_Block = root.find("Data_Block")		
	List_of_OSVs = Data_Block.find("List_of_OSVs")

	count = int(List_of_OSVs.attrib["count"])

	t = zeros(count)
	x = zeros(count)
	y = zeros(count)
	z = zeros(count)
	vx = zeros(count)
	vy = zeros(count)
	vz = zeros(count)

	set_ut1 = False
	for i, osv in enumerate(List_of_OSVs):
		UTC = osv.find("UTC").text[4:]
		
		x[i] = osv.find("X").text   #conversion de string a double es implicita
		y[i] = osv.find("Y").text
		z[i] = osv.find("Z").text
		vx[i] = osv.find("VX").text
		vy[i] = osv.find("VY").text
		vz[i] = osv.find("VZ").text

		if not set_ut1:
			ut1 = UTC
			set_ut1 = True

		t[i] = utc2time(UTC, ut1)

	return t, x, y, z, vx, vy, vz


tti = cron()

archivo = argv[1]
#archivo = "S1B_OPER_AUX_POEORB_OPOD_20200824T111152_V20200803T225942_20200805T005942.EOF"
at, ax, ay, az, avx, avy, avz = leer_eof(archivo)

Pred = archivo.replace(".EOF",".PRED")

hr = 3600. #s
km = 1000. #m

Omg = 7.2921*10**-5 #rad/s
Mtierra = 5.97219 * 10**24 #kg
Rtierra = 6371*km #m
G = 6.67430*10**-11 #N*m^2/kg^2
Atmosfera = Rtierra + 80.*km #metros
hi = 700.* km # altura inicial

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

    zp[3:6] = R.T@(Fg + FJ2 + FJ3 - (2 * (R_p@xp) + (R_p_p@x))) 
   
    return zp



t = at
x0 = Rtierra + hi
vt = 6820.*3.6 #m/s

z0 = np.array([ax[0], ay[0], az[0], avx[0], avy[0], avz[0]])


sol = odeint(satelite, z0, t )

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]


delta = at[-1]
t = np.linspace(0, delta, len(at))
t1=cron()

sol_odeint= odeint(satelite, z0, t)
t2=cron()


tiempo_odeint = t2-t1


x_odeint = sol_odeint[:,0]
y_odeint = sol_odeint[:,1]
z_odeint = sol_odeint[:,2]
vx_odeint = sol_odeint[:,3]
vy_odeint = sol_odeint[:,4]
vz_odeint = sol_odeint[:,5]


grad_x = np.gradient(vx_odeint,at)
grad_y = np.gradient(vy_odeint,at)
grad_z = np.gradient(vz_odeint,at)

grad_ax = np.gradient(avx,at)
grad_ay = np.gradient(avy,at)
grad_az = np.gradient(avz,at)

deriva = np.sqrt((x_odeint[-1] - ax[-1])**2 + (y_odeint[-1] - ay[-1])**2 + (z_odeint[-1] - az[-1])**2)

#print (f"Deriva: {deriva/km} km")


ttf = cron()
ttotal = ttf-tti

#print(f"Tiempo total de ejecución {ttotal} s")







with open(Pred,"w") as fout:
    Nt = len(t)
    fout.write("<?xml version=\"1.0\"?>\n")
    fout.write("<Earth_Explorer_File>\n")
    fout.write("  <Earth_Explorer_Header>\n")
    fout.write("    <Fixed_Header>\n")
    fout.write(f"      <File_Name>{archivo}</File_Name>\n")
    fout.write("      <File_Description>Precise Orbit Ephemerides (POE) Orbit File</File_Description>\n")
    fout.write("      <Notes></Notes>\n")
    fout.write("      <Mission>Sentinel-1B</Mission>\n")
    fout.write("      <File_Class>OPER</File_Class>\n")
    fout.write("      <File_Type>AUX_POEORB</File_Type>\n")
    fout.write("      <Validity_Period>\n")
    fout.write("        <Validity_Start>UTC=2020-08-03T22:59:42</Validity_Start>\n")
    fout.write("        <Validity_Stop>UTC=2020-08-05T00:59:42</Validity_Stop>\n")
    fout.write("      </Validity_Period>\n")
    fout.write("      <File_Version>0001</File_Version>\n")
    fout.write("      <Source>\n")
    fout.write("        <System>OPOD</System>\n")
    fout.write("        <Creator>OPOD</Creator>\n")
    fout.write("        <Creator_Version>0.0</Creator_Version>\n")
    fout.write("        <Creation_Date>UTC=2020-08-24T11:11:52</Creation_Date>\n")
    fout.write("      </Source>\n")
    fout.write("    </Fixed_Header>\n")
    fout.write("    <Variable_Header>\n")
    fout.write("      <Ref_Frame>EARTH_FIXED</Ref_Frame>\n")
    fout.write("      <Time_Reference>UTC</Time_Reference>\n")
    fout.write("    </Variable_Header>\n")
    fout.write("  </Earth_Explorer_Header>\n")
    fout.write("<Data_Block type=\"xml\">\n")
    fout.write(f"  <List_of_OSVs count=\"{len(at)}\">\n")
    for i in range(len(at)):
        obj = dt.datetime(2020,8,3,22,59,42,000000)
        fecha = (obj + dt.timedelta(seconds=t[i])).strftime("%Y-%m-%dT%H:%M:%S.%f")
        fout.write("    <OSV>\n")
        fout.write(f"      <UTC>UTC={fecha}</UTC>\n")
        fout.write("      <Absolute_Orbit>+22765</Absolute_Orbit>\n")
        fout.write(f"      <X unit=\"m\">{x_odeint[i]}</X>\n")
        fout.write(f"      <Y unit=\"m\">{y_odeint[i]}</Y>\n")
        fout.write(f"      <Z unit=\"m\">{z_odeint[i]}</Z>\n")
        fout.write(f"      <VX unit=\"m/s\">{vx_odeint[i]}</VX>\n")
        fout.write(f"      <VY unit=\"m/s\">{vy_odeint[i]}</VY>\n")
        fout.write(f"      <VZ unit=\"m/s\">{vz_odeint[i]}</VZ>\n")
        fout.write("      <Quality>NOMINAL</Quality>\n")
        fout.write("    </OSV>\n")
    fout.write("  </List_of_OSVs>\n")
    fout.write("</Data_Block>\n</")
    fout.write("Earth_Explorer_File>")       


