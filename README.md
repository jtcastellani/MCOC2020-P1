# MCOC2020-P1
# Integración de ecuaciones diferenciales
![Entrega1](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%201/Entrega_1.png)
# Primeras predicciones con la EDM básica del satélite
![G1](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%202/Grafico1.png) 
![G2](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%202/Grafico2.png) 
![G3](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%202/Grafico3.png) 
* Si se coloca el satélite a una altura de 700km (en la dirección x, ojo con el radio de la tierra) y se le da una velocidad inicial tangencial y'(0)= vt ¿Cuanto debe valer vt de modo que el satélite efectivamente orbite sin caer de vuelta dentro de la atmosfera (asuma que esta comienza a una altura de 80km) ?  
Se encontró de manera aproximada que el valor a la velocidad mínima para que no entre a la atmósfera terrestre es de vt=8.020 m/s  
* ¿Como encontró vt?  
Básicamente por tanteo al encontrar que la curva no toque el eje de la atmósfera.
# Mejoras al modelo y estudio de convergencia
* P1. (2pt) Grafíque, como arriba, la posición (x,y,z) en el tiempo del vector de estado de Sentinel 1A/B que le tocó. Para esto, descargue y utilice la función leer_eof.py (Enlaces a un sitio externo.) para poder trabajar con los archivos EOF. Sintaxis  
![P1](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/GraficoP1.png)  
* P2. (5pt) Usando la condición inicial (primer OSV) de su archivo, compare la solución entre odeint y eulerint. Use Nsubdiviciones=1. Grafíque la deriva en el tiempo como arriba ¿Cuánto deriva eulerint de odeint en este caso al final del tiempo? (Esta pregunta solo compara algoritmos, no se usa más que la condición inicial del archivo EOF). ¿Cuanto se demora odeint y eulerint respectivamente en producir los resultados?  
![P2Euler](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/Grafico2Eulerint.png)
![P2Odeint](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/Grafico2Odeint.png)
![P2EulerVsOdeint](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/GraficoP2DifEuleyOdeint.png)  
Luego de varias repeticiones los tiempos de demora fueron aproximadamente:  
Odeint: 0,3  
Eulerint: 0,8  
* P3. (3pt) ¿Cuantas subdivisiones hay que usar para que la predicción con eulerint al final del tiempo esté en menos de un 1% de error? Grafique la deriva en el tiempo. Comente con respecto del tiempo de ejecución de eulerint ahora.  
!.[P3](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/Grafico3Eulerint.png)
Luego de 20 minutos y 1.500 subdivisiones se logró un error de 10%, por lo que no se pudo lograr al pocentaje de error deseado, ya que al minimizar el error aumenta el tiempo de ejecución.  
* P4. (5pt) Implemente las correcciones J2 y J3. Repita los gráficos de arriba para su caso particular. ¿Cuánta deriva incurre al agregar las correcciones J2 y J3? ¿Cuanto se demora su código en correr?  
![P4XYZL2L3](https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/GraficoP4J2J3.png)
![P4OEulerJ2J3]https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/Grafico4EulerintJ2J3.png)
![P4OdeintJ2J3]https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/Grafico4OdeintJ2J3.png)
![P4VsJ2J3]https://raw.githubusercontent.com/jtcastellani/MCOC2020-P1/master/Entrega%205/GraficoP4DifEuleyOdeintJ2J3.png)
El código demora 4,34 s
