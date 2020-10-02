# ASC
 Repositorio para la resolución del proyecto para la asignatura ASC (Aplicaciones de Soft Computing) (2020)

## Objetivo del proyecto
Realizar una implementación de un algoritmo multiobjetivo basado en agregación, con posibilidad de manejo de restricciones y cuyas prestaciones sobre los problemas propuestos superen un umbral concreto.

## Breve descripción
Consiste en la creación de un algoritmo que trabaja descomponiendo el problema de optimización multi-objetivos en varios subproblemas de un único objetivo. Para ello, vamos a utilizar la formulación de Tchebychef.
El problema a solucionar se llama ZDT3, y consiste en obtener el frontal pareto-óptimo más aproximado al ideal.

[![N|solid](https://www.researchgate.net/profile/Luis_Gonzalez54/publication/27467964/figure/fig1/AS:339723550773251@1458007814832/Figure4-Pareto-optimal-front-for-ZDT3.png)](https://www.researchgate.net/figure/Figure4-Pareto-optimal-front-for-ZDT3_fig1_27467964)

La carpeta del proyecto contiene dos ficheros propios necesarios para las comparaciones:
- PF.dat: Frente paleto ideal del problema ZDT3.
- all _ popm _ seed9.out: Contiene las soluciones de todas las generaciones generados por el profesor.
## Ejecución
Para ejecutar el proyecto, en un entorno GNU/Linux, en este caso, Debian 10, necesitamos tener instalado GNUPlot. 

#### Instalación GNUPlot
```sh
$ sudo apt-get update
$ sudo apt-get install gnuplot -y
```

### Compilación del proyecto
Una vez instalado, vamos a compilar el proyecto, para ello, ejecutamos el siguiente comando sobre la carpeta donde se encuentra el archivo "main.c".

```sh
$ gcc main.c -lm 
$ ./a.out
```

Con la opción -lm, le indicamos al compilador que añada las librerias de matemáticas, necesarias para este proyecto.

## Resultados
Se obtienen dos gráficas resultantes, el frente paleto y todas las generaciones generadas. Cada una de ellas tiene una comparación entre los resultados generados por el algoritmo junto con los resultados obtenidos por el profesor.
En el frente, la amarilla son los puntos generados por el algoritmo y el morado por el profesor. En cambio, en la representación de todas las generaciones, en verde se ven los puntos generados por el algoritmo y en verde las del profesor.

En la carpeta del proyecto, se generan varios archivos .temp:
- vpesos.temp: Contiene todos los pesos de cada individuo.
- solucionesU.temp: Contiene los resultados de la última generación.
- soluciones.temp: Contiene todas las soluciones de cada generacion, generadas por el algoritmo.
- puntosReferencias.temp: Contiene el mejor punto de referencia de cada generación.
- individuos.temp: Contiene el espacio de búsqueda de cada individuo durante cada generación.

