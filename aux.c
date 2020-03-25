#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include "string.h"

/*
    Calculo de la distancia euclidea de dos puntos (x1, y1) y (x2,y2)
*/
float euclidea_distance(float p1[2], float p2[2]){
    return (float) sqrt(pow((p2[0] - p1[0]), 2.0) + pow((p2[1] - p1[1]), 2.0));
}

/*
    Funcion para comparacion y ordenacion de un array de elementos
*/
int comparation_function(const void * p1, const void * p2){
    return ( *(int *)p1 - *(int *) p2);
}