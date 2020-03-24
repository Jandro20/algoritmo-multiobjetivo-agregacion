#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include "string.h"

#include "general.c"


/* 
    Inicializamos vectores peso con vectores equiespacioados de forma euclidea
*/
void inicialiceAlphaVector(float vector[N][2]){

    float paso = (float) 1/(N - 1);

    vector[0][0] = 1.0;
    vector[0][1] = 0.0;

    for (int i = 1; i < N; i++)
    {
        vector[i][0] = vector[i-1][0] - paso; 
        vector[i][1] = vector[i-1][1] + paso; 
    }
    
}

/*
    Calculo de la distancia euclidea de dos puntos (x1, y1) y (x2,y2)
*/
float euclidea_distance(float p1[2], float p2[2]){

    return (float) sqrt(pow((p2[0] - p1[0]), 2) - pow((p2[1] - p1[1]), 2));
}

/*
    Funcion para comparacion y ordenacion de un array de elementos
*/
int comparation_function(const void * p1, const void * p2){
    return ( *(int *)p1 - *(int *) p2);
}

/*
    Realizamos la vecindad de cada vector, con sus T-1 vecinos (vector incluido)
*/
void inicialiceNeighbourVector(float vector[N][2], int neightbours[N][T]){
    
    float aux[N];
    float *vecinos = malloc(T * sizeof(float));

    // ELemento i -> Calcular distancia euclidea con todo valor j
    for (int i = 0; i < N; i++)
    {
        
        /* 
        Obtengo array con las distancias euclideas desde i a todos los demás
        */
        for (int j = 0; j < N; j++)
        {   
            aux[j] = euclidea_distance(vector[i], vector[j]);
            printf("valor %d : %f \n",j, aux[j]);
        }

        /*
            Ordeno el array de menos a mayor distancia
            y obtengo los T primeros elementos.
        */
        qsort(aux, N, sizeof(float), comparation_function);

        vecinos = aux;
        
        
    }
    
}

int main(){
    alpha_vector[N][2];
    neightbours[N][T];
    inicialiceAlphaVector(alpha_vector);
    inicialiceNeighbourVector(alpha_vector, neightbours);

    return 0;
}