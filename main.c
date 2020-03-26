#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Poblacion
int const N = 200;

// Vecindad (%) 15% de la poblacion total (200*0.15)
int const T = 30;

// Generaciones
int const G = 200;


/*  FUNCIONES AUXILIARES    */
int const PI = 3.14159265359;

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

/*  FUNCIONES PRINCIPALES   */

/*  INICIALIZACION  */

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
    Realizamos la vecindad de cada vector, con sus T-1 vecinos (vector incluido)
*/
void inicialiceNeighbourVector(float vector[N][2], float neightbours[N][T]){

    // Array para obtener los T primeros valores del array.
    float *vecinos = (float *)malloc(T * sizeof(float));

    // Array auxiliar para no machacar valores
    float aux[N];

    // ELemento i -> Calcular distancia euclidea con todo valor j
    for (int i = 0; i < N; i++)
    {
        
        /* 
        Obtengo array con las distancias euclideas desde i a todos los demás
        */
        for (int j = 0; j < N; j++)
        {   
            aux[j] = euclidea_distance(vector[i], vector[j]);
        }

        /*
            Array ordenado, obtengo los T primeros valores y lo añado a la lista de vecinos.
        */
        
        vecinos = aux;
        
        for (int v = 0; v < T; v++)
        {
            neightbours[i][v] = vecinos[v];
        }
         
    }
    
}

/*
    Inicializamos una poblacion de N individuos (cada uno con T dimensiones)
    en un rango de 0 <= xi <= range
*/
void inicialicePopulation(float population[N][T]){
    
    float range = 1.0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < T; j++)
        {
            population[i][j] = (float)rand()/(float)(RAND_MAX/range);
        }

    }
    
}

/*
    Realizamos la evaluacion de la poblacion mediante las funciones de ZDT3
*/
void evaluate_zdt3(float population[N][T], float evaluation[N][2], float pReference[2]){

    float tmp = 0.0;

    float bestX = 1000000.0;
    float bestY = 1000000.0;

    for (int individuo = 0; individuo < N; individuo++)
    {
        float f1 = population[individuo][0];

        //Limpio tmp para siguiente iteracion
        tmp = 0.0;

        //Realizo sumatorio
        for (int j = 1; j < T; j++)
        {
            tmp += population[individuo][j];
        }

        float g = 1+((9*tmp)/(T-1));
        float h = 1-sqrt(f1/g)-(f1/g)*sin(10*PI*f1);

        float f2 = g*h;

        // Guardamos la evaluacion para los individuos
        evaluation[individuo][0] = f1;
        evaluation[individuo][1] = f2;

        //Vamos comprobando los valores para el punto de referencia
        if(bestX > f1) bestX = f1;
        if(bestY > f2) bestY = f2;
    }

    // Almaceno los mejores valores de cada objetivo (funciones que minimizan)
    pReference[0] = bestX;
    pReference[1] = bestY;
    
}

/*  ACCIONES POR ITERACION  */

void iteraciones(float population[N][T], float neightbours[N][T]){
    
    // Realizamos G iteraciones x N subproblemas = G*N = 4000
    for (int iteracion = 0; iteracion < G; iteracion++)
    {
        for (int subproblema = 0; subproblema < N; subproblema++)
        {
            /*REPRODUCCION*/

            /*EVALUACION*/

            /*ACTUALIZACION_PUNTO_REFERENCIA*/

            /*ACTUALIZACION_VECINOS*/
        }
        
    }
    
}

/*  MAIN    */
int main(){
    // Vectores alpha
    float alpha_vector[N][2];

    // Vectores de los vecinos
    float neightbours[N][T];

    // Vector de N individuos
    float population[N][T];

    // Vector para evaluacion [individuos](f1,f2)
    float evaluation[N][2];

    //Vector para punto de referencia (x,y)
    float pReference[2];

    inicialiceAlphaVector(alpha_vector);
    inicialiceNeighbourVector(alpha_vector, neightbours);
    inicialicePopulation(population);
    evaluate_zdt3(population, evaluation, pReference);

    printf("%f", pReference[0]);

    return 0;
}