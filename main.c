#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

void swap(float *v1, float *v2){
    float temp = *v1;
    *v1=*v2;
    *v2=temp;
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
void inicialiceNeighbourVector(float vector[N][2], int neightbours[N][T]){

    // Array auxiliar para no machacar valores
    float aux[N][2];

    // ELemento i -> Calcular distancia euclidea con todo valor j
    for (int i = 0; i < N; i++)
    {
        /* 
            Obtengo array con las distancias euclideas desde i a todos los demás
        */
        for (int j = 0; j < N; j++)
        {   
            aux[j][0] = j; // Pnt al ser comparado
            aux[j][1] = euclidea_distance(vector[i], vector[j]);    // Distancia
        }
        printf("- x = %f - y = %f - \n", vector[i][0], vector[i][1]);

        /* 
            Aplicamos ordenacion por bubble sort a las distancias desde el pnt i al resto.
        */
        for (int k1 = 0; k1 < N-1; k1++)
        {
            for (int k2 = 0; k2 < N-k1-1; k2++)
            {
                if(aux[k2][1] > aux[k2+1][1])
                {
                    // Intercambio indices y valores
                    swap(&aux[k2][0], &aux[k2+1][0]);
                    swap(&aux[k2][1], &aux[k2+1][1]);
                }
            }
        }
        
        /*
            Almaceno los indices que forman el conjuntos de vecinos para el subproblema i
        */
        for (int pt = 0; pt < T; pt++)
        {
            neightbours[i][pt] = (int) aux[pt][0];
            int indice = neightbours[i][pt];
            printf("pt(%f , %f) - INDICE: %d -> %d = (%f , %f) \n",vector[i][0], vector[i][1], i, neightbours[i][pt], vector[indice][0], vector[indice][1]);

        }        
        
        printf("\n");

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

void iteraciones(float population[N][T], int neightbours[N][T]){
    
    //Operadores
    float F = 0.5;  //<- Mutacion
    float CR = 0.5; //<- Cruce


    // Realizamos G iteraciones x N subproblemas = G*N = 4000
    for (int iteracion = 0; iteracion < G; iteracion++)
    {
        for (int subproblema = 0; subproblema < N; subproblema++)
        {
            /* REPRODUCCION */
            /* MUTACION Y CRUCE */

            //Conjunto aleatorio
            float vectoresEscogidos[3][T];
            //Seleccionamos un conjuntos aleatorio de indiviuos 
                //Escogo los 3 vecinos de forma aleatoria
                int r1, r2, r3;

                do
                {
                    r1 = (int)rand() % T;
                } while( r1==subproblema );
                do
                {
                    r2 = (int)rand() % T;
                } while( r2==subproblema || r2==r1);
                do
                {
                    r3 = (int)rand() % T;
                } while( r3==subproblema || r3==r1 || r3==r2 );



                int aleatorio = 0;
                for (int nE = 0; nE < 3; nE++)
                {
                    aleatorio = rand() % T;
                    for (int sub = 0; sub < T; sub++)
                    {
                        vectoresEscogidos[nE][sub] = neightbours[subproblema][sub];
                    }
                }
                
            // Aplicacion de operadores evolutivos.
                //Operador de mutacion
                //Resultado de la mutacion
                float v_i[T];
                
                for (int mutation = 0; mutation < T; mutation++)
                {
                    v_i[T] = vectoresEscogidos[0][T] + F*(vectoresEscogidos[1][T] - vectoresEscogidos[2][T]);
                }
                


            /* EVALUACION */

            /* ACTUALIZACION_PUNTO_REFERENCIA */

            /* ACTUALIZACION_VECINOS */
        }
        
    }
    
}

/*  MAIN    */
int main(){
    // Vectores alpha
    float alpha_vector[N][2];

    // Vectores de los vecinos
    int neightbours[N][T];

    // Vector de N individuos
    float population[N][T];

    // Vector para evaluacion [individuos](f1,f2)
    float evaluation[N][2];

    //Vector para punto de referencia (x,y)
    float pReference[2];

    // Mejoramos la aleatoriedad de los rand()
    srand(time(NULL));

    /* INICIALIZACION */
    inicialiceAlphaVector(alpha_vector);
    inicialiceNeighbourVector(alpha_vector, neightbours);
    inicialicePopulation(population);
    evaluate_zdt3(population, evaluation, pReference);

    /* OPERACIONES */
    iteraciones(population, neightbours);


    return 0;
}