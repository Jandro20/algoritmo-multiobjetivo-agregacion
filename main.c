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

// Valor random
#define URAND	((double)rand()/((double)RAND_MAX + 1.0))

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

    // Array auxiliar { [0]->indices - [1]-> valores distancia }
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
        //printf("- x = %f - y = %f - \n", vector[i][0], vector[i][1]);

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
void evaluate_zdt3_all(float population[N][T], float evaluation[N][2], float pReference[2]){

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

/*
    Realizamos la evaluacion de un individuo mediante las funciones de ZDT3
*/
void evaluate_zdt3(float individuo[T], float evaluation[2], float pReference[2]){

    float tmp = 0.0;

    float f1 = individuo[0];

    //Realizo sumatorio
    for (int j = 1; j < T; j++)
    {
        tmp += individuo[j];
    }

    float g = 1+((9*tmp)/(T-1));
    float h = 1-sqrt(f1/g)-(f1/g)*sin(10*PI*f1);

    float f2 = g*h;

    // Guardamos la evaluacion para los individuos
    evaluation[0] = f1;
    evaluation[1] = f2;

    // Almaceno los mejores valores de cada objetivo (funciones que minimizan)
    pReference[0] = f1;
    pReference[1] = f2;
    
}

void evaluacionGTE(float y[T], float vPesos[2], float pReference[2], float res){
    
    float f[2];
    float r[2];

    float valor1, valor2;

    evaluate_zdt3(y, f, r);

    valor1 = vPesos[0] * abs( f[0] - pReference[0] );
    valor2 = vPesos[1] * abs( f[1] - pReference[1] );

    if(valor1 > valor2) res = valor1;
    else res = valor2;    
}

/*  ACCIONES POR ITERACION  */

void iteraciones(float population[N][T], int neightbours[N][T], float pReferenceGlobal[2], float pesos[N][2], float mejorPunto[T]){
    
    //Operadores
    float F = 0.5;  //<- Mutacion
    float CR = 0.5; //<- Cruce

    //Salida de la mutacion
    float v[N][T];

    //Salida de la evaluacion
    float evaluation[2];
    
    //Puntos de referencia por evaluacion
    float pReferenceLocal[2];

    float mejorValor =0.;

    // Realizamos G iteraciones x N subproblemas = G*N = 4000
    for (int iteracion = 0; iteracion < G; iteracion++)
    {
        for (int subproblema = 0; subproblema < N; subproblema++)
        {
            /* REPRODUCCION */
            /* MUTACION Y CRUCE */

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

                //printf("%d -> %d - %d - %d \n", subproblema, r1, r2, r3);

                //Inicializo los tres vectores con los que voy a trabajar
                float xr1[T], xr2[T], xr3[T];

                /*
                    OPERADOR DE CRUCE Y MUTACION
                    - Si el valor URAND es menor al de cruce o el subproblema coincide con un valor calculado al azar
                        - Aplico mutacion y añado a v
                    - Sino
                        - Añado el subproblema actual a v
                */

                int prand = (int) URAND*subproblema;                
                

                if (URAND < CR || subproblema == prand)
                {
                    for (int t = 0; t < T; t++) {
                        xr1[t] = population[r1][t];
                        xr2[t] = population[r2][t];
                        xr3[t] = population[r3][t];

                        v[subproblema][t] = xr1[t] + F*(xr2[t] - xr3[t]);

                        //Restrinjo los valores a los limites (0 < x < 1)
                        if(v[subproblema][t] > 1) v[subproblema][t]=1;
                        if(v[subproblema][t] < 0) v[subproblema][t]=0;

                        //printf("%f %f %f = %f\n", xr1[t], xr2[t], xr3[t], v[subproblema][t]);
                    }
                }
                else{
                    for (int t = 0; t < T; t++) {
                        v[subproblema][t] = population[subproblema][t];
                    }
                }

                //SE DA LA MUTACION EN TORNO AL 50% DE LA VECINDAD
           
                //MUTACION GAUSSIANA

                                                   

            /* EVALUACION */
                evaluate_zdt3(v[subproblema], evaluation, pReferenceLocal);
            /* ACTUALIZACION_PUNTO_REFERENCIA */

            if(pReferenceGlobal[0]>pReferenceLocal[0]) {
                pReferenceGlobal[0]=pReferenceLocal[0];
                printf("Global en subproblema %d>  x: %f (%f) \n", subproblema, pReferenceGlobal[0], pReferenceLocal[0]);
            }
            if(pReferenceGlobal[1]>pReferenceLocal[1]) {
                pReferenceGlobal[1]=pReferenceLocal[1];
                printf("Global en subproblema %d>  y: %f (%f) \n", subproblema, pReferenceGlobal[1], pReferenceLocal[1]);
            }
            
            //printf("Local>  x: %f , y: %f \n", pReferenceLocal[0], pReferenceLocal[1]);
            //printf("Global en subproblema %d>  x: %f , y: %f \n", subproblema, pReferenceGlobal[0], pReferenceGlobal[1]);
            
            /* ACTUALIZACION_VECINOS */

            float valorIntermedio = 0.;

            evaluacionGTE(v[subproblema], pesos[subproblema] , pReferenceGlobal, valorIntermedio);

            if(mejorValor <= valorIntermedio){
                mejorValor=valorIntermedio;
                for (int l = 0; l < T; l++)
                {
                    mejorPunto[l]=v[subproblema][l];
                }
                
            }
                
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
    
    //Resultado
    float mejorPunto[T];

    // Mejoramos la aleatoriedad de los rand()
    srand(time(NULL));

    /* INICIALIZACION */
    inicialiceAlphaVector(alpha_vector);
    inicialiceNeighbourVector(alpha_vector, neightbours);
    inicialicePopulation(population);
    evaluate_zdt3_all(population, evaluation, pReference);

    /* OPERACIONES */
    iteraciones(population, neightbours, pReference, alpha_vector, mejorPunto);

    printf("GLOBAL ---->  x: %f , y: %f \n", pReference[0], pReference[1]);
    printf("Mejores valores ");
    for (int p = 0; p < T; p++)
    {
        printf("%f \n", mejorPunto[p]);
    }


    return 0;
}