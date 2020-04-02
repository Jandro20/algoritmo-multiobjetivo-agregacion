#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gnu-versions.h>
# define GNUPLOT_COMMAND "gnuplot -persist"

// Poblacion
int const N = 200;

// Vecindad (%) 15% de la poblacion total (200*0.15)
int const T = 30;

//Espacio de busqueda
int const EB = 30;

// Generaciones
int const G = 50;

/*  FUNCIONES AUXILIARES    */
float const PI = 3.14159265359;

float const e =  2.7182;

// Valor random
#define URAND	((double)rand()/((double)RAND_MAX + 1.0))

/*
    Calculo de la distancia euclidea de dos puntos (x1, y1) y (x2,y2)
*/
float euclidea_distance(float p1[2], float p2[2])
{
    return (float) sqrt(pow((p2[0] - p1[0]), 2.0) + pow((p2[1] - p1[1]), 2.0));
}

/*
    Funcion para comparacion y ordenacion de un array de elementos
*/
int comparation_function(const void * p1, const void * p2)
{
    return ( *(int *)p1 - *(int *) p2);
}

void swap(float *v1, float *v2)
{
    float temp = *v1;
    *v1=*v2;
    *v2=temp;
}
//FUNCIONES AUX PARA RANDOM
    double seed;
    double oldrand[55];
    int jrand;

    /* Create next batch of 55 random numbers */
    void advance_random ()
    {
        int j1;
        double new_random;
        for(j1=0; j1<24; j1++)
        {
            new_random = oldrand[j1]-oldrand[j1+31];
            if(new_random<0.0)
            {
                new_random = new_random+1.0;
            }
            oldrand[j1] = new_random;
        }
        for(j1=24; j1<55; j1++)
        {
            new_random = oldrand[j1]-oldrand[j1-24];
            if(new_random<0.0)
            {
                new_random = new_random+1.0;
            }
            oldrand[j1] = new_random;
        }
    }

    /* Get randomize off and running */
    void warmup_random (double seed)
    {
        int j1, ii;
        double new_random, prev_random;
        oldrand[54] = seed;
        new_random = 0.000000001;
        prev_random = seed;
        for(j1=1; j1<=54; j1++)
        {
            ii = (21*j1)%54;
            oldrand[ii] = new_random;
            new_random = prev_random-new_random;
            if(new_random<0.0)
            {
                new_random += 1.0;
            }
            prev_random = oldrand[ii];
        }
        advance_random ();
        advance_random ();
        advance_random ();
        jrand = 0;
        return;
    }

    /* Get seed number for random and start it up */
    void randomize()
    {
        int j1;
        for(j1=0; j1<=54; j1++)
        {
            oldrand[j1] = 0.0;
        }
        jrand=0;
        warmup_random (seed);
        return;
    }

    

    /* Fetch a single random number between 0.0 and 1.0 */
    double randomperc()
    {
        jrand++;
        if(jrand>=55)
        {
            jrand = 1;
            advance_random();
        }
        return((double)oldrand[jrand]);
    }

    /* Fetch a single random integer between low and high including the bounds */
    int rnd (int low, int high)
    {
        int res;
        if (low >= high)
        {
            res = low;
        }
        else
        {
            res = low + (randomperc()*(high-low+1));
            if (res > high)
            {
                res = high;
            }
        }
        return (res);
    }

    /* Fetch a single random real number between low and high including the bounds */
    double rndreal (double low, double high)
    {
        return (low + (high-low)*randomperc());
    }



/*  FUNCIONES PRINCIPALES   */

/*  INICIALIZACION  */

/* 
    Inicializamos vectores peso con vectores equiespacioados de forma euclidea
*/
void inicialiceAlphaVector(float vector[N][2])
{

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
void inicialiceNeighbourVector(float vector[N][2], int neightbours[N][T])
{


    // ELemento i -> Calcular distancia euclidea con todo valor j
    for (int i = 0; i < N; i++)
    {
        // Array auxiliar { [0]->indices - [1]-> valores distancia }
        float aux[N][2];
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
void inicialicePopulation(float population[N][EB])
{
    
    float range = 1.0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < EB; j++)
        {
            population[i][j] = (float)rand()/(float)(RAND_MAX/range);
        }

    }
    
}

/*
    Realizamos la evaluacion de la poblacion mediante las funciones de ZDT3
*/
void evaluate_zdt3_all(float population[N][EB], float evaluation[N][2], float pReference[2])
{

    float tmp = 0.0;

    float bestX = 1000000.0;
    float bestY = 1000000.0;

    for (int individuo = 0; individuo < N; individuo++)
    {
        float f1 = population[individuo][0];

        //Limpio tmp para siguiente iteracion
        tmp = 0.0;

        //Realizo sumatorio
        for (int j = 1; j < EB; j++)
        {
            tmp += population[individuo][j];
        }

        float g = 1+((9*tmp)/(EB-1));
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
void evaluate_zdt3(float individuo[EB], float evaluation[2], float pReference[2])
{

    float tmp = 0.0;

    float f1 = individuo[0];

    //Realizo sumatorio
    for (int j = 1; j < EB; j++)
    {
        tmp += individuo[j];
    }

    float g = 1+((9*tmp)/(EB-1));
    float h = 1-sqrt(f1/g)-(f1/g)*sin(10*PI*f1);

    float f2 = g*h;

    // Guardamos la evaluacion para los individuos
    evaluation[0] = f1;
    evaluation[1] = f2;


    // Almaceno los mejores valores de cada objetivo (funciones que minimizan)
    pReference[0] = f1;
    pReference[1] = f2;
    
}

float evaluacionGTE(float individuo[EB], float vPesos[2], float pReference[2])
{
    
    float f[2];
    float r[2];

    float valor1, valor2, res;

    evaluate_zdt3(individuo, f, r);

    valor1 = vPesos[0] * abs( f[0] - pReference[0] );
    valor2 = vPesos[1] * abs( f[1] - pReference[1] );

    if(valor1 > valor2){
        res = valor1;
    
    }else{
        res = valor2;
    } 


    return res;
}

/*  ACCIONES POR ITERACION  */

//TORNEO Y SELECCION
    void cruce(float padre1[EB], float padre2[EB], float hijo[EB], float hijo2[EB])
    {

        double CR = 0.7;

        double prand =  URAND;
        double y1, y2, yl, yu;
        double rand;
        double c1, c2;
        double alpha, beta, betaq;
        // Indice de distribucion para el cruce
        int eta_c = 10;

        if(prand <= CR){
            for (int i = 0; i < EB; i++)
            {
                if (randomperc()<=0.5 )
                {
                    if (fabs(padre1[i]-padre2[i]) > 1.0e-14)
                    {
                        if(padre1[i] < padre2[i])
                        {
                            y1 = padre1[i];
                            y2 = padre2[i];
                        }else
                        {
                            y1 = padre2[i];
                            y2 = padre1[i];
                        }

                        yl = 1.0;
                        yu = 0.0;

                        rand = randomperc();
                        beta = 1.0 + (.0)*(y1-yl)/(y2-y1);
                        alpha = 2.0 - pow(beta,-(eta_c+1.0)); 
                        if (rand <= (1.0/alpha))
                        {
                            betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                        }
                        else{
                            betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                        }
                        c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                        if (c1<yl)
                            c1=yl;
                        if (c2<yl)
                            c2=yl;
                        if (c1>yu)
                            c1=yu;
                        if (c2>yu)
                            c2=yu;
                        if (randomperc()<=0.5){
                            hijo[i] = c1;
                            hijo2[i]= c2;
                        }   
                        else
                        {
                            hijo[i] = c2;
                            hijo2[i] = c1;
                        }
                        }
                    else
                    {
                        hijo[i] = padre1[i];
                        hijo2[i] = padre2[i];
                    }
                    }
                    else
                    {
                        hijo[i] = padre2[i];
                        hijo2[i] = padre1[i];
                    }
                }
            }
            else
            {
                for (int i=0; i<EB; i++)
                {
                    hijo[i] = padre1[i];
                    hijo2[i] = padre2[i];
                }
            }

    }
/* Routine for binary tournament */
    void tournament (float individuo1[EB], float individuo2[EB], float parent[EB])
    {
        if ((randomperc()) <= 0.5)
        {
            for (int i = 0; i < EB; i++)
            {
                parent[i] = individuo1[i];
            }
            
        }
        else
        {
            for (int i = 0; i < EB; i++)
            {
                parent[i] = individuo2[i];
            }
            
        }
    }

    void selection (float OldPopulation[N][EB], float NewPopulation[N][EB])
    {
        int popsize = N;
        int *a1, *a2;
        int temp;
        int i;
        int rand;
        float parent1[EB], parent2[EB];
        a1 = (int *)malloc(popsize*sizeof(int));
        a2 = (int *)malloc(popsize*sizeof(int));
        for (i=0; i<popsize; i++)
        {
            a1[i] = a2[i] = i;
        }
        for (i=0; i<popsize; i++)
        {
            rand = rnd (i, popsize-1);
            temp = a1[rand];
            a1[rand] = a1[i];
            a1[i] = temp;
            rand = rnd (i, popsize-1);
            temp = a2[rand];
            a2[rand] = a2[i];
            a2[i] = temp;
        }
        for (i=0; i<popsize; i+=4)
        {
            tournament (OldPopulation[i], OldPopulation[i+1], parent1);
            tournament (OldPopulation[i+2], OldPopulation[i+3], parent2);
            cruce (parent1, parent2, NewPopulation[i], NewPopulation[i+1]);
            tournament (OldPopulation[i], OldPopulation[i+1], parent1);
            tournament (OldPopulation[i+2], OldPopulation[i+3], parent2);
            cruce (parent1, parent2, NewPopulation[i+2], NewPopulation[i+3]);
        }
        free (a1);
        free (a2);
        return;
    }

/* Routine to copy an individual 'ind1' into another individual 'ind2' */
    void copy_ind (float individuo1[EB], float individuo2[EB])
    {
        int i, j;
        if (T !=0)
        {
            for (i=0; i<EB; i++)
            {
                individuo2[i] = individuo1[i];
            }
        }
        return;
    }

    /* Routine to merge two populations into one */
    void merge(float OldPopulation1[N][EB], float OldPopulation2[N][EB], float NewPopulation[N][EB])
    {
        int popsize = N;
        int i, k;
        for (i=0; i<popsize; i++)
        {
            copy_ind (OldPopulation1[i], NewPopulation[i]);
        }
        for (i=0, k=popsize; i<popsize; i++, k++)
        {
            copy_ind (OldPopulation2[i], NewPopulation[i]);
        }
        return;
    }


void iteraciones(float population[N][EB], float population2[N][EB], int neightbours[N][T], float pReferenceGlobal[2], float pesos[N][2], FILE *p, FILE *p2){
    
    // DEFINICION DE VARIABLES LOCALES
        //Operadores
        //ORIGINAL 0.5 ambos
        float F = 0.5;  //<- Mutacion
        float CR = 0.5; //<- Cruce

        //Salida de la mutacion
        float mutation[EB];

        //Salida de la evaluacion
        float evaluation[2];
        
        //Puntos de referencia por evaluacion
        float pReferenceLocal[2];

        float T_j = (30.-1.)/20.;
        float PR = 1./30.;

    // Realizamos G iteraciones x N subproblemas = G*N = 4000 or 10000
    for (int iteracion = 0; iteracion < G; iteracion++)
    {

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < EB; j++)
            {
                fprintf(p2, "%f ", population2[i][j]);
            }
            fprintf(p2, "\n");
        }
        fprintf(p2, "\n");
        

        selection(population, population2);
        
        for (int subproblema = 0; subproblema < N; subproblema++)
        {
            
            
            /* REPRODUCCION */
            /* MUTACION Y CRUCE */

                // CRUCE

                //  MUTACION
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
                float xr1[EB], xr2[EB], xr3[EB];

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
                    for (int t = 0; t < EB; t++) {
                        xr1[t] = population2[neightbours[subproblema][r1]][t];
                        xr2[t] = population2[neightbours[subproblema][r2]][t];
                        xr3[t] = population2[neightbours[subproblema][r3]][t];

                        mutation[t] = xr1[t] + F*(xr3[t] - xr2[t]);

                        //Restrinjo los valores a los limites (0 < x < 1)
                        if(mutation[t] > 1) mutation[t]=1;
                        if(mutation[t] < 0) mutation[t]=0;

                        //printf("%f %f %f = %f\n", xr1[t], xr2[t], xr3[t], v[subproblema][t]);
                    }
                }
                else{
                    for (int t = 0; t < T; t++) {
                        mutation[t] = population2[subproblema][t];
                    }
                }       

                //SE DA LA MUTACION EN TORNO AL 50% DE LA VECINDAD


                float ran = URAND;
                if(ran < PR){
                    for (int eb = 0; eb < EB; eb++)
                    {
                        mutation[eb] = mutation[eb] + (float)((1/(sqrt(2*PI)*T_j))*pow(e, -((pow(mutation[eb],2)/2*pow(T_j, 2)))));
                    }               
                }
                


            /* EVALUACION */
                evaluate_zdt3(mutation, evaluation, pReferenceLocal);
            /* ACTUALIZACION_PUNTO_REFERENCIA */
                
                if(pReferenceGlobal[0]>pReferenceLocal[0]) {
                    pReferenceGlobal[0]=pReferenceLocal[0];
                    //printf("Global en subproblema %d>  x: %f (%f) \n", subproblema, pReferenceGlobal[0], pReferenceLocal[0]);
                }
                if(pReferenceGlobal[1]>pReferenceLocal[1]) {
                    pReferenceGlobal[1]=pReferenceLocal[1];
                    //printf("Global en subproblema %d>  y: %f (%f) \n", subproblema, pReferenceGlobal[1], pReferenceLocal[1]);
                }
            
           
            /* ACTUALIZACION_VECINOS */

                //Valor del individuo actual
                float valorActual = evaluacionGTE(mutation, pesos[subproblema], pReferenceGlobal);
                float valorVecino = 0.;
                int mejor = 0;

                float vecinoActual[EB];
                //Calculamos valores para los vecinos del individuo
                for (int vecino = 0; vecino < EB; vecino++)
                {
                    for (int p = 0; p < EB; p++)
                    {
                        vecinoActual[p] = population2[neightbours[subproblema][vecino]][p];
                    }

                    valorVecino = evaluacionGTE(vecinoActual, pesos[subproblema], pReferenceGlobal);
                    
                    if (valorActual >= valorVecino){
                        mejor = 1;
                        //Sustituimos al vecino por la mejor solucion hasta el momento
                        
                        population2[subproblema][neightbours[subproblema][vecino]] = mutation[vecino];
                        
                    }
                }

                //printf("valor actual: %f - valor nuevo: %f - mejor? = %d\n", valorActual, valorVecino, mejor);

            fprintf(p, "%1.3f %1.3f %f \n", evaluation[0], evaluation[1], 0.0);
        }
    }
}

/*  MAIN    */
int main(){
    // Vectores alpha
    float alpha_vector[N][2];

    // Vectores de los vecinos
    int neightbours[N][EB];

    // Vector de N individuos
    float population[N][EB];
    float population2[N][T];

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
    inicialicePopulation(population2);
    
    evaluate_zdt3_all(population, evaluation, pReference);

    char * commandsForGnuplot[] = {"set title \"TITLE\"", "set yrange [-1:6]", "plot 'soluciones.temp'"};
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    FILE * temp = fopen("soluciones.temp", "w");
    FILE * temp2 = fopen("individuos.temp", "w");
    
    
    /* OPERACIONES */
    iteraciones(population,population2, neightbours, pReference, alpha_vector, temp, temp2);

    //fprintf(temp, "%1.3f %1.3f \n", pReference[0], pReference[1]);

    for (int i=0; i < 3; i++)
    {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }

    return 0;
}