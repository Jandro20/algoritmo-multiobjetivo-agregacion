#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gnu-versions.h>


# define GNUPLOT_COMMAND "gnuplot -persist"

typedef FILE * Fichero;

// Poblacion
int const N = 100;

// Vecindad (%) 20% de la poblacion total (200*0.15)
int const T = 20;

//Espacio de busqueda
int const EB = 30;

// Generaciones
int const G = 100;

/*  FUNCIONES AUXILIARES    */
float const PI = 3.14159265359;

float const e =  2.7182;

float const F = 0.5;
float const CR = 0.5;

//mandar abajo al terminar
char * commandsForGnuplot2[] = {"set title \"Todas las generaciones\"", "set yrange [-0.8:1]", "plot 'soluciones.temp', 'all_popm_seed9.out'"};

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




/* 
    Inicializamos vectores peso con vectores equiespacioados de forma euclidea
*/
void inicialiceAlphaVector(float vector[N][2])
{

    float paso = (float) 1/(N - 1);

    vector[0][0] = 1.0;
    vector[0][1] = 0.0;

    for (int i = 1; i < N-1; i++)
    {
        vector[i][0] = vector[i-1][0] - paso; 
        vector[i][1] = vector[i-1][1] + paso; 
    }

    vector[N-1][0] = 0.;
    vector[N-1][1] = 1.;
    
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
    float g,h,f1, f2;

    f1 = individuo[0];

    //Realizo sumatorio
    for (int eb = 1; eb < EB; eb++)
    {
        tmp += individuo[eb];
    }

    g = 1+((9*tmp)/(EB-1));
    h = 1-sqrt(f1/g)-(individuo[0]/g)*sin(10*PI*individuo[0]);

    f2 = g*h;

    /* if(isnan(f1) || isnan(f2)){
        for (int i = 0; i < EB; i++)
        {
            printf("%f \n", individuo[i]);
        }
        printf("\n");
    } */

    // Guardamos la evaluacion para el individuo
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

    valor1 = vPesos[0] * (float) fabsf( f[0] - pReference[0] );
    valor2 = vPesos[1] * (float) fabsf( f[1] - pReference[1] );

    //printf("v2: %f * (%f - %f) \n", vPesos[1], f[1], pReference[1]);

    //printf("gte: valor1 %f valor2 %f \n", valor1, valor2);

    if(valor1 >= valor2){
        res = valor1;
    }else{
        res = valor2;
    } 

    return res;
}

/*  ACCIONES POR ITERACION  */


void dominancia(Fichero fsoluciones){
    float variable[3];
    if(fsoluciones == NULL)
    {
        printf("error");
        exit(1);
    }
    else
    {
        //while (!feof(fsoluciones)) {
            fscanf(fsoluciones, "%f %f %f", &variable[0], &variable[1], &variable[2]);
            printf("%f %f %f \n", variable[0], variable[1], variable[2]);
        //}
    }
}

void fmutation(float population[N][EB], int subproblema, int neightbours[N][T], float mutation[EB], Fichero fprueba){
    
    int r1, r2, r3; 

    do
    {
        r1 = (int) (URAND*T);
    } while( r1==subproblema );
    do
    {
        r2 = (int) (URAND*T);
    } while( r2==subproblema || r2==r1);
    do
    {
        r3 = (int) (URAND*T);
    } while( r3==subproblema || r3==r1 || r3==r2 );
    
    float xr1[EB], xr2[EB], xr3[EB];

    for (int t = 0; t < EB; t++) {

        xr1[t] = population[neightbours[subproblema][r1]][t];
        xr2[t] = population[neightbours[subproblema][r2]][t];
        xr3[t] = population[neightbours[subproblema][r3]][t];

        mutation[t] = xr1[t] + F*(xr2[t] - xr3[t]);

        //printf("MUTATION: sub: %d vecinos: %f (%d) + 0.5*(%f (%d) -  %f (%d)) = %f \n", subproblema, xr1[t], r1, xr2[t], r2, xr3[t], r3, mutation[t] );
        if(mutation[t] > 1.){ 
            mutation[t]=1.;
        }
        if(mutation[t] < 0.) {
            mutation[t]=0.;
        }
        //printf("%f %f %f = %f\n", xr1[t], xr2[t], xr3[t], mutation[subproblema][t]);
        #ifdef degug
        fprintf(fprueba, "MUTATION: sub: %d vecinos: %f (%d) + 0.5*(%f (%d) -  %f (%d)) = %f \n", subproblema, xr1[t], r1, xr2[t], r2, xr3[t], r3, mutation[t] );
        #endif
    }
}


void iteraciones(float population[N][EB], int neightbours[N][T], float pReferenceGlobal[2], float pesos[N][2]
                    , Fichero fsoluciones, Fichero findividuos, Fichero fpRefencia, Fichero fUsoluciones, Fichero fprueba){
    
    // DEFINICION DE VARIABLES LOCALES

        //Salida de la mutacion
        float mutation[EB];

        //y
        float posibleIndividuo[EB];

        //Salida de la evaluacion
        float evaluation[2];
        float evaluationP[2];

        float evaluationN[N][2];

        evaluation[0] = 0.;
        evaluation[1] = 0.;

        //Puntos de referencia por evaluacion
        float pReferenceLocal[2];

        float T_j = (float) (1.-0.)/20.;
        float PR = (float) 1./(float)EB;

    // Realizamos G iteraciones x N subproblemas = G*N = 4000 or 10000
    for (int iteracion = 0; iteracion < G; iteracion++)
    {

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < EB; j++)
            {
                fprintf(findividuos, "%f ", population[i][j]);
            }
            fprintf(findividuos, "\n");
        }
        fprintf(findividuos, "\n");

        for (int subproblema = 0; subproblema < N; subproblema++)
        {
            //  MUTACION
            #ifdef degug
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < EB; j++)
                {
                    fprintf(fprueba, "%f ", population[i][j]);
                }
                fprintf(fprueba, "\n");
            }
            fprintf(fprueba, "\n");
            #endif
            fmutation(population, subproblema, neightbours, mutation, fprueba);

            //CRUCE:
            float rd = 0.;
            int delta = (int)rand() % EB;

            #ifdef debug
            fprintf(fprueba, "CRUCE: sub: %d \n", subproblema);
            #endif

            for (int eb = 0; eb < EB; eb++)
            {
                rd = URAND;
                if(rd <= CR || eb == delta){
                    posibleIndividuo[eb] = mutation[eb];
                    
                    #ifdef degug
                    fprintf(fprueba, "EB: %d - rd: %f, mutation: %f \n", eb, subproblema, rd, mutation[eb]);
                    #endif

                }else if(rd > CR || eb != delta){
                    posibleIndividuo[eb] = population[subproblema][eb];
                    
                    #ifdef degug
                    fprintf(fprueba, "EB: %d - rd: %f, population: %f \n", eb, subproblema, rd, population[subproblema][eb]);
                    #endif

                }
            }

            
            
            //Mutacion gausiana
            float dpG = 0.;
            float nA, nB;
            float ran;

            #ifdef degug
            fprintf(fprueba, "MGaussiana: sub: %d \n", subproblema);
            #endif

            for (int eb = 0; eb < EB; eb++)
            {
                ran = URAND;
                if(ran <= PR){
                    //Distribuciones uniformes en el rango [0,1]
                    nA=(float)(rand())/((float) RAND_MAX );
                    nB=(float)(rand())/((float) RAND_MAX );
                    
                    dpG = (float) sqrt((-2.)*log(nA))*cos(2.*PI*nB);
                    
                    posibleIndividuo[eb] = posibleIndividuo[eb] + dpG*T_j;

                    if(posibleIndividuo[eb] < 0.0){
                        posibleIndividuo[eb] = (float) 0.0;
                    }
                    if(posibleIndividuo[eb] > 1.0){
                        posibleIndividuo[eb] = (float) 1.0;
                    }
                    
                    #ifdef degug
                    fprintf(fprueba, "EB: %d - ran: %f, modificador: %f, posibleIndividuo: %f \n", eb, subproblema, ran, dpG*T_j, posibleIndividuo[eb]);
                    #endif

                }               

            }
            
           /*  for (int i = 0; i < EB; i++)
            {
                printf("%f ", posibleIndividuo[i]);
            }
            printf("\n"); */
            
            /* EVALUACION */
                evaluate_zdt3(posibleIndividuo, evaluationP, pReferenceLocal);             

                #ifdef degug
                fprintf(fprueba, "Evaluacion %f %f ,pReferenciaGlobal %f %f ,pReferenciaLocal %f %f \n", evaluation[0], evaluation[1], pReferenceGlobal[0], pReferenceGlobal[1], pReferenceLocal[0], pReferenceLocal[1]);
                #endif

            /* ACTUALIZACION_PUNTO_REFERENCIA */
                
                if(pReferenceGlobal[0]>pReferenceLocal[0]) {
                    pReferenceGlobal[0]=pReferenceLocal[0];
                }
                if(pReferenceGlobal[1]>pReferenceLocal[1]) {
                    pReferenceGlobal[1]=pReferenceLocal[1];
                }
            
            /* ACTUALIZACION_VECINOS */

                #ifdef degug
                fprintf(fprueba, "Actualizacion vecinos: sub: %d \n", subproblema);    
                #endif

                //Valor del individuo actual
                float valorActual = 0.;
                float valorVecino = 0.;
                int vecino = 0;
                float vecinoActual[EB];
                int mejor = 0;
                //Calculamos valores para los vecinos del individuo
                //printf("Subproblema %d \n", subproblema);
                for (int t = 0; t < T; t++)
                {
                    vecino = neightbours[subproblema][t];

                    valorActual = evaluacionGTE(posibleIndividuo, pesos[vecino], pReferenceGlobal);

                    for (int p = 0; p < EB; p++)
                    {
                        vecinoActual[p] = population[vecino][p];
                    }

                    valorVecino = evaluacionGTE(vecinoActual, pesos[vecino], pReferenceGlobal);
                    
                    /* printf("A: ");
                    for (int p = 0; p < EB; p++)
                    {
                        printf("%f ", population[vecino][p]);
                    }
                    printf("\n"); */

                    if (valorActual <= valorVecino){
                        //Sustituimos al vecino por la mejor solucion hasta el momento
                        //printf("B: ");
                        mejor = 1;

                        for (int eb = 0; eb < EB; eb++)
                        {
                            population[vecino][eb] = posibleIndividuo[eb];
                            //printf("%f ", population[vecino][eb]);
                        }
                        //printf("\n");
                    }


                    #ifdef degug
                    fprintf(fprueba, "valor actual %f, valor vecino: %f, mejor = %d \n", valorActual, valorVecino, mejor);
                    #endif

                    //printf("Subproblema %d - Vecino: %d - valorVecino: %f - ValorActual: %f \n", subproblema, vecino, valorVecino, valorActual);
                }

                if(mejor == 1){
                    evaluation[0] = evaluationP[0];
                    evaluation[1] = evaluationP[1];

                    for (int i = 0; i < 2; i++)
                    {
                        evaluationN[subproblema][i] = evaluation[i];
                    }
                    
                }
                
                #ifdef degug
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < EB; j++)
                    {
                        fprintf(fprueba, "%f ", population[i][j]);
                    }
                    fprintf(fprueba, "\n");
                }
                fprintf(fprueba, "\n");
                #endif
                //printf("\n");

            fprintf(fsoluciones, "%f %f %f \n", evaluation[0], evaluation[1], 0.0);
            fprintf(fpRefencia, "%f %f \n", pReferenceGlobal[0], pReferenceGlobal[1]);


            


        }

        float value1 = 0.;
        
        float value2 = 0.;
        int flag = 0;
        for (int i = 0; i < N; i++)
        {
            value1 = evaluationN[i][0];
            value2 = evaluationN[i][1];

            for(int n = 1; n < N; n++){
                if((value1 < evaluationN[n][0]) && (value2 < evaluationN[n][1])){
                    flag = 1;
                }else if((value1 < evaluationN[n][0]) && (value2 > evaluationN[n][1])){
                    flag = 2;
                }else if((value1 > evaluationN[n][0]) && (value2 < evaluationN[n][1])){
                    flag = 3;
                }else{
                    flag = 0;
                }
            }
            if(flag == 1){
                evaluationN[i][0]= value1;
                evaluationN[i][1]= value2;
            }
        }
            
            
        
    }
    for (int i = 0; i < N; i++)
    {
        fprintf(fUsoluciones, "%f %f %f \n", evaluationN[i][0], evaluationN[i][1], 0.0);
    }
    
}

/*  MAIN    */
int main(){

    //Ficheros
    char * commandsForGnuplot[] = {"set title \"Frente\"", "set yrange [-0.8:1]", "plot 'solucionesU.temp' linetype 4 lw 2, 'PF.dat' linetype 1 lw 0.5"};
    
    Fichero gnuplotPipe = popen ("gnuplot -persistent", "w");
    Fichero gnuplotPipe2 = popen ("gnuplot -persistent", "w");
    Fichero fsoluciones = fopen("soluciones.temp", "w");
    Fichero fUsoluciones = fopen("solucionesU.temp", "w");
    Fichero findividuos = fopen("individuos.temp", "w");
    Fichero fvpesos = fopen("vpesos.temp", "w");
    Fichero fpRefencia = fopen("puntosReferencia.temp", "w");

    Fichero fpruebas = fopen("pruebas.temp", "w");

    // Vectores alpha
    float alpha_vector[N][2];

    // Vectores de los vecinos
    int neightbours[N][T];

    // Vector de N individuos
    float population[N][EB];

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
    
    evaluate_zdt3_all(population, evaluation, pReference);

    for (int p = 0; p < N; p++)
    {
        //Almaceno las primeras soluciones del algorirmo
        fprintf(fsoluciones, "%f %f %f \n", evaluation[p][0], evaluation[p][1], 0.0);
        fprintf(fvpesos, "%f %f \n", alpha_vector[p][0], alpha_vector[p][1]);
    }  

    /* OPERACIONES */
    iteraciones(population, neightbours, pReference, alpha_vector, fsoluciones, findividuos, fpRefencia, fUsoluciones, fpruebas);

    //fprintf(temp, "%1.3f %1.3f \n", pReference[0], pReference[1]);

    for (int i=0; i < 3; i++)
    {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
    for (int i=0; i < 3; i++)
    {
        fprintf(gnuplotPipe2, "%s \n", commandsForGnuplot2[i]); //Send commands to gnuplot one by one.
    }

    fclose(fsoluciones);
    fclose(findividuos);
    fclose(fvpesos);
    fclose(fpRefencia);
    fclose(fUsoluciones);

    return 0;
}


