#include "general.c"
#include "aux.c"


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
    float *vecinos = malloc(T * sizeof(float));

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

void inicialicePopulation(){

}

int main(){
    alpha_vector[N][2];
    neightbours[N][T];
    inicialiceAlphaVector(alpha_vector);
    inicialiceNeighbourVector(alpha_vector, neightbours);

    return 0;
}