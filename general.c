// Poblacion
#define N 200

// Generaciones
#define G 200

// Vecindad (%) 15% de la poblacion total (200*0.15)
#define T 30

// Espacio de busqueda
int const Xinferior = 0;
int const Xsuperior = 1;
int const p = 30;

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
