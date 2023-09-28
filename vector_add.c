#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void Allocate_vectors(double** x_pp, double** y_pp, double** z_pp, int n);
void Generate_random_vector(double a[], int n);
void Print_first_last_elements(double b[], int n, char vec_name[]);
void Vector_sum(double x[], double y[], double z[], int n);

/*---------------------------------------------------------------------*/
int main(void) {
    int n = 100000; // Tamaño de los vectores
    double *x, *y, *z;
    clock_t start_time, end_time; // Variables para medir el tiempo

    Allocate_vectors(&x, &y, &z, n);

    Generate_random_vector(x, n);
    Generate_random_vector(y, n);

    printf("Vector x:\n");
    Print_first_last_elements(x, n, "x");

    printf("Vector y:\n");
    Print_first_last_elements(y, n, "y");

    start_time = clock(); // Inicio de la medición de tiempo
    Vector_sum(x, y, z, n);
    end_time = clock(); // Fin de la medición de tiempo

    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC; // Tiempo total en segundos

    printf("The sum (vector z):\n");
    Print_first_last_elements(z, n, "z");

    printf("Tiempo de ejecución: %.6f segundos\n", total_time); // Imprimir el tiempo de ejecución

    free(x);
    free(y);
    free(z);

    return 0;
}

/*---------------------------------------------------------------------
 * Function:  Allocate_vectors
 * Purpose:   Allocate storage for the vectors
 * In arg:    n:  the order of the vectors
 * Out args:  x_pp, y_pp, z_pp:  pointers to storage for the vectors
 *
 * Errors:    If one of the mallocs fails, the program terminates
 */
void Allocate_vectors(
      double**  x_pp  /* out */, 
      double**  y_pp  /* out */, 
      double**  z_pp  /* out */, 
      int       n     /* in  */) {
   *x_pp = malloc(n*sizeof(double));
   *y_pp = malloc(n*sizeof(double));
   *z_pp = malloc(n*sizeof(double));
   if (*x_pp == NULL || *y_pp == NULL || *z_pp == NULL) {
      fprintf(stderr, "Can't allocate vectors\n");
      exit(-1);
   }
}  /* Allocate_vectors */

/*---------------------------------------------------------------------
 * Function:  Generate_random_vector
 * Purpose:   Generate random values for a vector
 * In args:   n:  order of the vector
 * Out arg:   a:  the vector with random values
 */
void Generate_random_vector(
      double  a[]  /* out */, 
      int     n    /* in  */) {
   int i;

   srand(time(NULL)); // Seed the random number generator
   for (i = 0; i < n; i++) {
      a[i] = ((double)rand() / RAND_MAX) * 100.0; // Generate random values between 0 and 100
   }
}  /* Generate_random_vector */

/*---------------------------------------------------------------------
 * Function:  Print_first_last_elements
 * Purpose:   Print the first and last 10 elements of a vector
 * In args:   b:  the vector to be printed
 *            n:  the order of the vector
 *            vec_name:  name of the vector
 */
void Print_first_last_elements(
      double  b[]        /* in */, 
      int     n          /* in */, 
      char    vec_name[] /* in */) {
   int i;
   printf("%s (First and Last 10 elements):\n", vec_name);
   for (i = 0; i < 10; i++) {
      printf("%f ", b[i]);
   }
   printf("... ");
   for (i = n - 10; i < n; i++) {
      printf("%f ", b[i]);
   }
   printf("\n");
}  /* Print_first_last_elements */

/*---------------------------------------------------------------------
 * Function:  Vector_sum
 * Purpose:   Add two vectors
 * In args:   x:  the first vector to be added
 *            y:  the second vector to be added
 *            n:  the order of the vectors
 * Out arg:   z:  the sum vector
 */
void Vector_sum(
      double  x[]  /* in  */, 
      double  y[]  /* in  */, 
      double  z[]  /* out */, 
      int     n    /* in  */) {
   int i;

   for (i = 0; i < n; i++)
      z[i] = x[i] + y[i];
}  /* Vector_sum */
