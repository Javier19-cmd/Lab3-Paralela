/**
 * @file inciso1.c
 * @author Javier Valle, Mario de Leòn
 * @brief Programa que calcula el producto punto de dos vectores y el producto de un vector con un escalar.
 * @class inciso1
 * @c 20159, 19029
 * @version 0.1
 * @date 2023-09-28
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Prototipos de funciones.

void Check_for_error(int local_ok, char fname[], char message[], MPI_Comm comm);
void Read_n(int* n_p, int* local_n_p, int my_rank, int comm_sz, MPI_Comm comm);
void Allocate_vectors(double** local_x_pp, double** local_y_pp, double** local_z_pp, int local_n, MPI_Comm comm);
void Read_vector(double local_a[], int local_n, int n, char vec_name[], int my_rank, MPI_Comm comm);
void Print_vector(double local_b[], int local_n, int n, char title[], int my_rank, MPI_Comm comm);
void Parallel_vector_sum(double local_x[], double local_y[], double local_z[], int local_n);
double Scalar_product(double local_x[], double local_y[], int local_n);
void Scalar_vector_product(double local_x[], double scalar, double local_result[], int local_n);

int main(void) {

    // Variables a usar.

    int n, local_n;
    int comm_sz, my_rank;
    double *local_x, *local_y, *local_z;
    double scalar;
    MPI_Comm comm;

    // Inicializa MPI.

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    // Lee el orden de los vectores y lo distribuye entre los procesos.

    Read_n(&n, &local_n, my_rank, comm_sz, comm);

    // Asigna memoria para los vectores locales.

    Allocate_vectors(&local_x, &local_y, &local_z, local_n, comm);

    // Lee los vectores x & y.

    Read_vector(local_x, local_n, n, "x", my_rank, comm);
    Read_vector(local_y, local_n, n, "y", my_rank, comm);

    // Broadcast del escalar.
    if (my_rank == 0) {
        printf("Ingresa el vector escalar: ");
        scanf("%lf", &scalar);
    }
    MPI_Bcast(&scalar, 1, MPI_DOUBLE, 0, comm);

    // Haciendo los cálculos en paralelo.
    Parallel_vector_sum(local_x, local_y, local_z, local_n);
    double dot_product = Scalar_product(local_x, local_y, local_n);
    double scalar_result[local_n];
    Scalar_vector_product(local_x, scalar, scalar_result, local_n);

    // Imprimiendo los resultados.
    if (my_rank == 0) {
        printf("Producto punto de x & y: %lf\n", dot_product);
        printf("Vector escalar del producto de x y un escalar: ");
        for (int i = 0; i < n; i++) {
            printf("%lf ", scalar_result[i]);
        }
        printf("\n");
    }

    free(local_x);
    free(local_y);
    free(local_z);

    MPI_Finalize();

    return 0;
}

// Imprime un mensaje de error y termina el programa.
void Check_for_error(int local_ok, char fname[], char message[], MPI_Comm comm) {
    int ok;

    MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
    if (ok == 0) {
        int my_rank;
        MPI_Comm_rank(comm, &my_rank);
        if (my_rank == 0) {
            fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname, message);
            fflush(stderr);
        }
        MPI_Finalize();
        exit(-1);
    }
}

// Lee el orden de los vectores y lo distribuye entre los procesos.
void Read_n(int* n_p, int* local_n_p, int my_rank, int comm_sz, MPI_Comm comm) {
    int local_ok = 1;
    char *fname = "Read_n";

    if (my_rank == 0) {
        printf("¿Cuàl es el orden de los vectores?\n");
        scanf("%d", n_p);
    }
    MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
    if (*n_p <= 0 || *n_p % comm_sz != 0) local_ok = 0;
    Check_for_error(local_ok, fname, "n deberìa ser > 0 y divisible entre comm_sz", comm);
    *local_n_p = *n_p / comm_sz;
}

// Asigna memoria para los vectores locales.
void Allocate_vectors(double** local_x_pp, double** local_y_pp, double** local_z_pp, int local_n, MPI_Comm comm) {
    int local_ok = 1;
    char* fname = "Allocate_vectors";

    *local_x_pp = malloc(local_n * sizeof(double));
    *local_y_pp = malloc(local_n * sizeof(double));
    *local_z_pp = malloc(local_n * sizeof(double));

    if (*local_x_pp == NULL || *local_y_pp == NULL || *local_z_pp == NULL) local_ok = 0;
    Check_for_error(local_ok, fname, "No se pueden asignar vector(es) locales", comm);
}

// Lee un vector de enteros de la entrada estándar.
void Read_vector(double local_a[], int local_n, int n, char vec_name[], int my_rank, MPI_Comm comm) {
    double* a = NULL;
    int i;
    int local_ok = 1;
    char* fname = "Read_vector";

    if (my_rank == 0) {
        a = malloc(n * sizeof(double));
        if (a == NULL) local_ok = 0;
        Check_for_error(local_ok, fname, "No se pueden asginar vectores temporales", comm);
        printf("Enter the vector %s\n", vec_name);
        for (i = 0; i < n; i++)
            scanf("%lf", &a[i]);
        MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, 0, comm);
        free(a);
    } else {
        Check_for_error(local_ok, fname, "No se pueden asginar vectores temporales", comm);
        MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, 0, comm);
    }
}

// Imprime un vector distribuido a través de todos los procesos.
void Print_vector(double local_b[], int local_n, int n, char title[], int my_rank, MPI_Comm comm) {
    double* b = NULL;
    int i;
    int local_ok = 1;
    char* fname = "Print_vector";

    if (my_rank == 0) {
        b = malloc(n * sizeof(double));
        if (b == NULL) local_ok = 0;
        Check_for_error(local_ok, fname, "No se pueden asginar vectores temporales", comm);
        MPI_Gather(local_b, local_n, MPI_DOUBLE, b, local_n, MPI_DOUBLE, 0, comm);
        printf("%s\n", title);
        for (i = 0; i < n; i++)
            printf("%f ", b[i]);
        printf("\n");
        free(b);
    } else {
        Check_for_error(local_ok, fname, "No se pueden asginar vectores temporales", comm);
        MPI_Gather(local_b, local_n, MPI_DOUBLE, b, local_n, MPI_DOUBLE, 0, comm);
    }
}

// Suma paralela de los vectores locales.
void Parallel_vector_sum(double local_x[], double local_y[], double local_z[], int local_n) {
    for (int local_i = 0; local_i < local_n; local_i++)
        local_z[local_i] = local_x[local_i] + local_y[local_i];
}

// Producto escalar de los vectores.
double Scalar_product(double local_x[], double local_y[], int local_n) {
    double local_dot_product = 0.0;
    for (int i = 0; i < local_n; i++) {
        local_dot_product += local_x[i] * local_y[i];
    }

    double global_dot_product;
    MPI_Allreduce(&local_dot_product, &global_dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_dot_product;
}

// Producto escalar de un vector con un escalar.
void Scalar_vector_product(double local_x[], double scalar, double local_result[], int local_n) {
    for (int i = 0; i < local_n; i++) {
        local_result[i] = scalar * local_x[i];
    }
}
