#define VECSIZE     65536
#define NB_EXPE     10

#ifndef MNBLAS_H
#define MNBLAS_H
#include "mnblas.h"
#endif

#ifndef COMPEXE_H
#define COMPEXE_H
#include "complexe.h"
#endif

typedef enum {TYPE_FLOAT, TYPE_DOUBLE,TYPE_COMPLEXE_FLOAT,TYPE_COMPLEXE_DOUBLE} VTYPE;

void init_flop () ;
void calcul_flop (char *message, int nb_operations_flottantes, unsigned long long int cycles) ;
void free_vm (void* address);

complexe_float_t gen_complexe_float (const float real, const float imaginary);
complexe_double_t gen_complexe_double (const double real, const double imaginary);




/* ||||||||||||||||| VECTOR ||||||||||||||||||| */

/*
    VECTOR_INIT
*/

void void_vector_sinit (float* V, const float x, const register unsigned int len);

void void_vector_dinit (double* V, const double x, const register unsigned int len);

void void_vector_cinit (complexe_float_t* V, const complexe_float_t x, const register unsigned int len);

void void_vector_zinit (complexe_double_t* V, const complexe_double_t x, const register unsigned int len);

float* vector_sinit (const float x, const register unsigned int len);

double* vector_dinit (const double x, const register unsigned int len);

complexe_float_t* vector_cinit (const complexe_float_t x, const register unsigned int len);

complexe_double_t* vector_zinit (const complexe_double_t x, const register unsigned int len);

complexe_float_t* vector_cinit2 (const float real, const float imaginary, const register unsigned int len);

complexe_double_t* vector_zinit2 (const double real, const double imaginary, const register unsigned int len);

/*
    END VECTOR_INIT
*/


/*
    VECTOR_PRINT
*/

void vector_print (void* V, VTYPE type, const register unsigned int len);

/*
    END VECTOR_PRINT
*/



/* ||||||||||||||||| MATRIX ||||||||||||||||||| */

/*
    MATRIX_INIT
*/

void void_matrix_sinit (float* V, const MNCBLAS_LAYOUT layout, const float x, const int M, const int N);

void void_matrix_dinit (double* V, const MNCBLAS_LAYOUT layout, const double x, const int M, const int N);

void void_matrix_cinit (complexe_float_t* V, const MNCBLAS_LAYOUT layout, const complexe_float_t x, const int M, const int N);

void void_matrix_zinit (complexe_double_t* V, const MNCBLAS_LAYOUT layout, const complexe_double_t x, const int M, const int N);

void void_matrix_cinit2 (complexe_float_t* V, const MNCBLAS_LAYOUT layout, const float real, const float imaginary, const int M, const int N);

void void_matrix_zinit2 (complexe_double_t* V, const MNCBLAS_LAYOUT layout, const double real, const double imaginary, const int M, const int N);

float* matrix_sinit (const MNCBLAS_LAYOUT layout, const float x, const int M, const int N);

double* matrix_dinit (const MNCBLAS_LAYOUT layout, const double x, const int M, const int N);

complexe_float_t* matrix_cinit (const MNCBLAS_LAYOUT layout, const complexe_float_t x, const int M, const int N);

complexe_double_t* matrix_zinit (const MNCBLAS_LAYOUT layout, const complexe_double_t x, const int M, const int N);

complexe_float_t* matrix_cinit2 (const MNCBLAS_LAYOUT layout, const float real, const float imaginary, const int M, const int N);

complexe_double_t* matrix_zinit2 (const MNCBLAS_LAYOUT layout, const double real, const double imaginary, const int M, const int N);

/*
    END MATRIX_INIT
*/


/*
    MATRIX_PRINT
*/

void matrix_print (const void* V, const VTYPE type, const MNCBLAS_LAYOUT layout ,const int M ,const int N);

/*
    END MATRIX_PRINT
*/
