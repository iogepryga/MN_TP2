#include "mnblas.h"
#include "complexe.h"

float  mnblas_snrm2(const int N, const float *X, const int incX) {
    register float sum = 0;
    for (register unsigned int i = 0; i < N ; i += incX) {
        sum += X[i] * X[i];
    }
    return sqrtf(sum);
}

double mnblas_dnrm2(const int N, const double *X, const int incX) {
    register double sum = 0;
    for (register unsigned int i = 0; i < N ; i += incX) {
        sum += X[i] * X[i];
    }
    return sqrt(sum);
}

float  mnblas_scnrm2(const int N, const void *X, const int incX) {
    register float sum = 0;
    for (register unsigned int i = 0; i < N ; i += incX) {
        sum += ((complexe_float_t*)X+i)->real * ((complexe_float_t*)X+i)->real + ((complexe_float_t*)X+i)->imaginary * ((complexe_float_t*)X+i)->imaginary;
    }
    return sqrtf(sum);
}

double mnblas_dznrm2(const int N, const void *X, const int incX) {
    register double sum = 0;
    for (register unsigned int i = 0; i < N ; i += incX) {
        sum += ((complexe_double_t*)X+i)->real * ((complexe_double_t*)X+i)->real + ((complexe_double_t*)X+i)->imaginary * ((complexe_double_t*)X+i)->imaginary;
    }
    return sqrt(sum);
}