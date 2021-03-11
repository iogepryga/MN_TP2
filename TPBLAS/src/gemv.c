#include "mnblas.h"
#include "complexe.h"


void mncblas_sgemv(const MNCBLAS_LAYOUT layout, const MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const float alpha, const float *A, const int lda, const float *X, const int incX, const float beta, float *Y, const int incY) {
    if(layout == MNCblasRowMajor) {
        if(TransA == MNCblasNoTrans) {
            register float sum = 0;
            for(register int i = 0, j; i < M ; i+= incY) {
                for(j = 0; j < N ; j+= incX) {
                    sum += *(A+i*N+j)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register float sum = 0;
            for(register int i = 0, j; i < M ; i+= incY) {
                for(j = 0; j < N ; j+= incX) {
                    sum += *(A+j*N+i)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        }
    } else if (layout == MNCblasColMajor){
        if(TransA == MNCblasNoTrans) {
            register float sum = 0;
            for(register int i = 0, j; i < M ; i+= incY) {
                for(j = 0; j < N ; j+= incX) {
                    sum += *(A+j*M+i)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        } else if (TransA == MNCblasTrans) {
            register float sum = 0;
            for(register int i = 0, j; i < M ; i+= incY) {
                for(j = 0; j < N ; j+= incX) {
                    sum += *(A+i*M+j)*X[j];
                }
                Y[i] = alpha*sum + beta*Y[i];
            }
        }
    }

}

void mncblas_dgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const double alpha, const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY) {
    // A completer
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const void *alpha, const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY) {
    // A completer
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                     const int M, const int N,
                     const void *alpha, const void *A, const int lda, const void *X, const int incX, const void *beta, void *Y, const int incY){
    // A completer
}