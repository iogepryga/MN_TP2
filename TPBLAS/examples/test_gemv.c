#include <stdio.h>
#include <x86intrin.h>

#ifndef MNBLAS_H
#define MNBLAS_H
#include "mnblas.h"
#endif

#ifndef COMPEXE_H
#define COMPEXE_H
#include "complexe.h"
#endif

#ifndef TESTUILS_H
#define TESTUILS_H
#include "testutils.h"
#endif

#define M_RESULTAT 2
#define N_RESULTAT 4
#define RAND_MAXIMUM 10

int main (int argc, char **argv) {
    int MAX_SIZE = (M_RESULTAT < N_RESULTAT ? N_RESULTAT : M_RESULTAT);
    float* V1s = (float*)malloc(MAX_SIZE*sizeof(float)); float* V2s = (float*)malloc(MAX_SIZE*sizeof(float)); float* V2s_copy = (float*)malloc(MAX_SIZE*sizeof(float));
    double* V1d = (double*)malloc(MAX_SIZE*sizeof(double)); double* V2d = (double*)malloc(MAX_SIZE*sizeof(double)); double* V2d_copy = (double*)malloc(MAX_SIZE*sizeof(double));
    complexe_float_t* V1c = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t)); complexe_float_t* V2c = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t)); complexe_float_t* V2c_copy = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t));
    complexe_double_t* V1z = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t)); complexe_double_t* V2z = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t)); complexe_double_t* V2z_copy = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t));
    
    float* Ms = (float*)malloc(M_RESULTAT*N_RESULTAT*sizeof(float));
    double* Md = (double*)malloc(M_RESULTAT*N_RESULTAT*sizeof(double));
    complexe_float_t* Mc = (complexe_float_t*)malloc(M_RESULTAT*N_RESULTAT*sizeof(complexe_float_t));
    complexe_double_t* Mz = (complexe_double_t*)malloc(M_RESULTAT*N_RESULTAT*sizeof(complexe_double_t));
    
    
    
    complexe_float_t tmpc; tmpc.real = 2; tmpc.imaginary = 0;
    complexe_double_t tmpz; tmpz.real = 2; tmpz.imaginary = 0;

    


    printf("                   TEST_GEMV\n|||||||||||||||||||||||||||||||||||||||||||||||||||||\n                           1: TEST DE BON RESULTAT\n<------------------------------------------>\n                   float\n");
    printf("-------------------------- test 1 (V2s=2*Ms*V1s+2*V2s) (MNCblasRowMajor): \n");
    void_matrix_sinit(Ms,1,M_RESULTAT,N_RESULTAT);
    void_vector_sinit(V1s, 1,N_RESULTAT);
    void_vector_sinit(V2s, 2,M_RESULTAT);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,N_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    printf("-------------------------- test 2 (V2s=2*Ms*V1s+2*V2s) (MNCblasColMajor): \n");
    void_vector_sinit(V2s, 2,M_RESULTAT);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,N_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    mncblas_sgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);


    printf("-------------------------- test 3 (V2s=2*Ms*V1s+2*V2s) (MNCblasRowMajor): \n");
    void_matrix_sinit_rand(Ms,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,N_RESULTAT);
    void_vector_sinit_rand(V2s, RAND_MAXIMUM,M_RESULTAT);
    mncblas_scopy(M_RESULTAT,V2s,1,V2s_copy,1);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,N_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,M_RESULTAT);
    printf("-------------------------- test 4 (V2s=2*Ms*V1s+2*V2s) (MNCblasColMajor): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_FLOAT,Ms,M_RESULTAT,N_RESULTAT);
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,N_RESULTAT);
    printf("V2s : "); vector_print(V2s_copy,TYPE_FLOAT,M_RESULTAT);
    mncblas_sgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s_copy,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s_copy,TYPE_FLOAT,M_RESULTAT);

    // Trans :

    printf("-------------------------- test 5 (V2s=2*Ms*V1s+2*V2s) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_sinit(Ms,1,M_RESULTAT,N_RESULTAT);
    void_vector_sinit(V1s, 1,M_RESULTAT);
    void_vector_sinit(V2s, 2,N_RESULTAT);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,M_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    mncblas_sgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    printf("-------------------------- test 6 (V2s=2*Ms*V1s+2*V2s) (MNCblasColMajor,MNCblasTrans): \n");
    void_vector_sinit(V2s, 2,N_RESULTAT);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,M_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    mncblas_sgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);


    printf("-------------------------- test 7 (V2s=2*Ms*V1s+2*V2s) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_sinit_rand(Ms,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,M_RESULTAT);
    void_vector_sinit_rand(V2s, RAND_MAXIMUM,N_RESULTAT);
    mncblas_scopy(N_RESULTAT,V2s,1,V2s_copy,1);
    printf("---- Avant :\n");
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,M_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    mncblas_sgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,N_RESULTAT);
    printf("-------------------------- test 8 (V2s=2*Ms*V1s+2*V2s) (MNCblasColMajor,MNCblasTrans): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_FLOAT,Ms,M_RESULTAT,N_RESULTAT);
    printf("Ms : "); matrix_print(Ms,TYPE_FLOAT,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,M_RESULTAT);
    printf("V2s : "); vector_print(V2s_copy,TYPE_FLOAT,N_RESULTAT);
    mncblas_sgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Ms,1,V1s,1,2,V2s_copy,1);
    printf("---- Apres :\n");
    printf("V2s : "); vector_print(V2s_copy,TYPE_FLOAT,N_RESULTAT);

    printf("<------------------------------------------>\n                   double\n");

    printf("-------------------------- test 1 (V2d=2*Md*V1d+2*V2d) (MNCblasRowMajor): \n");
    void_matrix_dinit(Md,1,M_RESULTAT,N_RESULTAT);
    void_vector_dinit(V1d, 1,N_RESULTAT);
    void_vector_dinit(V2d, 2,M_RESULTAT);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,N_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    mncblas_dgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    printf("-------------------------- test 2 (V2d=2*Md*V1d+2*V2d) (MNCblasColMajor): \n");
    void_vector_dinit(V2d, 2,M_RESULTAT);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,N_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    mncblas_dgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);


    printf("-------------------------- test 3 (V2d=2*Md*V1d+2*V2d) (MNCblasRowMajor): \n");
    void_matrix_dinit_rand(Md,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,N_RESULTAT);
    void_vector_dinit_rand(V2d, RAND_MAXIMUM,M_RESULTAT);
    mncblas_dcopy(M_RESULTAT,V2d,1,V2d_copy,1);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,N_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    mncblas_dgemv(MNCblasRowMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,M_RESULTAT);
    printf("-------------------------- test 4 (V2d=2*Md*V1d+2*V2d) (MNCblasColMajor): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_DOUBLE,Md,M_RESULTAT,N_RESULTAT);
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,N_RESULTAT);
    printf("V2d : "); vector_print(V2d_copy,TYPE_DOUBLE,M_RESULTAT);
    mncblas_dgemv(MNCblasColMajor,MNCblasNoTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d_copy,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d_copy,TYPE_DOUBLE,M_RESULTAT);

    // Trans :

    printf("-------------------------- test 5 (V2d=2*Md*V1d+2*V2d) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_dinit(Md,1,M_RESULTAT,N_RESULTAT);
    void_vector_dinit(V1d, 1,M_RESULTAT);
    void_vector_dinit(V2d, 2,N_RESULTAT);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,M_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    mncblas_dgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    printf("-------------------------- test 6 (V2d=2*Md*V1d+2*V2d) (MNCblasColMajor,MNCblasTrans): \n");
    void_vector_dinit(V2d, 2,N_RESULTAT);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,M_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    mncblas_dgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);


    printf("-------------------------- test 7 (V2d=2*Md*V1d+2*V2d) (MNCblasRowMajor,MNCblasTrans): \n");
    void_matrix_dinit_rand(Md,RAND_MAXIMUM,M_RESULTAT,N_RESULTAT);
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,M_RESULTAT);
    void_vector_dinit_rand(V2d, RAND_MAXIMUM,N_RESULTAT);
    mncblas_dcopy(N_RESULTAT,V2d,1,V2d_copy,1);
    printf("---- Avant :\n");
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasRowMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,M_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    mncblas_dgemv(MNCblasRowMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,N_RESULTAT);
    printf("-------------------------- test 8 (V2d=2*Md*V1d+2*V2d) (MNCblasColMajor,MNCblasTrans): \n");
    printf("---- Avant :\n");
    row_to_col_major(TYPE_DOUBLE,Md,M_RESULTAT,N_RESULTAT);
    printf("Md : "); matrix_print(Md,TYPE_DOUBLE,MNCblasColMajor,M_RESULTAT,N_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,M_RESULTAT);
    printf("V2d : "); vector_print(V2d_copy,TYPE_DOUBLE,N_RESULTAT);
    mncblas_dgemv(MNCblasColMajor,MNCblasTrans,M_RESULTAT,N_RESULTAT,2,Md,1,V1d,1,2,V2d_copy,1);
    printf("---- Apres :\n");
    printf("V2d : "); vector_print(V2d_copy,TYPE_DOUBLE,N_RESULTAT);

/*

    printf("<------------------------------------------>\n                   complexe_float_t\n");

    printf("---------- test 1 (V2d=2*V1d+V2d): \n");
    void_vector_cinit(V1c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    void_vector_cinit2(V2c,2,0,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mnblas_caxpy(VECSIZE_RESULTAT,&tmpc,V1c,1,V2c,1);
    printf("---- Apres :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("---------- test 2 (V2d=2*V1d+V2d): \n");
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_cinit_rand(V2c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    V1c[VECSIZE_RESULTAT-1] = gen_complexe_float(0,0);
    V2c[VECSIZE_RESULTAT-1] = gen_complexe_float(1,0);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mnblas_caxpy(VECSIZE_RESULTAT,&tmpc,V1c,1,V2c,1);
    printf("---- Apres :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);

    printf("<------------------------------------------>\n                   complexe_double_t\n");

    printf("---------- test 1 (V2d=2*V1d+V2d): \n");
    void_vector_zinit(V1z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    void_vector_zinit2(V2z,2,0,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mnblas_zaxpy(VECSIZE_RESULTAT,&tmpz,V1z,1,V2z,1);
    printf("---- Apres :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("---------- test 2 (V2d=2*V1d+V2d): \n");
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_zinit_rand(V2z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    V1z[VECSIZE_RESULTAT-1] = gen_complexe_double(0,0);
    V2z[VECSIZE_RESULTAT-1] = gen_complexe_double(1,0);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mnblas_zaxpy(VECSIZE_RESULTAT,&tmpz,V1z,1,V2z,1);
    printf("---- Apres :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);

*/
    free_vm(V1s);
    free_vm(V1d);
    free_vm(V1c);
    free_vm(V1z);
    free_vm(V2s);
    free_vm(V2s_copy);
    free_vm(V2d);
    free_vm(V2d_copy);
    free_vm(V2c);
    free_vm(V2c_copy);
    free_vm(V2z);
    free_vm(V2z_copy);

    #define M_FLOPS    10000
    #define N_FLOPS    10000
    #define NB_EXPE     10
    #define NB_OPE_REEL 1 // ??
    #define NB_OPE_COMPLEXE 2 // ??

    MAX_SIZE = (M_RESULTAT < N_RESULTAT ? N_RESULTAT : M_RESULTAT);
    V1s = (float*)malloc(MAX_SIZE*sizeof(float)), V2s = (float*)malloc(MAX_SIZE*sizeof(float));
    V1d = (double*)malloc(MAX_SIZE*sizeof(double)), V2d = (double*)malloc(MAX_SIZE*sizeof(double));
    V1c = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t)), V2c = (complexe_float_t*)malloc(MAX_SIZE*sizeof(complexe_float_t));
    V1z = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t)), V2z = (complexe_double_t*)malloc(MAX_SIZE*sizeof(complexe_double_t));

    Ms = (float*)malloc(M_FLOPS*N_FLOPS*sizeof(float));
    Md = (double*)malloc(M_FLOPS*N_FLOPS*sizeof(double));
    Mc = (complexe_float_t*)malloc(M_FLOPS*N_FLOPS*sizeof(complexe_float_t));
    Mz = (complexe_double_t*)malloc(M_FLOPS*N_FLOPS*sizeof(complexe_double_t));


    printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n                  2 : FLOPS\n <-------------------------------------------------->\n                     float\n");
    unsigned long long int start, end ; 
    for(int i = 0; i < NB_EXPE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Ms,1,V1s,1,2,V2s,1);
        end = _rdtsc();
        calcul_flop("mncblas_sgemv : ", NB_OPE_REEL*M_FLOPS*N_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      float sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Ms,1,V1s,1,2,V2s,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_sgemv : ", NB_EXPE*NB_OPE_REEL*M_FLOPS*N_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      double\n");
    for(int i = 0; i < NB_EXPE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_dgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Md,1,V1d,1,2,V2d,1);
        end = _rdtsc();
        calcul_flop("mncblas_dgemv : ", NB_OPE_REEL*M_FLOPS*N_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      double sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_dgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Md,1,V1d,1,2,V2d,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_dgemv : ", NB_EXPE*NB_OPE_REEL*M_FLOPS*N_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_float_t\n");
    for(int i = 0; i < NB_EXPE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_cgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Mc,1,V1c,1,2,V2c,1);
        end = _rdtsc();
        calcul_flop("mncblas_cgemv : ", NB_OPE_COMPLEXE*M_FLOPS*N_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_float_t sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_cgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Mc,1,V1c,1,2,V2c,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_cgemv : ", NB_EXPE*NB_OPE_COMPLEXE*M_FLOPS*N_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_double_t\n");
    for(int i = 0; i < NB_EXPE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_zgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Mz,1,V1z,1,2,V2z,1);
        end = _rdtsc();
        calcul_flop("mncblas_zgemv : ", NB_OPE_COMPLEXE*M_FLOPS*N_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_double_t sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_zgemv(MNCblasRowMajor,MNCblasNoTrans,M_FLOPS,N_FLOPS,2,Mz,1,V1z,1,2,V2z,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_zgemv : ", NB_EXPE*NB_OPE_COMPLEXE*M_FLOPS*N_FLOPS ,end-start);


    free_vm(V1s);
    free_vm(V1d);
    free_vm(V1c);
    free_vm(V1z);
    free_vm(V2s);
    free_vm(V2d);
    free_vm(V2c);
    free_vm(V2z);
}