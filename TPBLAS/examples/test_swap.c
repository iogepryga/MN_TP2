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

#define VECSIZE_RESULTAT     20
#define RAND_MAXIMUM 10

int main (int argc, char **argv) {
    float* V1s = (float*)malloc(VECSIZE_RESULTAT*sizeof(float)); float* V2s = (float*)malloc(VECSIZE_RESULTAT*sizeof(float));
    double* V1d = (double*)malloc(VECSIZE_RESULTAT*sizeof(double)); double* V2d = (double*)malloc(VECSIZE_RESULTAT*sizeof(double));
    complexe_float_t* V1c = (complexe_float_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_float_t)); complexe_float_t* V2c = (complexe_float_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_float_t));
    complexe_double_t* V1z = (complexe_double_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_double_t)); complexe_double_t* V2z = (complexe_double_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_double_t));

    


    printf("                   TEST_SWAP\n|||||||||||||||||||||||||||||||||||||||||||||||||||||\n                           1: TEST DE BON RESULTAT\n<------------------------------------------>\n                   float\n");
    printf("---------- test 1 : \n");
    void_vector_sinit(V1s, 1,VECSIZE_RESULTAT);
    void_vector_sinit(V2s, 2,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    mncblas_sswap(VECSIZE_RESULTAT,V1s,1,V2s,1);
    printf("---- Apres :\n");
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("---------- test 2 : \n");
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_sinit_rand(V2s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    V1s[VECSIZE_RESULTAT-1] = 0;
    V2s[VECSIZE_RESULTAT-1] = 1;
    printf("---- Avant :\n");
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);
    mncblas_sswap(VECSIZE_RESULTAT,V1s,1,V2s,1);
    printf("---- Apres :\n");
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("V2s : "); vector_print(V2s,TYPE_FLOAT,VECSIZE_RESULTAT);

    printf("<------------------------------------------>\n                   double\n");

    printf("---------- test 1 : \n");
    void_vector_dinit(V1d, 1,VECSIZE_RESULTAT);
    void_vector_dinit(V2d, 2,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_dswap(VECSIZE_RESULTAT,V1d,1,V2d,1);
    printf("---- Apres :\n");
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("---------- test 2 : \n");
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_dinit_rand(V2d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    V1d[VECSIZE_RESULTAT-1] = 0;
    V2d[VECSIZE_RESULTAT-1] = 1;
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_dswap(VECSIZE_RESULTAT,V1d,1,V2d,1);
    printf("---- Apres :\n");
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2d : "); vector_print(V2d,TYPE_DOUBLE,VECSIZE_RESULTAT);

    printf("<------------------------------------------>\n                   complexe_float_t\n");

    printf("---------- test 1 : \n");
    void_vector_cinit(V1c, gen_complexe_float(1,0),VECSIZE_RESULTAT);
    void_vector_cinit2(V2c,2,0,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cswap(VECSIZE_RESULTAT,V1c,1,V2c,1);
    printf("---- Apres :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("---------- test 2 : \n");
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_cinit_rand(V2c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    V1c[VECSIZE_RESULTAT-1] = gen_complexe_float(0,0);
    V2c[VECSIZE_RESULTAT-1] = gen_complexe_float(1,0);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    mncblas_cswap(VECSIZE_RESULTAT,V1c,1,V2c,1);
    printf("---- Apres :\n");
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("V2c : "); vector_print(V2c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);

    printf("<------------------------------------------>\n                   complexe_double_t\n");

    printf("---------- test 1 : \n");
    void_vector_zinit(V1z, gen_complexe_double(1,0),VECSIZE_RESULTAT);
    void_vector_zinit2(V2z,2,0,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zswap(VECSIZE_RESULTAT,V1z,1,V2z,1);
    printf("---- Apres :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("---------- test 2 : \n");
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    void_vector_zinit_rand(V2z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("---- Avant :\n");
    V1z[VECSIZE_RESULTAT-1] = gen_complexe_double(0,0);
    V2z[VECSIZE_RESULTAT-1] = gen_complexe_double(1,0);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    mncblas_zswap(VECSIZE_RESULTAT,V1z,1,V2z,1);
    printf("---- Apres :\n");
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("V2z : "); vector_print(V2z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);




    #define VECSIZE_FLOPS    100000
    #define NB_EXPE     10
    #define NB_OPE_REEL 1 // ??
    #define NB_OPE_COMPLEXE 2 // ??

    V1s = (float*)malloc(VECSIZE_FLOPS*sizeof(float)), V2s = (float*)malloc(VECSIZE_FLOPS*sizeof(float));
    V1d = (double*)malloc(VECSIZE_FLOPS*sizeof(double)), V2d = (double*)malloc(VECSIZE_FLOPS*sizeof(double));
    V1c = (complexe_float_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_float_t)), V2c = (complexe_float_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_float_t));
    V1z = (complexe_double_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_double_t)), V2z = (complexe_double_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_double_t));


    printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n                  2 : FLOPS\n <-------------------------------------------------->\n                     float\n");
    unsigned long long int start, end ; 
    for(int i = 0; i < NB_EXPE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_sswap(VECSIZE_FLOPS,V1s,1,V2s,1);
        end = _rdtsc();
        calcul_flop("mncblas_sswap : ", NB_OPE_REEL*VECSIZE_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      float sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_sswap(VECSIZE_FLOPS,V1s,1,V2s,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_sswap : ", NB_EXPE*NB_OPE_REEL*VECSIZE_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      double\n");
    for(int i = 0; i < NB_EXPE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_dswap(VECSIZE_FLOPS,V1d,1,V2d,1);
        end = _rdtsc();
        calcul_flop("mncblas_dswap : ", NB_OPE_REEL*VECSIZE_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      double sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_dswap(VECSIZE_FLOPS,V1d,1,V2d,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_dswap : ", NB_EXPE*NB_OPE_REEL*VECSIZE_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_float_t\n");
    for(int i = 0; i < NB_EXPE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_cswap(VECSIZE_FLOPS,V1c,1,V2c,1);
        end = _rdtsc();
        calcul_flop("mncblas_cswap : ", NB_OPE_COMPLEXE*VECSIZE_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_float_t sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_cswap(VECSIZE_FLOPS,V1c,1,V2c,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_cswap : ", NB_EXPE*NB_OPE_COMPLEXE*VECSIZE_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_double_t\n");
    for(int i = 0; i < NB_EXPE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mncblas_zswap(VECSIZE_FLOPS,V1z,1,V2z,1);
        end = _rdtsc();
        calcul_flop("mncblas_zswap : ", NB_OPE_COMPLEXE*VECSIZE_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_double_t sur NB_EXPE\n");
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mncblas_zswap(VECSIZE_FLOPS,V1z,1,V2z,1);
    }
    end = _rdtsc();
    calcul_flop("mncblas_zswap : ", NB_EXPE*NB_OPE_COMPLEXE*VECSIZE_FLOPS ,end-start);
}