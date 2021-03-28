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
    float* V1s = (float*)malloc(VECSIZE_RESULTAT*sizeof(float));
    double* V1d = (double*)malloc(VECSIZE_RESULTAT*sizeof(double));
    complexe_float_t* V1c = (complexe_float_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_float_t));
    complexe_double_t* V1z = (complexe_double_t*)malloc(VECSIZE_RESULTAT*sizeof(complexe_double_t));
    complexe_float_t tmpc = gen_complexe_float(1,0);
    complexe_double_t tmpz = gen_complexe_double(1,0);
    


    printf("                   TEST_IAMAX\n|||||||||||||||||||||||||||||||||||||||||||||||||||||\n                           1: TEST DE BON RESULTAT\n<------------------------------------------>\n                   float\n");
    printf("---------- test 1 : \n");
    void_vector_sinit(V1s, 1,VECSIZE_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_isamax : %ld\n",mnblas_isamax(VECSIZE_RESULTAT,V1s,1));
    printf("---------- test 2 : \n");
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_isamax : %ld\n",mnblas_isamax(VECSIZE_RESULTAT,V1s,1));
    printf("---------- test 3 : \n");
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_isamax : %ld\n",mnblas_isamax(VECSIZE_RESULTAT,V1s,1));
    printf("---------- test 4 : \n");
    void_vector_sinit_rand(V1s, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1s : "); vector_print(V1s,TYPE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_isamax : %ld\n",mnblas_isamax(VECSIZE_RESULTAT,V1s,1));

    printf("<------------------------------------------>\n                   double\n");

    printf("---------- test 1 : \n");
    void_vector_dinit(V1d, 1,VECSIZE_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_idamax : %ld\n",mnblas_idamax(VECSIZE_RESULTAT,V1d,1));
    printf("---------- test 2 : \n");
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_idamax : %ld\n",mnblas_idamax(VECSIZE_RESULTAT,V1d,1));
    printf("---------- test 3 : \n");
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_idamax : %ld\n",mnblas_idamax(VECSIZE_RESULTAT,V1d,1));
    printf("---------- test 4 : \n");
    void_vector_dinit_rand(V1d, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1d : "); vector_print(V1d,TYPE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_idamax : %ld\n",mnblas_idamax(VECSIZE_RESULTAT,V1d,1));

    printf("<------------------------------------------>\n                   complexe_float_t\n");

    printf("---------- test 1 : \n");
    void_vector_cinit(V1c, tmpc,VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_icamax : %ld\n",mnblas_icamax(VECSIZE_RESULTAT,V1c,1));
    printf("---------- test 2 : \n");
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_icamax : %ld\n",mnblas_icamax(VECSIZE_RESULTAT,V1c,1));
    printf("---------- test 3 : \n");
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_icamax : %ld\n",mnblas_icamax(VECSIZE_RESULTAT,V1c,1));
    printf("---------- test 4 : \n");
    void_vector_cinit_rand(V1c, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1c : "); vector_print(V1c,TYPE_COMPLEXE_FLOAT,VECSIZE_RESULTAT);
    printf("mnblas_icamax : %ld\n",mnblas_icamax(VECSIZE_RESULTAT,V1c,1));

    printf("<------------------------------------------>\n                   complexe_double_t\n");

    printf("---------- test 1 : \n");
    void_vector_zinit(V1z, tmpz,VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_izamax : %ld\n",mnblas_izamax(VECSIZE_RESULTAT,V1z,1));
    printf("---------- test 2 : \n");
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_izamax : %ld\n",mnblas_izamax(VECSIZE_RESULTAT,V1z,1));
    printf("---------- test 3 : \n");
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_izamax : %ld\n",mnblas_izamax(VECSIZE_RESULTAT,V1z,1));
    printf("---------- test 4 : \n");
    void_vector_zinit_rand(V1z, RAND_MAXIMUM,VECSIZE_RESULTAT);
    printf("V1z : "); vector_print(V1z,TYPE_COMPLEXE_DOUBLE,VECSIZE_RESULTAT);
    printf("mnblas_izamax : %ld\n",mnblas_izamax(VECSIZE_RESULTAT,V1z,1));


    free_vm(V1s);
    free_vm(V1d);
    free_vm(V1c);
    free_vm(V1z);

    #define VECSIZE_FLOPS    100000
    #define NB_EXPE_VISIBLE     6
    #define NB_EXPE     1000
    #define NB_OPE_REEL 1
    #define NB_OPE_COMPLEXE 2

    V1s = (float*)malloc(VECSIZE_FLOPS*sizeof(float));
    V1d = (double*)malloc(VECSIZE_FLOPS*sizeof(double));
    V1c = (complexe_float_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_float_t));
    V1z = (complexe_double_t*)malloc(VECSIZE_FLOPS*sizeof(complexe_double_t));


    printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n                  2 : FLOPS\n <-------------------------------------------------->\n                     float\n");
    unsigned long long int start, end ; 
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mnblas_isamax(VECSIZE_FLOPS,V1s,1);
        end = _rdtsc();
        calcul_flop("mnblas_isamax : ", NB_OPE_REEL*VECSIZE_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      float sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mnblas_isamax(VECSIZE_FLOPS,V1s,1);
    }
    end = _rdtsc();
    calcul_flop("mnblas_isamax : ", NB_EXPE*NB_OPE_REEL*VECSIZE_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      double\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mnblas_idamax(VECSIZE_FLOPS,V1d,1);
        end = _rdtsc();
        calcul_flop("mnblas_idamax : ", NB_OPE_REEL*VECSIZE_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      double sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mnblas_idamax(VECSIZE_FLOPS,V1d,1);
    }
    end = _rdtsc();
    calcul_flop("mnblas_idamax : ", NB_EXPE*NB_OPE_REEL*VECSIZE_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_float_t\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mnblas_icamax(VECSIZE_FLOPS,V1c,1);
        end = _rdtsc();
        calcul_flop("mnblas_icamax : ", NB_OPE_COMPLEXE*VECSIZE_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_float_t sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mnblas_icamax(VECSIZE_FLOPS,V1c,1);
    }
    end = _rdtsc();
    calcul_flop("mnblas_icamax : ", NB_EXPE*NB_OPE_COMPLEXE*VECSIZE_FLOPS ,end-start);
    printf("<--------------------------------------------------------------->\n                      complexe_double_t\n");
    for(int i = 0; i < NB_EXPE_VISIBLE; i++) {
        printf("------------------------------------------------\n");
        start = _rdtsc();
        mnblas_izamax(VECSIZE_FLOPS,V1z,1);
        end = _rdtsc();
        calcul_flop("mnblas_izamax : ", NB_OPE_COMPLEXE*VECSIZE_FLOPS ,end-start);
    }
    printf("<--------------------------------------------------------------->\n                      complexe_double_t sur NB_EXPE (%d)\n",NB_EXPE);
    start = _rdtsc();
    for(int i = 0; i < NB_EXPE; i++) {
        mnblas_izamax(VECSIZE_FLOPS,V1z,1);
    }
    end = _rdtsc();
    calcul_flop("mnblas_izamax : ", NB_EXPE*NB_OPE_COMPLEXE*VECSIZE_FLOPS ,end-start);



    free_vm(V1s);
    free_vm(V1d);
    free_vm(V1c);
    free_vm(V1z);
}