#include <stdio.h>
#include <x86intrin.h>

#ifndef TESTUILS_H
#define TESTUILS_H
#include "testutils.h"
#endif

// La frequence du processeur est de 2.6 GHZ
static const float duree_cycle = (float) 1 / (float) 1.8 ;
// duree du cycle en nano seconde 10^-9

static unsigned long long int residu ;

void init_flop () {
  unsigned long long int start, end ;
  start = _rdtsc () ;
  end =_rdtsc () ;
  residu = end - start ;
}

void calcul_flop (char *message, int nb_operations_flottantes, unsigned long long int cycles) {
  printf ("%s %d operations %5.3f GFLOP/s\n", message, nb_operations_flottantes, ((float) nb_operations_flottantes) / (((float) (cycles - residu)) * duree_cycle)) ;
}

void free_vm (void* address) {
    free(address);
}

complexe_float_t gen_complexe_float (const float real, const float imaginary) {
    complexe_float_t tmp; tmp.real = real; tmp.imaginary = imaginary;
    return tmp;
}
complexe_double_t gen_complexe_double (const double real, const double imaginary) {
    complexe_double_t tmp; tmp.real = real; tmp.imaginary = imaginary;
    return tmp;
}



/*
    VECTOR_INIT
*/

void void_vector_sinit (float* V, const float x, const register unsigned int len) {
    for(register unsigned int i = 0; i < len; i++) {
        V[i] = x;
    }
}

void void_vector_dinit (double* V, const double x, const register unsigned int len) {
    for(register unsigned int i = 0; i < len; i++) {
        V[i] = x;
    }
}

void void_vector_cinit (complexe_float_t* V, const complexe_float_t x, const register unsigned int len) {
    for(register unsigned int i = 0; i < len; i++) {
        V[i] = x;
    }
}

void void_vector_zinit (complexe_double_t* V, const complexe_double_t x, const register unsigned int len) {
    for(register unsigned int i = 0; i < len; i++) {
        V[i] = x;
    }
}

void void_vector_cinit2 (complexe_float_t* V, const float real, const float imaginary, const register unsigned int len) {
    for(register unsigned int i = 0; i < len; i++) {
        V[i].real = real;
        V[i].imaginary = imaginary;
    }
}

void void_vector_zinit2 (complexe_double_t* V, const double real, const double imaginary, const register unsigned int len) {
    for(register unsigned int i = 0; i < len; i++) {
        V[i].real = real;
        V[i].imaginary = imaginary;
    }
}

float* vector_sinit (const float x, const register unsigned int len) {
    float* V = (float*)malloc(len*sizeof(float));
    for(register unsigned int i = 0; i < len; i++) {
        V[i] = x;
    }
    return V;
}

double* vector_dinit (const double x, const register unsigned int len) {
    double* V = (double*)malloc(len*sizeof(double));
    for(register unsigned int i = 0; i < len; i++) {
        V[i] = x;
    }
    return V;
}

complexe_float_t* vector_cinit (const complexe_float_t x, const register unsigned int len) {
    complexe_float_t* V = (complexe_float_t*)malloc(len*sizeof(complexe_float_t));
    for(register unsigned int i = 0; i < len; i++) {
        V[i] = x;
    }
    return V;
}

complexe_double_t* vector_zinit (const complexe_double_t x, const register unsigned int len) {
    complexe_double_t* V = (complexe_double_t*)malloc(len*sizeof(complexe_double_t));
    for(register unsigned int i = 0; i < len; i++) {
        V[i] = x;
    }
    return V;
}

complexe_float_t* vector_cinit2 (const float real, const float imaginary, const register unsigned int len) {
    complexe_float_t* V = (complexe_float_t*)malloc(len*sizeof(complexe_float_t));
    for(register unsigned int i = 0; i < len; i++) {
        V[i].real = real;
        V[i].imaginary = imaginary;
    }
    return V;
}

complexe_double_t* vector_zinit2 (const double real, const double imaginary, const register unsigned int len) {
    complexe_double_t* V = (complexe_double_t*)malloc(len*sizeof(complexe_double_t));
    for(register unsigned int i = 0; i < len; i++) {
        V[i].real = real;
        V[i].imaginary = imaginary;
    }
    return V;
}

/*
    END VECTOR_INIT
*/


/*
    VECTOR_PRINT
*/

void vector_print (void* V, VTYPE type, const register unsigned int len) {
    register unsigned int i;
    printf("( ");
    if(type == TYPE_FLOAT) {
        float* tmp = (float*) V;
        for(i = 0; i < len-1; i++)
            printf ("%.2f, ", tmp[i]);
        printf ("%.2f", tmp[len-1]);
    } else if(type == TYPE_DOUBLE){
        double* tmp = (double*) V;
        for(i = 0; i < len-1; i++)
            printf ("%.2f, ", tmp[i]);
        printf ("%.2f", tmp[len-1]);
    } else if (type == TYPE_COMPLEXE_FLOAT) {
        complexe_float_t* tmp = (complexe_float_t*) V;
        for (i = 0; i < len-1; i++)
            printf ("%.2f+%.2fi, ", tmp[i].real,tmp[i].imaginary);
        printf ("%.2f+%.2fi", tmp[len-1].real,tmp[len-1].imaginary);
    } else if (type == TYPE_COMPLEXE_DOUBLE) {
        complexe_double_t* tmp = (complexe_double_t*) V;
        for (i = 0; i < len-1; i++)
            printf ("%.2f+%.2fi, ", tmp[i].real,tmp[i].imaginary);
        printf ("%.2f+%.2fi", tmp[len-1].real,tmp[len-1].imaginary);
    }
    printf(" )\n");
}

/*
    END VECTOR_PRINT
*/



/* ||||||||||||||||| MATRIX ||||||||||||||||||| */

/*
    MATRIX_INIT
*/

void void_matrix_sinit (float* V, const MNCBLAS_LAYOUT layout, const float x, const int M, const int N) {
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[i*N+j] = x;
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[j*N+i] = x;
        }
    }
}

void void_matrix_dinit (double* V, const MNCBLAS_LAYOUT layout, const double x, const int M, const int N) {
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[i*N+j] = x;
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[j*N+i] = x;
        }
    }
}

void void_matrix_cinit (complexe_float_t* V, const MNCBLAS_LAYOUT layout, const complexe_float_t x, const int M, const int N) {
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[i*N+j] = x;
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[j*N+i] = x;
        }
    }
}

void void_matrix_zinit (complexe_double_t* V, const MNCBLAS_LAYOUT layout, const complexe_double_t x, const int M, const int N) {
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[i*N+j] = x;
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[j*N+i] = x;
        }
    }
}

void void_matrix_cinit2 (complexe_float_t* V, const MNCBLAS_LAYOUT layout, const float real, const float imaginary, const int M, const int N) {
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++) {
                V[i*N+j].real = real;
                V[i*N+j].imaginary = imaginary;
            }
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++) {
                V[j*N+i].real = real;
                V[j*N+i].imaginary = imaginary;
            }
        }
    }
}

void void_matrix_zinit2 (complexe_double_t* V, const MNCBLAS_LAYOUT layout, const double real, const double imaginary, const int M, const int N) {
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++) {
                V[i*N+j].real = real;
                V[i*N+j].imaginary = imaginary;
            }
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++) {
                V[j*N+i].real = real;
                V[j*N+i].imaginary = imaginary;
            }
        }
    }
}

float* matrix_sinit (const MNCBLAS_LAYOUT layout, const float x, const int M, const int N) {
    float* V = (float*)malloc(M*N*sizeof(float));
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[i*N+j] = x;
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[j*N+i] = x;
        }
    }
    return V;
}

double* matrix_dinit (const MNCBLAS_LAYOUT layout, const double x, const int M, const int N) {
    double* V = (double*)malloc(M*N*sizeof(double));
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[i*N+j] = x;
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[j*N+i] = x;
        }
    }
    return V;
}

complexe_float_t* matrix_cinit (const MNCBLAS_LAYOUT layout, const complexe_float_t x, const int M, const int N) {
    complexe_float_t* V = (complexe_float_t*)malloc(M*N*sizeof(complexe_float_t));
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[i*N+j] = x;
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[j*N+i] = x;
        }
    }
    return V;
}

complexe_double_t* matrix_zinit (const MNCBLAS_LAYOUT layout, const complexe_double_t x, const int M, const int N) {
    complexe_double_t* V = (complexe_double_t*)malloc(M*N*sizeof(complexe_double_t));
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[i*N+j] = x;
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++)
                V[j*N+i] = x;
        }
    }
    return V;
}

complexe_float_t* matrix_cinit2 (const MNCBLAS_LAYOUT layout, const float real, const float imaginary, const int M, const int N) {
    complexe_float_t* V = (complexe_float_t*)malloc(M*N*sizeof(complexe_float_t));
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++) {
                V[i*N+j].real = real;
                V[i*N+j].imaginary = imaginary;
            }
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++) {
                V[j*N+i].real = real;
                V[j*N+i].imaginary = imaginary;
            }
        }
    }
    return V;
}

complexe_double_t* matrix_zinit2 (const MNCBLAS_LAYOUT layout, const double real, const double imaginary, const int M, const int N) {
    complexe_double_t* V = (complexe_double_t*)malloc(M*N*sizeof(complexe_double_t));
    if(layout == MNCblasRowMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++) {
                V[i*N+j].real = real;
                V[i*N+j].imaginary = imaginary;
            }
        }
    } else if (layout == MNCblasColMajor) {
        for(register int i = 0, j; i < M; i++) {
            for(j = 0; j < M; j++) {
                V[j*N+i].real = real;
                V[j*N+i].imaginary = imaginary;
            }
        }
    }
    return V;
}

/*
    END MATRIX_INIT
*/


/*
    MATRIX_PRINT
*/

void matrix_print (const void* V, const VTYPE type, const MNCBLAS_LAYOUT layout ,const int M ,const int N) {
    register unsigned int i,j;
    if(type == TYPE_FLOAT) {
        float* tmp = (float*) V;
        if(layout == MNCblasRowMajor) {
            for(i = 0; i < M; i++) {
                printf("| "); for(j = 0 ; j < N ; j++) { printf ("%.2f ", tmp[i*N+j]); } printf("|\n");
            }
        } else if (layout == MNCblasColMajor) {
            for(i = 0; i < M; i++) {
                printf("| "); for(j = 0 ; j < N ; j++) { printf ("%.2f ", tmp[j*N+i]); } printf("|\n");
            }
        }
    } else if(type == TYPE_DOUBLE){
        double* tmp = (double*) V;
        if(layout == MNCblasRowMajor) {
            for(i = 0; i < M; i++) {
                printf("| "); for(j = 0 ; j < N ; j++) { printf ("%.2f ", tmp[i*N+j]); } printf("|\n");
            }
        } else if (layout == MNCblasColMajor) {
            for(i = 0; i < M; i++) {
                printf("| "); for(j = 0 ; j < N ; j++) { printf ("%.2f ", tmp[j*N+i]); } printf("|\n");
            }
        }
    } else if (type == TYPE_COMPLEXE_FLOAT) {
        complexe_float_t* tmp = (complexe_float_t*) V;
        if(layout == MNCblasRowMajor) {
            for(i = 0; i < M; i++) {
                printf("| "); for(j = 0 ; j < N ; j++) { printf ("%.2f+%.2fi ", tmp[i*N+j].real,tmp[i*N+j].imaginary); } printf("|\n");
            }
        } else if (layout == MNCblasColMajor) {
            for(i = 0; i < M; i++) {
                printf("| "); for(j = 0 ; j < N ; j++) { printf ("%.2f+%.2fi ", tmp[j*N+i].real, tmp[j*N+i].imaginary); } printf("|\n");
            }
        }
    } else if (type == TYPE_COMPLEXE_DOUBLE) {
        complexe_double_t* tmp = (complexe_double_t*) V;
        if(layout == MNCblasRowMajor) {
            for(i = 0; i < M; i++) {
                printf("| "); for(j = 0 ; j < N ; j++) { printf ("%.2f+%.2fi ", tmp[i*N+j].real,tmp[i*N+j].imaginary); } printf("|\n");
            }
        } else if (layout == MNCblasColMajor) {
            for(i = 0; i < M; i++) {
                printf("| "); for(j = 0 ; j < N ; j++) { printf ("%.2f+%.2fi ", tmp[j*N+i].real, tmp[j*N+i].imaginary); } printf("|\n");
            }
        }
    }
}

/*
    END MATRIX_PRINT
*/