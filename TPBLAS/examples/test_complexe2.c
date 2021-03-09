#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe2.h"

#define    NB_FOIS        4194304

#include "flop.h"

int main (int argc, char **argv)
{
 complexe_float_t c1= {1.0, 2.0} ;
 complexe_float_t c2= {3.0, 6.0} ;
 complexe_float_t c3 ;

 complexe_double_t cd1 ;
 complexe_double_t cd2 ;
 complexe_double_t cd3 ;

 unsigned long long int start, end ;
 int i ;

 init_flop () ;
 
 c3 = add_complexe_float (c1, c2) ;

 printf ("c3.r %f c3.i %f\n", c3.real, c3.imaginary) ;

 cd1 = (complexe_double_t) {10.0, 7.0} ;
 cd2 = (complexe_double_t) {25.0, 32.0} ;

 cd3 = add_complexe_double (cd1, cd2) ;

 printf ("cd3.r %f cd3.i %f\n", cd3.real, cd3.imaginary) ;

 start =_rdtsc () ;
 
 for (i = 0 ; i < NB_FOIS; i++)
   {
     cd3 = add_complexe_double (cd1, cd2) ;
   }

 end = _rdtsc () ;

  printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", cd3.real, cd3.imaginary, end-start) ;
  calcul_flop ("calcul complexe ", NB_FOIS*2, end-start) ;

  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Multiplication complexe float\n----------------------------------------------\n");
  c3 = mult_complexe_float(c1,c2);
  printf("c1 * c2 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c2,c1);
  printf("c2 * c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c1,c1);
  printf("c1 * c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = mult_complexe_float(c2,c2);
  printf("c2 * c2 = %f+%fi\n",c3.real,c3.imaginary);
  start =_rdtsc () ;
  for (i = 0 ; i < NB_FOIS; i++)
    {
      c3 = mult_complexe_float (c1, c2) ;
    }
  end = _rdtsc () ;
  printf("Calcul flop : \n");
  calcul_flop ("Mult complexe float", NB_FOIS*6, end-start) ;
  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Multiplication complexe double\n----------------------------------------------\n");
  cd3 = mult_complexe_double(cd1,cd2);
  printf("cd1 * cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd2,cd1);
  printf("cd2 * cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd1,cd1);
  printf("cd1 * cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = mult_complexe_double(cd2,cd2);
  printf("cd2 * cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  start =_rdtsc () ;
  for (i = 0 ; i < NB_FOIS; i++)
    {
      cd3 = mult_complexe_double (cd1, cd2) ;
    }
  end = _rdtsc () ;
  printf("Calcul flop : \n");
  calcul_flop ("Mult complexe double", NB_FOIS*6, end-start) ;

  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Division complexe float\n----------------------------------------------\n");
  c3 = div_complexe_float(c1,c2);
  printf("c1 / c2 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = div_complexe_float(c2,c1);
  printf("c2 / c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = div_complexe_float(c1,c1);
  printf("c1 / c1 = %f+%fi\n",c3.real,c3.imaginary);
  c3 = div_complexe_float(c2,c2);
  printf("c2 / c2 = %f+%fi\n",c3.real,c3.imaginary);
  start =_rdtsc () ;
  for (i = 0 ; i < NB_FOIS; i++)
    {
      c3 = div_complexe_float (c1, c2) ;
    }
  end = _rdtsc () ;
  printf("Calcul flop : \n");
  calcul_flop ("div complexe float", NB_FOIS*6, end-start) ;
  printf("||||||||||||||||||||||||||||||||||||||||||||||\n        Division complexe double\n----------------------------------------------\n");
  cd3 = div_complexe_double(cd1,cd2);
  printf("cd1 / cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = div_complexe_double(cd2,cd1);
  printf("cd2 / cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = div_complexe_double(cd1,cd1);
  printf("cd1 / cd1 = %f+%fi\n",cd3.real,cd3.imaginary);
  cd3 = div_complexe_double(cd2,cd2);
  printf("cd2 / cd2 = %f+%fi\n",cd3.real,cd3.imaginary);
  start =_rdtsc () ;
  for (i = 0 ; i < NB_FOIS; i++)
    {
      cd3 = div_complexe_double (cd1, cd2) ;
    }
  end = _rdtsc () ;
  printf("Calcul flop : \n");
  calcul_flop ("div complexe double", NB_FOIS*6, end-start) ;
  exit (0) ;
}


