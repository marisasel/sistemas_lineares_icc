#include <stdio.h>

int main ()
{

  char vetor_01[5], var_01;
  int  vetor_02[6], var_02;
  float vetor_03[7], var_03;
  double  vetor_04[6], var_04;
  long int var_05;
  long long int var_06;
  short int var_07;

  struct cls {
    int a;
    char b[5];
    int c[10];
  };

  printf("\n");
  printf("\n");

  printf("\n");
  printf(" ================================== \n");
  printf("|    TIPO       | TAMANHO EM BYTES |\n");
  printf("|==================================|\n");
  printf("| char          |       %2lu         |\n", sizeof(char));
  printf("| int           |       %2lu         |\n", sizeof(int));
  printf("| unsig int     |       %2lu         |\n", sizeof(unsigned int));
  printf("| short int     |       %2lu         |\n", sizeof(short int));
  printf("| long int      |       %2lu         |\n", sizeof(long int));
  printf("| long long int |       %2lu         |\n", sizeof(long long int));
  printf("| float         |       %2lu         |\n", sizeof(float));
  printf("| double        |       %2lu         |\n", sizeof(double));
  printf("| long double   |       %2lu         |\n", sizeof(long double));
  printf("| struct cls    |       %2lu         |\n", sizeof(struct cls));
  printf(" ==================================|\n");
  printf("| void *        |       %2lu         |\n", sizeof(void *));
  printf("| char *        |       %2lu         |\n", sizeof(char *));
  printf("| int  *        |       %2lu         |\n", sizeof(int *));
  printf("| short int *   |       %2lu         |\n", sizeof(short int *));
  printf("| long int  *   |       %2lu         |\n", sizeof(long int *));
  printf("| float     *   |       %2lu         |\n", sizeof(float *));
  printf("| double    *   |       %2lu         |\n", sizeof(double *));
  printf("| long double * |       %2lu         |\n", sizeof(long double *));
  printf(" ==================================  \n");

  printf("\n");
  printf("\n");

  printf(" ====================================== \n");
  printf("|      VARIAVEL     | TAMANHO EM BYTES |\n");
  printf("|======================================|\n");
  printf("| char   v[5]       |       %2lu         |\n", sizeof vetor_01 );
  printf("| int    v[6]       |       %2lu         |\n", sizeof vetor_02 );
  printf("| float  v[7]       |       %2lu         |\n", sizeof vetor_03 );
  printf("| double v[6]       |       %2lu         |\n", sizeof(vetor_04));
  printf("| char   var        |       %2lu         |\n", sizeof var_01 );
  printf("| int    var        |       %2lu         |\n", sizeof var_02 );
  printf("| short int var     |       %2lu         |\n", sizeof(var_07));
  printf("| long int var      |       %2lu         |\n", sizeof(var_05));
  printf("| long long int var |       %2lu         |\n", sizeof(var_06));
  printf("| float  var        |       %2lu         |\n", sizeof var_03 );
  printf("| double var        |       %2lu         |\n", sizeof(var_04));
  printf(" ====================================== \n");

  printf("\n");
  printf("\n");

}
