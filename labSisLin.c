//***************************************
//  Exercício: sistemas lineares        *
//  Disciplina: ICC                     *
//  Curso: Informática Biomédica - UFPR *
//  Aluna: Marisa Sel Franco            *
//***************************************

#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

int main () {    
  
  int tam = 9;
  int n[9] = {10, 30, 50, 128, 256, 512, 1000, 2000, 3000};
  char colunas[] = {"n     |   t_egp     |     t_gs     |  it_gs  | normaResiduo_gs    |     t_ref       |  it_ref | normaResiduo_ref\n"};

  printf("\n\n");  
  printf("Soluciona sistemas com matrizes diagonalmente dominantes alocados como 'pontPont':\n\n");
  printf("%s", colunas); 
  printf("=======================================================================================================================\n");

  for (int i = 0; i < tam; i++) {
    printf("%d\t", n[i]);
    // inicializa gerador de números aleatórios
    srand(202201); 
    int iteracoes = 0;
    real_t *x, tTotal;
    real_t normaResiduo_egp, normaResiduo_gs, normaResiduo_ref; 
    SistLinear_t *SL;                       // aloca ponteiro para sistema linear
    SL = alocaSisLin(n[i],pontPont);
    iniSisLin(SL, diagDominante, COEF_MAX); // inicializa sistema linear diagonalmente dominante
    x = (real_t*)malloc(n[i] * sizeof(real_t));
    
    // Método da Eliminação de Gauss
    iteracoes = eliminacaoGauss (SL, x, &tTotal);
    printf("%lf\t", tTotal);

    // Método de Gauss-Seidel
    for (int l = 0; l < SL->n; l++)             // zera o vetor de x para dar um chute inicial 
      x[l] = 0.0;
    iteracoes = gaussSeidel(SL, x, ERRO, &tTotal);
    printf("%lf\t", tTotal);
    printf("%d\t", iteracoes);
    normaResiduo_gs = normaL2Residuo(SL, x);
    printf("%e\t\t", normaResiduo_gs);

    // Método de Refinamento
    iteracoes = refinamento(SL, x, ERRO, &tTotal);
    printf("%lf\t", tTotal);
    printf("%d\t", iteracoes);
    normaResiduo_ref = normaL2Residuo(SL, x);
    printf("%e\n", normaResiduo_ref);

    liberaSisLin(SL); //LIBERA MEMÓRIA
    free(x);
  }

  printf("\n\n");
  printf("Soluciona sistemas com matrizes dominantes alocados como 'pontVet':\n\n");
  printf("%s", colunas); 
  printf("=======================================================================================================================\n");
  
  for (int i = 0; i < tam; i++) {
    printf("%d\t", n[i]);
    // inicializa gerador de números aleatórios
    srand(202201); 
    int iteracoes = 0;
    real_t *x, tTotal;
    real_t normaResiduo_egp, normaResiduo_gs, normaResiduo_ref; 
    SistLinear_t *SL;                       // aloca ponteiro para sistema linear
    SL = alocaSisLin(n[i],pontVet);
    iniSisLin(SL, diagDominante, COEF_MAX); // inicializa sistema linear diagonalmente dominante
    x = (real_t*)malloc(n[i] * sizeof(real_t));
    
    // Método da Eliminação de Gauss
    iteracoes = eliminacaoGauss (SL, x, &tTotal);
    printf("%lf\t", tTotal);

    // Método de Gauss-Seidel
    for (int l = 0; l < SL->n; l++)             // zera o vetor de x para dar um chute inicial 
      x[l] = 0.0;
    iteracoes = gaussSeidel(SL, x, ERRO, &tTotal);
    printf("%lf\t", tTotal);
    printf("%d\t", iteracoes);
    normaResiduo_gs = normaL2Residuo(SL, x);
    printf("%e\t\t", normaResiduo_gs);

    // Método de Refinamento
    iteracoes = refinamento(SL, x, ERRO, &tTotal);
    printf("%lf\t", tTotal);
    printf("%d\t", iteracoes);
    normaResiduo_ref = normaL2Residuo(SL, x);
    printf("%e\n", normaResiduo_ref);

    liberaSisLin(SL); //LIBERA MEMÓRIA
    free(x);
  }

  printf("\n\n");
  printf("Soluciona sistemas com matrizes de Hilbert alocados como 'pontPont':\n\n");
  printf("%s", colunas); 
  printf("=======================================================================================================================\n");

  for (int i = 0; i < tam; i++) {
    printf("%d\t", n[i]);
    // inicializa gerador de números aleatórios
    srand(202201); 
    int iteracoes = 0;
    real_t *x, tTotal;
    real_t normaResiduo_egp, normaResiduo_gs, normaResiduo_ref; 
    SistLinear_t *SL;                       // aloca ponteiro para sistema linear
    SL = alocaSisLin(n[i],pontPont);
    iniSisLin(SL, hilbert, COEF_MAX);       // inicializa sistema linear com matriz do tipo Hilbert
    x = (real_t*)malloc(n[i] * sizeof(real_t));
    
    // Método da Eliminação de Gauss
    iteracoes = eliminacaoGauss (SL, x, &tTotal);
    printf("%lf\t", tTotal);

    // Método de Gauss-Seidel
    for (int l = 0; l < SL->n; l++)             // zera o vetor de x para dar um chute inicial 
      x[l] = 0.0;
    iteracoes = gaussSeidel(SL, x, ERRO, &tTotal);
    printf("%lf\t", tTotal);
    printf("%d\t", iteracoes);
    normaResiduo_gs = normaL2Residuo(SL, x);
    printf("%e\t\t", normaResiduo_gs);

    // Método de Refinamento
    iteracoes = refinamento(SL, x, ERRO, &tTotal);
    printf("%lf\t", tTotal);
    printf("%d\t", iteracoes);
    normaResiduo_ref = normaL2Residuo(SL, x);
    printf("%e\n", normaResiduo_ref);

    liberaSisLin(SL); //LIBERA MEMÓRIA
    free(x);
  }

  printf("\n\n");
  printf("Soluciona sistemas com matrizes de Hilbert alocados como 'pontVet':\n\n");
  printf("%s", colunas); 
  printf("=======================================================================================================================\n");
  
  for (int i = 0; i < tam; i++) {
    printf("%d\t", n[i]);
    // inicializa gerador de números aleatórios
    srand(202201); 
    int iteracoes = 0;
    real_t *x, tTotal;
    real_t normaResiduo_egp, normaResiduo_gs, normaResiduo_ref; 
    SistLinear_t *SL;                       // aloca ponteiro para sistema linear
    SL = alocaSisLin(n[i],pontVet);
    iniSisLin(SL, hilbert, COEF_MAX);       // inicializa sistema linear com matriz do tipo Hilbert
    x = (real_t*)malloc(n[i] * sizeof(real_t));
    
    // Método da Eliminação de Gauss
    iteracoes = eliminacaoGauss (SL, x, &tTotal);
    printf("%lf\t", tTotal);

    // Método de Gauss-Seidel
    for (int l = 0; l < SL->n; l++)             // zera o vetor de x para dar um chute inicial 
      x[l] = 0.0;
    iteracoes = gaussSeidel(SL, x, ERRO, &tTotal);
    printf("%lf\t", tTotal);
    printf("%d\t", iteracoes);
    normaResiduo_gs = normaL2Residuo(SL, x);
    printf("%e\t\t", normaResiduo_gs);

    // Método de Refinamento
    iteracoes = refinamento(SL, x, ERRO, &tTotal);
    printf("%lf\t", tTotal);
    printf("%d\t", iteracoes);
    normaResiduo_ref = normaL2Residuo(SL, x);
    printf("%e\n", normaResiduo_ref);

    liberaSisLin(SL); //LIBERA MEMÓRIA
    free(x);
  }

  return 1;
}