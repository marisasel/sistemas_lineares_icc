#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

// Função que faz uma cópia idêntica do SL 
SistLinear_t *copiaSistemaLinear(SistLinear_t *SL) {

  SistLinear_t *SL_alterado;
  int tipo_alocacao = SL->tipoAloc_A;
  unsigned int n = SL->n;
  
  if (tipo_alocacao == 0) {
    SL_alterado = alocaSisLin(n, pontPont); // aloca sistema como matriz como vetor de N ponteiros para vetores de tamanho N
    for (int i = 0; i < n; ++i)
      memcpy(SL_alterado->A[i], SL->A[i], n * sizeof(real_t));
    memcpy(SL_alterado->b, SL->b, n * sizeof(real_t));
  }
  else {
    SL_alterado = alocaSisLin(n, pontVet);  // aloca sistema como matriz como vetor de N ponteiros para um único vetor de tamanho N*N
    memcpy(SL_alterado->A[0], SL->A[0], (n * n) * sizeof(real_t));
    memcpy(SL_alterado->b, SL->b, n * sizeof(real_t));
  }
  return SL_alterado;
}

// Função que encontra o maior valor da coluna para definir o pivô e devolve o índice da linha correspondente
int encontraMax(SistLinear_t *SL, int indice) {
  int indice_pivo = indice;
  unsigned int n = SL->n;
  //printf("A coluna em análise é %d.\n", indice_pivo);  //APAGAR
  for(int i = indice + 1; i < n; i++) {
    if(fabs(SL->A[i][indice]) > fabs(SL->A[indice_pivo][indice])) {
      indice_pivo = i;
    }
  }
  //printf("Terminou o laço e viu todas as linhas. O índice da linha pivô é %d. O valor é: %10g \n", indice_pivo, SL->A[indice_pivo][indice]);  // APAGA
  return indice_pivo; 
}

// Função que troca a linha pela linha com o pivô
void trocaLinha (SistLinear_t *SL, int i, int indice_pivo) {
  unsigned int n = SL->n;
  real_t *aux = (real_t*)malloc (n * sizeof(real_t)), aux2;     // vetores auxiliares

  if (!aux) {
    fprintf(stderr, "Erro ao alocar vetor auxiliar para trocar linhas.\n");  //TESTAR
    //return -1;
  }

  memcpy(aux, SL->A[i], n * sizeof(real_t));                    // copia a linha i para o vetor auxiliar
  memcpy(SL->A[i], SL->A[indice_pivo], n * sizeof(real_t));     // copia a linha com o pivô para a linha i
  memcpy(SL->A[indice_pivo], aux, n * sizeof(real_t));          // copia o vetor auxiliar para a linha onde estava o pivô

  aux2 = SL->b[i];                // copia o termo independente i para o vetor auxiliar
  SL->b[i] = SL->b[indice_pivo];  // copia o termo independente da linha com o pivô para a linha i
  SL->b[indice_pivo] = aux2;      // copia o termo independente no vetor auxiliar para a linha onde estava o pivô

  free(aux);
  return;
}

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal) {
  *tTotal = timestamp();
  unsigned int n = SL->n;

  SistLinear_t *SL_alterado = copiaSistemaLinear(SL);         // aloca ponteiro e copia o sistema linear original
  // ESCALONAMENTO: PIVOTEAMENTO PARCIAL DO SL COM TRIANGULAÇÃO RESULTANDO EM MATRIZ SUPERIOR
  for (int i = 0; i < n; ++i) {                               // para cada equação (linha) do sistema linear
    int indice_linha_pivo = encontraMax(SL_alterado, i);      // encontra o maior valor da coluna correspondente
    if (i != indice_linha_pivo)                               // se o índice da liha pivô for diferente do índice da equação atuak
      trocaLinha(SL_alterado, i, indice_linha_pivo);          // troca a linha da equação pela linha pivô
    for (int k = i + 1; k < n ; ++k) {                        // para cada coluna/ componente subsequente
      double m = SL_alterado->A[k][i] / SL_alterado->A[i][i]; // calcula o valor do multiplicador m
      if (isfinite(m) == 0) {
        fprintf(stderr, "O cálculo do multiplicador na EGP deu infinty ou NaN.\n");
        return -1;
      }
      SL_alterado->A[k][i] = 0.0;                             // atribui ao 1º elemento o valor de 0.0
      for (int j = i + 1; j < n; ++j) {                       // para o 2º elemento da linha em diante
        SL_alterado->A[k][j] -= SL_alterado->A[i][j] * m;     // calcula o novo valor do elemento
        if (isfinite(SL_alterado->A[k][j]) == 0) {
          fprintf(stderr, "O cálculo da novo valor do elemento na EGP deu underflow ou overflow.\n");
          return -2;
        }
      }
      SL_alterado->b[k] -= SL_alterado->b[i] * m;             // atualiza o valor do termo independente da equação
      if (isfinite(SL_alterado->b[k]) == 0) {
        fprintf(stderr, "O cálculo da termo independente na EGP deu underflow ou overflow.\n");
        return -3;
      }
    }  
  }
  //printf("Imprime a cópia do sistema dentro da função da Elimição de Gauss:\n"); // APAGAR
  //prnSisLin(SL_alterado); //APAGAR
  
  // RETROSUBSTITUIÇÃO
  real_t dividendo = 0.0;
  int i = n - 1;
  
  x[i] = SL_alterado->b[i] / SL_alterado->A[i][i];            // calcula o iésimo x
  if (isfinite(x[i]) == 0) {
    fprintf(stderr, "O cálculo do iésimo x na EGP deu infinty ou NaN.\n");
    return -4;
  }
  //printf("O valor do x %d é %g\n", i, x[i]); // APAGAR
  for (i = n - 2; i >= 0; i--) {                              // para cada equação
    dividendo = SL_alterado->b[i];                            // insere no dividendo o valor do termo independente da equação
    for (int j = i + 1; j < n; j++) {                         // para os demais elementos da equação
      dividendo -= (SL_alterado->A[i][j] * x[j]);             // substrai o valor do cálculo de x do dividendo
      if (isfinite(dividendo) == 0) {
        fprintf(stderr, "O cálculo do dividendo na EGP deu underflow ou overflow.\n");
        return -5;
      }
    }
    x[i] = dividendo / SL_alterado->A[i][i];                  // ao final das subtrações, divide o dividendo encontrado pelo multiplicador do iésimo x na equação
    if (isfinite(x[i]) == 0) {
      fprintf(stderr, "O cálculo do iésimo x na EGP deu infinty ou NaN.\n");
      return -4;
    }
    //printf("O valor do x %d é %g\n", i, x[i]); // APAGAR
  }

  liberaSisLin(SL_alterado); //LIBERA MEMÓRIA
  *tTotal = timestamp() - *tTotal;
  return 0;
}


/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/

real_t normaL2Residuo(SistLinear_t *SL, real_t *x) {
  real_t linha, norma = 0.0;
  unsigned int n = SL->n;

  for (int i = 0; i < n; i++) {
    linha = 0.0;
    for (int j = 0; j < n; j ++) 
      linha += SL->A[i][j] * x[j];
    norma += (SL->b[i] - linha) * (SL->b[i] - linha);
  }
  if (isfinite(norma) == 0) {
    fprintf(stderr, "O cálculo da norma do resíduo deu underflow ou overflow.\n");
    return -1;
  }
  return sqrt(norma);
}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal) {
  *tTotal = timestamp();
  int k, i, j;
  real_t soma, diff, xk;
  real_t norma = 1.0 + erro;                                  // valor atribuído à norma apenas para entrar no laço
  unsigned int n = SL->n;

  // até que sejam atingidos os critérios de parada: 
  // norma máxima do erro absoluto aproximado em X e o número máximo de iterações
  for (k = 0; norma > erro && k < MAXIT; ++k) {
    norma = 0.0;                                              // zera o valor da norma
    for (i = 0; i < n; ++i) {                                 // para cada equação
      for (soma = 0, j = 0; j < i; ++j)                       // para cada coluna/ componente da equação anterior a i
        soma += SL->A[i][j] * x[j];                           // faz o somatória
      for (j = i + 1; j < n; ++j)                             // para cada coluna/ componente da equação posterior a i
        soma += SL->A[i][j] * x[j];                           // faz o somatória
      if (isfinite(soma) == 0) {
        fprintf(stderr, "O cálculo do somatório no GS deu underflow ou overflow.\n");
        return (-k - 1);
      }
      xk = (SL->b[i] - soma) / SL->A[i][i];                   // calcula o valor de x para a iteração i
      if (isfinite(xk) == 0) {
        fprintf(stderr, "O cálculo do iésimo x no GS deu infinty ou NaN.\n");
        return (-k - 1);
      }
      // faz o cálculo da norma máxima do erro absoluto aproximado em x
      diff = fabs(xk - x[i]);
      x[i] = xk;
      if (diff > norma)
        norma = diff;
    }
  }
  if (k >= (MAXIT)){
    fprintf(stderr, "O método de Gauss-Seidel não convergiu porque a matriz não é diagonalmente dominante.\n");
    *tTotal = timestamp() - *tTotal;
    return -k;
  }
  else {
    *tTotal = timestamp() - *tTotal;
    return k;
  }    
}

int calculaResiduo(SistLinear_t *SL, real_t *x, real_t *residuo) {
  real_t linha;
  unsigned int n = SL->n;

  for (int j = 0; j < n; j++) {
    linha = 0.0;
    for (int k = 0; k < n; k++) { 
      linha += SL->A[j][k] * x[k];
    }
    residuo[j] = (SL->b[j] - linha);
  }
  return 0;
}

/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int refinamento (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal) {
  *tTotal = timestamp();
  int iteracoes = 0, i, res;
  unsigned int n = SL->n;
  real_t tGauss, *residuo, *w_vetor_diferencas;
  real_t norma = 1.0 + erro;                                  // valor atribuído à norma apenas para entrar no laço

  SistLinear_t *SL_alterado = copiaSistemaLinear(SL);         // aloca ponteiro e copia o sistema linear original

  residuo = (real_t*)malloc(n * sizeof(real_t));
  if (!residuo) {
    fprintf(stderr, "Erro no refinamento: falha ao alocar vetor para resíduos.\n");
    return -1;
  }
  w_vetor_diferencas = (real_t*)malloc(n * sizeof(real_t));
  if (!w_vetor_diferencas) {
    fprintf(stderr, "Erro no refinamento: falha ao alocar vetor para diferenças (w).\n");
    return -1;
  }

  // 1º passo: obtém solução inicial x(0) resolvendo Ax = b e inicializar i = 0
  // -> feito fora desta função, por GS

  // Enquanto não atingir os critérios de parada
  for (i = 0; norma > erro && i < MAXIT; ++i) {
    // 2º passo: calcula o resíduo r = b - Ax^(i)
    res = calculaResiduo(SL, x, residuo);
    if (res != 0)
      return -2;
    // 3º passo: resolve Aw = r
    memcpy(SL_alterado->b, residuo, n * sizeof(real_t));     //copia o resíduo para o vetor b do sistema
    iteracoes = eliminacaoGauss (SL_alterado, w_vetor_diferencas, &tGauss);
    // 4º passo: obter nova solução x´ = x + w e testar critério de parada
    for (int j = 0; j < n; j++)
      x[j] = x[j] + w_vetor_diferencas[j];
    norma = normaL2Residuo(SL, x);
    if (norma == -1)
      return (-i - 1);
  }
  if (i >= (MAXIT)){
    fprintf(stderr, "O método de Refinamento não convergiu porque a matriz não é diagonalmente dominante.\n");
    *tTotal = timestamp() - *tTotal;
    return -i;
  }
  else {
    *tTotal = timestamp() - *tTotal;
    return i;
  }  
}