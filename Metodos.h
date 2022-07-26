#ifndef __METODOS_H__
#define __METODOS_H__

// Parâmetros para teste de convergência
#define MAXIT 100           // número máximo de iterações para métodos iterativos
#define ERRO 1.0e-6         // tolerância para critérios de parada em métodos iterativos
//#define P_INF (1.0/0.0)     // infinto positivo 
//#define N_INF (-1.0/0.0)    // infinito negativo
//#define NaN (0.0/0.0)       // not a number

// Calcula a normaL2 do resíduo
real_t normaL2Residuo(SistLinear_t *SL, real_t *x);

// Método da Eliminação de Gauss
int eliminacaoGauss(SistLinear_t *SL, real_t *x, double *tTotal);

// Método de Refinamento
int refinamento(SistLinear_t *SL, real_t *x, real_t erro, double *tTotal);

// Método de Gauss-Seidel
int gaussSeidel(SistLinear_t *SL, real_t *x, real_t erro, double *tTotal);

// Função que copia o sistema linear
SistLinear_t *copiaSistemaLinear(SistLinear_t *SL);

// Função que encontra o maior valor da coluna para definir o pivô
// e devolve o índice da linha correspondente
int encontraMax(SistLinear_t *SL, int indice);

// Função que troca a linha pela linha com o pivô
void trocaLinha(SistLinear_t *SL, int i, int indice_pivo);

// Calcula a normaL2 do resíduo
int calculaResiduo(SistLinear_t *SL, real_t *x, real_t *residuo);

#endif // __METODOS_H__

