<h1> Biblioteca para resolução de sistemas lineares </h1>
Implementação de uma biblioteca, em C, para resolução de sistemas lineares pelos métodos de Eliminação de Gauss, Gauss-Seidel e Refinamento. Trabalho desenvolvido para a disciplina de Introdução à Computação Científica, do curso de bacharelado em Informática Biomédica da UFPR.

A especificação do trabalho está disponível no arquivo: especificacao_CI1164_2022_1_sistemas_lineares.pdf

## Funções utilitárias 

Foram fornecidas pelo prof. Dr. Armando Delgado para uso neste uso livre neste exercício um conjunto de funções utilitárias, disponíveis via Gitlab C3SL no link https://gitlab.c3sl.ufpr.br/nicolui/ci1164-utils ou via linha de comando: git clone git@gitlab.c3sl.ufpr.br:nicolui/ci1164-utils.git

O módulo sislin.* contém funções para definir um sistema linear. O módulo utils.* contém a definição da função timestamp() usadas na implementação final. 

## Como compilar
No terminal, execute: 

```
make
```
## Como executar
No terminal, execute: 

```
./labSisLin
```

## Observações sobre os resultados obtidos (tempos, número de iterações, tamanho do SL)

A partir dos resultados obtidos, é possível estimar que o método de EGP é mais vantajoso que o método de Gauss-Seidel para resolução de sistemas lineares quando os sistemas são mal-condicionados, como os que possuem matrizes do tipo Hilbert. Isso porque o método de Gauss-Siedel não converge em sistemas que não são diagonalmente dominantes.

Já quanto ao tipo de alocação de matriz, a alocação do tipo pontVet - que aloca uma matriz como vetor de N ponteiros para um único vetor com N*N elementos - é muito mais eficiente em termos de tempo de execução do que a alocação do tipo pontPont - que aloca uma matriz  como  vetor de  N  ponteiros.
