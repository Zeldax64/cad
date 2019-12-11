# CAD - Computação de Alto Desempenho

Este repositório contém as atividades da disciplina de Computação de Alto Desempenho (CAD) ministrada na Universidade Federal do Ceará (UFC). O repositório contém duas atividades, a primeira corresponde a multiplicação de matrizes utilizando múltiplos processos e MPI (mpi-mm). A segunda atividade é um algoritmo para executar o ordenamento em um vetor usando (mpi-bucket).

Cada pasta contém um makefile para compilar a aplicação. Duas flags de compilação são utilizadas para melhorar a performance da aplicação mpi-mm, são elas:

* -march=native: compila o programa utilizando a arquitetura do processador do computador como alvo permitindo um melhor tempo de execução. 
* -funroll-loops: permite que o GCC otimize loops melhorando o tempo de execução mas, em contrapartida, aumentando o tamanho do ELF.

Melhores resultados podem ser obtidos utilizando a flag `-ffast-math`. Esta opção encontra-se comentada no Makefile pois ela melhora a performance a custo de uma menor precisão de operações com float.