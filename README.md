<p align="center">
L2C - Soluções em Computação Científica <br><br>
Curso:<br>
Do Cálculo à Simulação Computacional - Fundamentos de Métodos Numéricos com Aplicações<br><br><br>
Trabalho de fim de curso
</p>
<br>
<br>
Resumo: 
<br>
<br>
O objetivo deste trabalho foi a implementação em C++ dos algorítmos para solução iterativa de sistemas lineares, usando os seguintes métodos:<br>
Gauss-Siedel;<br>
Gradiente Conjugado (GC);<br>
Gradiente Conjugado pré-condicionado (PGC);<br>
Gradiente Bi-Conjugado (GBiC);<br>
Gradiente Bi-Conjugado pré-condicionado (PGBiC);<br>
Gradiente Bi-Conjugado pré-condicionado estabilizado (PGBiC-Stab).
<br>
<br>
Componentes:
<br>
<br>
1 - fvTestCase.cpp <br><br>
Este programa gera os vetores "A" e "B", do sistema Ax=B, a partir das equações discretizadas do caso ilustrado abaixo.
<br>
<br>

![Alt text](images/mesh.png)

A discretização foi feita utilizando o método dos volumes finitos (FVM), com esquema de interpolação linear (diferenças centrais), tanto para o termo advectivo, quanto para o termo difusivo. A figura abaixo ilustra o processo.
<br>
<br>

![Alt text](images/eq1.png)

<br>
<br>
2 - gaussSiedel.cpp <br><br>
Este programa resolve o sistema Ax=B, usando o método Gauss-Siedel. Os dados de entrada são os vetores "A.dat" e "B.dat", gerados pelo usuário, ou pelo programa fvTestCase.cpp. O vetor "X0" deve necessáriamente ser fornecido pelo usuário, tomando o cuidado de ser compatível com as matrizes "A" e "B". O programa fornece a matriz solução no arquivo "X.dat".
