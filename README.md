<p align="center">
L2C - Soluções em Computação Científica <br>
Do cálculo à simulação computacional<br>
Fundamentos de Métodos Numéricos com Aplicações<br><br>
Trabalho de fim de curso
</p>
<br>
<br>
Componentes:
<br>
<br>
1 - fvTestCase.cpp <br><br>
Gera os vetores "A" e "B", do sistema Ax=B, a partir das equações discretizadas do caso ilustrado abaixo.
<br>
<br>

![Alt text](images/mesh.png)

A discretização foi feita utilizando o método dos volumes finitos (FVM), com esquema de interpolação linear (diferenças centrais), tanto para o termo advectivo, quanto para o termo difusivo.
<br>
<br>
2 - gaussSiedel.cpp <br><br>
Resolve o sistema Ax=B, usando o método Gauss-Siedel. Os dados de entrada são os vetores "A.dat" e "B.dat", gerados pelo usuário, ou pelo programa fvTestCase.cpp. O vetor "X0" deve necessáriamente ser fornecido pelo usuário, tomando o cuidado de ser compatível com as matrizes "A" e "B". O programa fornece a matriz solução no arquivo "X.dat".
