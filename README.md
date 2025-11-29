[mesh.html](https://github.com/user-attachments/files/23835759/mesh.html)Trabalho de final de curso

Componentes:

1 - fvTestCase.cpp

Gera os vetores "A" e "B", do sistema Ax=B, a partir das equações discretizadas de um caso, obtidas pelo método de volumes finitos. Descrever o caso, usando figuras para mostrar a geometria, malha 1D, equações na forma diferencial e equações discretizadas.

<img width="904" height="149" alt="image" src="https://github.com/user-attachments/assets/d892e4b4-09e4-4c71-999a-6020b3febca2" />

2 - gaussSiedel.cpp

Resolve o sistema Ax=B, usando o método Gauss-Siedel. Os dados de entrada são os vetores "A.dat" e "B.dat", gerados pelo usuário, ou pelo programa fvTestCase.cpp. O vetor "X0" deve necessáriamente ser fornecido pelo usuário, tomando o cuidado de ser compatível com as matrizes "A" e "B". O programa fornece a matriz solução no arquivo "X.dat".
