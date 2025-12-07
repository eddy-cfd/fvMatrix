<p align="center">
L2C - Soluções em Computação Científica <br><br>
Curso:<br>
Do Cálculo à Simulação Computacional - Fundamentos de Métodos Numéricos com Aplicações<br><br>
Professor Rafael Gabler Gontijo<br><br><br>
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
<br>
<p align="center">
Revisão do Método dos Volumes Finitos (FVM)
</p>
<br>
Considere o escoamento de ar através de um duto unidimensional (1D), de comprimento L [m], em regime estacionário, com velocidade U [m/s] e temperaturas fixadas na entrada e saída do duto, Ti [K] e To [K], respectivamente, ilustrado abaixo. Determinar a termperatura do ar ao longo do duto.
<br>
<br>

![Alt text](images/figura1.png)

O perfil de temperatura ao longo do duto pode ser determinado analiticamente usando-se a equação 4. Pode-se também obter uma solução numérica, usando-se o método dos volumes finitos (FVM), através da discretização da geometria do duto em volumes de controle, onde se realiza a integração da equação 2 em cada um dos volumes discretizados.
<br>
<br>

![Alt text](images/figura2.png)

A integração da equação 2 é feita através do Teorema de Gauss, que converte uma integral de volume em uma integral de superfície (a superfície em torno do volume, que delimita o volume). Como neste caso a malha é 1D, a integral de superfície deve ser resolvida considerando-se a superfície esquerda (subindice "e") e a superfície direita (subíndice "d"), de cada volume de controle. A integração pelo método da quadratura de Gauss permite expressar a integral de superfície (superfície de cada face) através de um valor situado no centróide de cada face, possibilitando que a integral seja convertida em um somatório das faces. Observa-se, no entanto, que a operação de integração do volume de controle resultou em uma expressão para o cálculo dos valores da variável de interesse nas faces dos volumes de controle. Porém, necessita-se de uma expressão para cálculo da variável de interesse no centróide do volume de controle. Obtem-se a expressão na forma desejada, equação 6,  utilizando esquemas de interpolação convenientes  para escrever os valores no centróide dos volumes de controle, em função dos valores nas suas faces. Neste caso, utilizou-se a interpolação linear (diferenças centrais), tanto para o termo advectivo, quanto para o termo difusivo. Convém ressaltar que a forma de interpolação não é única, podendo ser adotadas várias formas, em função, principalmente, do tipo de equação diferencial sendo resolvida.
<br>
<br>

![Alt text](images/figura3.png)

<br>
<br>
<p align="center">
Componentes de software
</p>
<br>
1 - fvTestCase.cpp <br><br>
Este programa gera os vetores "A" e "B", do sistema Ax=B, a partir das equações discretizadas do caso ilustrado acima.
<br>
<br>
2 - gaussSiedel.cpp <br><br>
Este programa resolve o sistema Ax=B, usando o método Gauss-Siedel. A convergência do método é garatida se, e somente se, a matriz "A" for diagonal dominante. A matriz "A" não precisa ser simétrica. Os dados de entrada são os vetores "A.dat", "B.dat" e X0.dat, gerados pelo usuário, ou pelo programa fvTestCase.cpp. O programa fornece o vetor solução no arquivo "X.dat".
<br>
<br>
3 - gradConjugado.cpp <br><br>
Este programa resolve o sistema Ax=B, usando o método do gradiente conjungado. Para tanto, é necessário que a matriz "A" seja simétrica e positivo definida. Os dados de entrada são os vetores "A.dat", "B.dat" e X0.dat, gerados pelo usuário, ou pelo programa fvTestCase.cpp. O programa fornece o vetor solução no arquivo "X.dat". 
