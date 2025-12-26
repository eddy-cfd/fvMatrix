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
O objetivo deste trabalho foi a implementação em C++ dos algorítmos para solução iterativa de sistemas lineares, usando os seguintes métodos:<br><br>
Gauss-Siedel;<br>
Gradiente Conjugado (GC);<br>
Gradiente Conjugado pré-condicionado (PGC);<br>
Gradiente Bi-Conjugado (GBiC);<br>
Gradiente Bi-Conjugado pré-condicionado (PGBiC);<br>
Gradiente Bi-Conjugado pré-condicionado estabilizado (PGBiC-Stab).
<br>
<br>
Para complementar o trabalho, foi desenvolvido um programa para gerar um sistema de equações, na forma A.X = B, baseado na discretização de um problema físico, utilizando o Método dos Volumes Finitos (FVM). Para se entender o funcionamento do programa, será feito primeiramente uma revisão do FVM, já no contexto do enunciado do problema físico que se pretende resolver.
<br>
<br>
<br>
<br>
<p align="center">
Revisão do Método dos Volumes Finitos (FVM)
</p>
<br>
Considere o escoamento de ar através de um duto unidimensional (1D), de comprimento L [m], em regime estacionário, com velocidade U [m/s], temperaturas fixadas na entrada e saída do duto, Ti [K] e To [K], respectivamente, com geração interna de calor, conforme ilustrado abaixo. Determinar a termperatura do ar ao longo do duto.
<br>
<br>
O perfil de temperatura ao longo do duto pode ser determinado analiticamente usando-se a equação 4, válida para o caso em que não há geração interna de calor no duto.
<br>
<br>

![Alt text](images/figura1.png)

Pode-se também obter uma solução numérica, usando-se o método dos volumes finitos (FVM), através da discretização da geometria do duto em volumes de controle, onde se resolve a equação 3 em cada um dos volumes discretizados.
<br>
<br>

![Alt text](images/figura2.png)

Para obter a solução numérica do problema, primeiramente aplica-se o Teorema de Gauss na equação 3 (ou equação 5, já escrita utilizando a variável genérica phi), que converte uma integral de volume em uma integral de superfície (a superfície em torno do volume, que delimita o volume). Como neste caso a malha é 1D, a integral de superfície deve ser resolvida considerando-se a superfície esquerda (subíndice "e") e a superfície direita (subíndice "d"), de cada volume de controle. A integração pelo método da quadratura de Gauss permite expressar a integral de superfície (superfície de cada face) através de um valor situado no centróide de cada face, possibilitando que a integral seja convertida em um somatório das faces. Observa-se, no entanto, que a operação de integração do volume de controle resultou em uma expressão para o cálculo dos valores da variável de interesse nas faces dos volumes de controle. Porém, necessita-se de uma expressão para cálculo da variável de interesse que utilize os valores de centróide do volume de controle. Obtem-se a expressão na forma desejada, equação 6,  utilizando esquemas de interpolação convenientes  para escrever os valores nas faces dos volumes de controle, em função dos valores de centróide dos volumes de controle. Neste caso, utilizou-se a interpolação linear (diferenças centrais), tanto para o termo advectivo, quanto para o termo difusivo. Convém ressaltar que a forma de interpolação não é única, podendo ser adotadas várias formas, em função, principalmente, do tipo de equação diferencial sendo resolvida.
<br>
<br>

![Alt text](images/figura3.png)

A dedução acima vale para volumes de controle no interior do domínio. No caso de volumes de controle na fronteira do dominio, ou seja, volumes de controle cuja uma, ou mais faces, formam o contorno do domínio, a equação 6 toma formas diferentes. No caso do volume de controle do lado esquerdo do duto, os valores das condições de contorno são aplicados diretamente à face esquerda do volume, cujo afastamento em relação ao centróide do volume de controle é igual a metade do comprimento do volume de controle. Analogamente, no volume de controle do lado direito do duto, os valores das condições de contorno são aplicados diretamente à face direita do volume de controle.
<br>
<br>

![Alt text](images/figura4.png)

![Alt text](images/figura5.png)

Basicamente, cada volume de controle tem uma equação algébrica associada com subindice P, que o relaciona com os volumes vizinhos, correspondendo a uma linha de um sistema de equações algébricas lineares. O volume de controle em questão tem o subindice "P, o volume à esquerda, subindice "E" e o volume a direita, subindice "D". As equações 6, 7 e 8 foram deduzidas e já colocadas na forma padrão (equação 9), o que facilita o reconhecimento dos coeficientes para montagem das matrizes. 
<br>
<br>

![Alt text](images/figura6.png)

![Alt text](images/figura7.png)
<br>
<br>
Depois de obter as equações algébricas para o cálculo da variável de interesse do escoamento, é necessário arranjar as equações na forma de um sistema linear, no formato Ax=B, montar os vetores A e B, para que então o sistema possa ser resolvido utilizando algum método numérico implementado neste repositório.
<br>
<br>
![Alt text](images/figura8.png)

<br>
<br>
<br>
<p align="center">
Componentes de software
</p>
<br>
1 - fvTestCase.cpp <br><br>
Sintaxe:<br>
"fvTestCase [quantidade de volumes de controle] [velocidade do escoamento]"<br><br>
Exemplo:<br>
"fvTestCase 100 0.8" --> 100 elementos de malha (volumes de controle) e velocidade do escoamento de 0.8 m/s.
<br>
<br>
Este programa gera os vetores "A" e "B", do sistema A.X = B e o vetor "X0" de solução inicial, a partir das equações discretizadas do caso ilustrado acima. Esses vetores são usados pelos programas de solução de sistemas lineares deste repositório. Além disso, é gerada a solução analítica, válida para o caso sem geração interna de calor.
<br><br>
Propriedades do fluido e condições de contorno são informadas no código fonte. Quando há alteração, é necessário compilar novamente.
<br><br>
Cada vez que o fvTestCase é executado, são gerados/editados os seguintes arquivos:
<br><br>
-->  A.dat - Vetor A (matriz de coeficientes).<br>
-->  B.dat - Vetor B (termos fonte).<br>
-->  X0.dat - Vetor X0 (solução inicial).<br>
-->  r.dat - Coordenadas (afastamentos) dos centróides dos volumes de controle, ao longo do duto.<br>
-->  solAnalitica.dat - Solução analítica do caso. Está no formato apropriado para ser usado no GNUPlot.
<br>
<br>
<br>
<br>
2 - gaussSiedel.cpp <br><br>
Este programa resolve o sistema A.X = B, usando o método Gauss-Siedel. A convergência do método é garatida se, e somente se, a matriz "A" for diagonal dominante. A matriz "A" não precisa ser simétrica. Os dados de entrada são os vetores "A.dat", "B.dat" e X0.dat, gerados pelo usuário, ou pelo programa fvTestCase.cpp. O programa fornece o vetor solução no arquivo "X.dat".
<br>
<br>
<br>
<br>
3 - gradConjugado.cpp <br><br>
Este programa resolve o sistema A.X = B, usando o método do gradiente conjungado. Para tanto, é necessário que a matriz "A" seja simétrica e positivo definida. Os dados de entrada são os vetores "A.dat", "B.dat" e X0.dat, gerados pelo usuário, ou pelo programa fvTestCase.cpp. O programa fornece o vetor solução no arquivo "X.dat". 
