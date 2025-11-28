#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

int main() {

  //declaração e inicialização de variáveis
  //
  //
  int N = 10;            //numero de celulas
  double A = 0.1;        //área da seção transversal
  double d = 1.0;        //distancia entre centróides
  double u = 0.01;        //velocidade do escoamento
  double Tin = 100.0;      //temperatura inlet
  double Tout = 200.0;     //temperatura outlet
  double sBar = 1000.0;  //taxa volumetrica de geracao de calor
  double k = 100.0;      //condutividade térmica
  double rho = 1.0;      //densidade
  double cp = 1000.0;    //calor especifico
 
  double D = 0.0; D=k*A/d;
  double F = 0.0; F=rho*cp*u*A;
  double V = 0.0; V=A*d;

  double i_aL = 0.0; i_aL = D + 0.5*F;
  double i_aR = 0.0; i_aR = D - 0.5*F;
  double i_Sp = 0.0;
  double i_aP = 0.0; i_aP = i_aL + i_aR + F - F - i_Sp;
  double i_Su = 0.0; i_Su = sBar * V;
  
  double l_aL = 0.0;
  double l_aR = 0.0; l_aR = D - 0.5*F;
  double l_Sp = 0.0; l_Sp = -2 * D - F;
  double l_aP = 0.0; l_aP = l_aL + l_aR + F - F - l_Sp;
  double l_Su = 0.0; l_Su = Tin * (2 * D + F) + sBar * V;

  double r_aL = 0.0; r_aL = D + 0.5*F;
  double r_aR = 0.0;
  double r_Sp = 0.0; r_Sp = -2 * D + F;
  double r_aP = 0.0; r_aP = r_aL + r_aR + F - F - r_Sp;
  double r_Su = 0.0; r_Su = Tout * (2 * D - F) + sBar * V;

  //-------------------------------------------------------------------------------------------------
  //Criar o vetore A (coeficientes) e o arquivo A.dat
  //-------------------------------------------------------------------------------------------------
  vector <double> coef_A(N*N,0);
  coef_A[0] = l_aP;
  coef_A[1] = -l_aR;
  coef_A[coef_A.size()-2] = -r_aL;
  coef_A[coef_A.size()-1] = r_aP;

  ofstream outputFileA("A.dat"); //Opção de construtor que cria o objeto e o arquivo na mesma instrução
  
  if (outputFileA.is_open()) {
    for (const double& num : coef_A) {outputFileA << num << "\n";} //Escreve cada número seguido por um newline
    outputFileA.close(); //Fecha o stream
    cout << "Arquivo A.dat (termos fonte) criado/atualizado com sucesso." << "\n";
  } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}
  

  //--------------------------------------------------------------------------------------------------
  //Criar o vetor B (termos fonte) e o arquivo B.dat
  //--------------------------------------------------------------------------------------------------
  vector <double> coef_B(N);
  coef_B[0] = l_Su;  //primeiro elemento do vetor (elemento de face)
  coef_B[N-1] = r_Su;  //ultimo elemento do vetor (elemento de face)
  for (size_t i=1; i <= N-2; i++) {coef_B[i] = i_Su;} //elementos internos

  ofstream outputFileB("B.dat"); //Opção de construtor que cria o objeto e o arquivo na mesma instrução
  
  if (outputFileB.is_open()) {
    for (const double& num : coef_B) {outputFileB << num << "\n";} //Escreve cada número seguido por um newline
    outputFileB.close(); //Fecha o stream
    cout << "Arquivo B.dat (termos fonte) criado/atualizado com sucesso." << "\n";
  } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}

 return 0;
}
