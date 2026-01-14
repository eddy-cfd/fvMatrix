#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

int main(int argc, char* argv[]) {
  
  //declaração e inicialização de variáveis
  //
  double N = 10.0;                      //numero de células
  double u = 0.0;                       //velocidade do escoamento  
  double qdot = 1000.0;                 //taxa volumetrica de geracao de calor
  
  cout << "\n---------------------------------- fvMatrix ----------------------------------";
  cout << "\n-------------------- Gerador de Sistemas de Equações FVM ---------------------"; 
  cout << "\n---------------------- Duto de ar / Regime estacinário -----------------------\n";

  //Altera o número de elementos e a velocidade do escoamento, caso usuário forneça os valores na execução do programa
  if (argc > 3 && argc < 5) {
    N = stod(argv[1]);
    u = stod(argv[2]);
    qdot = stod(argv[3]);
    cout << "\nN = " << N << " elementos" << "\n";   
    cout << "u = " << u << " m/s" << "\n";
    cout << "qdot = " << qdot << " W/m3" << "\n\n"; 
  } 
  else{
    cout << "\nParametros default\n"; 
    cout << "N = " << N << " elementos" << "\n";   
    cout << "u = " << u << " m/s" << "\n";
    cout << "qdot = " << qdot << " W/m^3" << "\n\n";
  }

  double L = 5.0;                       //comprimento do duto
  double d = 0.0; d=L/N;                //tamanho do elemento e distancia entre centróides
  double Tin = 100.0;                   //temperatura inlet
  double Tout = 200.0;                  //temperatura outlet
  double k = 100.0;                     //condutividade térmica
  double rho = 1.0;                     //densidade
  double cp = 1000.0;                   //calor especifico
  double alfa = 0.0; alfa=k/(rho*cp);   //difusividade térmica

  double i_aE = 0.0; i_aE = (alfa/d) + 0.5*u;
  double i_aD = 0.0; i_aD = (alfa/d) - 0.5*u;
  double i_Fp = 0.0;
  double i_cFp = 0.0;
  double i_Fphi = 0.0; i_Fphi = (qdot*d)/(rho*cp);
  double i_aP = 0.0; i_aP = i_aE + i_aD + i_Fp;
  double i_bp = 0.0; i_bp = i_Fphi + (i_cFp * i_Fp);
   
  double e_aE = 0.0;
  double e_aD = 0.0; e_aD = (alfa / d) - (0.5 * u);
  double e_Fp = 0.0; e_Fp = u + (2 * alfa / d);
  double e_cFp = 0.0; e_cFp = Tin;
  double e_Fphi = 0.0; e_Fphi = (qdot*d)/(rho*cp);
  double e_aP = 0.0; e_aP = e_aE + e_aD + e_Fp;
  double e_bp = 0.0; e_bp = e_Fphi + (e_cFp * e_Fp);
 
  double d_aE = 0.0; d_aE = (alfa / d) + (0.5 * u);
  double d_aD = 0.0;
  double d_Fp = 0.0; d_Fp = (2 * alfa / d) - u;
  double d_cFp = 0.0; d_cFp = Tout;
  double d_Fphi = 0.0; d_Fphi = (qdot*d)/(rho*cp);
  double d_aP = 0.0; d_aP = d_aE + d_aD + d_Fp;
  double d_bp = 0.0; d_bp = d_Fphi + (d_cFp * d_Fp);


  //declaração de objetos para I/O
  ofstream outputFile;

  //--------------------------------------------------------------------------------------------------
  //Criar o vetor r (coordenadas dos centróides dos volumes de controle) e o arquivo r.dat
  //--------------------------------------------------------------------------------------------------
  vector <double> r(N);
  r[0] = d/2;  //primeiro elemento do vetor (elemento de face)
  for (size_t i=1; i < N; i++) {r[i] = r[i-1]+d;} //elementos internos

  outputFile.open ("r.dat"); //Abre o stream
  
  if (outputFile.is_open()) {
    for (const double& num : r) {outputFile << num << "\n";} //Escreve cada número seguido por um newline
    outputFile.close(); //Fecha o stream
    cout << "Arquivo r.dat (coordenadas dos centróides) criado/atualizado com sucesso." << "\n";
  } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}


  //-------------------------------------------------------------------------------------------------
  //Criar o vetor A (coeficientes) e o arquivo A.dat
  //-------------------------------------------------------------------------------------------------
  vector <double> coef_A(N*N,0);

  //Alocação dos termos da dianonal principal
  coef_A[0] = e_aP; //elemento da esquerda (inlet)
  int j = 1;
  for (int i = 1; i < N-1; i++) {coef_A[i*N+j] = i_aP; j++;} //elementos internos
  coef_A[coef_A.size()-1] = d_aP; //elelento da direita (outlet)
  
  //Alocação dos termos da diagonal superior
  j = 2;
  coef_A[1] = -e_aD; //elemento a direita do elemento de inlet
  for (int i = 1; i < N-1; i++) {coef_A[i*N+j] = -i_aD; j++;} //elementos internos

  //Alocação dos termos da diagonal inferior
  j = 0;
  for (int i = 1; i < N-1; i++) {coef_A[i*N+j] = -i_aE; j++;} //elementos internos
  coef_A[coef_A.size()-2] = -d_aE; //elemento a esquerda do elemento de outlet

  outputFile.open("A.dat"); //Abre o stream
  
  if (outputFile.is_open()) {
    for (const double& num : coef_A) {outputFile << num << "\n";} //Escreve cada número seguido por um newline
    outputFile.close(); //Fecha o stream
    cout << "Arquivo A.dat (coeficientes) criado/atualizado com sucesso." << "\n";
  } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}
  

  //--------------------------------------------------------------------------------------------------
  //Criar o vetor B (termos fonte) e o arquivo B.dat
  //--------------------------------------------------------------------------------------------------
  vector <double> coef_B(N);
  coef_B[0] = e_bp;  //primeiro elemento do vetor (elemento de face)
  coef_B[N-1] = d_bp;  //ultimo elemento do vetor (elemento de face)
  for (size_t i=1; i <= N-2; i++) {coef_B[i] = i_bp;} //elementos internos

  outputFile.open("B.dat"); //Abre o stream 
  
  if (outputFile.is_open()) {
    for (const double& num : coef_B) {outputFile << num << "\n";} //Escreve cada número seguido por um newline
    outputFile.close(); //Fecha o stream
    cout << "Arquivo B.dat (termos fonte) criado/atualizado com sucesso." << "\n";
  } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}

  //--------------------------------------------------------------------------------------------------
  //Criar o vetor X0 (estimativa inicial) e o arquivo X0.dat
  //--------------------------------------------------------------------------------------------------
  vector <double> coef_X0(N);
  for (size_t i=0; i <= N-1; i++) {coef_X0[i] = 0.5*(Tin+Tout);} //elementos internos

  outputFile.open("X0.dat"); //Abre o stream
  
  if (outputFile.is_open()) {
    for (const double& num : coef_X0) {outputFile << num << "\n";} //Escreve cada número seguido por um newline
    outputFile.close(); //Fecha o stream
    cout << "Arquivo X0.dat (estimativa inicial) criado/atualizado com sucesso." << "\n";
  } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}

  //--------------------------------------------------------------------------------------------------
  //Criar arquivo com a solução analítica
  //--------------------------------------------------------------------------------------------------
  double x = 0.0;
  double T = 0.0;

  if(qdot==0 && u != 0){
    outputFile.open("solAnalitica.dat"); //Abre o stream
    if (outputFile.is_open()){
      for (x = 0; x <= L; x+=0.01) { 
        T = ((Tout-Tin) * (exp(u*x/alfa)-1)/((exp(u*L/alfa)-1)))+Tin;
        outputFile << x << " " << T << "\n";
      }
      outputFile.close(); //Fecha o stream
      cout << "Arquivo solAnalitica.dat (solução analítica) criado/atualizado com sucesso." << "\n";
    } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}
  } 
  else if(u == 0) {
     outputFile.open("solAnalitica.dat"); //Abre o stream
    if (outputFile.is_open()){
      for (x = 0; x <= L; x+=0.01) { 
        T = Tin + (x/L)*(Tout-Tin)+(0.5*x*qdot/k)*(L-x);
        outputFile << x << " " << T << "\n";
      }
      outputFile.close(); //Fecha o stream
      cout << "Arquivo solAnalitica.dat (solução analítica) criado/atualizado com sucesso." << "\n";
    } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}
  } 
  else cout << "Não existe solução analítica para esta combinação de dados de entrada." << "\n";    

 return 0;
}
