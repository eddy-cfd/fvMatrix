#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

//prototipos
void residuo(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
void multMatrizVetor(vector<double>&, vector<double>&, vector<double>&);
double prodEscVetor(vector<double>&, vector<double>&);
double normaVetor(vector<double>&);
void multEscVetor(double, vector<double>&, vector<double>&);
void somaVetor(vector<double>&, vector<double>&, vector<double>&);


int main(){
  //variáveis de dados de entrada
  vector <double> A;
  vector <double> B;
  vector <double> r;  
  vector <double> X;
  ifstream dataFromFile;
  ofstream dataToFile;
  double numero;

  //carregar vetor A
  dataFromFile.open("A.dat", ios::in);
  if(dataFromFile.is_open())
    while (dataFromFile >> numero) A.push_back(numero);
  else 
    {cout << "Impossível abrir arquivos de dados\n\n"; return 1;}
  dataFromFile.close();

  //carregar vetor B
  dataFromFile.open("B.dat", ios::in);
  if(dataFromFile.is_open())
    while (dataFromFile >> numero) B.push_back(numero);
  else 
    {cout << "Impossível abrir arquivos de dados\n\n"; return 1;}
  dataFromFile.close();

  //carregar vetor X com os valores iniciais do arquivo X0.dat
  dataFromFile.open("X0.dat", ios::in);
  if(dataFromFile.is_open())
    while (dataFromFile >> numero) X.push_back(numero);
  else 
    {cout << "Impossível abrir arquivos de dados\n\n"; return 1;}  


  //verificação da consistência dos dados
  cout << "Método grandiente conjugado\n\n";
  if(X.size() == B.size() || A.size() == pow(B.size(),2)) 
    cout << "Dados carregados......[OK]\nDados consistentes....[OK]\n\n";
  else  
  {cout << "Problema com arquivos de dados ou inconsistência de dados\n\n"; return 1;}

  //variáveis internas do algoritmo
  int N = X.size();
  int iter = 0;
  double normaR = 0.0;
  double ALFA = 0.0;
  double ALFANum = 0.0;
  double ALFADen = 1.0;
  double BETA = 0.0;
  double BETANum = 0.0;
  double BETADen = 1.0;
  vector <double> R(N,0);
  vector <double> R1(N,0);
  vector <double> D(N,0);
  vector <double> ADn(N,0);
  vector <double> ALFAnDn(N,0);
  vector <double> BETAnDn(N,0);
  vector <double> Xn1(N,0); 

  //inicialização do vetor resíduo(R) de direção(D)  para inicio das iterações
  residuo(A,X,B,R);   
  D=R;
  
  //loop 
  do{
    iter = iter+1;
    //cálculo do alfa
    ALFAnNum = prodEscVetor(D, R);
    multMatrizVetor(A, D, ADn);
    ALFAnDen = prodEscVetor(D, ADn);  
    ALFAn = ALFAnNum / ALFAnDen;

    //cálculo do Xn+1
    multEscVetor(ALFAn, D, ALFAnDn);
    somaVetor(X, ALFAnDn, Xn1);
  
    //cálculo do resíduo de Xn+1
    residuo(A,Xn1,B,R1);
    normaR1 = normaVetor(Rn1);
    
    //cálculo do beta
    BETAnNum = prodEscVetor(R1, Rn1);
    BETAnDen = prodEscVetor(Rn, Rn);
    BETAn = BETAnNum / BETAnDen;

    //atualiza valor  do D;
    multEscVetor(BETAn, D, BETAnDn);
    somaVetor(R1, BETAnDn, D);

    //atualiza valores para próximo loop
    X=Xn1;
    Rn=R1;

    } while (normaR1 > 0.001);

    cout << "Solução do sistema de equações:\n";
    //for(int i=0; i < N; i++)    cout << "X[" << i << "] = " << Xn[i] << "\n";      
    cout << "Norma do resíduo = " << normaR1 << "\n";
    cout << "Número de iterações = " << iter << "\n"; 
 
    //--------------------------------------------------------------------------------------------------
  //Criar o arquivo RX.dat com a solução do sistema para plotagem
  //--------------------------------------------------------------------------------------------------
 
  dataToFile.open ("RX.dat");
  if (dataToFile.is_open()) {
    for (int i=0; i<N; i++){  
      dataToFile << R[i] << " ";
      dataToFile << X[i] << "\n"; //Escreve cada número seguido por um newline
    }
    dataToFile.close(); //Fecha o stream
    cout << "Arquivo RX.dat (R:posição, X:solução do sistema) criado/atualizado com sucesso." << "\n";
  } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}
    
    //--------------------------------------------------------------------------------------------------
    //Criar o arquivo X.dat com a solução do sistema
    //--------------------------------------------------------------------------------------------------
  ofstream outputFile("X.dat");

  if (outputFile.is_open()) {
    for (const double& num : X) {outputFile << num << "\n";} //Escreve cada número seguido por um newline
    outputFile.close(); //Fecha o stream
    cout << "Arquivo X.dat (solução do sistema) criado/atualizado com sucesso." << "\n";
  } else {cerr << "Erro: Impossível criar/abrir arquivo.\n";}
    return 0;
}

void residuo(vector<double>& mtx_A, vector<double>& vet_X, vector<double>& vet_B, vector<double>& vet_R){
//--------------------------------------------------------------------------------
//-------------------- resolve a equação R=B-AX ----------------------------------
//--------------------------------------------------------------------------------
  int N = vet_X.size();

  for(int i = 0; i < N; i++){
    double sum = 0;
    for(int j = 0; j < N; j++)  sum += mtx_A[i * N + j] * vet_X[j];
    vet_R[i] = vet_B[i] - sum;
  }
}

void multMatrizVetor(vector<double>& mtx_A, vector<double>& vet_X, vector<double>& vet_AX){
//--------------------------------------------------------------------------------
//-------------------- resolve a equação R=AX ------------------------------------
//--------------------------------------------------------------------------------
  int N = vet_X.size();

  for(int i = 0; i < N; i++){
    double sum = 0;
    for(int j = 0; j < N; j++)  sum += mtx_A[i * N + j] * vet_X[j];
    vet_AX[i] = sum;
  }
}

double prodEscVetor(vector<double>& vet_U, vector<double>& vet_V){
//--------------------------------------------------------------------------------
//-------------------- resolve a equação U.V  ------------------------------------
//--------------------------------------------------------------------------------
  int N = vet_U.size();
  double sum = 0;

  for(int i = 0; i < N; i++){
    sum += vet_U[i]*vet_V[i];
  }
  return sum;
}
double normaVetor(vector<double>& vet_U){
//--------------------------------------------------------------------------------
//-------------------- calcula norma de um vetor  --------------------------------
//--------------------------------------------------------------------------------
  int N = vet_U.size();
  double sum = 0.0;
  double norm = 0.0;

  for(int i = 0; i < N; i++){
    sum += vet_U[i]*vet_U[i];
  }
  norm=sqrt(sum);
  return norm;
}

void multEscVetor(double esc_a, vector<double>& vet_V, vector<double>& vet_aV){
//--------------------------------------------------------------------------------
//-------------------- resolve a equação aV  ------------------------------------
//--------------------------------------------------------------------------------
  int N = vet_V.size();

  for(int i = 0; i < N; i++)  vet_aV[i] = esc_a * vet_V[i];
}

void somaVetor(vector<double>& vet_U, vector<double>& vet_V, vector<double>& vet_UmaisV){
//--------------------------------------------------------------------------------
//-------------------- resolve a equação U+V  ------------------------------------
//--------------------------------------------------------------------------------
  int N = vet_V.size();

  for(int i = 0; i < N; i++)  vet_UmaisV[i] = vet_U[i] + vet_V[i];
}
