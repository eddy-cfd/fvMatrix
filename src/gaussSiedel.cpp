#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

//prototipos
void residuo(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
double normaVetor(vector<double>&);


int main(){
  //variáveis de dados de entrada
  vector <double> A;
  vector <double> B;
  vector <double> Xn;
  ifstream dataFromFile;
  double numero;

  //popular vetor A
  dataFromFile.open("A.dat", ios::in);
  if(dataFromFile.is_open())
    while (dataFromFile >> numero) A.push_back(numero);
  else 
    {cout << "Impossível abrir arquivos de dados\n\n"; return 1;}
  dataFromFile.close();

  //popular vetor B
  dataFromFile.open("B.dat", ios::in);
  if(dataFromFile.is_open())
    while (dataFromFile >> numero) B.push_back(numero);
  else 
    {cout << "Impossível abrir arquivos de dados\n\n"; return 1;}
  dataFromFile.close();

  //popular vetor Xn com os valores iniciais do arquivo X0.dat
  dataFromFile.open("X0.dat", ios::in);
  if(dataFromFile.is_open())
    while (dataFromFile >> numero) Xn.push_back(numero);
  else 
    {cout << "Impossível abrir arquivos de dados\n\n"; return 1;}
  dataFromFile.close();


  //verificação da consistência dos dados
  cout << "Método Gauss Siedel\n\n";
  if(Xn.size() == B.size() || A.size() == pow(B.size(),2)) 
    cout << "Dados carregados......[OK]\nDados consistentes....[OK]\n\n";
  else  
  {cout << "Problema com arquivos de dados ou inconsistência de dados\n\n"; return 1;}

  //verificação se matriz A é diagonal dominante
  //-----------FALTA IMPLEMENTAR---------------------------
 
  //variáveis internas do algoritmo
  int N = Xn.size();
  int iter = 0;
  double somaN = 0.0;
  double somaN_ = 0.0;
  double normaRn = 0.0;
  vector <double> Rn(N,0.0);
 
  //algoritmo de solução do sistema
  do{
    iter = iter+1;
    for(int i=0; i < N; i++){
      somaN=0.0;
      somaN_=0.0; 
      if (i-1 >= 0) for(int j=0; j <= i-1; j++) somaN += A[i*N+j]*Xn[j];
      for(int j=i+1; j < N; j++) somaN_ += A[i*N+j]*Xn[j]; 
      Xn[i]=(1/A[i*(N+1)])*(B[i]-somaN_-somaN);
    }
    residuo(A,Xn,B,Rn);
    normaRn = normaVetor(Rn);
  } while (normaRn > 0.001);
  
  //escrever na tela a solução do sistema
  cout << "Solução do sistema:" << "\n";
  //for(int i=0; i < N; i++)    cout << "X[" << i << "] = " << Xn[i] << "\n";      
  cout << "Norma do resíduo = " << normaRn << "\n";
  cout << "Número de iterações = " << iter << "\n";

  //--------------------------------------------------------------------------------------------------
  //Criar o arquivo X.dat com a solução do sistema
  //--------------------------------------------------------------------------------------------------
  ofstream outputFile("X.dat");

  if (outputFile.is_open()) {
    for (const double& num : Xn) {outputFile << num << "\n";} //Escreve cada número seguido por um newline
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


