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
  vector <double> X0;
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

  //popular vetor X0
  dataFromFile.open("X0.dat", ios::in);
  if(dataFromFile.is_open())
    while (dataFromFile >> numero) X0.push_back(numero);
  else 
    {cout << "Impossível abrir arquivos de dados\n\n"; return 1;}
  dataFromFile.close();


  //verificação da consistência dos dados
  cout << "Método grandiente conjugado\n\n";
  if(X0.size() == B.size() || A.size() == pow(B.size(),2)) 
    cout << "Dados carregados......[OK]\nDados consistentes....[OK]\n\n";
  else  
  {cout << "Problema com arquivos de dados ou inconsistência de dados\n\n"; return 1;}

  //variáveis internas do algoritmo
  int N = X0.size();
  int iter = 0;
  double normaRn1 = 0.0;
  double ALFAn = 0.0;
  double ALFAnNum = 0.0;
  double ALFAnDen = 1.0;
  double BETAn = 0.0;
  double BETAnNum = 0.0;
  double BETAnDen = 1.0;
  vector <double> Rn(N,0);
  vector <double> Rn1(N,0);
  vector <double> Dn(N,0);
  vector <double> ADn(N,0);
  vector <double> ALFAnDn(N,0);
  vector <double> BETAnDn(N,0);
  vector <double> Xn=X0;
  vector <double> Xn1(N,0); 

  //inicialização do vetor resíduo(Rn) de direção(Dn)  para inicio das iterações
  residuo(A,X0,B,Rn);   
  Dn=Rn;
  
  //loop 
  do{
    iter = iter+1;
    //cálculo do alfa
    ALFAnNum = prodEscVetor(Dn, Rn);
    multMatrizVetor(A, Dn, ADn);
    ALFAnDen = prodEscVetor(Dn, ADn);  
    ALFAn = ALFAnNum / ALFAnDen;

    //cálculo do Xn+1
    multEscVetor(ALFAn, Dn, ALFAnDn);
    somaVetor(Xn, ALFAnDn, Xn1);
  
    //cálculo do resíduo de Xn+1
    residuo(A,Xn1,B,Rn1);
    normaRn1 = normaVetor(Rn1);
    
    //cálculo do beta
    BETAnNum = prodEscVetor(Rn1, Rn1);
    BETAnDen = prodEscVetor(Rn, Rn);
    BETAn = BETAnNum / BETAnDen;

    //atualiza valor  do Dn;
    multEscVetor(BETAn, Dn, BETAnDn);
    somaVetor(Rn1, BETAnDn, Dn);

    //atualiza valores para próximo loop
    Xn=Xn1;
    Rn=Rn1;

    } while (normaRn1 > 0.001);

    cout << "Solução do sistema de equações:\n";
    for(int i=0; i < N; i++)    cout << "X[" << i << "] = " << Xn[i] << "\n";      
    cout << "Norma do resíduo = " << normaRn1 << "\n";
    cout << "Número de iterações = " << iter << "\n"; 
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
