#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

//prototipos
void residuo(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
void multMatrizVetor(vector<double>&, vector<double>&, vector<double>&);
double prodEscVetor(vector<double>&, vector<double>&);
double normaVetor(vector<double>&);
void multEscVetor(double, vector<double>&, vector<double>&);
void somaVetor(vector<double>&, vector<double>&, vector<double>&);

int main() {
  //variáveis de dados de entrada
  vector<double> A;
  vector<double> B;
  vector<double> r;
  vector<double> X;
  ifstream dataFromFile;
  ofstream dataToFile;
  double numero;
  double ti, tf, t_total;

  cout << "\n------------------------------ Método Gradiente Conjugado  ------------------------------\n";

  //carregar vetor A
  dataFromFile.open("A.dat", ios::in);
  if (dataFromFile.is_open())
    while (dataFromFile >> numero) A.push_back(numero);
  else {
    cout << "Impossível abrir arquivos de dados\n\n";
    return 1;
  }
  dataFromFile.close();

  //carregar vetor B
  dataFromFile.open("B.dat", ios::in);
  if (dataFromFile.is_open())
    while (dataFromFile >> numero) B.push_back(numero);
  else {
    cout << "Impossível abrir arquivos de dados\n\n";
    return 1;
  }
  dataFromFile.close();

  //carregar vetor X com os valores iniciais do arquivo X0.dat
  dataFromFile.open("X0.dat", ios::in);
  if (dataFromFile.is_open())
    while (dataFromFile >> numero) X.push_back(numero);
  else {
    cout << "Impossível abrir arquivos de dados\n\n";
    return 1;
  }
  dataFromFile.close();

  //carregar vetor r com coordenadas r.dat
  //usado quando está resolvendo caso gerado pelo fvMatrix
  dataFromFile.open("r.dat", ios::in);
  if (dataFromFile.is_open())
    while (dataFromFile >> numero) r.push_back(numero);
  else {
    cout << "Impossível abrir arquivos de dados\n\n";
    return 1;
  }
  dataFromFile.close();

  //verificação da consistência dos dados
  if (X.size() == B.size() || A.size() == pow(B.size(), 2))
    cout << "Dados carregados......[OK]\nDados consistentes....[OK]\n\n";
  else {
    cout << "Problema com arquivos de dados ou inconsistência de dados\n\n";
    return 1;
  }

  //variáveis internas do algoritmo
  int N = X.size();
  int n = 0;
  double normaR = 0.0;
  double ALFA = 0.0;
  double ALFANum = 0.0;
  double ALFADen = 1.0;
  double BETA = 0.0;
  double BETANum = 0.0;
  double BETADen = 1.0;
  vector<double> R(N, 0);
  vector<double> R_np1(N, 0);
  vector<double> D(N, 0);
  vector<double> AxD(N, 0);
  vector<double> ALFAxD(N, 0);
  vector<double> BETAxD(N, 0);
  vector<double> X_np1(N, 0);  //np1 denota n+1

  // ----------------- cpu_time(ti) -----------------
  ti = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;

  //inicialização do vetor resíduo(R) de direção(D)  para inicio das nações
  residuo(A, X, B, R);
  D = R;

  //loop
  do {
    n++;
    //cálculo do alfa
    ALFANum = prodEscVetor(D, R);
    multMatrizVetor(A, D, AxD);
    ALFADen = prodEscVetor(D, AxD);
    ALFA = ALFANum / ALFADen;

    //cálculo do Xn+1
    multEscVetor(ALFA, D, ALFAxD);
    somaVetor(X, ALFAxD, X_np1);

    //cálculo do resíduo de Xn+1
    residuo(A, X_np1, B, R_np1);
    normaR = normaVetor(R_np1);

    //cálculo do beta
    BETANum = prodEscVetor(R_np1, R_np1);
    BETADen = prodEscVetor(R, R);
    BETA = BETANum / BETADen;

    //atualiza valor  do D;
    multEscVetor(BETA, D, BETAxD);
    somaVetor(R_np1, BETAxD, D);

    //atualiza valores para próximo loop
    X = X_np1;
    R = R_np1;

  } while (normaR > 0.001);

  // ----------------- cpu_time(tf) -----------------
  tf = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
  t_total = tf - ti;

  cout << "Solução do sistema de equações:\n";
  cout << "Norma do resíduo = " << normaR << "\n";
  cout << "Número de iterações = " << n << "\n";
  cout << "Tempo de CPU = " << t_total << "s\n";

  //--------------------------------------------------------------------------------------------------
  //Criar o arquivo X.dat com a solução do sistema
  //--------------------------------------------------------------------------------------------------
  dataToFile.open("X.dat");
  if (dataToFile.is_open()) {
    for (const double& num : X) { dataToFile << num << "\n"; }  //Escreve cada número seguido por um newline
    dataToFile.close();                                         //Fecha o stream
    cout << "Arquivo X.dat (solução do sistema) criado/atualizado com sucesso." << "\n";
  } else {
    cerr << "Erro: Impossível criar/abrir arquivo.\n";
  }

  //--------------------------------------------------------------------------------------------------
  //Criar o arquivo solNumerica.dat com a solução do sistema para plotagem
  //--------------------------------------------------------------------------------------------------
  dataToFile.open("solNumerica.dat");
  if (dataToFile.is_open()) {
    for (int i = 0; i < N; i++) {
      dataToFile << r[i] << " ";
      dataToFile << X[i] << "\n";  //Escreve cada número seguido por um newline
    }
    dataToFile.close();  //Fecha o stream
    cout << "Arquivo solNumerica.dat (r:posição, X:solução do sistema) criado/atualizado com sucesso." << "\n\n";
  } else {
    cerr << "Erro: Impossível criar/abrir arquivo.\n";
  }

  return 0;
}

void residuo(vector<double>& mtx_A, vector<double>& vet_X, vector<double>& vet_B, vector<double>& vet_R) {
  //--------------------------------------------------------------------------------
  //-------------------- resolve a equação R=B-AX ----------------------------------
  //--------------------------------------------------------------------------------
  int N = vet_X.size();

  for (int i = 0; i < N; i++) {
    double sum = 0;
    for (int j = 0; j < N; j++) sum += mtx_A[i * N + j] * vet_X[j];
    vet_R[i] = vet_B[i] - sum;
  }
}

void multMatrizVetor(vector<double>& mtx_A, vector<double>& vet_X, vector<double>& vet_AX) {
  //--------------------------------------------------------------------------------
  //-------------------- resolve a equação R=AX ------------------------------------
  //--------------------------------------------------------------------------------
  int N = vet_X.size();

  for (int i = 0; i < N; i++) {
    double sum = 0;
    for (int j = 0; j < N; j++) sum += mtx_A[i * N + j] * vet_X[j];
    vet_AX[i] = sum;
  }
}

double prodEscVetor(vector<double>& vet_U, vector<double>& vet_V) {
  //--------------------------------------------------------------------------------
  //-------------------- resolve a equação U.V  ------------------------------------
  //--------------------------------------------------------------------------------
  int N = vet_U.size();
  double sum = 0;

  for (int i = 0; i < N; i++) { sum += vet_U[i] * vet_V[i]; }
  return sum;
}
double normaVetor(vector<double>& vet_U) {
  //--------------------------------------------------------------------------------
  //-------------------- calcula norma de um vetor  --------------------------------
  //--------------------------------------------------------------------------------
  int N = vet_U.size();
  double sum = 0.0;
  double norm = 0.0;

  for (int i = 0; i < N; i++) { sum += vet_U[i] * vet_U[i]; }
  norm = sqrt(sum);
  return norm;
}

void multEscVetor(double esc_a, vector<double>& vet_V, vector<double>& vet_aV) {
  //--------------------------------------------------------------------------------
  //-------------------- resolve a equação aV  ------------------------------------
  //--------------------------------------------------------------------------------
  int N = vet_V.size();

  for (int i = 0; i < N; i++) vet_aV[i] = esc_a * vet_V[i];
}

void somaVetor(vector<double>& vet_U, vector<double>& vet_V, vector<double>& vet_UmaisV) {
  //--------------------------------------------------------------------------------
  //-------------------- resolve a equação U+V  ------------------------------------
  //--------------------------------------------------------------------------------
  int N = vet_V.size();

  for (int i = 0; i < N; i++) vet_UmaisV[i] = vet_U[i] + vet_V[i];
}
