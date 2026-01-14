#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
using namespace std;

//prototipos
void residuo(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
void multMatrizVetor(vector<double>&, vector<double>&, vector<double>&);
double prodEscVetor(vector<double>&, vector<double>&);
double normaVetor(vector<double>&);
void multEscVetor(double, vector<double>&, vector<double>&);
void somaVetor(vector<double>&, vector<double>&, vector<double>&);
void matrizTransposta(vector<double>&, vector<double>&);

int main() {
  //variáveis de dados de entrada
  vector<double> A;  //vetor (matriz) de coeficientes
  vector<double> B;  //termos fonte
  vector<double> r;  //afastamento dos centroides em relação a origem.
  vector<double> X;  //vetor solução
  ifstream dataFromFile;
  ofstream dataToFile;
  double numero;
  double ti, tf, t_total;

  cout << "\n------------------------------ Método Gradiente Biconjugado  ------------------------------\n";

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
  double BETA_np1 = 0.0;
  double BETANum_np1 = 0.0;
  double BETADen_np1 = 1.0;
  vector<double> At(N * N, 0);
  vector<double> R(N, 0);
  vector<double> Rhat(N, 0);
  vector<double> R_np1(N, 0);
  vector<double> Rhat_np1(N, 0);
  vector<double> D(N, 0);
  vector<double> Dhat(N, 0);
  vector<double> AxD(N, 0);
  vector<double> AtxDhat(N, 0);
  vector<double> ALFAxD(N, 0);
  vector<double> negALFAxAxD(N, 0);
  vector<double> negALFAxAtxDhat(N, 0);
  vector<double> ALFAxDhat(N, 0);
  vector<double> BETA_np1xD(N, 0);
  vector<double> BETA_np1xDhat(N, 0);
  vector<double> X_np1(N, 0);  //np1 denota n+1

// ----------------- cpu_time(ti) -----------------
  ti = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
  
  //inicialização do vetor resíduo(R) de direção(D)  para inicio das nações
  residuo(A, X, B, R);
  Rhat = R;
  D = R;
  Dhat = R;
  matrizTransposta(A, At);

  //loop
  do {
    n++;
    //cálculo do alfa
    ALFANum = prodEscVetor(Rhat, R);
    multMatrizVetor(A, D, AxD);
    ALFADen = prodEscVetor(Dhat, AxD);
    ALFA = ALFANum / ALFADen;

    //cálculo do Xn+1
    multEscVetor(ALFA, D, ALFAxD);
    somaVetor(X, ALFAxD, X_np1);

    //cálculo do novo resíduo R_np1;
    multEscVetor(-ALFA, AxD, negALFAxAxD);
    somaVetor(R, negALFAxAxD, R_np1);

    //cálculo do novo resíduo Rhat_np1;
    multMatrizVetor(At, Dhat, AtxDhat);
    multEscVetor(-ALFA, AtxDhat, negALFAxAtxDhat);
    somaVetor(Rhat, negALFAxAtxDhat, Rhat_np1);

    //cálculo do beta
    BETANum_np1 = prodEscVetor(Rhat_np1, R_np1);
    BETADen_np1 = prodEscVetor(Rhat, R);
    BETA_np1 = BETANum_np1 / BETADen_np1;

    //atualiza valor  do D;
    multEscVetor(BETA_np1, D, BETA_np1xD);
    somaVetor(R_np1, BETA_np1xD, D);

    //atualiza valor  do Dhat;
    multEscVetor(BETA_np1, Dhat, BETA_np1xDhat);
    somaVetor(Rhat_np1, BETA_np1xDhat, Dhat);

    //atualiza valores para próximo loop
    X = X_np1;
    R = R_np1;
    Rhat = Rhat_np1;
    normaR = normaVetor(R);

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

void matrizTransposta(vector<double>& mtx_A, vector<double>& mtx_At) {
  //--------------------------------------------------------------------------------
  //---------------------Calcula a transposta de uma matrix-------------------------
  //--------------------------------------------------------------------------------
  int N = sqrt(mtx_A.size());

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i != j) {
        mtx_At[i * N + j] = mtx_A[j * N + i];
      } else {
        mtx_At[i * N + j] = mtx_A[i * N + j];
      }
    }
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
  //-------------------- resolve a equação aV  -------------------------------------
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
