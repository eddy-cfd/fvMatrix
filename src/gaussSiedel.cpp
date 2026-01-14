#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

//prototipos
void residuo(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
double normaVetor(vector<double>&);

int main(int argc, char* argv[]) {
  //variáveis de dados de entrada/saída
  vector<double> A;
  vector<double> B;
  vector<double> X;
  vector<double> r;
  double lambda = 1;
  ifstream dataFromFile;
  ofstream dataToFile;
  double numero;
  double ti, tf, t_total;

  cout << "\n------------------------------ Método Gauss-Siedel ------------------------------\n";

  if (argc > 1) {
    lambda = stod(argv[1]);
    cout << "\n"
         << "lambda = " << lambda << "\n";
  } else {
    cout << "\n"
         << "lambda (default) = " << lambda << "\n";
  }

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

  //variáveis internas do algoritmo
  int N = X.size();
  int k = 0;                  //conta o número de iterações
  double soma_k = 0.0;        //soma relativa a solução na iteração corrente
  double soma_kmenos1 = 0.0;  //soma relativa a solução na iteração anterior
  double somaTermosND = 0.0;  //soma termos não-diagonais
  vector<double> R(N);        //resíduo
  double normaR = 0.0;        //norma do resíduo, usada como critério de parada do algoritmo gaussSiedel

  //verificação da consistência dos dados
  if (X.size() == B.size() || A.size() == pow(B.size(), 2)) {
    cout << "Dados carregados......[OK]\nDados consistentes....[OK]\n\n";
  } else {
    cout << "Problema com arquivos de dados ou inconsistência de dados\n\n";
    return 1;
  }

  //verificação da dominancia diagonal
  for (int i = 0; i < N; i++) {
    somaTermosND = 0.0;
    for (int j = 0; j < N; j++) { somaTermosND += abs(A[i * N + j]); }
    somaTermosND = somaTermosND - A[i * (N + 1)];
    if (somaTermosND > A[i * (N + 1)]) {
      cout << "Sistema não é diagonal dominante\n";
      return 0;
    }
  }

  // ----------------- cpu_time(ti) -----------------
  ti = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;

  //algoritmo de Gauss-Siedel
  do {
    k++;
    for (int i = 0; i < N; i++) {
      soma_k = 0.0;
      soma_kmenos1 = 0.0;
      if (i - 1 >= 0) {
        for (int j = 0; j <= i - 1; j++) { soma_k += A[i * N + j] * X[j]; }
      }
      for (int j = i + 1; j < N; j++) { soma_kmenos1 += A[i * N + j] * X[j]; }
      X[i] = (lambda * (1 / A[i * (N + 1)]) * (B[i] - soma_kmenos1 - soma_k)) + ((1 - lambda) * X[i]);
    }
    residuo(A, X, B, R);
    normaR = normaVetor(R);
  } while (normaR > 0.001);
  tf = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
  t_total = tf - ti;

  //escrever na tela a solução do sistema
  cout << "Solução do sistema:" << "\n";
  cout << "Norma do resíduo = " << normaR << "\n";
  cout << "Número de iterações = " << k << "\n";
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
