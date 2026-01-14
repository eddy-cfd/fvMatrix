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
  ifstream dataFromFile;
  ofstream dataToFile;
  double numero;
  double ti, tf, t_total;

  cout << "\n------------------------------ Método de Thomas ------------------------------\n";

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

  //carregar vetor X com os valores iniciais do arquivo X0.dat:q
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
  if (X.size() == B.size() || A.size() == pow(B.size(), 2)) {
    cout << "Dados carregados......[OK]\nDados consistentes....[OK]\n\n";
  } else {
    cout << "Problema com arquivos de dados ou inconsistência de dados\n\n";
    return 1;
  }

  //variáveis internas do algoritmo
  int N = X.size();
  int k = 0;
  vector<double> E;
  vector<double> F;
  vector<double> G;

  // ----------------- cpu_time(ti) -----------------
  ti = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;

  //Separar os elementos tridiagonais em vetores
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) {
        F.push_back(A[i * N + j]);
      } else if (j == i + 1 && i < N - 1) {
        G.push_back(A[i * N + j]);
      } else if (j == i - 1 && i > 0) {
        E.push_back(A[i * N + j]);
      }
    }
  }

  //Decomposição LU
  k = 0;
  do {
    E[k] = E[k] / F[k];
    F[k + 1] = F[k + 1] - E[k] * G[k];
    k++;
  } while (k < N - 1);

  //Substitução progressiva - cálculo de {d}
  k = 1;
  do {
    B[k] = B[k] - E[k - 1] * B[k - 1];
    k++;
  } while (k < N);

  //Substituição regressiva - cálculo de {x}
  X[N - 1] = B[N - 1] / F[N - 1];
  k = N - 2;
  do {
    X[k] = (B[k] - G[k] * X[k + 1]) / F[k];
    k--;
  } while (k >= 0);

  // ----------------- cpu_time(tf) -----------------
  tf = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
  t_total = tf - ti;
  //
  //escrever na tela a solução do sistema
  cout << "Solução do sistema:" << "\n";
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
