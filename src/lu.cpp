#include <cmath>
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
  if (X.size() == B.size() || A.size() == pow(B.size(), 2)) {
    cout << "Dados carregados......[OK]\nDados consistentes....[OK]\n\n";
  } else {
    cout << "Problema com arquivos de dados ou inconsistência de dados\n\n";
    return 1;
  }

  //variáveis internas do algoritmo
  int N = X.size();
  int k = 0;
  vector<double> L(N * N, 0);  //vetor que armazena os fatores
  vector<double> U(N * N, 0);  //vetor em que se opera a eliminação de Gauss

  //Antes de iniciar a eliminação de Gauss, montar L de forma a ficar uma matriz identidade
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) {
        L[i * N + j] = 1;
      } else {
        L[i * N + j] = 0;
      }
    }
  }

  //montar L com os fatores
  for (int j = 0; j < N; j++) {
    for (int i = 1; i < N; i++)
      if (j < i) { L[i * N + j] = A[i * N + j] / A[j * N + j]; }
  }

  //Aplicar a eleiminação de Gauss no vetor A

  /*
  do {
    cout << L[k] << "\n";
    k++;
  } while (k < N * N);
*/
  //escrever na tela a solução do sistema
  cout << "Solução do sistema:" << "\n";

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
