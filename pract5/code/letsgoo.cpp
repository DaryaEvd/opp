#include <fstream>
#include <iostream>

struct MyMatrix {
  double *data = nullptr;
  size_t colmns;
  size_t rows;

  MyMatrix(size_t rows, size_t colmns) {
    data = new double[colmns * rows];
    this->rows = rows;
    this->colmns = colmns;
  }
};

void generateGlider(MyMatrix matrix) {
  matrix.data[0 * matrix.colmns + 1] = 1;
  matrix.data[1 * matrix.colmns + 2] = 1;
  matrix.data[2 * matrix.colmns + 0] = 1;
  matrix.data[2 * matrix.colmns + 1] = 1;
  matrix.data[2 * matrix.colmns + 2] = 1;
}

void initMatrixWithZeroes(MyMatrix matrix) {
  for (size_t i = 0; i < matrix.rows; ++i) {
    for (size_t j = 0; j < matrix.colmns; ++j) {
      matrix.data[i * matrix.colmns + j] = 0;
    }
  }
}

void freeMatrix(MyMatrix matrix) { delete[] matrix.data; }

void printMatrixToFile(MyMatrix matrix, std::fstream &file) {
  for (size_t i = 0; i < matrix.rows; ++i) {
    for (size_t j = 0; j < matrix.colmns; ++j) {
      file << matrix.data[i * matrix.colmns + j] << " ";
    }
    file << "\n";
  }
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cout
        << "Bad amount of arguments!\n"
        << "Enter rows amount, then enter columns amount.\n"
        << "Also enter mode of work 'g' - glider or 'r' - random."
        << std::endl;
    return 0;
  }

  const int rowsAmount = atoi(argv[1]);
  const int columnsAmount = atoi(argv[2]);
  const char modeToWork = argv[3][0];

  std::fstream inputFile;
  inputFile.open("start.txt");
  if (!inputFile) {
    std::cout << "Can't open input file\n";
    return 0;
  }

  MyMatrix matrixStart = MyMatrix(rowsAmount, columnsAmount);
  initMatrixWithZeroes(matrixStart);

  if (modeToWork == 'r') {
    std::cout << "You're in random mode" << std::endl;

  }

  if (modeToWork == 'g') {
    std::cout << "You're in Glider mode" << std::endl;

    generateGlider(matrixStart);

  } 
  
  printMatrixToFile(matrixStart, inputFile);

  inputFile.close();

  freeMatrix(matrixStart);

  return 0;
}
