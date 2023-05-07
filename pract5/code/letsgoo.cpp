#include <fstream>
#include <iostream>

struct MyMatrix {
  int *data = nullptr;
  size_t colmns;
  size_t rows;

  MyMatrix(size_t rows, size_t colmns) {
    data = new int[colmns * rows];
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

void copyMatrix(MyMatrix oldMatr, MyMatrix newMatr) {
  for (size_t i = 0; i < oldMatr.rows; ++i) {
    for (size_t j = 0; j < oldMatr.colmns; ++j) {
      newMatr.data[i * newMatr.colmns + j] =
          oldMatr.data[i * oldMatr.colmns + j];
    }
  }
}

int countNeighbors(MyMatrix matrix, int x, int y) {
  int sum = 0;

  for (size_t i = -1; i < 2; ++i) {
    for (size_t j = -1; j < 2; ++j) {
  
      int column = (x + j + matrix.colmns) % matrix.colmns;
      int row = (y + i + matrix.rows) % matrix.rows;

      sum += matrix.data[row * matrix.colmns + column];
    }
  }

  sum -= matrix.data[x * matrix.colmns + y];

  return sum;
}

void computeNextGeneration(MyMatrix matrix, MyMatrix nextMatrix) {
  for (size_t i = 0; i < matrix.rows; ++i) {
    for (size_t j = 0; j < matrix.colmns; ++j) {
      
      int state = matrix.data[i * matrix.colmns + j];

      int neighborsAmount = countNeighbors(matrix, i, j);

      if(state == 0 && neighborsAmount == 3) {
        nextMatrix.data[i * nextMatrix.colmns + j] = 1;
      }
      else if(state == 1 && (neighborsAmount < 2 || neighborsAmount > 3)) {
        nextMatrix.data[i * nextMatrix.colmns + j] = 0;
      }
      else {
        nextMatrix.data[i * nextMatrix.colmns + j] =
            matrix.data[i * nextMatrix.colmns + j] = 0;
      }


    }
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

  MyMatrix matrixCopy = MyMatrix(rowsAmount, columnsAmount);
  copyMatrix(matrixStart, matrixCopy);

  std::fstream interFile;
  interFile.open("inter.txt");
  if (!interFile) {
    std::cout << "Can't open inter file\n";
    return 0;
  }

  computeNextGeneration(matrixStart);

  inputFile.close();

  freeMatrix(matrixStart);

  return 0;
}
