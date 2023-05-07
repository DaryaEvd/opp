#include <fstream>
#include <iostream>

struct MyMatrix {
  int *data = nullptr;
  size_t colmns;
  size_t rows;

  MyMatrix(size_t rows, size_t colmns) {
    data = new int[colmns * rows]();
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

// x - rows, y - columns
int countNeighbors(MyMatrix matrix, int xMatr, int yMatr) {
  int sum = 0;

  for (int i = -1; i < 2; ++i) {   // rows
    for (int j = -1; j < 2; ++j) { // columns

      int currRow = (xMatr + i + matrix.rows) % matrix.rows;
      int currColumn = (yMatr + j + matrix.colmns) % matrix.colmns;

      sum += matrix.data[currRow * matrix.colmns + currColumn];
    }
  }

  sum -= matrix.data[xMatr * matrix.colmns + yMatr];
  std::cout << "SUMMA is: " << sum << std::endl << std::endl;

  return sum;
}

void computeNextGeneration(MyMatrix matrix, MyMatrix nextMatrix) {
  int clown = 0;
  for (int i = 0; i < matrix.rows; ++i) {
    for (int j = 0; j < matrix.colmns; ++j) {

      int state = matrix.data[i * matrix.colmns + j];

      std::cout << "CLOWN: " << clown << ": " << state << std::endl;

      int neighborsAmount = countNeighbors(matrix, i, j);

      if (state == 0 && neighborsAmount == 3) {
        nextMatrix.data[i * nextMatrix.colmns + j] = 1;
      } else if (state == 1 &&
                 (neighborsAmount < 2 || neighborsAmount > 3)) {
        nextMatrix.data[i * nextMatrix.colmns + j] = 0;
      } else {
        nextMatrix.data[i * nextMatrix.colmns + j] =
            matrix.data[i * nextMatrix.colmns + j] = state;
      }

      clown++;
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

  std::fstream outputFile;
  outputFile.open("end.txt");
  if (!outputFile) {
    std::cout << "Can't open output file\n";
    freeMatrix(matrixStart);
    inputFile.close();
    return 0;
  }

  MyMatrix nextGen = MyMatrix(rowsAmount, columnsAmount);
  computeNextGeneration(matrixStart, nextGen);

  printMatrixToFile(nextGen, outputFile);

  inputFile.close();
  outputFile.close();

  freeMatrix(matrixStart);

  return 0;
}
