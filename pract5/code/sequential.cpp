#include <fstream>
#include <iostream>

void generateGlider(int *data, int rows, int columns) {
  data[0 * columns + 1] = 1;
  data[1 * columns + 2] = 1;
  data[2 * columns + 0] = 1;
  data[2 * columns + 1] = 1;
  data[2 * columns + 2] = 1;
}

void initMatrixWithZeroes(int *data, int rows, int columns) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      data[i * columns + j] = 0;
    }
  }
}

void printMatrixToFile(int *data, int rows, int columns,
                       std::fstream &file) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      file << data[i * columns + j] << " ";
    }
    file << "\n";
  }
}

void copyMatrix(int *dest, int *src, int rows, int columns) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      dest[i * columns + j] = src[i * columns + j];
    }
  }
}

// x - rows, y - columns
int countNeighbors(int *data, int rows, int columns, int xMatr,
                   int yMatr) {
  int sum = 0;

  for (int i = -1; i < 2; ++i) {   // rows
    for (int j = -1; j < 2; ++j) { // columns

      int currRow = (xMatr + i + rows) % rows;
      int currColumn = (yMatr + j + columns) % columns;

      sum += data[currRow * columns + currColumn];
    }
  }

  sum -= data[xMatr * columns + yMatr]; // cell itself
  // std::cout << "SUMMA is: " << sum << std::endl << std::endl;

  return sum;
}

void computeNextGeneration(int *oldData, int *nextData, int rows,
                           int columns) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {

      int state = oldData[i * columns + j];

      int neighborsAmount =
          countNeighbors(oldData, rows, columns, i, j);

      if (state == 0 && neighborsAmount == 3) {
        nextData[i * columns + j] = 1;
      } else if (state == 1 &&
                 (neighborsAmount < 2 || neighborsAmount > 3)) {
        nextData[i * columns + j] = 0;
      } else {
        nextData[i * columns + j] = oldData[i * columns + j] = state;
      }
    }
  }
}

bool equalsToPrevEvolution(int *prevMatr, int *currMatrix, int rows,
                           int columns) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      if (currMatrix[i * columns + j] != prevMatr[i * columns + j]) {
        return false;
      }
    }
  }
  return true;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Bad amount of arguments!\n"
              << "Enter rows amount, then enter columns amount.\n"
              << std::endl;
    return 0;
  }

  const int rowsAmount = atoi(argv[1]);
  const int columnsAmount = atoi(argv[2]);

  std::fstream inputFile;
  inputFile.open("start.txt",
                 std::ios::in | std::ios::out | std::ios::trunc);
  if (!inputFile) {
    std::cout << "Can't open input file\n";
    return 0;
  }

  int *currentGen = new int[rowsAmount * columnsAmount]();

  generateGlider(currentGen, rowsAmount, columnsAmount);

  printMatrixToFile(currentGen, rowsAmount, columnsAmount, inputFile);

  std::fstream outputFile;
  outputFile.open("end.txt", std::ios::out | std::ios::trunc);
  if (!outputFile) {
    std::cout << "Can't open output file\n";
    delete[] currentGen;
    inputFile.close();
    return 0;
  }

  outputFile << "start matrix -------------" << std::endl;
  printMatrixToFile(currentGen, rowsAmount, columnsAmount,
                    outputFile);

  int *nextGen = new int[rowsAmount * columnsAmount]();

  const long maxIterations = 1000; 

  int **prevEvolution = new int *[maxIterations]();

  int iterCurr = 0;
  bool repeated = false;

  while (iterCurr < maxIterations && !repeated) {
    prevEvolution[iterCurr] = new int[rowsAmount * columnsAmount]();

    copyMatrix(prevEvolution[iterCurr], currentGen, rowsAmount,
               columnsAmount);

    computeNextGeneration(currentGen, nextGen, rowsAmount,
                          columnsAmount);
    outputFile << "after iter: " << iterCurr <<
    "-----------------"
               << std::endl;
    printMatrixToFile(nextGen, rowsAmount, columnsAmount,
    outputFile);

    copyMatrix(currentGen, nextGen, rowsAmount, columnsAmount);

    for (int i = iterCurr; i > -1; --i) {
      if (equalsToPrevEvolution(prevEvolution[i], nextGen, rowsAmount,
                                columnsAmount)) {

        // std::cout << "iter " << iterCurr << ": equals" <<
        // std::endl; // in real: iter + 1
        repeated = true;
        break;
      }
    }

    iterCurr++;
    // std::cout << "curr iteration: " << iterCurr << std::endl;
  }

  if(repeated) {
    std::cout << "First repeat of cell automat stage after: "
      << iterCurr << " iterCurr" << std::endl;
  }
  else {
    std::cout << "Finished after iteration: " << iterCurr << std::endl;
  }

  inputFile.close();
  outputFile.close();

  delete[] currentGen;
  delete[] nextGen;

  for (int i = 0; i < maxIterations; i++) {
    delete[] prevEvolution[i];
  }
  delete[] prevEvolution;

  return 0;
}