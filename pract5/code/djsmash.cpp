#include <fstream>
#include <iostream>

void generateGlider(int *data, int rows, int columns) {
  data[0 * columns + 1] = 1;
  data[1 * columns + 2] = 1;
  data[2 * columns + 0] = 1;
  data[2 * columns + 1] = 1;
  data[2 * columns + 2] = 1;

  // data[3 * columns + 4] = 1;
  // data[0 * columns + 4] = 1;
}

void generateBlinker(int *data, int rows, int columns) {
  data[1 * columns + 2] = 1;
  data[2 * columns + 2] = 1;
  data[3 * columns + 2] = 1;
}

void generateBlock(int *data, int rows, int columns) {
  data[1 * columns + 2] = 1;
  data[1 * columns + 3] = 1;

  data[2 * columns + 2] = 1;
  data[2 * columns + 3] = 1;
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
  inputFile.open("start.txt",
                 std::ios::in | std::ios::out | std::ios::trunc);
  if (!inputFile) {
    std::cout << "Can't open input file\n";
    return 0;
  }
  
  int *currentGen = new int[rowsAmount * columnsAmount]();

  if (modeToWork == 'r') {
    std::cout << "You're in random mode" << std::endl;
  }

  else if (modeToWork == 'g') {
    std::cout << "You're in Glider mode" << std::endl;
    generateGlider(currentGen, rowsAmount, columnsAmount);
  }

  // else if (modeToWork == 'b') {
  //   std::cout << "You're in Blinker mode" << std::endl;
  //   generateBlinker(currentGen, rowsAmount, columnsAmount);
  // }

  // else if (modeToWork == 's') {
  //   std::cout << "You're in Sqare (Block) mode" << std::endl;
  //   generateBlock(currentGen, rowsAmount, columnsAmount);
  // } else {
  //   std::cout << "bad input of mode" << std::endl;
  //   return 0;
  // }

  printMatrixToFile(currentGen, rowsAmount, columnsAmount, inputFile);
 
  // for(int i = 0; i < rowsAmount; i++) {
  //   for(int j = 0; j < columnsAmount; j++) {
  //     std::cout << "cell " << i * columnsAmount + j << " has ";
  //     int neighbCount = countNeighbors(currentGen, rowsAmount, columnsAmount, i, j);
  //     std::cout << neighbCount << std::endl;
  //   }
  // }

  std::fstream outputFile;
  outputFile.open("end.txt", std::ios::out | std::ios::trunc);
  if (!outputFile) {
    std::cout << "Can't open output file\n";
    delete[] currentGen;
    inputFile.close();
    return 0;
  }
  
  outputFile << "start matrix -------------" << std::endl;
  printMatrixToFile(currentGen, rowsAmount, columnsAmount, outputFile);

  int *nextGen = new int[rowsAmount * columnsAmount]();
  // printMatrixToFile(nextGen, rowsAmount, columnsAmount, outputFile);

  const int maxIterations = 30; // change it :)

  int **prevEvolution = new int *[maxIterations]();
  // for (int i = 0; i < maxIterations; i++) {
  //   prevEvolution[i] = new int[rowsAmount * columnsAmount]();
  // }

  int iterations = 0;

  while (iterations < maxIterations) {
    prevEvolution[iterations] = new int[rowsAmount * columnsAmount]();

    copyMatrix(prevEvolution[iterations], currentGen, rowsAmount,
               columnsAmount);

    computeNextGeneration(currentGen, nextGen, rowsAmount,
                          columnsAmount);
    outputFile << "after iter: " << iterations << "-----------------"
               << std::endl;
    printMatrixToFile(nextGen, rowsAmount, columnsAmount, outputFile);

    copyMatrix(currentGen, nextGen,rowsAmount, columnsAmount);

    for (int i = iterations; i >= 0; --i) {
      if (equalsToPrevEvolution(prevEvolution[i], nextGen, rowsAmount,
                                columnsAmount)) {
        if(i == 0) {
          continue;
        }
        std::cout << "iter " << iterations << ": equals" << std::endl;
        break;
      }
    }

    iterations++;
    // std::cout << "curr iteration: " << iterations << std::endl;
    
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
