#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>

void generateGlider(int *data, int rows, int columns) {
  data[0 * columns + 1] = 1;
  data[1 * columns + 2] = 1;
  data[2 * columns + 0] = 1;
  data[2 * columns + 1] = 1;
  data[2 * columns + 2] = 1;
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

int *countRowsForEachProc(int rows, int columns, int amountOfProcs) {
  int *arrayAmountOfRowsPerOneProc = new int[amountOfProcs];
  int minAmountOfRowsPerOneProc = rows / amountOfProcs;

  int remainder = rows % amountOfProcs;

  for (int currProc = 0; currProc < amountOfProcs; currProc++) {
    arrayAmountOfRowsPerOneProc[currProc] = minAmountOfRowsPerOneProc;

    if (currProc < remainder) {
      arrayAmountOfRowsPerOneProc[currProc]++;
    }
  }
  return arrayAmountOfRowsPerOneProc;
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

  MPI_Init(&argc, &argv);

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);


   int *currentGen = new int[rowsAmount * columnsAmount]();

  if (rankOfCurrProc == 0) {
    std::fstream inputFile;
    inputFile.open("begin.txt",
                   std::ios::in | std::ios::out | std::ios::trunc);
    if (!inputFile) {
      std::cout << "Can't open input file\n";
      return 0;
    }

    generateGlider(currentGen, rowsAmount, columnsAmount);

    printMatrixToFile(currentGen, rowsAmount, columnsAmount,
                      inputFile);
  }

  int *realAmountOfRowsPerOneProc =
      countRowsForEachProc(rowsAmount, columnsAmount, amountOfProcs);

  if (rankOfCurrProc == 0) {
    for (int i = 0; i < amountOfProcs; i++) {
      std::cout << i << ": " << realAmountOfRowsPerOneProc[i]
                << std::endl;
    }
  }

  MPI_Finalize();
  return 0;
}
