#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>

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

  std::cout << "rows amount: " << rowsAmount << std::endl;

  MPI_Init(&argc, &argv);

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  // std::cout << "procs amount: " << amountOfProcs << std::endl;

  int minAmountOfRowsPerOneProc;

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
