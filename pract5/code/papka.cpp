#include <fstream>
#include <iostream>
#include <cmath>
#include <mpi.h>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Bad amount of arguments!\n"
              << "Enter rows amount, then enter columns amount.\n"
              << std::endl;
    return 0;
  }

  const int rowsAmount = atoi(argv[1]);
  const int columnsAmount = atoi(argv[2]);

  std::cout << "rows amount: " <<  rowsAmount << std::endl;

  MPI_Init(&argc, &argv);

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  // std::cout << "procs amount: " << amountOfProcs << std::endl;

  int minAmountOfRowsPerOneProc;

  int realAmountOfRowsPerOneProc[amountOfProcs];

  if (rowsAmount % amountOfProcs == 0) {
    minAmountOfRowsPerOneProc = rowsAmount / amountOfProcs;

    if (rankOfCurrProc == 0) {
      std::cout << "min am: " << minAmountOfRowsPerOneProc
                << std::endl;
    }
  } else {
    minAmountOfRowsPerOneProc = rowsAmount / amountOfProcs ;

    // std::cout << "cloooown" << std::endl;

    for(int i = 0; i < amountOfProcs; i++) {
      realAmountOfRowsPerOneProc[i] = minAmountOfRowsPerOneProc;

      if(i < (rowsAmount % amountOfProcs)) {
        realAmountOfRowsPerOneProc[i]++;
      }
    }

    if (rankOfCurrProc == 0) {
      for (int i = 0; i < amountOfProcs; i++) {
        std::cout << "proc " << i
                  << ", realAmount: " << realAmountOfRowsPerOneProc[i]
                  << std::endl;
      }
    }
  }

  // std::cout << minAmountOfRowsPerOneProc << std::endl;

  MPI_Finalize();
  return 0;
}
