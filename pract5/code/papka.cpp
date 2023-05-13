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

int *countElemsNumInEachProc(int amountOfProcs, int rows,
                             int columns) {
  int *elemsNum = new int[amountOfProcs];

  int basicRowsCount = rows / amountOfProcs;
  int restRowsCount = rows % amountOfProcs;

  for (int i = 0; i < amountOfProcs; i++) {
    elemsNum[i] = basicRowsCount * columns;
    if (restRowsCount > 0) {
      elemsNum[i] += columns;
      --restRowsCount;
    }
  }
  return elemsNum;
}

int *countRowsInEachProcess(const int *elementsNumberArr,
                            int amountOfProcs, int columns) {
  int *rowsNum = new int[amountOfProcs];
  for (int i = 0; i < amountOfProcs; i++) {
    rowsNum[i] = elementsNumberArr[i] / columns;
  }
  return rowsNum;
}

int *createElemsOffsetArr(const int *elemsNum, int amountOfProcs) {
  int *elementsOffsetArray = new int[amountOfProcs];
  int elementsOffset = 0;
  for (int i = 0; i < amountOfProcs; i++) {
    elementsOffsetArray[i] = elementsOffset;
    elementsOffset += elemsNum[i];
  }
  return elementsOffsetArray;
}

int *createRowsOffsetArr(const int *elementsOffsetArray,
                         int amountOfProcs, int rows) {
  int *rowsOffsetArray = new int[amountOfProcs];
  for (int i = 0; i < amountOfProcs; i++) {
    rowsOffsetArray[i] = elementsOffsetArray[i] / rows;
  }

  return rowsOffsetArray;
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

  // amount of elems, handling each process
  int *elemsNum = countElemsNumInEachProc(amountOfProcs, rowsAmount,
                                          columnsAmount);

  // amount of rows, which each process handles
  int *rowsNum =
      countRowsInEachProcess(elemsNum, amountOfProcs, columnsAmount);

  // count, from which element data sends to each process
  int *elementsOffsetArray =
      createElemsOffsetArr(elemsNum, amountOfProcs);

  // count, from which row data sends to each process
  int *rowsOffsetArray = createRowsOffsetArr(
      elementsOffsetArray, amountOfProcs, rowsAmount);

  if (rankOfCurrProc == 0) {
    for (int i = 0; i < amountOfProcs; i++) {
      std::cout << i << ": " << rowsNum[i] << std::endl;
    }
  }

  int *partMatrGen = new int[elemsNum[rankOfCurrProc]];

  MPI_Scatterv(currentGen, elemsNum,
               elementsOffsetArray, MPI_INT,
               partMatrGen, elemsNum[rankOfCurrProc],
               MPI_INT, 0, MPI_COMM_WORLD);

  

  MPI_Finalize();
  return 0;
}
