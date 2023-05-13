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
  int *elemsNumArray = new int[amountOfProcs];

  int basicRowsCount = rows / amountOfProcs;
  int restRowsCount = rows % amountOfProcs;

  for (int i = 0; i < amountOfProcs; i++) {
    elemsNumArray[i] = basicRowsCount * columns;
    if (restRowsCount > 0) {
      elemsNumArray[i] += columns;
      --restRowsCount;
    }
  }
  return elemsNumArray;
}

int *countRowsInEachProcess(const int *elementsNumberArr,
                            int amountOfProcs, int columns) {
  int *rowsNumArray = new int[amountOfProcs];
  for (int i = 0; i < amountOfProcs; i++) {
    rowsNumArray[i] = elementsNumberArr[i] / columns;
  }
  return rowsNumArray;
}

int *createElemsOffsetArr(const int *elemsNumArray,
                          int amountOfProcs) {
  int *elementsOffsetArray = new int[amountOfProcs];
  int elementsOffset = 0;
  for (int i = 0; i < amountOfProcs; i++) {
    elementsOffsetArray[i] = elementsOffset;
    elementsOffset += elemsNumArray[i];
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

// bool *countStopVector(int lengthVector) {
//   bool *vector = new bool[lengthVector]();
// }

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
  int *elemsNumArray = countElemsNumInEachProc(
      amountOfProcs, rowsAmount, columnsAmount);

  // amount of rows, which each process handles
  int *rowsNumArray = countRowsInEachProcess(
      elemsNumArray, amountOfProcs, columnsAmount);

  // count, from which element data sends to each process
  int *elementsOffsetArray =
      createElemsOffsetArr(elemsNumArray, amountOfProcs);

  // count, from which row data sends to each process
  int *rowsOffsetArray = createRowsOffsetArr(
      elementsOffsetArray, amountOfProcs, rowsAmount);

  if (rankOfCurrProc == 0) {
    for (int i = 0; i < amountOfProcs; i++) {
      std::cout << i << ": " << rowsNumArray[i] << std::endl;
    }
  }

  int *partCurrGen = new int[elemsNumArray[rankOfCurrProc]]();

  MPI_Scatterv(currentGen, elemsNumArray, elementsOffsetArray,
               MPI_INT, partCurrGen, elemsNumArray[rankOfCurrProc],
               MPI_INT, 0, MPI_COMM_WORLD);

  int *upperLine = new int[columnsAmount]();
  int *lowerLine = new int[columnsAmount]();

  const long maxIterations = 1000;

  int **historyOfEvolution =
      new int *[maxIterations](); // TODO: don't forget to init mem

  int iterCurr = 0;
  bool repeated = false;

  // determine lower and upper nrighbour

  int rankPrev = rankOfCurrProc - 1;
  int rankNext = rankOfCurrProc + 1;

  int *nextGen = new int[rowsAmount * columnsAmount]();

  int tagFirstLine = 0;
  int tagLastLine = 1;

  MPI_Request requestFirstLineSend;
  MPI_Request requestLastLineSend;
  MPI_Request requestGetLastLine;
  MPI_Request requestGetFirstLine;

  while (iterCurr < maxIterations && !repeated) {

    // initiation of sending first line to the prev core
    MPI_Isend(partCurrGen, elemsNumArray[rankOfCurrProc], MPI_INT,
              rankPrev, tagFirstLine, MPI_COMM_WORLD,
              &requestFirstLineSend);

    // initiation of sending last line the to next core
    MPI_Isend(
        // partCurrGen[elemsNumArray[rankOfCurrProc] - columnsAmount],
        &partCurrGen[rowsOffsetArray[rankOfCurrProc] -
                     1], // mb not -1, but -2, idk yet ...
        elemsNumArray[rankOfCurrProc], MPI_INT, rankNext, tagLastLine,
        MPI_COMM_WORLD, &requestLastLineSend);

    // initiation of receiving last line from the previous core
    MPI_Irecv(upperLine, elemsNumArray[rankOfCurrProc], MPI_INT,
              rankPrev, tagLastLine, MPI_COMM_WORLD,
              &requestGetLastLine);

    // initiation of receiving first line from the next core
    MPI_Irecv(lowerLine, elemsNumArray[rankOfCurrProc], MPI_INT,
              rankNext, tagFirstLine, MPI_COMM_WORLD,
              &requestGetFirstLine);

    bool *vectorStopFlag = new bool[iterCurr - 1];
    for (int core = 0; core < amountOfProcs; core++) {
      // equalsToPrevEvolution()
    }
  }

  MPI_Finalize();
  return 0;
}
