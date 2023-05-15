#include "mpi.h"
#include <iostream>
#include <iterator>
#include <vector>

void generateGlider(bool *extendedPartMatr, int columnsAmount) {
  extendedPartMatr[2 * columnsAmount + 0] = true;
  extendedPartMatr[2 * columnsAmount + 1] = true;
  extendedPartMatr[2 * columnsAmount + 2] = true;
  extendedPartMatr[1 * columnsAmount + 2] = true;
  extendedPartMatr[0 * columnsAmount + 1] = true;
}

bool CompareMatrices(const bool *first, const bool *second,
                     size_t from, size_t to) {
  for (size_t i = from; i < to; ++i) {
    if (first[i] != second[i]) {
      return false;
    }
  }
  return true;
}

void CalculateStopVector(std::vector<bool *> historyOfEvolution,
                         bool *stop_vector, bool *extendedPartMatr,
                         int rowsAmount, int columnsAmount) {
  size_t vector_size = historyOfEvolution.size() - 1;
  auto it = historyOfEvolution.begin();
  for (int i = 0; i < vector_size; ++i) {
    stop_vector[i] =
        CompareMatrices(*it, extendedPartMatr, columnsAmount,
                        columnsAmount * (rowsAmount + 1));
    it++;
  }
}

bool CheckIsEnd(int rowsAmount, int columnsAmount,
                const bool *stop_matrix) {
  for (int i = 0; i < columnsAmount; ++i) {
    bool stop = true;
    for (int j = 0; j < rowsAmount; ++j) {
      stop &= stop_matrix[j * columnsAmount + i];
    }
    if (stop)
      return true;
  }
  return false;
}

// bool CheckIsEnd(bool *stop, int num, int size) {
//   for (int i = 0; i < num; i++) {
//     int res = 0;
//     for (int j = 0; j < size; j++)
//       res += stop[i + j * num];
//     if (res == size)
//       return true;
//   }

//   return false;
// }

bool MakeDecision(bool prev, int cnt) {
  if (prev) {
    if (cnt < 2 || cnt > 3)
      return false;
  } else {
    if (cnt == 3)
      return true;
  }
  return prev;
}

void CalcNext(int rowsAmount, int columnsAmount, const bool *prev,
              bool *next) {
  for (int i = 1; i < rowsAmount - 1; ++i) {
    for (int j = 0; j < columnsAmount; ++j) {
      int cnt = 0;
      if (prev[i * columnsAmount + (j + 1) % columnsAmount])
        cnt++;
      if (prev[i * columnsAmount +
               (j + columnsAmount - 1) % columnsAmount])
        cnt++;
      if (prev[(i + 1) * columnsAmount + (j + 1) % columnsAmount])
        cnt++;
      if (prev[(i + 1) * columnsAmount +
               (j + columnsAmount - 1) % columnsAmount])
        cnt++;
      if (prev[(i - 1) * columnsAmount + (j + 1) % columnsAmount])
        cnt++;
      if (prev[(i - 1) * columnsAmount +
               (j + columnsAmount - 1) % columnsAmount])
        cnt++;
      if (prev[(i + 1) * columnsAmount + j])
        cnt++;
      if (prev[(i - 1) * columnsAmount + j])
        cnt++;
      next[i * columnsAmount + j] =
          MakeDecision(prev[i * columnsAmount + j], cnt);
    }
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
  int *rowsNumArr = new int[amountOfProcs];
  for (int i = 0; i < amountOfProcs; i++) {
    rowsNumArr[i] = elementsNumberArr[i] / columns;
  }
  return rowsNumArr;
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

void startLife(int amountOfProcs, int rankOfCurrProc,
               bool *startMatrix, int rowsAmount, int columnsAmount) {
  // amount of elems, handling each process
  int *counts = countElemsNumInEachProc(amountOfProcs, rowsAmount,
                                        columnsAmount);

  // amount of rows, which each process handles
  int *rowsNumArr =
      countRowsInEachProcess(counts, amountOfProcs, columnsAmount);

  // count, from which element data sends to each process
  int *offsets = createElemsOffsetArr(counts, amountOfProcs);

  // count, from which row data sends to each process
  int *rowsOffsetArray =
      createRowsOffsetArr(offsets, amountOfProcs, rowsAmount);

  // if (rankOfCurrProc == 0) {
  //   for (int i = 0; i < amountOfProcs; i++) {
  //     std::cout << i << ": " << rowsNumArr[i] << std::endl;
  //   }
  // }

  bool *extendedPartMatr =
      new bool[counts[rankOfCurrProc] + columnsAmount * 2];
  bool *basePartMatr = extendedPartMatr + columnsAmount;

  MPI_Scatterv(startMatrix, counts, offsets, MPI_C_BOOL, basePartMatr,
               counts[rankOfCurrProc], MPI_C_BOOL, 0, MPI_COMM_WORLD);

  int prevRank = (rankOfCurrProc + amountOfProcs - 1) % amountOfProcs;
  int nextRank = (rankOfCurrProc + 1) % amountOfProcs;

  std::vector<bool *> historyOfEvolution;

  int iteration = 0;
  bool stop = false;
  while (!stop) {

    // Save iteration extendedPartMatr
    bool *extendedNextPartMatr =
        new bool[counts[rankOfCurrProc] + 2 * columnsAmount];
    bool *baseNextPartMatr = extendedNextPartMatr + columnsAmount;
    historyOfEvolution.push_back(extendedPartMatr);

    // Start iteration
    iteration++;
    MPI_Request requestSendFirstLine, requestSendLastLine;
    MPI_Request requestGetLastLine, requestGetFirstLine;

    MPI_Isend(basePartMatr, columnsAmount, MPI_C_BOOL, prevRank, 1,
              MPI_COMM_WORLD, &requestSendFirstLine);
    MPI_Isend(basePartMatr + counts[rankOfCurrProc] - columnsAmount,
              columnsAmount, MPI_C_BOOL, nextRank, 0, MPI_COMM_WORLD,
              &requestSendLastLine);
    MPI_Irecv(extendedPartMatr, columnsAmount, MPI_C_BOOL, prevRank,
              0, MPI_COMM_WORLD, &requestGetLastLine);
    MPI_Irecv(basePartMatr + counts[rankOfCurrProc], columnsAmount,
              MPI_C_BOOL, nextRank, 1, MPI_COMM_WORLD,
              &requestGetFirstLine);

    // Check
    MPI_Request flagsReq;
    bool *stop_vector;
    bool *stop_matrix;
    int vector_size = historyOfEvolution.size() - 1;
    if (vector_size > 1) {
      stop_vector = new bool[vector_size];
      CalculateStopVector(historyOfEvolution, stop_vector,
                          extendedPartMatr,
                          rowsNumArr[rankOfCurrProc], columnsAmount);
      stop_matrix = new bool[vector_size * amountOfProcs];
      MPI_Iallgather(stop_vector, (int)vector_size, MPI_C_BOOL,
                     stop_matrix, (int)vector_size, MPI_C_BOOL,
                     MPI_COMM_WORLD, &flagsReq);

      // MPI_Ialltoall(stop_vector, vector_size, MPI_C_BOOL,
      // stop_matrix,
      //               vector_size, MPI_C_BOOL, MPI_COMM_WORLD,
      // &flagsReq);

      // MPI_Wait(&flagsReq, MPI_STATUS_IGNORE);

      // stop = CheckIsEnd(amountOfProcs, (int)vector_size,
      // stop_matrix);
      // // stop = CheckIsEnd(stop_vector, iteration, amountOfProcs);

      // delete[] stop_vector;
      // delete[] stop_matrix;
    }

    // if (stop)
    //   break;

    CalcNext(rowsNumArr[rankOfCurrProc], columnsAmount, basePartMatr,
             baseNextPartMatr);
    MPI_Wait(&requestSendFirstLine, MPI_STATUS_IGNORE);
    MPI_Wait(&requestGetLastLine, MPI_STATUS_IGNORE);
    CalcNext(3, columnsAmount, extendedPartMatr,
             extendedNextPartMatr);
    MPI_Wait(&requestSendLastLine, MPI_STATUS_IGNORE);
    MPI_Wait(&requestGetFirstLine, MPI_STATUS_IGNORE);
    CalcNext(3, columnsAmount,
             basePartMatr +
                 (rowsNumArr[rankOfCurrProc] - 2) * columnsAmount,
             baseNextPartMatr +
                 (rowsNumArr[rankOfCurrProc] - 2) * columnsAmount);

    if (vector_size > 1) {
      MPI_Wait(&flagsReq, MPI_STATUS_IGNORE);

      stop = CheckIsEnd(amountOfProcs, (int)vector_size, stop_matrix);
      // stop = CheckIsEnd(stop_vector, iteration, amountOfProcs);

      delete[] stop_vector;
      delete[] stop_matrix;
    }

    if (stop)
      break;

    // Switch main extendedPartMatr -- next iteration
    extendedPartMatr = extendedNextPartMatr;
    basePartMatr = baseNextPartMatr;
  }
  if (rankOfCurrProc == 0) {
    std::cout << "Repeat after: " << iteration - 1 << " iterations"
              << std::endl;
  }

  for (bool *matrixDump : historyOfEvolution) {
    delete[] (matrixDump);
  }

  delete[] offsets;
  delete[] counts;
}

int main(int argc, char *argv[]) {

  if (argc != 3) {
    std::cout
        << "Bad amount of arguments!\n"
        << "Enter rowsAmount amount, then enter columns amount.\n"
        << std::endl;
    return 0;
  }

  const int rowsAmount = std::atoi(argv[1]);
  const int columnsAmount = std::atoi(argv[2]);

  MPI_Init(&argc, &argv);

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  bool *startMatrix = nullptr;
  if (rankOfCurrProc == 0) {
    startMatrix = new bool[rowsAmount * columnsAmount];
    std::fill(startMatrix, startMatrix + columnsAmount * rowsAmount,
              false);
    generateGlider(startMatrix, columnsAmount);
  }

  double startt = MPI_Wtime();
  startLife(amountOfProcs, rankOfCurrProc, startMatrix, rowsAmount,
            columnsAmount);
  double endt = MPI_Wtime();

  if (rankOfCurrProc == 0) {
    std::cout << "Time taken: " << endt - startt << " sec."
              << std::endl;
  }

  if (rankOfCurrProc == 0)
    delete[] startMatrix;
  MPI_Finalize();
  return 0;
}