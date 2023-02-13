#include <float.h> //for DBL_MAX
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void fillVector(unsigned long long *vector1,
                unsigned long long *vector2, size_t sizeVect) {
  // srand(time(NULL));
  // for (size_t i = 0; i < sizeVect; i++) {
  //   vector1[i] = rand() % 100;
  //   vector2[i] = rand() % 100;
  // }

  /* for checking the correct computations*/
  for (size_t i = 0; i < sizeVect; i++) {
    vector1[i] = 1;
    vector2[i] = 2;
  }
}

// void printVector(unsigned long long *vector, size_t sizeVect) {
//   for (size_t i = 0; i < sizeVect; i++) {
//     printf("%lld ", vector[i]);
//   }
//   printf("\n");
// }

unsigned long long countScalarMult(unsigned long long *vector1,
                                   unsigned long long *vector2,
                                   size_t sizeVect1,
                                   size_t sizeVect2) {
  unsigned long long sum = 0;
  for (size_t i = 0; i < sizeVect1; i++) {
    for (size_t j = 0; j < sizeVect2; j++) {
      sum += vector1[i] * vector2[j];
    }
  }
  return sum;
}

int main(int argc, char *argv[]) {
  // size_t sizeVect = 32 * 32; //for checking the correct
  // computations
  size_t sizeVect = 512 * 24 * 10;
  printf("Size of vector is: %ld\n", sizeVect);

  MPI_Init(&argc, &argv);

  int amountOfAvailableProc;
  MPI_Comm_size(MPI_COMM_WORLD, &amountOfAvailableProc);
  printf("Amount of processes: %d\n", amountOfAvailableProc);

  int rankOfCurrentProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrentProc);

  const size_t lengthOfPortion = sizeVect / amountOfAvailableProc;

  double minTime = DBL_MAX;

  if (rankOfCurrentProc == 0) {
    unsigned long long *vector1 = (unsigned long long *)calloc(
        sizeVect, sizeof(unsigned long long));
    unsigned long long *vector2 = (unsigned long long *)calloc(
        sizeVect, sizeof(unsigned long long));

    double startTime = MPI_Wtime();
    fillVector(vector1, vector2, sizeVect);

    // printVector(vector1, sizeVect);
    // printVector(vector2, sizeVect);

    for (size_t i = 1; i < amountOfAvailableProc; i++) {
      int dest = i;
      int msgTag1 = 1;
      int msgTag2 = 2;

      // sendBuff, countOfElemsInBuff, datatype, destNumberOfReciever,
      // tagMsg, comm
      MPI_Send(&vector1[i * lengthOfPortion], lengthOfPortion,
               MPI_UNSIGNED_LONG_LONG, dest, msgTag1, MPI_COMM_WORLD);
      MPI_Send(vector2, sizeVect, MPI_UNSIGNED_LONG_LONG, dest,
               msgTag2, MPI_COMM_WORLD);
    }

    unsigned long long result =
        countScalarMult(vector1, vector2, lengthOfPortion, sizeVect);

    for (size_t i = 1; i < amountOfAvailableProc; i++) {
      unsigned long long bufferOfPartialResults = 0;
      // data will save here
      int source = i;
      int tag = i;

      // recieveBuff, countOfRecvBuffer, datatype, sourceNumofSender,
      // tagMsg, comm, status
      MPI_Recv(&bufferOfPartialResults, 1, MPI_UNSIGNED_LONG_LONG,
               source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      result += bufferOfPartialResults;
    }

    double endTime = MPI_Wtime();
    double timeForThisRepeat = endTime - startTime;
    if (timeForThisRepeat < minTime) {
      minTime = timeForThisRepeat;
    }

    printf("Result is: %llu\n", result);
    printf("Time taken: %f sec\n", minTime);

    free(vector1);
    free(vector2);
  }

  else {
    unsigned long long *vector1 = (unsigned long long *)calloc(
        sizeVect, sizeof(unsigned long long));
    unsigned long long *vector2 = (unsigned long long *)calloc(
        sizeVect, sizeof(unsigned long long));

    fillVector(vector1, vector2, sizeVect);

    // recieveBuff, countOfRecvBuffer, datatype, sourceNumofSender,
    // tagMsg, comm, status
    MPI_Recv(vector1, lengthOfPortion, MPI_UNSIGNED_LONG_LONG, 0, 1,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(vector2, sizeVect, MPI_UNSIGNED_LONG_LONG, 0, 2,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    unsigned long long result =
        countScalarMult(vector1, vector2, lengthOfPortion, sizeVect);

    MPI_Send(&result, 1, MPI_UNSIGNED_LONG_LONG, 0, rankOfCurrentProc,
             MPI_COMM_WORLD);

    free(vector1);
    free(vector2);
  }

  MPI_Finalize();
}
