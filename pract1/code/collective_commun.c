#include <float.h> //for DBL_MAX
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

const int root = 0;

void fillVector(unsigned long long *vector1, unsigned long long *vector2,
                size_t sizeVect) {
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

unsigned long long countScalarMult(unsigned long long *vector1, unsigned long long *vector2,
                         size_t sizeVect1, size_t sizeVect2) {
  unsigned long long sum = 0;
  for (size_t i = 0; i < sizeVect1; i++) {
    for (size_t j = 0; j < sizeVect2; j++) {
      sum += vector1[i] * vector2[j];
    }
  }
  return sum;
}

int main(int argc, char *argv[]) {
  // size_t sizeVect = 32 * 32; // for checking the correct computations
  size_t sizeVect = 512 * 24 * 10;
  printf("Size of vector is: %ld\n", sizeVect);

  MPI_Init(&argc, &argv);

  int amountOfAvailableProc;
  MPI_Comm_size(MPI_COMM_WORLD, &amountOfAvailableProc);
  printf("Amount of processes: %d\n", amountOfAvailableProc);

  int rankOfCurrentProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrentProc);

  unsigned long long *baseVector1 = NULL;
  unsigned long long *baseVector2 = NULL;

  const size_t lengthOfPortion = sizeVect / amountOfAvailableProc;
  unsigned long long *bufferedVector1 =
      (unsigned long long *)calloc(lengthOfPortion, sizeof(unsigned long long ));
  unsigned long long *bufferedVector2 =
      (unsigned long long *)calloc(sizeVect, sizeof(unsigned long long ));

  double startTime = 0;

  if (rankOfCurrentProc == 0) {
    baseVector1 = (unsigned long long *)calloc(sizeVect, sizeof(unsigned long long ));
    baseVector2 = (unsigned long long *)calloc(sizeVect, sizeof(unsigned long long ));
    fillVector(baseVector1, baseVector2, sizeVect);

    // printVector(baseVector1, sizeVect);
    // printVector(baseVector2, sizeVect);

    memcpy(bufferedVector2, baseVector2, sizeVect * sizeof(unsigned long long ));
    startTime = MPI_Wtime();
  }

  // 1st vector uses Scatter
  // sendbuf, sendcount, sendtype, recievebuf, recievetype, root, comm
  MPI_Scatter(baseVector1, lengthOfPortion, MPI_UNSIGNED_LONG_LONG,
              bufferedVector1, lengthOfPortion, MPI_UNSIGNED_LONG_LONG, root,
              MPI_COMM_WORLD);

  // 2nd vector uses Broadcast
  // inout buffer, countOfSendingElems, datatype, root, comm
  MPI_Bcast(bufferedVector2, sizeVect, MPI_UNSIGNED_LONG_LONG, root,
            MPI_COMM_WORLD);

  unsigned long long tmpRes = countScalarMult(bufferedVector1, bufferedVector2,
                                    lengthOfPortion, sizeVect);

  unsigned long long finalRes = 0;

  // Broadcast + Reduce for the 2nd vector

  // sendbuff, recieverbuff, countOfElems in sendBuff, datatype,
  // operation, root, comm
  MPI_Reduce(&tmpRes, &finalRes, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, root,
             MPI_COMM_WORLD);

  double minTime = DBL_MAX;
  if (rankOfCurrentProc == 0) {
    double endTime = MPI_Wtime();
    double timeForThisRepeat = endTime - startTime;
    if (timeForThisRepeat < minTime) {
      minTime = timeForThisRepeat;
    }
    printf("Result is: %llu\n", finalRes);
    printf("Time taken: %f sec\n", minTime);
  }

  MPI_Finalize();

  free(baseVector1);
  free(baseVector2);
  free(bufferedVector1);
  free(bufferedVector2);
}
