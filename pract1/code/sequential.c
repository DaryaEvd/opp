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

int main() {
  // size_t sizeVect = 32 * 32; //for checking the correct computations
  size_t sizeVect = 512 * 24 * 10;
  printf("Size of vector is: %ld\n", sizeVect);

  unsigned long long *vector1 = (unsigned long long *)calloc(
      sizeVect, sizeof(unsigned long long  ));

  unsigned long long *vector2 = (unsigned long long *)calloc(
      sizeVect, sizeof(unsigned long long  ));

  struct timespec start, end;
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  fillVector(vector1, vector2, sizeVect);

  // printVector(vector1, sizeVect);
  // printVector(vector2, sizeVect);

  unsigned long long result =
      countScalarMult(vector1, vector2, sizeVect, sizeVect);
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);

  printf("Result is: %llu\n", result);
  printf("Time taken: %lf sec\n",
         end.tv_sec - start.tv_sec +
             0.000000001 * (end.tv_nsec - start.tv_nsec));
  free(vector1);
  free(vector2);
}
