#include <iostream>
#include <mpi.h>
#include <thread.h>

struct DataThread {
  pthread_t threadID;
  int threadNumber;
  int amountOfThreads;
  int proccessID;
}

void doTask() {

}

int main(int argc, char **argv) {
  if (argc != 1) {
    std::cout << "Bad amount of arguments!\n"
              // << "Enter process amount. \n"
              << std::endl;
    return 0;
  }

  MPI_Init(&argc, &argv);

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  int amountOfThreadsInEachProc = 2; // ?????????????

  struct DataThread *tData =
      malloc(sizeof(*tData) * amountOfThreadsInEachProc);
  if (tData == NULL) {
    fprintf(stderr, "No enough mem \n");
    exit(1);
  }

  for(int i = 1; i < amountOfThreadsInEachProc; i++) {
    tData[i].threadNumber = i;
    tData[i].amountOfThreads = amountOfThreadsInEachProc;

    if(pthread_create(&tData[i].threadID, NULL, doTask, &tData[i]) != 0) {
      fprintf(stderr, "Can't create thread \n");
      exit(1);
    }
  }

  tData[0].threadNumber = 0;
  tData[1].amountOfThreads = amountOfThreadsInEachProc;
  doTask(&tData[0]);

  return 0;
}