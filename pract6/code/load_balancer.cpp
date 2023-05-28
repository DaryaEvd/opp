#include <iostream>
#include <mpi.h>
#include <thread.h>

#define WORK_THREAD 0
#define RECV_THREAD 1

double globalRes = 0.0;

// arg - элемент структуры
void workThread(void *arg) {
  struct DataThread *data = (struct DataThread *)arg;
}

void receiveThread(void *arg) {}

int main(int argc, char **argv) {
  // if (argc != 1) {
  //   std::cout << "Bad amount of arguments!\n"
  //             // << "Enter process amount. \n"
  //             << std::endl;
  //   return 0;
  // }

  MPI_Init_thread(int *argc, char ***argv, int required,
                  int *provided);

  double startt = MPI_Wtime();

  int reqiredLevel =
      MPI_THREAD_MULTIPLE; // we want this level of supporting threads
  int providedLevel;       // real level of supporting threads

  MPI_Init_thread(&argc, &argv, reqiredLevel, &providedLevel);
  if (providedLevel != reqiredLevel) {
    MPI_Finalize();
    perror("Can't load reqired level");
    return 0;
  }

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  pthread_attr_t attrs;
  if (pthread_attr_init(&attrs) != 0) {
    MPI_Finalize();
    perror("Can't init attrs");
    abort();
  }

  if (pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE)) {
    perror("Error in setting attrs");
    abort();
  }

  int amountOfThreadsInEachProc = 2;
  pthread_t myThreads[amountOfThreadsInEachProc];

  if (pthread_create(&myThreads[WORK_THREAD], &attrs,
                     workThread, &)) {
  }

  if (pthread_create(&myThreads[RECV_THREAD], attrs,
                     receiveThread, &) != 0) {
    perror("Can't create thread");
    abort();
  }

  pthread_attr_destroy(&attrs);

  for (int i = 0; i < amountOfThreadsInEachProc; i++) {
    if (pthreads_join(myThreads[i], NULL) != 0) {
      perror("Can not join a thread");
      abort();
    }
  }

  int amountOfTasksInList = 100 * amountOfProcs;

  double endt = MPI_Wtime();

  double time = endt - startt;

  if (rankOfCurrProc == 0) {
    std::cout << "Time: " << time << std::endl;
  }

  MPI_Finalize();
  return 0;
}
