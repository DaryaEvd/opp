#include <iostream>
#include <math.h>
#include <mpi.h>
#include <pthread.h>

#define AMOUNT_OF_LISTS 5
#define WEIGHT_COEFFICIENT 100

#define MIN_AMOUNT_OF_TASKS_TO_SHARE 20
#define TASKS_PER_PROCESS 400

#define TAG_REQUEST 0
#define TAG_REPLY 1

double RES_PER_ITERATION = 0.0;
double GLOABAL_RESULT_SIN = 0.0;

int rankOfCurrProc, amountOfProcs;

int *tasks;
int tasksInRemain;
int amountOfTasksAlreadyExecuted;

pthread_mutex_t mutexTasks;
pthread_mutex_t mutexTasksInRemain;

pthread_t recvThread;

void initTasksWeight() {
  pthread_mutex_lock(&mutexTasks);
  for (int i = 0; i < TASKS_PER_PROCESS; ++i) {
    tasks[i] = abs(50 - i % 100) *
               abs(rankOfCurrProc - (TASKS_PER_PROCESS % amountOfProcs)) *
               WEIGHT_COEFFICIENT;
  }
  pthread_mutex_unlock(&mutexTasks);
}

void calculateTask() {
  pthread_mutex_lock(&mutexTasksInRemain);

  for (int i = 0; tasksInRemain; ++i, tasksInRemain--) {
    pthread_mutex_unlock(&mutexTasksInRemain);

    // когда отдается часть заданий, функция не должна выоплнять часть
    // заздания которая возможна будет прередана
    pthread_mutex_lock(&mutexTasks);
    int task_weight = tasks[i];
    pthread_mutex_unlock(&mutexTasks);

    for (int j = 0; j < task_weight; ++j) {
      RES_PER_ITERATION += sin(j);
    }

    ++amountOfTasksAlreadyExecuted;

    pthread_mutex_lock(&mutexTasksInRemain);
  }
  pthread_mutex_unlock(&mutexTasksInRemain);
}

void *receiverThreadGO(void *args) {
  int tasksToSend;
  int rankRequestedTasks;

  while (true) {
    // receiving process rank that requests tasks
    MPI_Recv(&rankRequestedTasks, 1, MPI_INT, MPI_ANY_SOURCE,
             TAG_REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (rankRequestedTasks == rankOfCurrProc)
      break;

    pthread_mutex_lock(&mutexTasksInRemain);
    if (tasksInRemain >= MIN_AMOUNT_OF_TASKS_TO_SHARE) {
      tasksToSend = tasksInRemain / 2;
      tasksInRemain -= tasksToSend;

      // sending number of tasks that shares
      MPI_Send(&tasksToSend, 1, MPI_INT, rankRequestedTasks,
               TAG_REPLY, MPI_COMM_WORLD);

      pthread_mutex_lock(&mutexTasks);
      // sending tasks
      MPI_Send(tasks + amountOfTasksAlreadyExecuted + tasksInRemain -
                   1,
               tasksToSend, MPI_INT, rankRequestedTasks, TAG_REPLY,
               MPI_COMM_WORLD);
      pthread_mutex_unlock(&mutexTasks);
    } else {
      tasksToSend = 0;

      MPI_Send(&tasksToSend, 1, MPI_INT, rankRequestedTasks,
               TAG_REPLY, MPI_COMM_WORLD);
    }
    pthread_mutex_unlock(&mutexTasksInRemain);
  }
  return NULL;
}

void *workerThreadStart(void *args) {
  tasks = new int[TASKS_PER_PROCESS];

  double startt;
  double minTime, maxTime;

  for (int iterCounter = 0; iterCounter < AMOUNT_OF_LISTS;
       ++iterCounter) {
    initTasksWeight();

    pthread_mutex_lock(&mutexTasksInRemain);
    tasksInRemain = TASKS_PER_PROCESS;
    pthread_mutex_unlock(&mutexTasksInRemain);
    amountOfTasksAlreadyExecuted = 0;
    int amountOfAdditionalasks;

    startt = MPI_Wtime();

    calculateTask();

    //   requesting tasks from other processes
    for (int currentProc = 0; currentProc < amountOfProcs;
         ++currentProc) {
      if (currentProc == rankOfCurrProc)
        continue;

      // Send to other process that cur proc is free
      // and can compute additional tasks
      MPI_Send(&rankOfCurrProc, 1, MPI_INT, currentProc, TAG_REQUEST,
               MPI_COMM_WORLD);

      // Receive number of additional tasks
      MPI_Recv(&amountOfAdditionalasks, 1, MPI_INT, currentProc,
               TAG_REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if (amountOfAdditionalasks > 0) {
        // Receive tasks
        MPI_Recv(tasks, amountOfAdditionalasks, MPI_INT, currentProc,
                 TAG_REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        pthread_mutex_lock(&mutexTasksInRemain);
        tasksInRemain = amountOfAdditionalasks;
        pthread_mutex_unlock(&mutexTasksInRemain);

        calculateTask();
      }
    }

    double endt = MPI_Wtime();
    double resTime = endt - startt;

    MPI_Allreduce(&resTime, &minTime, 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);

    MPI_Allreduce(&resTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    if (rankOfCurrProc == 0) {

      std::cout << "================================================="
                << std::endl;

      std::cout << "Iteration numer: " << iterCounter << std::endl;
      std::cout << "Disbalance time: " << maxTime - minTime
                << std::endl;
      std::cout << "Disbalance percentage: "
                << (maxTime - minTime) / maxTime * 100 << std::endl;
      std::cout << "----------------------------------------------"
                << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int currentProc = 0; currentProc < amountOfProcs;
         currentProc++) {
      if (rankOfCurrProc == currentProc) {
        std::cout << "\t\tCurrent proc is: " << rankOfCurrProc
                  << std::endl;
        std::cout << "Amount of executed tasks: "
                  << amountOfTasksAlreadyExecuted << std::endl;
        std::cout << "Result of calculating is: " << RES_PER_ITERATION
                  << std::endl;
        std::cout << "Time per iteration: " << resTime << " seconds"
                  << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // Terminate Thread 'Receiver'
  // recv назходится в режиме ожидания сообшения, все процесс к
  // которому ты происодеинен закончил работу (говорим рисеиверу)
  MPI_Send(&rankOfCurrProc, 1, MPI_INT, rankOfCurrProc, 0,
           MPI_COMM_WORLD);

  MPI_Allreduce(&RES_PER_ITERATION, &GLOABAL_RESULT_SIN, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);

  delete tasks;

  return NULL;
}

void createAndStartThreads() {
  pthread_mutex_init(&mutexTasks, NULL);
  pthread_mutex_init(&mutexTasksInRemain, NULL);

  pthread_attr_t attributes;
  if (pthread_attr_init(&attributes) != 0) {
    MPI_Finalize();
    perror("Can't init attributes");
    abort();
  }

  if (pthread_attr_setdetachstate(&attributes,
                                  PTHREAD_CREATE_JOINABLE) != 0) {
    MPI_Finalize();
    perror("Error in setting attributes");
    abort();
  }

  if (pthread_create(&recvThread, &attributes, receiverThreadGO,
                     NULL) != 0) {
    MPI_Finalize();
    perror("Can't create thread");
    abort();
  }

  pthread_attr_destroy(&attributes);

  workerThreadStart(NULL); // it's main thread

  // main thread is waitng for recieiver thread to finish
  pthread_join(recvThread, NULL);

  pthread_mutex_destroy(&mutexTasks);
  pthread_mutex_destroy(&mutexTasksInRemain);
}

int main(int argc, char **argv) {
  int reqiredLevel =
      MPI_THREAD_MULTIPLE; // we want this level of supporting threads
  int providedLevel;       // real level of supporting threads

  MPI_Init_thread(&argc, &argv, reqiredLevel, &providedLevel);
  if (providedLevel != reqiredLevel) {
    MPI_Finalize();
    perror("Can't load reqired level");
    return 0;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  double startt = MPI_Wtime();
  createAndStartThreads();
  double endt = MPI_Wtime();

  double resTime = endt - startt;

  if (rankOfCurrProc == 0) {
    std::cout << "======================================="
              << std::endl;
    std::cout << "Time for all lists: " << resTime << "seconds"
              << std::endl;
  }

  MPI_Finalize();
  return 0;
}
