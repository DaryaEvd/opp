#include <iostream>
#include <mpi.h>

/* Axises are
|-----Y
|
|
X
*/

static int dimOfAnyGrid = 2;
static int X_AXIS = 0;
static int Y_AXIS = 1;

struct MyMatrix {
  size_t row;
  size_t column;
  double *data = nullptr;

  MyMatrix(size_t row, size_t column) {
    this->row = row;
    this->column = column;
    data = new double[row * column]();
  }
};

void freeMatrix(MyMatrix matrix) { delete[] matrix.data; }

void initMatrix(MyMatrix matrix) {
  for (size_t i = 0; i < matrix.row; ++i) {
    for (size_t j = 0; j < matrix.column; ++j) {
      matrix.data[i * matrix.column + j] = rand() % 100 + 15;
    }
  }
}

void multimplyMtrices(MyMatrix m1, MyMatrix m2, MyMatrix mRes) {
  for (size_t i = 0; i < m1.row; i++) {
    for (size_t j = 0; j < m2.column; j++) {
      for (size_t k = 0; k < m1.column; k++) {
        mRes.data[i * m2.column + j] +=
            m1.data[i * m1.column + k] * m2.data[k * m2.column + j];
      }
    }
  }
}

void printMat(MyMatrix matrix) {
  for (size_t i = 0; i < matrix.row; i++) {
    for (size_t j = 0; j < matrix.column; j++) {
      std::cout << matrix.data[i * matrix.column + j] << " ";
    }
    std::cout << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc != 6) {
    std::cout
        << "Error! Enter rowsInC and columns amount for 2 matrixes!"
        << std::endl;
    return 0;
  }

  // srand(time(nullptr)); // just for random vales each time
  srand(0); // for testint the same values in one task

  const size_t dim1 = atoi(argv[1]);
  const size_t dim2 = atoi(argv[2]);
  const size_t dim3 = atoi(argv[3]);

  const int rowsGrid = atoi(argv[4]);    // height of 2d grid
  const int columnsGrid = atoi(argv[5]); // weight of 2d grid

  MPI_Init(&argc, &argv);
  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  if (amountOfProcs != rowsGrid * columnsGrid) {
    if (rankOfCurrProc == 0) {
      std::cout << "Bad input! Amount of processes should be equal "
                   "to rowsGrid * columnsGrid"
                << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  int dimSize[] = {rowsGrid, columnsGrid};
  int periods[] = {0};

  MPI_Comm commGrid;

  double startt = MPI_Wtime();

  /**** MPI_Cart_create
    * @param  MPI_Comm_oldcomm - исходный коммуникатор, для процессов
  которого определяется декартова топология
    * @param int ndims - число размерности создаваемой декартоковй
  решетки
    * @param  int *dims - целочисленный массив, каждый элемент
  которого опредеяет размер по каждому измерению
    * @param int *periods - целочисленный массив флагов, определяющих
  периодичнось каждого измерения.
  Периодичность(циличность). При пересылке данных для решения
  некоторых задач удобно, чтобы послдний процесс пересылал данные в
  первый. Это можно автоматически обеспечить, задав этот параметр
  единицей. Если флаг = 1, то по этому измерению будет иметь место
  периодичность. В основном исопльзуется для cart_shift
    * @param int reorder - целочисленный флаг, определяющий, можно ли
  среде MPI автоматически менять порядок нумерации процессов
      Этот пареметр дает возможность самой системе MPI поменять
  порядок процессов в исходном коммуникаторе в процессе определения
  нового коммуникатора. То есть процессы будут тем же самым
  множеством, они будут организованы в виде решетки, но их ранги в
  новом коммуникаторе не будут совпадать с рангами исходного
  коммуникатора.
  Зачем это было предусмотрено? MPI может знать конфигурации
  компа, на котором будет запущена программа. И бывает удобно эту
  конфигурацию учесть.
  Если в самой декартовой топологии расположение процессов
  соответствует расположению физических устройств, на которых эти
  процессы запущены, то можно оптимизировать процесс передечи данных
  между процессами.
  Для нашей задачи нужно, чтобы ранги процессов в точности соппадали
  с рангами процессов в исходном коммуникаторе. Поэтому reorder = 0.
    * @param  MPI_Comm *cartcomm - результирующий коммуникатор с
  декартовой топологией
  */
  MPI_Cart_create(MPI_COMM_WORLD, dimOfAnyGrid, dimSize, periods, 0,
                  &commGrid);

  int coordsOfCurrProc[dimOfAnyGrid];

  /*MPI_Cart_coords
    Для каждого процесса в нашей решетке определим его координаты.
    Из любого процесса можем получить координаты любого другого
    процесса
  */
  MPI_Cart_coords(commGrid, rankOfCurrProc, dimOfAnyGrid,
                  coordsOfCurrProc);

  /*MPI_Cart_sub
    Хотим получить новые коммуникаторы, в каждом из которых будут
    располагаться процессы, входящие в одну и ту же строку.
    Их порядок будет соответствовать порядку в исходной двумерной
    решетке.
    Чтобы создать расщепелние, надо сказать, по каким
    координатам проводить объединение. Для создания коммуникатора
    строки хотим объединять по первой координате, которая потом
    пропадет в результате объединения. А вторая останется. И эта
    оставшаяся координата будет единственной в оставшихся одномерных
    расщепленных коммуникаторов. То есть объединение производится по
    координатам, которые мы убираем.
  */
  MPI_Comm commRow;
  int remainX[] = {X_AXIS, Y_AXIS}; // 1я координата пропадет, а
                                    // 2я координата останется
  MPI_Cart_sub(commGrid, remainX, &commRow);

  MPI_Comm commColumn;
  int remainY[] = {Y_AXIS, X_AXIS}; // 1я координата останется, а
                                    // 2я координата пропадет
  MPI_Cart_sub(commGrid, remainY, &commColumn);

  MyMatrix A(dim1, dim2);
  MyMatrix B(dim2, dim3);
  MyMatrix C(dim1, dim3);

  if (rankOfCurrProc == 0) {
    initMatrix(A);
    initMatrix(B);

    // std::cout << "A[" << A.row << " x " << A.column << "]"
    //           << std::endl;
    // printMat(A);
    // std::cout << "B[" << B.row << " x " << B.column << "]"
    //           << std::endl;
    // printMat(B);
  }

  MyMatrix partA(dim1 / rowsGrid, dim2);

  /*
    Раздача матрицы А по горизонтальным полосам на вертикальную
    линейку процессов (0;0), (1;0), (2;0), …, (p1 - 1; 0) при помощи
    MPI_Scatter
  */
  if (coordsOfCurrProc[Y_AXIS] == 0) {
    MPI_Scatter(A.data, partA.row * partA.column, MPI_DOUBLE,
                partA.data, partA.row * partA.column, MPI_DOUBLE, 0,
                commColumn);
  }

  /*
    Каждый из процессов в левой вертикальной колонке
     ( (1;0), (2;0), …, (p1 - 1; 0) ) при помощи MPI_Bcast раздает
    свою полосу матрицы A всем процессам своей горизонтали.
    Т.е. процесс (1;0) раздает свою полосу процессам (1;1), (1;2), ...
  */
  MPI_Bcast(partA.data, partA.row * partA.column, MPI_DOUBLE, 0,
            commRow);

  MyMatrix partB(dim2, dim3 / columnsGrid);

  /*
    Определение нового производного типа данных для выбора из
    матрицы B вертикальных полос
  */
  MPI_Datatype bSendType; // type of columns in matrix B

  /*** MPI_Type_vector
   * @param count -  число блоков
   * @param blocklength - число элементов в каждом блоке
   * @param stride - число элементов между началами каждого блока
   * @param oldtype - стрый тип данных
   * @param newtype - новый тип данных
   */
  MPI_Type_vector(partB.row, partB.column, dim3, MPI_DOUBLE,
                  &bSendType);

  /*MPI_type_commit
    Описанные типы являются локальными, т.е. могу вызываться только в
    некотороых процессах параллельного приложения - в тех, которые в
    дальшейшем принимают новые типы (для определения других новых
    типов или посылки/приёма данных)

    Если новый тип хотим
    использовать при пересылке сообщений между процессами, то его
    необходимо дополнительно "зрарегистрировать", вызвав для него
    функцию MPI_type_commit
  */
  MPI_Type_commit(&bSendType);

  /*
    Создаем новый тип, который определяет только один (!!!) столбец
    матрицы B
  */
  MPI_Datatype bRecvType; // type of only ONE (!!!) column
  MPI_Type_vector(1, partB.row * partB.column, 0, MPI_DOUBLE,
                  &bRecvType);
  MPI_Type_commit(&bRecvType);

  /*
    Выделяем данные в матрице В для дальнейшей раздачи

    ** Раздачу производит нулевой процесс
    Вообще, можно было бы написать 'if rankOfCurrProc == 0', но
    лучше утончить координаты
  */
  if (coordsOfCurrProc[X_AXIS] == 0 &&
      coordsOfCurrProc[Y_AXIS] == 0) {
    // цикл от 1, потому что нулевой процесс не может сам себе
    // высылать данные
    for (int i = 1; i < columnsGrid; i++) {
      MPI_Send(B.data + partB.column * i, 1, bSendType, i, 11,
               commRow);
    }
    /*
      Нулевой процесс отправляет свои данные.
      Делаем отдельно, тк нулевой процесс не может сам себе сендить
    */
    for (int i = 0; i < B.row; i++) {
      for (int j = 0; j < partB.column; j++) {
        partB.data[i * partB.column + j] = B.data[i * B.column + j];
      }
    }
  }

  /*
    Раздача матрицы B по вертикальным полосам на горизонтальную
    линейку процессов (0;0), (0;1), (0;2), …, (0; p2 – 1). Здесь
    происходит прием данных на нулевой столбец решётки
  */
  else if (coordsOfCurrProc[X_AXIS] == 0) {
    MPI_Recv(partB.data, 1, bRecvType, 0, 11, commRow,
             MPI_STATUS_IGNORE);
  }

  /*
    Раздаем матрицу B по вертикальным столбцам решётки
  */
  MPI_Bcast(partB.data, dim2 * partB.column, MPI_DOUBLE, 0,
            commColumn);

  MyMatrix partC(dim1 / rowsGrid, dim3 / columnsGrid);
  /*
    Перемножаем строку матрицы А на столбец матрицы В
  */
  multimplyMtrices(partA, partB, partC);

  MPI_Datatype cRecvType;
  MPI_Type_vector(partC.row, partC.column, dim3, MPI_DOUBLE,
                  &cRecvType);
  MPI_Type_commit(&cRecvType);

  int offset[amountOfProcs];
  for (int procRank = 0; procRank < amountOfProcs; procRank++) {
    MPI_Cart_coords(commGrid, procRank, dimOfAnyGrid,
                    coordsOfCurrProc);

    /*
      Отсчитываем смещение для каждого процесса.
      То есть вот мы отдельно посчитали каждый "квадратик" (partC),
      а теперь нам надо определить расположение этого квадратика
      относительно начала буфера матрицы С - левого верхнего угла,
      условно говоря

      *Можно воспринимать это размещение в духе:
      вправо столько-то клеточек, вниз - столько-то
      (ну типа того, поскольку тут посложнее малёк)
    */

    offset[procRank] =
        coordsOfCurrProc[X_AXIS] * partC.row * C.column +
        coordsOfCurrProc[Y_AXIS] * partC.column;
  }

  /*
    Собираем всю матрицу С на на процессе (0; 0)
  */
  if (coordsOfCurrProc[X_AXIS] == 0 &&
      coordsOfCurrProc[Y_AXIS] == 0) {
    /*
      Нулевой процесс отправляет (ну тип копирует) данные из части
      матрицы в целую матрицы, потому что отправлять самому себе
      нельзя))
    */
    for (int i = 0; i < partC.row; i++) {
      for (int j = 0; j < partC.column; j++) {
        C.data[i * C.column + j] = partC.data[i * partC.column + j];
      }
    }
    for (int i = 1; i < amountOfProcs; i++) {
      MPI_Recv(C.data + offset[i], 1, cRecvType, i, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  } else {
    // Отправляем полученные данные со всех процессов в нулевой
    // процесс
    MPI_Send(partC.data, partC.row * partC.column, MPI_DOUBLE, 0, 0,
             MPI_COMM_WORLD);
  }

  MPI_Type_free(&bSendType);
  MPI_Type_free(&cRecvType);

  MPI_Comm_free(&commGrid);
  MPI_Comm_free(&commColumn);
  MPI_Comm_free(&commRow);

  double endt = MPI_Wtime();

  if (rankOfCurrProc == 0) {
    // std::cout << "C[" << C.row << " x " << C.column << "]"
    //           << std::endl;
    // printMat(C);

    std::cout << "Time taken: " << endt - startt << " sec"
              << std::endl;
  }

  freeMatrix(A);
  freeMatrix(B);
  freeMatrix(C);

  freeMatrix(partA);
  freeMatrix(partB);
  freeMatrix(partC);

  MPI_Finalize();
  return 0;
}
