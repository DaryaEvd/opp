#include <iostream>
#include <fstream>

void initMatrix(int *matrix, int rowsAmount, int columnsAmount) {
  for(int i = 0; i < rowsAmount; i++) {
    for(int j = 0; j < columnsAmount; j++) {
      matrix[i * columnsAmount + j] = i;
    }
  }
}

int main(int argc, char **argv) {
  if(argc != 3) {
    printf("Bad amount of arguments!\n" \
    "Enter rows amount, then enter columns amount.\n");
    return 0;
  }
 
  const int rowsAmount = atoi(argv[1]);
  const int columnsAmount = atoi(argv[2]);

  int *matrix = new int[rowsAmount * columnsAmount];

  initMatrix(matrix, rowsAmount, columnsAmount);

  std::ofstream inputFile;
  inputFile.open("start.txt");
  if(!inputFile) {
    printf("Can't open input file\n");
    return 0;
  }

  for(size_t row = 0; row < rowsAmount; row++) {
    for(size_t column = 0; column < columnsAmount; column++) {
      inputFile << matrix[row * columnsAmount + column] << " ";
    }
    inputFile << "\n";
  }

  inputFile.close();

  delete [] matrix;

  return 0;
}
