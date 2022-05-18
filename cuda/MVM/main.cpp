// matrix vector multiplecation
#include <cstdlib>
#include <iostream>
#include <vector>

std::vector<float> cpuMV(std::vector<float>& matrix, std::vector<float>& vector,
                         size_t m, size_t n) {
  std::vector<float> ans;
  // the position of (0,0) is at the top left corner
  float tempv = 0;

  for (size_t i = 0; i < m; i++) {
    tempv = 0;
    for (size_t j = 0; j < n; j++) {
      tempv = tempv+(matrix[i * n + j] * vector[j]);
      printf("i %ld j %ld m %f v %f tempv %f\n", i, j, matrix[i * n + j],
             vector[j], tempv);
    }
    ans.push_back(tempv);
  }
  return ans;
}

int main(int argc, char** argv) {
  // int matrix size m*n
  size_t m = 10;
  size_t n = 20;

  // matrix (it is unnecessary to use the insert vector here)
  // row m column n
  std::vector<float> matrix(m * n);

  // vector is n*1
  std::vector<float> vector(n);

  // init matrix and vector
  srand(static_cast<unsigned>(time(0)));

  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      matrix[i * n + j] =
          static_cast<float>(rand()) / static_cast<float>(RAND_MAX);

      printf("i %ld j %ld v %lf\n", i, j, matrix[i * n + j]);
    }
  }

  // init vector
  for (size_t j = 0; j < n; j++) {
    vector[j] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    printf("vector j %ld, v %lf\n", j, vector[j]);
  }

  // ans is m*1
  std::vector<float> ansCPU(m);

  ansCPU = cpuMV(matrix, vector, m, n);

  for (size_t i = 0; i < m; i++) {
    std::cout << ansCPU[i] << ",";
  }

  return 0;
}