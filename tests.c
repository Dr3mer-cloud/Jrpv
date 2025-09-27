#include "header.h"
#include "solve.h"
#include "tests.h"

void InitializeIdentityMatrix(double *p, const size_t n) {
  size_t i = 0, j = 0;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i == j) {
        p[i * n + j] = 1.;
      }
      else {
        p[i * n + j] = 0.;
      }
    }
  }
}

double* GetMatrixProduct(const double *a, const double *b, const size_t n) {
  size_t k = 0, i = 0, j = 0;
  double *res;
  res = (double*)malloc(n * n * sizeof(double));
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      res[i * n + j] = 0.;
    }
  }
  for (i = 0; i < n; ++i) {
    for(j = 0; j < n; ++j) {
      for (k = 0; k < n; ++k) {
        res[i * n + j] += a[i * n + k] * b[k * n + j];
      }
    }
  }
  return res;
}

double EvaluateAccuracy(const double *a, const double *b, const size_t n) {
  double *e, *t;
  double sum = 1e10, tmp = 0.;
  size_t i = 0, j = 0;
  e = (double*)malloc(n * n * sizeof(double));
  InitializeIdentityMatrix(e, n);
  t = GetMatrixProduct(a, b, n);
    /*
    fprintf(stdout, "A*A^-1\n");
    PrintArray(t, n);
    */
  for (j = 0; j < n; ++j) { // будет медленно скакать по буферу
    for (i = 0; i < n; ++i) {
      tmp += fabs(t[i * n + j] - e[i * n + j]);
    }
    if (sum > tmp) {
      sum = tmp;
    }
    tmp = 0.;
  }
  free(e);
  free(t);
  return sum;
}

void GenerateYakobi(double *a, double c, double d, const size_t n) {
  size_t i = 0, j = 0;
  double tmp = 0.;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i == j) {
        tmp = d;
      }
      else if ((j == i + 1) || (i == j + 1)) {
        tmp = c;
      }
      else {
        tmp = 0.;
      }
      a[i * n + j] = tmp;
    }
  }
}

void PrintTime(const double *time1, const double *time2) {
  printf("%f\n", *time1);
  printf("%f\n", *time2);
  printf("time2 / time1: %f\n", (*time2)/(*time1));
}

void PrintArray(const double *arr, const size_t n, const long int how_many_print) {
  long int i = 0, j = 0;
  for (i = 0; i < how_many_print; ++i) {
    for (j = 0; j < how_many_print; ++j) {
      fprintf(stdout, "%e ", arr[i * n + j]);
    }
    fprintf(stdout, "\n");
  }
}
