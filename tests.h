#ifndef _TESTS_H_
#define _TESTS_H_

void InitializeIdentityMatrix(double *p, const size_t n);
double* GetMatrixProduct(const double *a, const double *b, const size_t n);
double EvaluateAccuracy(const double *a, const double *b, const size_t n);
void GenerateYakobi(double *a, double c, double d, const size_t n);
void PrintTime(const double *time1, const double *time2);
void PrintArray(const double *arr, const size_t n, const long int how_many_print);

#endif
