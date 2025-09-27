#ifndef _SOLVE_H_
#define _SOLVE_H_

#define SUCCESS       0
#define ERROR        (-1)
#define ERROR_CMDLN  (-2)
#define ERROR_MEMORY (-3)
#define ERROR_OPEN   (-4)
#define ERROR_READ   (-5)
#define ERROR_WRITE  (-6)
#define ERROR_STRING (-7)
#define ERROR_MATRIX (-8)

#define EPS 1e-16

void PrintBefore(void *a, void *b, size_t n, size_t how_many_print);
void PrintAfter(void *a, void *b, size_t n, size_t how_many_print);
void PrintAfterRearrangement(void *a, void *b, size_t n, size_t how_many_print);
int CreateMatrices(double **a, double **b, double **c, size_t size);
void free_3(void *a, void *b, void *c);
int PrintError(int res, char* info);
double ColumnNorm(const double *a, const size_t n);
int Jordan(double *a, double *b, const size_t n, const double norm);
int JordanWithPivot(double *a, double *b, const size_t n, const double norm);
int RearrangeStrings(double *a, double *b, const size_t n);
#endif
