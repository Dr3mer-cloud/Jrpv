#include "header.h"
#include "solve.h"
#include "tests.h"

int PrintError(int res, char* info) {
  if (res >= 0) return res;
  if (!info) {
    info = (char*)"";
  }
  switch (res) {
    case ERROR_OPEN: printf("Cannot open file %s\n", info); break;
    case ERROR_READ: printf("Cannot read file %s\n", info); break;
    case ERROR_WRITE: printf("Cannot write to file %s\n", info); break;
    case ERROR_MEMORY: printf("Not enough memory\n"); break;
    case ERROR_STRING: printf("Incorrect string %s\n", info); break;
    case ERROR_MATRIX: printf("Incorrect matrix\n"); break;
    default: printf("Unknown error %d\n",res); break;
  }
  return res;
}

int CreateMatrices(double **a, double **b, double **c, size_t size) {
  *a = (double*)malloc(size * size * sizeof(double));  // работяга
  *b = (double*)malloc(size * size * sizeof(double)); // хранит исходную матрицу
  *c = (double*)malloc(size * size * sizeof(double));  // в нем будет обратная
  if (a == NULL || b == NULL || c == NULL) {
    fprintf(stderr, "Failed to provide memory\n");
    free_3(a, b, c);
    return ERROR_MEMORY;
  }
  return SUCCESS;
}


void PrintBefore(void *a, void *b, size_t n, size_t how_many_print) {
  printf("\nBefore calling Jordan\n");
  fprintf(stdout, "A:\n");
  PrintArray(a, n, how_many_print);
  fprintf(stdout, "\nB ( = A^(-1)) :\n");
  PrintArray(b, n, how_many_print);
}

void PrintAfter(void *a, void *b, size_t n, size_t how_many_print) {
  printf("\n\nAfter calling Jordan\n");
  fprintf(stdout, "\nA:\n");
  PrintArray(a, n, how_many_print);
  fprintf(stdout, "\nB ( = A^(-1)) :\n");
  PrintArray(b, n, how_many_print);
}

void PrintAfterRearrangement(void *a, void *b, size_t n, size_t how_many_print) {
  printf("\n\n\nAfter rearrangement\n");
  fprintf(stdout, "\nA:\n");
  PrintArray(a, n, how_many_print);
  fprintf(stdout, "\nB ( = A^(-1)) :\n");
  PrintArray(b, n, how_many_print);
}

void free_3(void *a, void *b, void *c) {
  if (a) free(a);
  if (b) free(b);
  if (c) free(c);
}

double ColumnNorm(const double *a, const size_t n) {
  size_t i = 0, j = 0;
  double max_sum = 0., sum = 0.;
  for (j = 0; j < n; ++j) {
    sum = 0.;
    for (i = 0; i < n; ++i) {
      sum += fabs(a[i * n + j]);
    }
    if (sum > max_sum) max_sum = sum;
  }
  return max_sum;
}

int JordanWithPivot(double *a, double *b, const size_t n, const double norm) {
  double tmp = 0., eps = 1e-16 * norm, dif = 0.;
  size_t k = 0, i = 0, j = 0, max_row = 0, count_zeros = 0; // if matrix is singular, i would get string with only zeros
  double max_in_column = 0.;
  int *used_strings = (int*)malloc(n * sizeof(int));

  for (i = 0; i < n; ++i) { used_strings[i] = 0; }

  for (k = 0; k < n; ++k) { 
    max_in_column = 0.;
    // find row with max element, swap rows 
    for (i = 0; i < n; ++i) { // in submatrix cycle through rows THERE WAS i = k, I CHANGED IT 
                                  // SO NOW IT GETS ZEROS EXCEPT ONE 1 IN THE COLUMN, HAVE TO PROPERLY REVERSE THE MATRIX
      if (((used_strings[i] == 0) && fabs(a[i * n + k]) > max_in_column)) {
        max_in_column = fabs(a[i * n + k]); 
        max_row = i;
      }
    }
    used_strings[max_row] = 1;
    /*
    printf("\nk: %lu, max_row: %ld, max_in_column: %.2f\n", k, max_row, max_in_column);
    printf("used_strings: "); for (j = 0; j < n; ++j) printf("%d ", used_strings[j]); printf("\n\n");
    PrintArr(a, n);
    */

    if (max_in_column < eps) {
      fprintf(stderr, "max_in_column = %.2e < eps, matrix is singular\n", max_in_column);
      free(used_strings);
	    return -1; // column with zeros - matrix is singular
    }
    tmp = a[max_row * n + k]; // divide by this element to get 1 at first position
    //printf("-------------------------------------------a[max_row * n + k] = %.2f\n", tmp);
    for (j = 0; j < n; ++j) { // divide current k-th string by this diagonal element
                                  // maybe i can divide starting from k ?
      a[max_row * n + j] /= tmp;
      b[max_row * n + j] /= tmp;
    }
    /*
    printf("\nAfter division in row by this element and before adding to the other strings:\n");
    PrintArr(a, n);
    */
    for (i = 0; i < n; ++i) {
      // i < k -- reset to zero elements above current k-th string
      // i > k -- reset to zero elements below current k-th string
      count_zeros = 0;
      if (i == max_row) { continue; } // aii = akk, i turning elements to zero with akk, so should keep it 
      // in each row save element in k-th column (first element)
      dif = a[i * n + k]; // first non-zero element in that current (means string that will be subtracted from others) string
                          // rewriting this position will break everything
                          // to reset first element in current string to zero, i need to add k-th string
                          // multiplied by this dif (so a[i * n + k] --> 0)
      for (j = k; j < n; ++j) { // now through columns in current i-th row, starting from k - previous elements should be equal to zero
        a[i * n + j] -= dif * a[max_row * n + j]; // akj is element at the same column in k-th string
        if (fabs(tmp) < eps) { ++count_zeros; } //tmp == 0 
      }
      // and in this string the same for matrix b
      for (j = 0; j < n; ++j) {
        b[i * n + j] -= dif * b[max_row * n + j]; // there was b[k * n + j]
      }
      if (count_zeros == n) { return -1; }
    }
    /*
    printf("\nAfter adding to the other strings:\n");
    PrintArr(a, n);
    */
  }
  free(used_strings);
  return 0;
}


int Jordan(double *a, double *b, const size_t n, const double norm) {
  double tmp = 0., eps = 1e-16 * norm, dif = 0.;
  size_t k = 0, i = 0, j = 0, count_zeros = 0; // if matrix is singular, i would get string with only zeros
  for (k = 0; k < n; ++k) { 
    tmp = a[k * n + k];
    if (fabs(tmp) < eps) { // diagonal element = 0 - nothing to do (in simple gauss elimination)
      continue;
    }
    for (j = 0; j < n; ++j) { // divide current k-th string by this diagonal element
      a[k * n + j] /= tmp;
      b[k * n + j] /= tmp;
    }
    for (i = 0; i < n; ++i) {
      // i < k -- reset to zero elements above current k-th string
      // i > k -- reset to zero elements below current k-th string
      count_zeros = 0;
      if (i == k) { continue; } // aii = akk, i turning elements to zero with akk, so should keep it 
      // in each row save element in k-th column (first element)
      dif = a[i * n + k]; // first non-zero element in that current string, (if algorithm can be used)
                          // rewriting this position will break everything
                          // to reset first element in current string to zero, i need to add k-th string
                          // multiplied by this dif (so a[i * n + k] --> 0)
      for (j = k; j < n; ++j) { // now through columns in current i-th row, starting from k - previous elements should be equal to zero
        a[i * n + j] -= dif * a[k * n + j]; // akj is element at the same column in k-th string
        if (fabs(tmp) < eps) { ++count_zeros; } //tmp == 0 
      }
      // and in this string the same for matrix b
      for (j = 0; j < n; ++j) {
        b[i * n + j] -= dif * b[k * n + j];
      }
      if (count_zeros == n) { 
        fprintf(stderr, "%lu string with only zeros\n", i);
        return -1;
      }
    }
  }
  /*
  fprintf(stdout, "A\n");
  PrintArray(a, n);
  fprintf(stdout, "B\n");
  PrintArray(b, n);
  fprintf(stdout, "\n");
  */
  return 0;
}

int RearrangeStrings(double *a, double *b, const size_t n) {
  double *arr_a, *arr_b, *temp_a, *temp_b;
  double temp = 0.;
  size_t i = 0, j = 0, units_counter = 0, zero_counter = 0, pos = 0;
  long int string_number = -1;

  arr_a = (double*)malloc(n * sizeof(double));
  arr_b = (double*)malloc(n * sizeof(double));
  temp_a = (double*)malloc(n * n * sizeof(double)); // to store initial state
  temp_b = (double*)malloc(n * n * sizeof(double));
  if (arr_a == NULL || arr_b == NULL || temp_a == NULL || temp_b == NULL) {
    fprintf(stderr, "\nFailed to provide memory\n");
    return -1;
  }

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      pos = i * n + j;
      temp_a[pos] = a[pos];
      temp_b[pos] = b[pos];
    }
  }

  for (i = 0; i < n; ++i) { // find value = 1. in current i-th string
    string_number = -1;     // for check if it hasn't been updated - error or detA = 0;          
    zero_counter = 0; units_counter = 0; // for new string new counts
    for (j = 0; j < n; ++j) {
      temp = fabs(temp_a[i * n + j]); // element in i-th string
      if (temp > 1e-14) { // deviation from 0
        if (fabs(temp - 1) > 1e-14) { // deviation from 1
          fprintf(stderr, "\nA matrix has element with deviation from zero and one > 1e-14");
          return -1; 
        } else { // fabs(temp - 1) < 1e-14 i.e. it's close to 1
          ++units_counter;
          string_number = j; // number of j defines number of string where 1 on j-th j-th position in matrix
          //printf("\ni = %lu string_number = %lu\n", i, string_number);
        }
      }
      if (temp < 1e-14) { 
        ++zero_counter;
      }
    }
    
    if (string_number < 0) { 
      fprintf(stderr, "\nno elements equal 1 in matrix\n");
      return -1;
    }
    if (zero_counter != (n - 1)) {
      fprintf(stderr, "\nzero_counter in %lu string != (n - 1)", i);
      fprintf(stderr, "\nzero_counter = %lu, (n - 1) = %lu", zero_counter, (n-1));
      return -1;
    }
    if (units_counter != 1) {
      fprintf(stderr, "\nunits_counter in %lu string != 1", i);
      fprintf(stderr, "\nstring = %lu units_counter = %lu", i, units_counter);
      return -1;
    }

    for (j = 0; j < n; ++j) { // copying
      pos = i * n + j; // i because i'm copying current string
      arr_a[j] = temp_a[pos]; //////////////////////// copying that element in temporary array arr_a
      arr_b[j] = temp_b[pos];
    }


    // now i know string number and i-th string is correct - copying it from *arr_a to the *a - result matrix
    for (j = 0; j < n; ++j) {
      pos = string_number * n + j; // string number because it's proper index in matrix A for that string with 1 in position number_string
      a[pos] = arr_a[j];
      b[pos] = arr_b[j];
    }
  }
  free(arr_a);
  free(arr_b);
  free(temp_a);
  free(temp_b);
  return 0;
}

/*
int Jordan_with_pivot(double *a, double *b, const size_t n) {
  double tmp = 0., eps = 1e-14, dif = 0.;
  size_t k = 0, i = 0, j = 0, max_row = 0, count_zeros = 0; // if matrix is singular, i would get string with only zeros
  double max_in_column = 0.;
  long int *used_strings = (long int*)malloc(n * sizeof(long int));
  for (i = 0; i < n; ++i) { used_strings[i] = 0; }
  for (k = 0; k < n; ++k) { 
    // find row with max element, swap rows 
    for (i = 0; i < n; ++i) { // in submatrix cycle through rows THERE WAS i = k, I CHANGED IT 
                                  // SO NOW IT GETS ZEROS EXCEPT ONE 1 IN THE COLUMN, HAVE TO PROPERLY REVERSE THE MATRIX
      if (fabs(a[i * n + k]) > max_in_column) {
        max_in_column = fabs(a[i * n + k]); 
        max_row = i;
      }
    }
    printf("\nk: %lu, max_row: %ld, max_in_column: %.2f\n", k, max_row, max_in_column);
    PrintArr(a, n);

    if (max_in_column < eps) return -1; // column with zeros - matrix is singular
    tmp = a[max_row * n + k]; // divide by this element to get 1 at first position
    printf("-------------------------------------------a[max_row * n + k] = %.2f\n", tmp);
    if (fabs(tmp) < eps) { // diagonal element = 0 - nothing to do (in simple gauss elimination)
      continue;
    }
    for (j = 0; j < n; ++j) { // divide current k-th string by this diagonal element
                                  // maybe i can divide starting from k ?
      a[max_row * n + j] /= tmp;
      b[max_row * n + j] /= tmp;
    }
    printf("\nAfter division in row by this element and before adding to the other strings:\n");
    PrintArr(a, n);
    for (i = 0; i < n; ++i) {
      // i < k -- reset to zero elements above current k-th string
      // i > k -- reset to zero elements below current k-th string
      if (i == max_row) { continue; } // aii = akk, i turning elements to zero with akk, so should keep it 
      // in each row save element in k-th column (first element)
      dif = a[i * n + k]; // first non-zero element in that current (means string that will be subtracted from others) string
                          // rewriting this position will break everything
                          // to reset first element in current string to zero, i need to add k-th string
                          // multiplied by this dif (so a[i * n + k] --> 0)
      for (j = k; j < n; ++j) { // now through columns in current i-th row, starting from k - previous elements should be equal to zero
        a[i * n + j] -= dif * a[max_row * n + j]; // akj is element at the same column in k-th string
        if (fabs(tmp) < eps) { ++count_zeros; } //tmp == 0 
      }
      // and in this string the same for matrix b
      for (j = 0; j < n; ++j) {
        b[i * n + j] -= dif * b[max_row * n + j]; // there was b[k * n + j]
      }
      if (count_zeros == n) { return -1; }
      count_zeros = 0;
    }
    printf("\nAfter adding to the other strings:\n");
    PrintArr(a, n);
    max_in_column = 0.; // only if i use all declarations before starting implementation, else just define after (k=0; ...) string
  }
  free(used_strings);
  return 0;
}
*/

