#include "header.h"
#include "tests.h"
#include "solve.h"

maybe_t FieldToDouble(const char *str);
maybe_t FieldToInteger(const char *str);
int ProvideMemoryAndCallJordan(const size_t n, const long int how_many_print, const int formula_number, const int is_file, 
                                const int Jordan_with_pivot);

FILE *fin;

maybe_t FieldToDouble(const char *str) {
  double d;
  char *e; 
  maybe_t result;
  errno = 0;
  d = strtod(str, &e);
  if (!errno && (*e == '\0')) {
    result.has_value = 1;
    result.val.d = d;
    return result;
  }
  result.has_value = 0;
  return result;
}

maybe_t FieldToInteger(const char *str) {
  int i = 0;
  char *e;
  maybe_t result;
  errno = 0;
  i = strtol(str, &e, 10);
  if ((!errno) && (*e == '\0')) { 
    if ((INT_MIN <= i) && (i <= INT_MAX)) {
      result.val.l = i;
      result.has_value = 1;
    }
    else {
      result.has_value = -1;
    }
    return result;
  }
  result.has_value = 0;
  return result;
}

int ProvideMemoryAndCallJordan(const size_t n, const long int how_many_print, const int formula_number, const int is_file, 
                                const int Jordan_with_pivot) {
  // Jordan_with_pivot == 1 means i call default Jordan
  size_t i = 0, j = 0, pos = 0; // pos for not to recount (i * n + j)
  double *a = NULL, *b = NULL, *a0 = NULL;
  double tmp = 0., time = 0., norm = 0.;

  clock_t start, end;
  int status = INT_MIN;

  status = CreateMatrices(&a, &b, &a0, n); if (status < 0) return -1;

  InitializeIdentityMatrix(b, n);
  
  switch(is_file) {
    case 1: // input has file
      // fill matrices a and a0
      for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
          fscanf(fin, "%lf", &tmp);
          pos = i * n + j;
          a[pos] = tmp;
          a0[pos] = tmp;
        }
      }
      break;
    default: // input doesn't have file
      for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
          pos = i * n + j; // save for not to recount
          switch (formula_number) {
            case 1:
              a[pos] = n - (i > j ? i : j) + 1; break;
            case 2:
              a[pos] = (i > j ? i : j); break;
            case 3:
              a[pos] = (i > j ? (i - j) : (j - i)); break;
            case 4:
              a[pos] = 1/(i + j - 1); break;
            default:
              fprintf(stderr, "formula_number = %d; < 1 or > 4\n", formula_number);
              return -1;
          }
          a0[pos] = a[pos]; ///////// fill a0 matrix
        }
      }

  }

  if (Jordan_with_pivot == 0) printf ("\nThis is default Jordan//////////////////////////////////////////////////////////////////\n");
  else printf("\nJordan with pivot//////////////////////////////////////////////////////////////////\n");

  PrintBefore(a, b, n, how_many_print);
  
  norm = ColumnNorm(a, n); // calculate here, doesn't have to take Jordan function runtime
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  start = clock();
  status = (Jordan_with_pivot == 1 ? JordanWithPivot(a, b, n, norm) : Jordan(a, b, n, norm));
  end = clock();
  time = (double)(end - start) / CLOCKS_PER_SEC;
  if (status == -1) {
    fprintf(stderr, "detA = 0\n");
    free_3(a, a0, b);
    if (fin) { rewind(fin); }
    return -1;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  PrintAfter(a, b, n, how_many_print);

  fprintf(stdout, "\nTime: %.15e\n", time);
  
  if (Jordan_with_pivot == 1) { // call that only Jordan_with_pivot was called
    RearrangeStrings(a, b, n);
    PrintAfterRearrangement(a, b, n, how_many_print);
  }

  tmp = EvaluateAccuracy(a0, b, n);

  fprintf(stdout, "\n||A*A^-1 - E||: %.15e\n", tmp);
  free_3(a, a0, b);
  if (fin) { rewind(fin); }
  
  return SUCCESS;
}

int main(int argc, char **argv) {
  double tmp = 0;
  maybe_t check;
  size_t n = 0, cnt = 0;
  long int temp_n, how_many_print = 0, formula_number = 0;
  switch (argc) {
  case 1:
    fprintf(stdout, "Usage: ./program_name matrix_size how_many_print formula_number filename\n");
    return 0;

  case 5:
    check = FieldToInteger(argv[1]);
    if (check.has_value == -1 || check.has_value == 0) {
      fprintf(stderr, "n is not integer\n");
      return -1;
    }
    n = (size_t)(check.val.l);

    check = FieldToInteger(argv[2]);
    if (check.has_value == -1 || check.has_value == 0) {
      fprintf(stderr, "how_many_print is not integer\n");
      return -1;
    }
    how_many_print = check.val.l;

    check = FieldToInteger(argv[3]);
    if (check.has_value == -1 || check.has_value == 0) {
      fprintf(stderr, "formula number is not integer\n");
      return -1;
    }
    formula_number = check.val.l;

    if (formula_number != 0) {
      fprintf(stderr, "filename exists, but formula_number isn't zero\n");
      return -1;
    }

    fin = fopen(argv[4], "r");
    if (fin == NULL) {
      fprintf(stderr, "File wasn't open\n");
      return -1;
    }
    while (!feof(fin)) { // feof checks end-of-file indicator and returns 0 if end of file has been reached
      temp_n = fscanf(fin, "%lf", &tmp);
      if (temp_n == 1) { // 1 element has been scanned and writed to tmp
        ++cnt;
      }
      else if (temp_n != EOF) { // EOF is macros that expands to -1 and indicates end of file condition or error condition
        fprintf(stderr, "File messed up, read indicator variable = %ld (!= EOF)\nwere read %ld variables\nstore variable = %e\n", temp_n, cnt, tmp);
        return -1;
      }
    }
    temp_n = sqrt(cnt);
    if ((long int)(cnt) != temp_n * temp_n) {
      fprintf(stderr, "Matrix wasn't formed properly\n");
      return -1;
    }
    rewind(fin); // sets file pointer to the beginning

    //int ProvideMemoryAndCallJordan(const int n, const int how_many_print, const int formula_number, const int is_file, 
    //                            const int JordanWithPivot);
    ProvideMemoryAndCallJordan(n, how_many_print, formula_number, 1, 1);
    ProvideMemoryAndCallJordan(n, how_many_print, formula_number, 1, 0);

    fclose(fin);
    return SUCCESS;

  case 4:

    check = FieldToInteger(argv[1]);
    if (check.has_value == -1 || check.has_value == 0) {
      fprintf(stderr, "n is not integer\n");
      return -1;
    }
    n = check.val.l;

    check = FieldToInteger(argv[2]);
    if (check.has_value == -1 || check.has_value == 0) {
      fprintf(stderr, "how_many_print is not integer\n");
      return -1;
    }
    how_many_print = check.val.l;

    check = FieldToInteger(argv[3]);
    if (check.has_value == -1 || check.has_value == 0) {
      fprintf(stderr, "formula number is not integer\n");
      return -1;
    }
    formula_number = check.val.l;
    //int ProvideMemoryAndCallJordan(const int n, const int how_many_print, const int formula_number, const int is_file, 
    //                            const int JordanWithPivot);
    // first Jordan_with_pivot
    // 2 formula gives 0 in position (0, 0)
    ProvideMemoryAndCallJordan(n, how_many_print, formula_number, 0, 1);
    n *= 2;
    ProvideMemoryAndCallJordan(n, how_many_print, formula_number, 0, 1);
    //ProvideMemoryAndCallJordan(n, how_many_print, formula_number, 0, 0); 
    
    return SUCCESS;

    default:
      fprintf(stderr, "Extra or not enough input data, terminated\n");
      return -1;
    }

}

