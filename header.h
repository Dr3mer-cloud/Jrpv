#ifndef _MATRIX_HEADER_H_
#define _MATRIX_HEADER_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <time.h>

typedef struct {
	union {
		long int l;
		double d;
	} val;
	int has_value;
} maybe_t;

maybe_t FieldToInteger(const char *str);
maybe_t FieldToDouble(const char *str);

#endif
