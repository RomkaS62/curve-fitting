#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "curve-fit/gradient.h"

static double mulsum_d(
		const double *restrict x,
		const double *restrict y,
		size_t length)
{
	size_t i;
	double ret = 0.0;

	for (i = 0; i < length; i++) {
		ret += x[i] * y[i];
	}

	return ret;
}

static void power_series(double *ret, const double x0, size_t len)
{
	size_t i;

	if (len == 0) {
		return;
	}

	ret[0] = 1.0;

	if (len <= 1) {
		return;
	}

	ret[1] = x0;

	for (i = 2; i < len; i++) {
		ret[i] = ret[i - 1] * x0;
	}
}

struct polynomial_s {
	double *terms;
	int count;
	int dimensions;
};

static double polynomial_eval_1d(
		const double *terms,
		const double at,
		double *restrict pow_series_buf,
		const size_t length)
{
	power_series(pow_series_buf, at, length);
	return mulsum_d(pow_series_buf, terms, length);
}

static double polynomial_eval(
		const struct polynomial_s *p,
		double *restrict pow_series_buf,
		const double *restrict at)
{
	double ret = 0.0;
	double *dimension_terms;
	int i;

	for (i = 0; i < p->dimensions; i++) {
		dimension_terms = p->terms + p->count * i;
		ret += polynomial_eval_1d(dimension_terms, at[i], pow_series_buf, p->count);
	}

	return ret;
}

void print_da(FILE *out, const double *arr, size_t len)
{
	size_t i;

	for (i = 0; i < len - 1; i++) {
		fprintf(out, "%g ", arr[i]);
	}

	if (len) {
		fprintf(out, "%g\n", arr[i]);
	}
}

struct polynomial_e_s {
	struct polynomial_s p;
	double *psb;
};

static double polynomial_eval_e(struct polynomial_e_s *p, const double *at)
{
	return polynomial_eval(&p->p, p->psb, at);
}

int main(const int argc, const char **argv)
{
	double terms[] = {
		0, 1, 1,
		0, 2, 3
	};

	double at[] = { 2, 4 };
	double min_at[2] = { 0, 0 };
	double psb[3];
	struct polynomial_e_s p = { { terms, 3, 2 }, psb };

	gr_minimize(polynomial_eval_e, &p, at, min_at, 2, 0.3, 0.2, 50);
	print_da(stdout, min_at, 2);

	return 0;
}
