#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "curve-fit/gradient.h"

static void suba_d(double *a, const double *b, size_t len)
{
	size_t i;

	for (i = 0; i < len; i++) {
		a[i] -= b[i];
	}
}

static void mula_d(double *restrict x, const double y, size_t len)
{
	size_t i;

	for (i = 0; i < len; i++) {
		x[i] *= y;
	}
}

static double vec_len(const double *v, size_t len)
{
	double ret = 0.0;
	size_t i;

	for (i = 0; i < len; i++) {
		ret += v[i] * v[i];
	}

	return sqrt(ret);
}

static void calc_derivative(
		double (*fn)(void *, const double *),
		void *state,
		double *at,
		double *step_buf,
		double *ret,
		double step,
		int argn)
{
	double val_at;
	int i;

	val_at = fn(state, at);
	memcpy(step_buf, at, sizeof(*at) * argn);

	for (i = 0; i < argn; i++) {
		if (step == 0.0) {
			ret[i] = 0.0;
			continue;
		}

		step_buf[i] += step;
		ret[i] = (fn(state, step_buf) - val_at) / step;
		step_buf[i] = at[i];
	}
}

static void powm1_d(double *restrict x, size_t len)
{
	size_t i = 0;

	for (i = 0; i < len; i++) {
		x[i] = 1.0 / x[i];
	}
}

void gr_minimize(
		double (*fn)(void *, const double *),
		void *state,
		const double *start_vector,
		double *ret_vector,
		const int argn,
		const double start_step,
		const double step_mul,
		const int iterations)
{
	double step;
	double *derivative = NULL;
	double *step_buf = NULL;
	int i;

	step = start_step * step_mul;
	memcpy(ret_vector, start_vector, sizeof(*ret_vector) * argn);
	step_buf = malloc(sizeof(*step_buf) * argn);
	derivative = calloc(sizeof(*derivative), argn);

	for (i = 0; i < iterations && step != 0.0; i++) {
		calc_derivative(fn, state, ret_vector, step_buf, derivative, step, argn);
		mula_d(derivative, step_mul, argn);
		suba_d(ret_vector, derivative, argn);
		step = vec_len(derivative, argn) * step_mul;
	}

	free(derivative);
	free(step_buf);

	return;
}
