#ifndef GRADIENT_DESCENT_GRADIENT_H
#define GRADIENT_DESCENT_GRADIENT_H

void gr_minimize(
		double (*fn)(void *, const double *),
		void *state,
		const double *start_vector,
		double *ret_vector,
		const int argn,
		const double start_step,
		const double step_mul,
		const int iterations);

#endif /* GRADIENT_DESCENT_GRADIENT_H */
