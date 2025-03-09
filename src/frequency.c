/**
 * We want to compute integrals
 *  S(f) = \int s(t) e^{2pi i f t} dt
 * for frequencies f.
 *
 * Discretizing s on [0, 2pi) to an array of n points gives
 * X[k] = s(h * k)
 * for h = 2 pi / n
 *  t = 0, h, 2h, ..., (n - 1) h
 * 
 * S(f) = \sum_k s(h * k) e^{2pi i f (h * k)} * h
 *      = \sum_k X[k] e^{2pi i f k / n} / n
 *
 * Discretizing on general intervals [a, b) gives different exponents, but this
 * is unnecessary. We can transform t -> s(t) to t -> s((t - a) / (b - a))
 * to get an interval [0, 2 pi), and then use time and frequency scaling.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>

bool prob_zero(double complex x, double eps)
{
    double r = creal(x);
    double c = cimag(x);
    return (fabs(r) < eps && fabs(c) < eps);
}

double s(double x, int f)
{
    /**
     * Frequency f and -f are 0.5i, -0.5i, the rest is 0 because
     * 0.5i e^(-2 pi f x) + 0.5i e^(-2 pi(-f)x) = 
     * 0.5i (cos(-2pifx) + i sin(-2pifx)) - 0.5i (cos(2pifx) + i sin(2pifx))
     * = 0.5 i^2 * 2 * - sin(2 pi f x) = sin(2pifx)
     *
     * ?! I guess we get n - f because we discretize on [0, 2 pi) instead of
     * [-pi, pi). What about the 2pi?
     **/
    return sin(f * x);
}

void init(int n, double h, int f, double complex *X)
{
    for (int k = 0; k < n; k++) {
        X[k] = s(k * h, f) + 0 * I;
    }
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        printf("Usage: %s <N> <FREQUENCY>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int n    = atoi(argv[1]);
    int f    = atoi(argv[2]);
    double h = 2.0 * M_PI / n;
    double complex *S = fftw_malloc(n * sizeof(double complex));
    fftw_plan plan = fftw_plan_dft_1d(n, S, S, FFTW_FORWARD, FFTW_ESTIMATE);

    init(n, h, f, S);

    fftw_execute(plan);

    for (int k = 0; k < n; k++) {
        if (!prob_zero(S[k], 1e-10)) {
            printf("Frequency %d: %e + %e i\n", 
                   k,
                   creal(S[k]) / n,
                   cimag(S[k]) / n);
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(S);

    return EXIT_SUCCESS;
}
