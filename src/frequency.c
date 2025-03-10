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
 *      = \sum_k X[k] e^{(2pi i / n) * f * k} / n
 *
 * S(-f) = \sum_k X[k] e^{(2pi i / n) * -f * k} / n
 *       = \sum_k X[k] e^{(2pi i / n) * -f * k + 2pi i * n / n} / n
 *       = \sum_k X[k] e^{(2pi i / n) * (n - f) * k} / n
 *       = S(n - f)
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

bool prob_zero(double complex z, double eps)
{
    return (cabs(z) < eps);
}

double s(double x, int f)
{
    /**
     * Frequency f and -f are 0.5i, -0.5i, the rest is 0 because
     * 0.5i e^(-2 pi f x) - 0.5i e^(-2 pi(-f)x) = 
     * 0.5i (cos(-2pifx) + i sin(-2pifx)) - 0.5i (cos(2pifx) + i sin(2pifx))
     * = -0.5 i^2 * 2 * sin(2 pi f x) = sin(2pifx)
     *
     * Hm, what about the 2 pi?
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
        double complex ck = S[k] / n;
        if (!prob_zero(ck, 1e-15)) {
            printf("Frequency %d: %e + %e i\n", 
                   k,
                   creal(ck),
                   cimag(ck));
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(S);

    return EXIT_SUCCESS;
}
