#ifndef WORKDIR_HELP_METHODS_H
#define WORKDIR_HELP_METHODS_H

struct Options;
struct Options {
    double tolerance;
    int maxiter;
    double* v;
    double c;
};
double norm(double* a1, int N);
void subtract(double* a1, double* a2, int N, double * diff);

#endif //WORKDIR_HELP_METHODS_H