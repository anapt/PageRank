#include "help_methods.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double norm(double* a1, int N){
    double sum = 0;
    for (int i = 0; i < N; ++i) {

        double a = a1[i] * a1[i];
        sum = sum + a;
    }
    return (sqrt(sum));
}

void subtract(double* a1, double* a2, int N, double * diff){

    for (int i = 0; i < N; ++i) {
        diff[i] = a1[i] - a2[i];
    }
}

