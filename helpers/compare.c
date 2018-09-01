#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PARALLEL_FILE "./output/parallel.txt"
#define SERIAL_FILE "./output/serial.txt"

#define N 6012

#define TOLERANCE 0.001

int main() {

    FILE *f_parallel, *f_serial;

    f_parallel = fopen(PARALLEL_FILE,"r");
    f_serial = fopen(SERIAL_FILE,"r");
    if (f_serial == 0 || f_parallel == 0) {
        perror("[ERROR]: Failed to open input files\n");
        exit(EXIT_FAILURE);
    }

    int i=0;
    double gval, pval, diff, sum=0;

    for(i=0; i<N; i++){
        fscanf(f_parallel, "%lf", &pval);
        fscanf(f_serial, "%lf", &gval);

        diff = pval - gval;
        sum = sum + diff * diff;
    }
    fclose(f_parallel);
    fclose(f_serial);

    double norm = sqrt(sum);

    if (norm > TOLERANCE){
        printf("[INFO]: Output files have significant differences. Norm is %lf\n", norm);
    }
    else{
        printf("[INFO]: Output files are identical. Norm is %lf\n", norm);
    }

}
