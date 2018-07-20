#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "help_methods.h"

#define DAT_FILE "hollins.dat"
#define MATRIX_FILE "matrix.dat"


struct timeval startwtime, endwtime;
double seq_time;

int main() {

    int N; // number of nodes
    int num_of_lines;
    int i, j, k;

    gettimeofday (&startwtime, NULL);
    /** open the nodes file to obtain the number of nodes **/
    FILE *fnodes;
    fnodes = fopen(DAT_FILE,"r");
    if (fnodes == 0) {
        perror("failed to open input file\n");
        exit(EXIT_FAILURE);
    }
    fscanf(fnodes,"%d",&N);
    printf("Number of nodes: %d \n", N);
    fclose(fnodes);

    /** allocate memory for the adjacency matrix E **/
    double **E = (double **)malloc(N*sizeof(double *));
    for (i = 0; i < N; i ++){
        E[i] = (double *)malloc(N*sizeof(double));
        for (j = 0; j < N; j ++){
            E[i][j] = 0;
        }
    }

    /** open the file again to complete the adjacency matrix **/
    FILE *flist;
    flist = fopen(DAT_FILE,"r");
    if (flist == 0) {
        perror("failed to open input file\n");
        exit(EXIT_FAILURE);
    }
    // excluding the first two elements (number of nodes and number of neighbors)
    fseek(flist, sizeof(int), SEEK_SET);
    fscanf(flist, "%d", &num_of_lines);
    for (i = 0; i < num_of_lines; i ++){
        fscanf(flist,"%d %d",&j, &k);
        if ((j >= 0 && j < N) && (k >= 0 && k < N)){
            E[j-1][k-1] = 1;
        }
    }
    fclose(flist);

    gettimeofday (&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);
    printf("Reading data wall clock time = %f\n", seq_time);

//    /*** print the adjacency matrix ***/
//    FILE *fmatrix;
//    fmatrix = fopen(MATRIX_FILE,"w");
//    if (fmatrix == 0) {
//        perror("failed to open output file\n");
//        exit(EXIT_FAILURE);
//    }
//    for (i = 0; i < N; i ++){
//        for (j = 0; j < N; j ++){
//            fprintf(fmatrix,"%f ", E[i][j]);
//        }
//        fprintf(fmatrix,"\n");
//    }
//    fclose(fmatrix);


    struct Options options;
    options.tolerance = 1.0e-7;
    options.maxiter = 500;
    options.c = 0.85;
    /** construct personalization vector v **/
    double* temp;
    temp = (double*)malloc(N* sizeof(double));
    for (i = 0; i < N; i++) {
        temp[i] = (double)(1) / (double)(N);
    }
    options.v = temp;

    /** d is the n-dimensional column vector identifying the nodes with outdegree 0 **/
    double* d = (double*)malloc(N* sizeof(double));
    /** normalize matrix && identify rows with outdegree 0 **/
    for (i = 0; i < N; i++){
        double sum = 0;
        for (j = 0; j < N; j++) {
            sum = sum + E[i][j];
        }
        if (sum > 0){
            for (j = 0; j < N; j++){
                E[i][j] = (E[i][j] / sum);
            }
        }else if (sum == 0){
            d[i] = 1;
        }
    }

    /**
     * construct P' (transition possibility matrix) in the following manner
     * D = d*options.v_transpose
     * P' = P + D (here P = E)
     */


    /** allocate memory for P and D matrices **/
    double **P = (double **)malloc(N*sizeof(double *));
    for (i = 0; i < N; i ++){
        P[i] = (double *)malloc(N*sizeof(double));
        for (j = 0; j < N; j ++){
            P[i][j] = 0;
        }
    }
    double **D = (double **)malloc(N*sizeof(double *));
    for (i = 0; i < N; i ++){
        D[i] = (double *)malloc(N*sizeof(double));
        for (j = 0; j < N; j ++){
            D[i][j] = 0;
        }
    }
    for (i = 0; i < N; i ++){
        for (j = 0; j < N; j ++){
            D[i][j] = d[i]*options.v[j];
            P[i][j] = E[i][j] + D[i][j];
        }
    }


    int flag = 0;
    for (i = 0; i < N; i ++) {
        double sum = 0;
        for (j = 0; j < N; j++) {
            sum = sum + P[i][j];
        }
        if (sum == 0){
            flag =1;
        }
    }
    if (flag == 1){
        printf("ERROR \n");
    }

    /** construct irreducible Markov matrix P"
     * as follows
     * E = e * options.v_transpose (e is a N size vector with ones.)
     * P" = options.c * P' + (1 - options.c)*E
     * (here P' = P)
     **/

    /** allocate memory for the teleportation matrix E_t **/
    double **E_t = (double **)malloc(N*sizeof(double *));
    for (i = 0; i < N; i ++){
        E_t[i] = (double *)malloc(N*sizeof(double));
        for (j = 0; j < N; j ++){
            E_t[i][j] = 1 * options.v[j];
        }
    }
    for (i = 0; i < N; i ++){
        for (j = 0; j < N; j ++){
            P[i][j] = options.c * P[i][j] + (1 - options.c)*E_t[i][j];
        }
    }

    /**
     * for the Gauss-Seidel method we are solving the following equation
     * ([I]NxN - (options.c)*P.transpose) * x = ((1-options.c)/N )* [1]N
     * here : ([I]NxN - (options.c)*P.transpose) = K
     */
    double **K = (double **)malloc(N*sizeof(double *));
    for (i = 0; i < N; i ++){
        K[i] = (double *)malloc(N*sizeof(double));
        for (j = 0; j < N; j ++){
            K[i][j] = 0;
        }
    }
    for (i = 0; i < N; i ++) {
        for (j = 0; j < N; j++) {
            if (i==j){
                K[i][j] = 1 - (options.c * P[j][i]);
            }
            else{
                K[i][j] = (-1)* (options.c * P[j][i]);
            }

        }
    }

    /** GAUSS - SEIDEL IMPLEMENTATION **/
    double* x = (double*)malloc(N* sizeof(double));
    // initializing x as the personalization vector
    for (i=0; i<N; i++){
        x[i] = options.v[i];
    }

    double* x_old = (double*)malloc(N* sizeof(double));

    // initialize delta
    double delta = 1;
    int iter = 0;

    gettimeofday (&startwtime, NULL);

    double* diff = (double*)malloc(N* sizeof(double));

    double sigma;
    double t;
    while (delta > options.tolerance && iter < options.maxiter){

        /** gauss - seidel iteration **/

        for (j=0;j<N;j++){
            x_old[j] = x[j];
        }


        for (i = 0; i < N; i++){
            double sigma = 0;
            for (j=0; j < i-1; j++){
                sigma = sigma + (K[i][j]*x[j]);
            }
            for (j=i+1; j < N; j++){
                sigma = sigma + (K[i][j]*x_old[j]);
            }
            x[i] = (double)((((1- options.c)/N) - sigma));

        }
        subtract(x, x_old, N, diff);
        delta = norm(diff, N);
        printf("%f \n", delta); // nan
        iter = iter + 1;


    }
    gettimeofday (&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);
    printf("Serial PageRank wall clock time = %f\n", seq_time);
    double sum1 = 0;
    for (i=0;i<N;i++){
        sum1 = sum1 + x[i];
        printf("%f\n", x[i]);
    }
//    printf("sum %f \n", sum1);


    if (delta > options.tolerance || iter == options.maxiter){
        printf("Algorithm did not converged after %d iterations \n", iter);
    }else{
        printf("Converged after %d iterations. \n", iter);
    }


}