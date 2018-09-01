#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "../helpers/help_methods.h"

#define DAT_FILE "./dataset/hollins.dat"


struct timeval startwtime, endwtime;
double seq_time;

int main() {

    int NUM_OF_NODES;
    int NUM_OF_LINES;
    int i, j;


    gettimeofday (&startwtime, NULL);

    /** open the nodes file to obtain the number of nodes **/
    FILE *f;
    f = fopen(DAT_FILE,"r");
    if (f == 0) {
        perror("failed to open input file\n");
        exit(EXIT_FAILURE);
    }

    fscanf(f,"%d",&NUM_OF_NODES);
    fscanf(f, "%d", &NUM_OF_LINES);
    printf("Number of nodes: %d \n", NUM_OF_NODES);
    printf("Number of lines: %d \n", NUM_OF_LINES);

    /** allocate memory for the adjacency matrix E **/
    double **E = (double **) malloc(NUM_OF_NODES * sizeof(double *));
    for (i = 0; i < NUM_OF_NODES; i++){
        E[i] = (double *)malloc(NUM_OF_NODES * sizeof(double));
        for (j = 0; j < NUM_OF_NODES; j++){
            E[i][j] = 0;
        }
    }
    
    int src, dst;
    for (i = 0; i < NUM_OF_LINES; i++){
        fscanf(f,"%d %d",&src, &dst);
        if ((src >= 0 && src < NUM_OF_NODES) && (dst >= 0 && dst < NUM_OF_NODES)){
            E[src-1][dst-1] = 1;
        }
    }
    fclose(f);

    gettimeofday (&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);
    printf("Reading data wall clock time = %f\n", seq_time);




    /** construct personalization vector v **/
    double *v = (double*) malloc(NUM_OF_NODES * sizeof(double));
    for (i = 0; i < NUM_OF_NODES; i++) {
        v[i] = (double)(1) / (double)(NUM_OF_NODES);
    }

    /** define options struct **/
    struct Options options;
    options.tolerance = 1.0e-7;
    options.maxiter = 500;
    options.c = 0.85;
    options.v = v;


    /** d is the n-dimensional column vector identifying the nodes with outdegree 0 **/
    double* d = (double*) malloc(NUM_OF_NODES * sizeof(double));

    /** normalize matrix && identify rows with outdegree 0 **/
    for (i = 0; i < NUM_OF_NODES; i++){
        double sum = 0;

        for (j = 0; j < NUM_OF_NODES; j++) {
            sum += E[i][j];
        }

        if (sum > 0){
            //normalize E matrix
            for (j = 0; j < NUM_OF_NODES; j++){
                E[i][j] /= sum;
            }
            d[i] = 0;
        }
        else if (sum == 0){
            // here node i has outdegree 0
            d[i] = 1;
        }
    }


    /**
     * construct P' (transition possibility matrix) in the following manner
     * D = d*options.v_transpose
     * P' = P + D (here P = E)
     */


    /** allocate memory for P and D matrices **/
    double **D = (double **) malloc(NUM_OF_NODES * sizeof(double *));
    double **P = (double **) malloc(NUM_OF_NODES * sizeof(double *));
    
    for (i = 0; i < NUM_OF_NODES; i++){
        D[i] = (double *) malloc(NUM_OF_NODES * sizeof(double));
        P[i] = (double *) malloc(NUM_OF_NODES * sizeof(double));
        for (j = 0; j < NUM_OF_NODES; j++){
            // D[i][j] = 0;
            // P[i][j] = 0;
            D[i][j] = d[i] * options.v[j];
            P[i][j] = E[i][j] + D[i][j];
        }
    }


    //check that sum of P lines is not zero
    int flag = 0;
    for (i = 0; i < NUM_OF_NODES; i++) {
        double sum = 0;
        for (j = 0; j < NUM_OF_NODES; j++) {
            sum = sum + P[i][j];
        }
        if (sum == 0){
            flag = 1;
        }
    }
    if (flag == 1){
        printf("ERROR \n");
    }

    /** construct irreducible Markov matrix P"
     * as follows
     * E = e * options.v_transpose (e is a NUM_OF_NODES size vector with ones.)
     * P" = options.c * P' + (1 - options.c)*E
     * (here P' = P)
     **/

    for (i = 0; i < NUM_OF_NODES; i++){
        for (j = 0; j < NUM_OF_NODES; j ++){
            P[i][j] = options.c * P[i][j] + (1 - options.c) * options.v[j];            
        }
    }


    /**
     *  The rest of the computations are done using the transpose
     *  of array P.
     *  For simplicity and memory management we are using the same matrix
     *  and changing the sequence of indices
     */



    /************************************ 
    *   GAUSS - SEIDEL IMPLEMENTATION
    *************************************/

    // initializing x as the personalization vector
    double* x = (double*) malloc(NUM_OF_NODES * sizeof(double));
    for (i=0; i<NUM_OF_NODES; i++){
        x[i] = options.v[i];
    }
    
    double* x_old = (double*) malloc(NUM_OF_NODES * sizeof(double));

    // initialize delta
    double delta = 1;
    int iter = 0;

    gettimeofday (&startwtime, NULL);

    double* diff = (double*) malloc(NUM_OF_NODES * sizeof(double));
    while (delta > options.tolerance && iter < options.maxiter){

        /** GAUSS-SEIDEL iteration **/
        for (i=0; i<NUM_OF_NODES; i++){
            x_old[i] = x[i];
        }

//        for (i = 0; i < NUM_OF_NODES; i++){
//            double sigma = 0;
//            for (j=0; j < i-1; j++){
//                sigma = sigma + P[j][i]*x[j];
//            }
//            for (j=i+1; j < NUM_OF_NODES; j++){
//                sigma = sigma + P[j][i]*x_old[j];
//            }
//
////            x[i] = ((1- (options.c))/N)-((sigma));
//            x[i] = (((1-options.c)/N) + (options.c * sigma));
//
////            printf("%f \n", 1/P[i][i]);
////            printf("%f \n", ((1-options.c)/b[i] -(sigma)));
//        }
        /** classic Pagerank **/
        for (i=0; i<NUM_OF_NODES; i++){
            double sum=0;
            for (j=0; j < NUM_OF_NODES; j++){
                sum = sum + P[j][i] * x_old[j];
            }
            x[i] = sum;
        }



        subtract(x, x_old, NUM_OF_NODES, diff);
        delta = norm(diff, NUM_OF_NODES);
        iter++;
    }

    gettimeofday (&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);
    printf("Serial PageRank wall clock time = %f\n", seq_time);


    double sum_final = 0;
    for (i=0; i<NUM_OF_NODES; i++){
        sum_final = sum_final + x[i];
    }
    // printf("sum %f \n", sum_final);
    
    // double y = 1 / sum_final;
    // for (i=0; i<NUM_OF_NODES; i++){
    //     printf("%f\n", x[i] * y);
    // }

    // Write final personalization vector x into .txt file
    f = fopen("./output/classic.txt", "w");
    for (i = 0; i < NUM_OF_NODES; i++) {
      fprintf(f, "%f\n", x[i]);
    }
    fclose(f);

    // for (i=0; i<10; i++){
    //     printf("%f \n", x[i]);
    // }

    if (delta > options.tolerance || iter == options.maxiter){
        printf("Algorithm did not converged after %d iterations", iter);
    }else{
        printf("Converged after %d iterations. \n", iter);
    }
}
