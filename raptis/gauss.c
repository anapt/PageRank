//gcc -Wall gauss-rap.c help_methods.c -o gauss -lm

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "help_methods.h"

#define DAT_FILE "hollins.dat"


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


    /** ========================================== 
    *   Allocate memory for the adjacency matrix E  
    *   ==========================================
    **/
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




    /** ================================== 
    *   Construct personalization vector v 
    *   ==================================
    **/
    double *v = (double*) malloc(NUM_OF_NODES * sizeof(double));
    for (i = 0; i < NUM_OF_NODES; i++) {
        v[i] = (double)(1) / (double)(NUM_OF_NODES);
    }


    /** =====================
    *   Define options struct  
    *   =====================
    **/
    struct Options options;
    options.tolerance = 1.0e-7;
    options.maxiter = 500;
    options.c = 0.85;
    options.v = v;


    /** =========================================
    *   Construct vector d and normalize matrix E 
    *   =========================================
    **/

    /** d is the n-dimensional column vector identifying the nodes with outdegree 0 **/
    double* d= (double*) malloc(NUM_OF_NODES * sizeof(double));

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


    /** ==========================================
    *   Construct transition possibility matrix P'
    *   ==========================================
    * 
    *   Formula:
    *
    *       D = d * options.v_transpose
    *       
    *       P' = P + D (here P = E) 
    */

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


    /** ======================================
    *   Construct irreducible Markov matrix P'
    *   ======================================
    * 
    *   Formula:
    *
    *       E = e * options.v_transpose     with e = ones(NUM_OF_NODES)
    *
    *       P" = options.c * P' + (1 - options.c) * E       (here P' = P)
    */

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


    /**
     * We solve the following equation:
     * ([I]NxN - (options.c)*P.transpose) * x = ((1-options.c)/N )* [1]N
     * here : ([I]NxN - (options.c)*P.transpose) = K
     */
    double **K = (double **) malloc (NUM_OF_NODES * sizeof(double *));
    for (i=0; i < NUM_OF_NODES; i++) {
        K[i] = (double *) malloc(NUM_OF_NODES * sizeof(double));
        for (j=0; j < NUM_OF_NODES; j++){
            K[i][j] = 0;

            if (i==j) {
                K[i][j] = 1 - (options.c * P[j][i]);
            }
            else {
                K[i][j] = (-1)* (options.c * P[j][i]);
            }
        }
    }

    /************************************ 
    *   GAUSS - SEIDEL IMPLEMENTATION
    *************************************/

    // initializing x as the personalization vector
    double* x = (double*) malloc(NUM_OF_NODES * sizeof(double));
    for (i=0; i<NUM_OF_NODES; i++){
        x[i] = options.v[i];
    }
    
    double* x_old = (double*) malloc(NUM_OF_NODES * sizeof(double));
    double* diff = (double*) malloc(NUM_OF_NODES * sizeof(double));

    // initialize delta
    double delta = 1;
    int iter = 0;

    gettimeofday (&startwtime, NULL);

    while (delta > options.tolerance && iter < options.maxiter){

        for (i=0; i<NUM_OF_NODES; i++){
            x_old[i] = x[i];
        }

        for (i = 0; i < NUM_OF_NODES; i++){
            double sigma = 0;

            
            for (j=0; j < i-1; j++) {
                sigma += K[i][j]*x[j];
            }
        
            for (j=i+1; j < NUM_OF_NODES; j++) {
                sigma += K[i][j]*x_old[j];
            }

            x[i] = (double)((((1 - options.c) / NUM_OF_NODES) - sigma));
        }

        subtract(x, x_old, NUM_OF_NODES, diff);
        delta = norm(diff, NUM_OF_NODES);
        iter++;
    }


    gettimeofday (&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);
    printf("Serial Gauss-Seidel PageRank wall clock time = %f\n", seq_time);


    double sum_final = 0;
    for (i=0; i<NUM_OF_NODES; i++){
        sum_final = sum_final + x[i];
    }


    // Write final personalization vector x into .txt file
    f = fopen("./output/gaussOutput.txt", "w");
    printf("Writing output file..\n");
    for (i = 0; i < NUM_OF_NODES; i++) {
      fprintf(f, "%f\n", x[i]);
    }
    fclose(f);


    if (delta > options.tolerance || iter == options.maxiter){
        printf("Algorithm did not converged after %d iterations", iter);
    }
    else {
        printf("Converged after %d iterations. \n", iter);
    }

}