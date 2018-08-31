#include <math.h>
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

    /** =======================================================
    *   Open dataset file to read NUM_OF LINES and NUM_OF_NODES  
    *   =======================================================
    **/
    FILE *f;
    f = fopen(DAT_FILE,"r");
    if (f == 0) {
        perror("[ERROR]: Failed to open input file\n");
        exit(EXIT_FAILURE);
    }

    fscanf(f,"%d",&NUM_OF_NODES);
    fscanf(f, "%d", &NUM_OF_LINES);
    printf("[INFO]: Number of nodes: %d \n", NUM_OF_NODES);
    printf("[INFO]: Number of lines: %d \n", NUM_OF_LINES);


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
    printf("[INFO]: Reading data wall clock time = %f\n", seq_time);




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
            // normalize E matrix
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
            D[i][j] = d[i] * options.v[j];
            P[i][j] = E[i][j] + D[i][j];
        }
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


    /** =================================================================
     *  The rest of the computations are done using the transpose
     *  of array P.
     *  For simplicity and memory management we are using the same matrix
     *  and changing the sequence of indices
     *  =================================================================
     */


    /** =======================================================================
     *  We solve the following equation:
     *
     *      ([I]NxN - options.c * P.transpose) * x = ((1-options.c) / N) * [1]N
     *      
     *   => K = ([I]NxN - (options.c) * P.transpose) = K
     *  =======================================================================
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

    /** =============================
    *   GAUSS - SEIDEL IMPLEMENTATION
    *   =============================
    */

    // initialize x as the personalization vector
    double* x = (double*) malloc(NUM_OF_NODES * sizeof(double));
    for (i=0; i<NUM_OF_NODES; i++){
        x[i] = options.v[i];
    }
    
    double* x_old = (double*) malloc(NUM_OF_NODES * sizeof(double));

    // initialize delta
    double delta = 1;
    int iter = 0;

    gettimeofday (&startwtime, NULL);

    while (delta > options.tolerance && iter < options.maxiter){

        for (i=0; i<NUM_OF_NODES; i++){
            x_old[i] = x[i];
        }
        
        double sigma = 0;

        for (i = 0; i < NUM_OF_NODES; i++){
            sigma = 0;

            for (j=0; j < i-1; j++) {
                sigma += K[i][j]*x[j];
            }

            for (j=i+1; j < NUM_OF_NODES; j++) {
                sigma += K[i][j]*x_old[j];
            }

            x[i] = (double)((((1 - options.c) / NUM_OF_NODES) - sigma));            
        }   

        // calculate norm of diff = x-x_old vector
        double sum = 0, diff = 0;
        for (int i = 0; i < NUM_OF_NODES; ++i) {
            diff = x[i] - x_old[i];
            sum = sum + diff * diff;
        }

        delta = sqrt(sum);
        iter++;
    }

    gettimeofday (&endwtime, NULL);
    seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);
    printf("[INFO]: Serial PageRank Gauss-Seidel wall clock time = %f\n", seq_time);

    /* ===================================================
    *  Write final personalization vector x into .txt file
    *  ===================================================
    */
    f = fopen("./output/gauss.txt", "w");

    printf("[INFO]: Writing output file..\n");
    for (i = 0; i < NUM_OF_NODES; i++) {
      fprintf(f, "%f\n", x[i]);
    }
    
    fclose(f);


    if (delta > options.tolerance || iter == options.maxiter){
        printf("[INFO]: Algorithm did not converged after %d iterations", iter);
    }
    else {
        printf("[INFO]: Converged after %d iterations. \n", iter);
    }

}

