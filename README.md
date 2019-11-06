# PageRank
PageRank (PR) is an algorithm used by Google Search to rank web pages in their search engine results. According to Google:
PageRank works by counting the number and quality of links to a page to determine a rough estimate of how important the website is. The underlying assumption is that more important websites are likely to receive more links from other websites.

## Repository

This repository provides a serial implementation of the algorithm in C language based on Asynchronous Gauss-Seidel Iterations, as well as a parallel implementation using OpenMP. The project was undertaken as part of the "Parallel and distributed systems" course of AUTH university.

## Dependencies

For the serial algorithm only a compiler is needed (e.g. gcc).

To compile the parallel versions, please, install OpenMP to your system's requirements.

## Usage

Run the code with the commands:
```
// compile and run classic PageRank algorithm
$ make classic
$ ./classic

// compile and run serial Gauss-Seidel PageRank algorithm
$ make serial
$ ./serial

// compile and run parallel Gauss-Seidel PageRank algorithm
$ make parallel
$ ./parallel

// compile and run serial and parallel output compare script
$ make compare
$ ./compare
```
## References
For extra reading material on PageRank: 

- [Adaptive methods for the computation of PageRank](https://www.sciencedirect.com/science/article/pii/S0024379504000023)
- [A PageRank Algorithm based on Asynchronous Gauss-Seidel Iterations](https://ieeexplore.ieee.org/document/8431212)
