# Matrix Multiplication Using MPI

Highly researched topic in computer science and linear algebra.  
Implement a parallel matrix algorithm using OpenMP.  
  
The code is written in C and uses OpenMP and optimizes the multiplication kernel using Strassen's algorithm.  

I've defined 2 5000x5000 matrices to test the algorithm. These are further subdivided into 2x2 matrices and Strassen's algorithm reduces the number of multiplications for a 2x2 matrix to 7 multiplications.

I'll be providing a comparison on the time taken for the algorithm to run using different number of tasks, nodes and cpus per task.

|     Nodes     |     Tasks     |  Cores/Task   |     Time(s)   | 
| ------------- | ------------- | ------------- | ------------- |
|      5        |       50      |       2       |      38.923   |
|      5        |       30      |       2       |      51.105   |
|      5        |       15      |       2       |     102.286   |
|      5        |       10      |       2       |     220.365   |
|      5        |       5       |       2       |     356.715   |
