# Environment
The experiments are conducted on a Linux server with an AMD EPYC 7443 24-Core Processor and 1024G memory. All algorithms are implemented in C++. The code is compiled with g++ 8.5 under O3 optimization. All algorithms are evaluated using a single thread.

# Input graph 
Using 'pre.cpp' to preprocess input graph as a simple connected graph with no self-loop and no multi-edge, also rename nodes from 0 to n-1.

# Run

1. Preprocess the input dataset.
```shell
g++ -o pre pre_process.cpp  -std=c++11 -O3
./pre [dataset]
```

2. Compile 'QCB-[method].cpp' to generate cycle basis using different approaches. 'main.cpp' is the same as 'QCB-ST.cpp' without output. Using check.cpp to check the correctness of cycle basis when having output.

```shell
g++ -o main main.cpp check.cpp -std=c++11 -O3
./main [mc_dataset]
```

3. Compile 'st_query.cpp' to query arbitrary cycles containing a certain edge using different cycle bases.

```shell
g++ -o st st_query.cpp  -std=c++11 -O3
./st [dataset] [cycle_basis_dirname]
```

# Example
 ```shell
 g++ -o QS QCB-ST.cpp check.cpp -std=c++11 -O3
 ./QS mc_test 1
 g++ -o st st_query.cpp  -std=c++11 -O3
 ./st test ST
 ```

