

Marcelina Nagales
Homework 10

> g++ -fopenmp nagales_hw10.cpp
> ./a.out

// on steve 
> mpic++ nagales_hw10.cpp
> ./a.out

(also had to change omp_get_wtime to MPI::Wtime)

to install on linux, I used: 

> sudo apt-get install libopenmpi-dev
> sudo apt-get install openmpi-bin


(I am still getting the issue that the time is increasing as more cores are
used, and I believe this to be incorrect)
