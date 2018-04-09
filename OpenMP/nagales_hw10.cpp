/*
 * write a program that does numerical integration
 * i.e. the integral of (sinx)dx [0, pi] = 2
*/

// Header files
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

double f(double x);
double trap(double upper_bound, double lower_bound, long n, int threadnum);
double simp(double upper_bound, double lower_bound, long n, int threadnum);

double trap(double upper_bound, double lower_bound, long n, int threadnum){

    double sum1 = 0;

    double x_n = lower_bound;
    int i;

    #pragma omp parallel for \
        shared( upper_bound, lower_bound, n) private(i, x_n) \
        reduction(+: sum1) num_threads(threadnum)
    

    for (i =0;i<=n; i++){
        if(i ==0 || i ==  n){
            sum1 = sum1 + f(x_n);
        }
        sum1 = sum1 + 2*f(x_n);
        x_n = lower_bound + (i+1)*(upper_bound - lower_bound)/n;
    }

    return sum1;
}

double simp(double upper_bound, double lower_bound, long n, int threadnum){
    double sum1 = 0;
    double x_n = lower_bound;
    int j; 
    #pragma omp parallel for \
        shared( upper_bound, lower_bound, n) private(j, x_n) \
        reduction(+: sum1) num_threads(threadnum)

    for(j = 0; j<=n; j++){
        if(j == 0 || j == n){
            sum1 = sum1 + f(x_n);
        }
        else if(j%2 == 0){
            sum1 = sum1+ 2.0*f(x_n);
        }
        else{
            sum1 = sum1+ 4.0*f(x_n);
        }
        x_n = lower_bound + (j+1)*(upper_bound - lower_bound)/n;
    }

    return sum1;
}

int main (){
// Variables:
    double upper_bound = M_PI; 
    double lower_bound = 0;
    double trapezoidal, simpson, s_error, t_error;
    float sum1, sum3 = 0;
    double sum2, sum4 = 0;
    double x_n = lower_bound;
    long n=100000000; 
    double exact_integral = 2;
    int threadlist[] = {1, 2, 4, 8};
    int k =0;
    float wtime, wtime2, wtime3;
/*
 * the process of calculating the integral given the original
 * function, the upper bound and the lower bound, using the 
 * trapezoidal function
 *
*/

    while(k < sizeof(threadlist)/sizeof(threadlist[0]))
    {

        wtime = MPI::Wtime();
        sum1 = trap(upper_bound, lower_bound, n, threadlist[k]);
        wtime = MPI::Wtime() - wtime;

        if(k == 0)
            sum3 = sum1;

        wtime2 = MPI::Wtime(); 
        sum2 = simp(upper_bound, lower_bound, n, threadlist[k]);

        wtime2 = MPI::Wtime() - wtime2;
        if(k == 0)
            sum4 = sum2;
        wtime3 = wtime + wtime2;
        printf("Work took %12f seconds for %d thread(s)\n", wtime3, threadlist[k]);
        k++;

    }
    /*################ DELETE ME: iVE BEEN MOVED ########################*/

    trapezoidal = ((upper_bound - lower_bound)/(2.0*n))*sum3;
    t_error = (exact_integral -trapezoidal)/exact_integral;

    simpson = ((upper_bound - lower_bound)/(3.0*n))*sum4;
    s_error = (exact_integral -simpson)/exact_integral;

    /*################ END DELETE ME ########################*/
//    }
    printf("The solution for the Trapezoidal Rule approximate integral of sin(x) for [%f, %f] is %e, with the error %e\n",lower_bound, upper_bound, trapezoidal, t_error);
    
    printf("The solution for the Simpson's Rule approximate integral of sin(x) for [%f, %f] is %e, with the error %e\n",lower_bound, upper_bound, simpson, s_error);
    
}

//first function is f, which in this case is sin(x)
double f(double x){
    return sin(x);
}