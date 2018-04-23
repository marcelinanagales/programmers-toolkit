/*
|================================|
|       Lab 3: Stochastic        |
|      Techniques in Global      |
|         Optimization           |
|--------------------------------|
|    Name: Marcelina Nagales     |
|      Due Date: 02/02/18        |
|================================|
*/

#include <iostream>                 // print to csv
#include <armadillo>                // matrices and vec (row/col vectors)
#include <vector>                   // smallCloud (random_shuffle)
#include <math.h>                   // sin, cos, sqrt, exp
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>                //std::max_element, random_shuffle
#include <fstream>
#include <string>                   // to_string()
    
using namespace std;
using namespace arma;


double f(double x, double y){
    
    return 1 + sin(x)*sin(x) + sin(y)*sin(y) - 0.1*exp(-x*x-y*y);
}

double gradf(double x, double y, int a){
    double result;
    if (a == 1)
        result = 0.2*exp(-x*x - y*y)*x + 2*cos(x)*sin(x);
    else
        result = 0.2*exp(-x*x - y*y)*y + 2*cos(y)*sin(y);
    
    return result;
}

double gradff(double x, double y, int a){
    double result;

    if (a == 1)
        result = exp(-x*x - y*y)*(-0.4*x*x + 2*exp(-x*x - y*y)*cos(2*x) + 0.2);
    else if (a == 2)
        result = -0.4*exp(-x*x - y*y)*x*y;
    else
        result = exp(-x*x - y*y)*(-0.4*y*y + 2*exp(-x*x - y*y)*cos(2*y) + 0.2);

    return result;
}

int main(){
    // all initial variables
    int N = 50;                     // N = number of points in sample space
    mat A(N, 2);                    // A = sample space
    int lower = -10;                // lower bound of range for x,y
    int upper = 10;                 // upper bound of range for x,y
    double range = upper-lower;     // range for x, y values
    vec fA(N);                      // column vector: f(x, y) for all x,y in A
    vector<int> smallCloud;         // for permutation of indices <randomization>
    int n_iterations = 700;         // number of iterations in for loop
    double max, min, average;       // f(x, y) calculations for all iterations
    int xm = -1;                    // xm retains the index of max f(x,y) in A
    int x_min = -1;                 // x_min retains the index of min f(x,y) in A
    int y1, y2, y3;                 // 3 random points in the sample space A   
    vector<int> iter_cloud;         // the iterations that need to plot the 
    iter_cloud.push_back(0);        //  sample space A (set it up) 
    iter_cloud.push_back(99);           
    iter_cloud.push_back(199);           
    iter_cloud.push_back(699); 
    int k = 0;                      // iteration counter
    string filename;                // for the iterator counter for iter_cloud

    // variables for the centroid and x_star values
    double midpointx, midpointy, dy, dx, x_starX, x_starY;
    // conditional values for xm 
    bool conditionX, conditionY, conditionZ;
    srand(unsigned (time(0)));

    // setting up the plot information to send to a csv
    ofstream lab3File;
    lab3File.open("lab3plots.csv");
    lab3File << "Iteration; Average; Max; Min \n"; 
    
    // Chooses N random points in range -10 < x, y < 10
    // calculate all f(x, y) for sample points in A
    for (int i = 0; i < N; i++){
        A(i, 0) = lower + (rand() / (RAND_MAX / range));
        A(i, 1) = lower + (rand() / (RAND_MAX / range));
        fA(i) = f(A(i, 0), A(i, 1));
        smallCloud.push_back(i);
    }

    // Start loop
    while (k < n_iterations){

        // this loop is to determine the max f(x, y) and its respective index
        max = -100;
        min = 100;
        average = 0;
        for (int j = 0; j < N; j++){

            // if the current f(x,y) is greater than the current max
            if (fA(j) > max){
                max = fA(j);
                xm = j;    
            }
            else if ( fA(j) < min ){
                min = fA(j);
                x_min = j;
            }
            average = average + fA(j);
        }

        // determine average 
        average = average/N; 

        // outputs to csv file
        lab3File << k << ';' << average << ';' << max << ';' << min << endl;

        // random permutation of indexes for 3 random points
        random_shuffle(smallCloud.begin(), smallCloud.end());

        // first three numbers in random permutation
        y1 = smallCloud[0];
        y2 = smallCloud[1];
        y3 = smallCloud[2];

        // calculating centroid (x, y) from y1 and y2
        midpointx = (A(y1, 0) + A(y2, 0))/2.0;
        midpointy = (A(y1, 1) + A(y2, 1))/2.0;

        // dy, dx are the distances in x and y from y3 to centroid
        dy = midpointy - A(y3, 1);
        dx = midpointx - A(y3, 0);

        // reflected y3 across the midpoint: x_star
        // calculates the dy and dx for less if statements
        x_starX = midpointx + dx;
        x_starY = midpointy + dy;

        // booleans to check that x_star is within the bounds
        // and so that the final if statement is simplified
        conditionX = x_starX >= -10 && x_starX <= 10;
        conditionY = x_starY >= -10 && x_starY <= 10;
        conditionZ = conditionX && conditionY;

        // simplified if statement (is x_star in bounds and is f(x*) smaller than max)
        if (conditionZ && f(x_starX, x_starY) < max){
            
            // if so, we change fA(xm), xm(x), and xm(y) 
            fA(xm) = f(x_starX, x_starY);
            A(xm, 0) = x_starX;
            A(xm, 1) = x_starY;
        }

        if (find(iter_cloud.begin(), iter_cloud.end(), k) != iter_cloud.end()){
            filename = "A_" + to_string(k+1) + ".txt";
            A.save(filename, raw_ascii);
        }
        k++;
    }// End loop
    cout << "Min of f(x, y) at x= " << A(x_min, 0) << " and y= " << A(x_min, 1) << endl;
    cout << "Min of f(" << A(x_min, 0) << ", " << A(x_min, 1)  << ")=" << min << endl;

    lab3File.close();

    // *** Newton's Method *** //
    double tol = 10e-10;                 // tolerance for while loop
    mat x(2, 2);                        // matrix of initial x values
    x << 0.1 << 0.1 << endr
      << 5.0 << 5.0 << endr;
    double gradfx, gradfy, gradfxx, gradfyy;
    double gradfxy, gradfyx;
    double denomin, p0, p1; 
    double f1, f2;
    double alpha;
    double alpha0 = 1;
    double rho = 0.5;
    double gamma = 0.5;
    int itercount = 0;
    double stepchange =1; 
    vec xold(2);
    mat xoriginal(2, 2);
    xoriginal << 0.1 << 0.1 << endr
      << 5.0 << 5.0 << endr;

    // have to do this for each set of x, y in vec x
    for (int m = 0; m < 2; m++ ){
        xold(0) = x(m, 0);
        xold(1) = x(m, 1);
        itercount = 0;

        stepchange = 1;
        // while statement: tolerance of step size and iteration count
        while (stepchange >= tol && itercount < 2000){
            // Newtons determine's p
            // calculate hessian
            gradfxx = gradff(x(m, 0), x(m, 1), 1); 
            gradfxy = gradff(x(m, 0), x(m, 1), 2);
            gradfyx = gradff(x(m, 0), x(m, 1), 2);
            gradfyy = gradff(x(m, 0), x(m, 1), 0);

            //calculate the gradient
            gradfx = gradf(x(m, 0), x(m, 1), 1);
            gradfy = gradf(x(m, 0), x(m, 1), 0);
            
            // denominator for inverse hessian 
            denomin = gradfxx*gradfyy - gradfxy*gradfyx;
            
            // calculate p with the inverse hessian and the gradient
            p0 = -(gradfx*gradfyy - gradfy*gradfxy)/denomin;
            p1 = -(-gradfxy*gradfx + gradfxx*gradfy)/denomin;

            // reset alpha
            alpha = alpha0;

            // while statement values
            f1 = f(x(m, 0)+alpha*p0, x(m, 1)+alpha*p1);
            f2 = f(x(m, 0), x(m, 0)) + gamma*alpha*(gradfx * p0 + gradfy * p1);

            // backtrack algorithm (makes alpha small)
            while (f1 > f2){
                alpha = alpha*rho;
                f1 = f(x(m, 0)+(alpha*p0), x(m, 1)+(alpha*p1));
                f2 = f(x(m, 0), x(m, 1)) + gamma*alpha*(gradfx * p0 + gradfy * p1);
            }// end of backtrack/alpha while loop

            // save the old x for stepchange 
            xold(0) = x(m, 0);
            xold(1) = x(m, 1);

            // set the new x
            x(m, 0) = x(m, 0)+ alpha*p0;
            x(m, 1) = x(m, 1)+ alpha*p1;
            
            // increment 
            itercount++;

            // step change conditional ~ tolerance
            stepchange = sqrt((x(m, 0)-xold(0))*(x(m, 0)-xold(0)) + (x(m, 1) - xold(1))*(x(m, 1) - xold(1)));

        }// end of while loop

        cout << "(Newton) Min for (" << xoriginal(m, 0) << " , " << xoriginal(m, 1) << ")= " << f(x(m, 0), x(m, 1)) 
        << " at iteration number " << itercount << endl;
    }// end of for loop

}