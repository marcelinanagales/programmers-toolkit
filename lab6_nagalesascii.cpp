# include <iostream>			// cout 
# include <armadillo>			// mat, vec
# include <math.h>				// exp
# include <random>				// normal distribution
# include <fstream>				// cout output

using namespace std;
using namespace arma;


double f(double x){
	// an input for weighted_jacobi_step
	return -4* exp(2*(x));
}

double f_exact(double x){
	// task 1
	return exp(2*x);
}

void weighted_jacobi_step(vec u0, int n, vec& u1, double h, vec x){
	// variable declarations
	// omega is the weight
	double omega = 2.0/3;
	// vector of f values 
	vec f_xy(n +1);
	//calculating f values
	for (int i = 0; i < n+1; i++){

		f_xy(i) = f(x(i));
	}

	// the weird for loop bounds are to avoid the end points 
	// and focus only on the inner nodes
	for (int k = 1; k < n; k++){
		// finite difference approximation
		u1(k) = omega*(u0(k-1)+ u0(k+1)+ pow(h, 2)*f_xy(k))/2 + 
					(1-omega)*u0(k);  

	} // end for k

	// return u1
} // end of weighted Jacobi

void function_restriction(vec u_h, vec& u_2h){

	// boundary nodes
	// all fine boundary nodes = all coarse boundary nodes
	u_2h(0) = u_h(0);
	u_2h(u_2h.n_rows-1) = u_h(u_h.n_rows - 1); 
	//for loop (interior nodes)
	for (int i = 1; i < u_2h.n_rows-1; i++){

		//interior nodes
		u_2h(i) = 0.25*(u_h(2*i-1) + 2*u_h(i) + u_h(2*i+1));
	}// end of for loop

	//return restricted function values
} // end of function restriction

void function_prolongation(vec& u_h, vec u_2h){
	//set the last value in uh to the last value of u2h 
	u_h(u_h.n_rows-1) = u_2h(u_2h.n_rows - 1);
	//for loop: all nodes (except for the last one)
	for(int j = 0; j < u_h.n_rows-1; j++){

		// if i is even (mod2 == 0)
		if (j % 2 == 0){
			u_h(j) = u_2h(j/2);
		}
		else{
			u_h(j) = (u_2h((j-1)/2) + u_2h((j+1)/2))/2.0;
		}
	}
	//return fine mesh function values 
}

void thomas_solver(vec a, vec b, vec c, vec d, vec& x){
	// number of values in a, b, c, d and x
	int n = b.n_rows;

	// initialize c_prime and d_prime
	vec c_prime(n);
	vec d_prime(n);

	// set first values of c_prime and d_prime according to wiki
	c_prime(0) = c(0)/b(0);
	d_prime(0) = d(0)/b(0);

	// for loop increment forwards to set the rest of c_prime and d_prime
	for(int i = 1; i < n; i++ ){
		c_prime(i) = c(i) / (b(i) - a(i)*c_prime(i-1));
		d_prime(i) = (d(i) - a(i)*d_prime(i-1)) / (b(i) - a(i)*c_prime(i-1));
	}
	// set the last x to last value in d_prime
	x(n - 1) = d_prime(n-1);	
	int j = n-2;
	// iterate backwards through to set x values 
	while (j >= 0){
		x(j) = d_prime(j) - c_prime(j)*x(j+1);
		j--;
	}

	// return x
}

int multigrid(double h5, int num_cycles, vec& u0h_d6, bool flag7, double tol){
	// Multigrid method
 	int n5 = 1/h5;									// number of intervals 
 	int iter = 0;									// initial iteration
 	double sum;										// summation for norm
 	// weighted jacobi variables
	vec u1h_temp(n5 + 1);							// first jacobi result
	vec u1h_d6(n5+1);								// second jacobi result
	vec x_d5(n5 + 1);								// x vector for f(x) calc

	// set up u1_temp and u1_h with boundary conditions
	u1h_d6(0) = 1;									
	u1h_d6(n5) = exp(2);
	
	// temp vectors
	u1h_temp(0) = 1;
	u1h_temp(n5) = exp(2);

	// set up x values
	for (int m = 0; m < n5; m ++){
		x_d5(m) = m*h5; 
	}
	// calculate rh variables
	vec fh_d6(n5+1);
	vec rh_d6(n5+1);
	double sum1;

	for (int t = 0; t < n5+1; t++){
		fh_d6(t) = f(x_d5(t));
	}

	// restrict residual variables
	vec r2h_d6(n5/2 + 1);
	r2h_d6.fill(0.0);

	vec r2h_interior(n5/2 -1); 

	// calculate e2h
	vec a_d6(n5/2 - 1);
	vec b_d6(n5/2 - 1);
	vec c_d6(n5/2 - 1);
 	vec e2h_temp(n5/2 - 1);
 	vec e_2h(n5/2 + 1);
 	int counter2;
 	double frob_norm;
 	double frob_norm2;

 	a_d6.fill(pow((n5/2), 2));
	a_d6(0) = 0;
	
	b_d6.fill(-2*pow((n5/2), 2));

	c_d6.fill(pow((n5/2), 2));
	c_d6(n5/2-2) = 0;

 	e2h_temp.fill(0);
 	e_2h.fill(0);

 	// prolong e_2h -> e_h variables
 	vec eh_d6(n5 + 1);
	eh_d6.fill(0.0);

// while loop
	while(iter < num_cycles){
		// u1_h: apply 2 steps of the weighted jacobi 
		weighted_jacobi_step(u0h_d6, n5, u1h_temp, h5, x_d5);
		weighted_jacobi_step(u1h_temp, n5, u1h_d6, h5, x_d5);

		// Compute residual on the fine mesh
		rh_d6(0) = 0;
		rh_d6(rh_d6.n_rows - 1) = 0;
		for (int p= 1; p < n5; p++){
			rh_d6(p) = - fh_d6(p) - (u1h_d6(p-1) -2*u1h_d6(p)+
				u1h_d6(p+1))/(h5*h5);
		}

		// Restrict residual to a coarser mesh
		function_restriction( rh_d6, r2h_d6);;
		
		for (int j = 1; j < n5/2; j++ ){
			r2h_interior(j-1) = r2h_d6(j);
		}

		// Solve for e_2h on the interior nodes using the thomas algorithm
		thomas_solver(a_d6, b_d6, c_d6, r2h_interior, e2h_temp);

		counter2 = 0;
	 	for (int k = 1; k < n5/2; k++ ){
			e_2h(k) = e2h_temp(counter2);
			counter2++;
		}

		// Prolong e_2h to e_h
		function_prolongation(eh_d6, e_2h);

		if(flag7){
			// calculate l2 norm of residual rh (on the fine mesh)
			sum = 0;

			for (int k = 1; k < n5; k++){		
			// in calculation of residual
				sum = sum + pow(abs(eh_d6(k)), 2); 

			} // end for k
			frob_norm = sqrt(sum);
			cout << "Multigrid: " << frob_norm << endl;
			// if less than provided tolerance
			if (frob_norm < tol)
			// return number of iterations
				return iter;
			else
				iter++;
		}
		else{
			// correct u1
			// set u0 = u1
			for (int m = 0; m < n5+1; m++ ){
				u0h_d6(m) = u1h_d6(m) + eh_d6(m);
			}
			iter++;
		}
	}

	return iter;
}

int main(){


	//set up variables
	// h for the calculation of u0_h
	double h1 = 1.0/128;

	// n1 is the number of points in x and u
	int n1 = 1/h1 ;
	
	// for random number *NOISE* generator of u0_d1
	default_random_engine generator;
	normal_distribution<double> distribution(0.0, 2.0);

	// generate u0 for Deliverable 1
	vec u0_d1(n1+1);
	u0_d1.fill(0.0);
	for(int i = 0; i < n1+1; i++){
		u0_d1(i)= u0_d1(i) + distribution(generator); 
	}

	// generate empty u1 for deliverable
	vec u1_d1(n1+1);
	u1_d1.fill(0.0);
	u1_d1(0) = u0_d1(0);
	u1_d1(n1-1) = u0_d1(n1-1);

	// generate empty u2 for deliverable
	vec u2_d1(n1+1);
	u2_d1.fill(0.0);
	u2_d1(0) = u0_d1(0);
	u2_d1(n1-1) = u0_d1(n1-1);
	
	// generate x for 
	vec x_d1(n1+1);
	for (int j =0; j < n1 + 1; j++){
		x_d1(j) = j*h1;
	}

	// Deliverable #1: weighted Jacobi (2x)
	// weighted_jacobi_step(vec u0, int n, vec& u1, double h, vec x)

	weighted_jacobi_step(u0_d1, n1, u1_d1, h1, x_d1);
	weighted_jacobi_step(u1_d1, n1, u2_d1, h1, x_d1);

	u0_d1.save("u0_d1.txt", raw_ascii);
	u2_d1.save("u1_d1.txt", raw_ascii);

	// Deliverable #2: function restriction
	// using the u0 given for Deliverable #1

	// generate empty u2 for deliverable
	vec u2h_d2(n1 /2 + 1);
	u2h_d2.fill(0.0);
	function_restriction( u0_d1, u2h_d2);

	u2h_d2.save("u2h_d2.txt", raw_ascii);

	// Deliverable #3: function prolongation
	// uh -> u2h

	int n3 = 32;
	// generate empty u1 for deliverable
	vec u2h_d3(n3 + 1);
	for (int k = 0; k < n3+1 ; k++){
		u2h_d3(k) = distribution(generator);
	}

	// generate u2h for deliverable
	vec uh_d3(2*n3+1);
	uh_d3.fill(0.0);

	function_prolongation( uh_d3, u2h_d3);


	u2h_d3.save("u2h_d3.txt", raw_ascii);
	uh_d3.save("uh_d3.txt", raw_ascii);

	// Deliverable #4: Thomas Algorithm
	int n4 = 10;
	vec a(n4);
	a.fill(-1);
	a(0) = 0;
	
	vec b(n4);
	b.fill(3);

	vec c(n4);
	c.fill(-1);
	c(n4-1) = 0;

	vec d(n4);
	d.fill(1);
	d(0) = 2;
	d(n4-1) = 2;
 	
 	vec x(n4);
 	x.fill(0);

 	thomas_solver(a, b, c, d, x);

 	x.print("x: ");

 	// Task 5: weighted Jacobi method

 	// set up 
 	double loop_tol = 100;
	double frob_norm;
	int counter = 0;
	double sum = 0;
	double tol = 0.00001;

	double norm; 
	double sum2;

 	double h5 = 1.0/256;
 	int n5 = 1/h5;
 	
 	vec u0_d5(n5 + 1);
 	u0_d5.fill(0);
 	u0_d5(0) = 1;
 	u0_d5(n5) = exp(2);

	vec u1_d5(n5 +1);
 	u1_d5.fill(0);
 	u1_d5(0) = 1;
 	u1_d5(n5) = exp(2);

	vec x_d5(n5 + 1);
	for (int m = 0; m < n5; m ++){
		x_d5(m) = m*h5; 
	}

 	//ofstream stuff
	ofstream lab6File;

	lab6File.open("jacobi_1.csv");
	lab6File << "Iteration; L2 norm \n";
	while ( counter < 2000){
		sum = 0;
		sum2 = 0;

 		weighted_jacobi_step(u0_d5, n5, u1_d5, h5, x_d5);

		for (int k = 1; k < n5; k++){		
			// in calculation of residual
			sum = sum + pow(abs(u1_d5(k) - u0_d5(k)), 2); 
			sum2 = sum2 + pow(abs(u1_d5(k) - f_exact(x_d5(k))), 2); 

		} // end for k

		// Residual
		frob_norm = sqrt(sum);
		norm = sqrt(sum2);
		lab6File << counter << ';' << norm << endl;
		loop_tol = frob_norm;

		if ( frob_norm > tol){
			
			for (int p = 1; p < n5 ; p++){
					// finite difference approximation
					u0_d5(p) = u1_d5(p);  
			} // end for p
		}

		counter++;
	}
	lab6File.close();

	u1_d5.save("u1h_jacobi1.txt", raw_ascii);

 	// Task 6: multigrid
 
 	// Min of 1000 multigrid cycles
 	// Task 6: multigrid

	h5 = 1/256.0;
	n5 = 1/h5;
 	int num_cycles = 1000;
 	int iterations;

	vec u0h(n5+1);								// initial u0
	u0h.fill(0.0);
	// u0_h: Initial guess on the fine mesh
	// set up boundary conditions
	u0h(0) = 1;
	u0h(n5) = exp(2);
	
	iterations = multigrid(h5, num_cycles, u0h, false, 0);
	u0h.save("u1h_multigrid1.txt", raw_ascii);

	// Task 7: Computational Cost

	// loop through different stopping tolerances
	vec tol_tcc(5);
	for (int a1 = 0; a1 < tol_tcc.n_rows; a1++){
		tol_tcc(a1) = pow(10, -(a1+1));
	}
	// Jacobi Method
	// set up 
 	loop_tol = 100;
	counter = 0;
	sum2 = 0;

 	double h7 = 1.0/64;
 	int n7 = 1/h7;
 	
 	vec u0_d7(n7 + 1);
 	u0_d7.fill(0);
 	u0_d7(0) = 1;
 	u0_d7(n7) = exp(2);

	vec u1_d7(n7 +1);
 	u1_d7.fill(0);
 	u1_d7(0) = 1;
 	u1_d7(n7) = exp(2);

	vec x_d7(n7 + 1);
	for (int m = 0; m < n7; m ++){
		x_d7(m) = m*h7; 
	}

	cout << "Task 7: Computational Costs" << endl << endl;
	// multigrid 
	// set up
	num_cycles = 10000;

	for (int b1 = 0; b1 < tol_tcc.n_rows; b1++){
	 	u0_d7.fill(0);
	 	u0_d7(0) = 1;
	 	u0_d7(n7) = exp(2);

	 	u1_d7.fill(0);
	 	u1_d7(0) = 1;
	 	u1_d7(n7) = exp(2);
	 	counter = 0;
	 	loop_tol = 100;
		while ( counter < 10000 && loop_tol > tol_tcc(b1)){
			sum = 0;

	 		weighted_jacobi_step(u0_d7, n7, u1_d7, h7, x_d7);

			for (int k = 1; k < n7; k++){		
				// in calculation of residual
				sum = sum + pow(abs(u1_d7(k) - u0_d7(x_d7(k))), 2); 
			} // end for k

			// Residual
			frob_norm = sqrt(sum);
			loop_tol = frob_norm;
		//	cout << "Jacobi: " << frob_norm << endl;

			if ( frob_norm > tol_tcc(b1)){
				
				for (int p = 1; p < n7 ; p++){
						// finite difference approximation
						u0_d7(p) = u1_d7(p);  
				} // end for p
			}
			else
				break;

			counter++;
		}
	/**/cout << tol_tcc(b1) << " Jacobi: " << counter << endl;

	}//end for loop (b1)
}