# include <iostream>			// cout 
# include <stdlib.h>			// rand
# include <fstream>				// cout output
# include <math.h>
# include <time.h>				// time
# include <mpi.h>				// OpenMPI Example 4

using namespace std;

// function for exercise 1 and 2
double monteHall(int num, bool flag){

//	srand(time(NULL));				// for the random chosen and picked

	double winPercent = 0;
	int chosen;  
	int choice;
	int empty;
	int new_choice = 0;

	for(int i = 0; i < num; i++){
		chosen = rand() % 3;  
		choice = rand() % 3;
		empty = chosen;

		// reveal loss door //
		while (empty == chosen || empty == choice){
			empty = rand() % 3;
		}

		// switch door //
		if (flag){
			new_choice = 0;
			while (new_choice == empty || new_choice == choice){
				new_choice++;
			}
			choice = new_choice;		
		}

		if (choice == chosen){
			winPercent++;
		}

	}

	return winPercent/num;

}

// function for exercise 3
double adjust_monteHall(int n, int m, int num){
	double winPercent = 0.0;
	int car;  
	int choice;

	int emptyIndex;
	bool reveal[m] = {false};		// the indexes are revealed doors 
									// 1s are revealed, others are not (0)
	int b = 0;

	int newChoiceOptions[m-n-1] = {0}; 
	int newChoiceIter = 0;
	int new_choice = 0;
	int door, reveals;
	double isame = 0;
	for(int i = 0; i < num; i++){
		b = 0;
		while (b < m){
			reveal[b] = false;
			b++;
		}
		newChoiceIter = 0;
		// list of new choices
		while (newChoiceIter < m-n-1){
			newChoiceOptions[newChoiceIter] = 0;
			newChoiceIter++; 
		}


		// select car and choice by user
		car = rand() % m;  
		choice = rand() % m;
		if (choice == car){
			isame ++;
		}

		// reveal loss door //
		reveals = 0; // should count all reveals
		door = 0; // should go thru all doors
		while(reveals < n){
			// picks a door
			emptyIndex = door;

			//if door is not choice or car, reveal
			if (emptyIndex != choice && emptyIndex != car){

				reveal[emptyIndex] = true;
				//cout << "revealed " << emptyIndex << endl;
				reveals++;
			}
			door++; // move on to next door
		}
		door = 0;
		newChoiceIter = 0;
		// list of new choices
		while (newChoiceIter < m-n-1){
			// if iter is not choice or revealed
			if (door != choice && !reveal[door]){
				// add to new choice options
				newChoiceOptions[newChoiceIter] = door;
				newChoiceIter++; 
			}
			// move on to next door
			door++; 
		}
		// switch door //

		new_choice = newChoiceOptions[rand() % (m-n-1)];

		if (new_choice == car){
			winPercent++;
		}

	}


	return winPercent/num;
}

// function for exercise 4
double openMPI_monteHall(int num, int cores){

//	srand(time(NULL));				// for the random chosen and picked

	double winPercent = 0;
	int chosen;  
	int choice;
	int empty;
	int new_choice = 0;

    #pragma omp parallel for \
    	shared(num) private( i, chosen, choice , empty) \
    	reduction(+: winPercent) num_threads(cores)

	for(int i = 0; i < num; i++){
		chosen = rand() % 3;  
		choice = rand() % 3;
		empty = chosen;

		// reveal loss door //
		while (empty == chosen || empty == choice){
			empty = rand() % 3;
		}

		if (choice == chosen){
			winPercent++;
		}

	}

	return winPercent/num;

}
int main(){

	int N = 10000;					// number of trials
 	ofstream outfile;
   	outfile.open("afile.txt");
   	float wtime;
   	float wtime2;
   	double winPercent;	
   	srand(time(NULL));				// for the random chosen and picked
	
	// Exercise 1
	/*	
		make 3 doors
		select chosen door
		pick random door
		reveal 1 lose door
		**keep selected door**
		if win, add to percentage
		repeat for N trials
	*/ 
	cout << "Exercise 1 Win Percentage: " << monteHall(N, false) << endl;
	// Exercise 2
	/*	
		make 3 doors
		select chosen door
		pick random door
		reveal 1 lose door
		**change selected door**
		if win, add to percentage
		repeat for N trials
	*/ 
	cout << "Exercise 2 Win Percentage: " << monteHall(N, true) << endl;
	// Exercise 3
	/*	
		make m doors
		select chosen door
		pick random door
		reveal m lose doors
		**change selected door**
		if win, add to percentage
		repeat for m < 100, n < m-2 
		output m, n, winPercentage
	*/ 
	for(int m = 3; m < 100; m++){
		for (int n = 1; n <= m-2; n++){
		outfile << m << " " << n << " " << adjust_monteHall(n, m, N) 
			<< endl;
		}
	}

	outfile.close();

	// Exercise 4
	/*	
		make 3 doors
		select chosen door
		pick random door
		reveal 1 lose door
		**keep selected door**
		if win, add to percentage
		repeat for N trials
		> open mpi
	*/
	// Initialize the MPI environment
    MPI_Init(NULL, NULL);

	for (int i = 1; i < 5; i++ ){
		if (i != 3){
			wtime = MPI::Wtime ( );
		    winPercent = openMPI_monteHall ( N , i);
		    wtime = wtime - MPI::Wtime ( ) ;
		    cout << "Time " << wtime << " seconds" << " for " 
		    << i << " threads."<< endl;
		}
	}

	// Finalize the MPI environment.
    MPI_Finalize();

}