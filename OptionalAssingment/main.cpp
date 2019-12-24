/* All necessary includes */
#include <iostream>
#include <vector>
#include <random>
#include <bitset>

using namespace std;

/* constant expressions */
constexpr auto N = 32;
constexpr auto popsize = 1000;		// population size
constexpr auto p_mut = 0.01;        // mutation rate
constexpr auto update_limit = 200;	// max number of iterations without updating fittest individual before halting

/* Everything related to random number generation */
// uniformly-distributed integer random number generator that produces non-deterministic random numbers
random_device rndgen("/dev/urandom");
// random number generator of doubles according to the uniform distribution between 0 and 1
uniform_real_distribution<double> unif(0.f, 1.0f);

/* user-defined headers */
#include "ugly_function.hpp"
#include "class_individuals.hpp"
#include "class_population.hpp"

int main(int argc, char* argv[]) {
    
// current generation
	population candidates(popsize);
// fittest individual so far
	individual* solution = new individual;
    *solution = *(candidates.fittest);
	int the_iter = 0;
// maximum number of iterations that are allowed without updating the fittest individual
	int max_iter = update_limit;
        
    cout << "population of " << popsize << " members" << endl;
    cout << "mutation rate: " << p_mut << endl;
	cout << "update limit: " << update_limit << " iterations" << endl;
    
    for (int i = 1; i <= max_iter; i++) {
        candidates.new_generation();			// make a new generation
	// update fittest individual
		if (candidates.fittest->fitness > solution->fitness) {
			*solution = *(candidates.fittest);	// update the solution
			the_iter = i;						// the iteration that has found the fittest individual so far
			max_iter = i + update_limit;		// push the maximum number of iterations
		}
	}
	
	cout << "fittest individual found during the " << the_iter << "-th iteration" << endl << endl;
	solution->print_individual();
    
    delete solution;

	return 0;
    
}
