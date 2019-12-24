#pragma once

/* Class population represents a set of candidate solutions */
class population {
public:
// class attributes
	vector <individual> members;
	individual* fittest;

// parametrized constructor: takes in parameter of population size
	population(int popsize)
	{
		for (int i = 0; i < popsize; i++) {		// generate popsize new random individuals
			individual new_ind;
			new_ind.random_individual();
			members.push_back(new_ind);
		}
		fittest = find_fittest();				// find the fittest individual
	}

// member functions
	// find the fittest individual
	individual* find_fittest() {
		individual* fittest = &members[0];
		double max_fitness = fittest->fitness;
		for (int i = 0; i < popsize; i++) {
			if (members[i].fitness > max_fitness) {
				fittest = &members[i];
				max_fitness = fittest->fitness;
			}
		}
		return fittest;
	}

	// select with replacement through a 1 vs 1 tournament
	individual* tournament() {
		int contestant1 = (int)(popsize * unif(rndgen));	// select two different contestants randomly
		int contestant2 = (int)(popsize * unif(rndgen));
		while (contestant2 == contestant1) contestant2 = (int)(popsize * unif(rndgen));
		// fight
		return ((members[contestant1].fitness > members[contestant2].fitness) ? &(members[contestant1]) : &(members[contestant2]));
	}

	// create a new generation and subtitute the old one
	void new_generation() {
		vector<individual> Q;
		for (int i = 0; i < popsize / 2; i++) {
			individual* parent1 = tournament();
			individual* parent2 = tournament();
			while (parent2 == parent1) parent2 = tournament();
			vector <individual> offspring = one_point_crossover(parent1, parent2);
			Q.insert(Q.end(), offspring.begin(), offspring.end());
		}
		members = Q;
		fittest = find_fittest();
	}
};
