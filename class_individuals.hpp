#pragma once

/* Class individual represents a candidate solution */
class individual {
public:
// data attributes
	// phenotype: 4 variables
	unsigned long x;
	unsigned long y;
	unsigned long z;
	unsigned long t;
	// fitness score: f(x, y, z, t)
	double fitness;

// member functions
	// random initialization of a individual
	void random_individual()
	{
		x = rndgen();
		y = rndgen();
		z = rndgen();
		t = rndgen();
		find_fitness();
	}

	// find the fitness of an individual
	void find_fitness()
	{
		fitness = f((double)x, (double)y, (double)z, (double)t);
	}
	
	// mutate an individual with a small probability
	void mutate() {
		// rewrite this with for auto on a vector of pointers to the chromosomes
		for (int i = 0; i < 4; i++) {
			double chance = unif(rndgen);
			if (chance < p_mut) {									// if an individual is allowed to mutate in that chromosome
				int num_mut = (int)expgen(rndgen);  				// compute how many mutations should be introduced
				if (num_mut == 0) num_mut = 1;						// always at least one mutation
				unsigned long* chr = NULL;							// pointer to the chromosome to be mutated
				switch (i) {
                    case 0: chr = &x; break;
                    case 1: chr = &y; break;
                    case 2: chr = &z; break;
                    case 3: chr = &t; break;
				}
				int j = 0;
				unsigned char pos;									// position that will be mutated
				while (j < num_mut) {
					pos = (unsigned char)(N * unif(rndgen));		// random position to introduce mutation
					*chr = (*chr)^(1U << pos);
					j++;
				}
			}
		}
		find_fitness();
	}

	// print individual (for debugging)
	void print_individual() {
		cout << "genotype:" << endl;
		cout << "\tchX: " << bitset<N>(x) << endl;
		cout << "\tchY: " << bitset<N>(y) << endl;
		cout << "\tchZ: " << bitset<N>(z) << endl;
		cout << "\tchT: " << bitset<N>(t) << endl;
		cout << endl << "phenotype:" << endl;
		cout << "\tx: " << x << endl;
		cout << "\ty: " << y << endl;
		cout << "\tz: " << z << endl;
		cout << "\tt: " << t << endl;
		cout << endl << "fitness: " << fitness << endl;
	}
};

/* Non-member function */
// return a vector two of new individuals (offspring) by passing the pointers to the parent
vector<individual> one_point_crossover(individual* parent1, individual* parent2) {
	vector<individual> offspring;					// vector of offspring individuals
	individual child1;
	individual child2;  
    int d;
    for (int i = 0; i < 4; i++) {
        d = (int)(unif(rndgen) * N);				// random crossover point
        unsigned int mask = 0xFFFFFFFFU << d;
        switch (i) {
            case 0: child1.x = (parent1->x & mask) | (parent2->x & ~mask); child2.x = (parent2->x & mask) | (parent1->x & ~mask); break;
            case 1: child1.y = (parent1->y & mask) | (parent2->y & ~mask); child2.y = (parent2->y & mask) | (parent1->y & ~mask); break;
            case 2: child1.z = (parent1->z & mask) | (parent2->z & ~mask); child2.z = (parent2->z & mask) | (parent1->z & ~mask); break;
            case 3: child1.t = (parent1->t & mask) | (parent2->t & ~mask); child2.t = (parent2->t & mask) | (parent1->t & ~mask); break;
        }
    }
	child1.mutate(); child2.mutate();				// mutate the children
	offspring.push_back(child1);
	offspring.push_back(child2);
	return offspring;
}






