#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

#define popsize 40
#define max_update 5

typedef struct {
  unsigned char Cij[(n_par-1)*d_par];
  double fitness;
} individual;


double uniform() {
  int randomData = open("/dev/urandom", O_RDONLY);
  if (randomData < 0) ExitError("Could not open /dev/urandom", 1);
  unsigned char num;
  ssize_t result = read(randomData, &num, 1);
  if (!result) ExitError("could not read from urandom", 1);
  double randNum = (double)num / 255.0;
  close(randomData);
  return randNum;
}


individual* find_fittest(individual* population) {
  individual* fittest = population;
  for (int i = 1; i < popsize; i++) if (population[i].fitness < fittest->fitness) fittest = &population[i];
  return fittest;
}


void Mutation(individual* candidate) {
  int overallPos = (4*(n_par-1)*d_par - 1) * uniform();
  int i = overallPos/4;
  int j = overallPos%4;
  candidate->Cij[i] = (candidate->Cij[i])^(1U << (7-j));
}


void printIndividual(individual* ind) {
  for (int i = 0; i < n_par-1; i++) {
    for (int j = 0; j < d_par; j++) printf("%2d ", ind->Cij[i*d_par + j]>>4);
    printf("\n");
  }
  printf("Fitness: %f\n\n", ind->fitness);
}


void OnePointCrossover(individual* parent1, individual* parent2, individual* child1, individual* child2) {
  unsigned short bitsPerInd = 4 * (n_par-1) * d_par;
  unsigned short overallPos = uniform() * (bitsPerInd - 1);
  unsigned short c = overallPos/4;
  unsigned short mask = 0xFFFFFFFFU << (8 - overallPos%4);
  for (int i = 0; i < c; i++) { child1->Cij[i] = parent1->Cij[i]; child2->Cij[i] = parent2->Cij[i]; }
  for (int i = c+1; i < (n_par-1)*d_par; i++) { child1->Cij[i] = parent2->Cij[i]; child2->Cij[i] = parent1->Cij[i]; }
  unsigned char p1 = parent1->Cij[c], p2 = parent2->Cij[c];
  child1->Cij[c] = (p1 & mask) | (p2 & ~mask);
  child2->Cij[c] = (p2 & mask) | (p1 & ~mask);
}


individual* TournamentSelection(individual* population) {
  int cont1 = uniform() * (popsize - 1), cont2 = uniform() * (popsize - 1);
  while (cont2 == cont1) cont2 = uniform() * (popsize - 1);
  if (population[cont1].fitness < population[cont2].fitness) return (population + cont1);
  return (population + cont2);
}


void writeResult(unsigned char* Cij) {
  register unsigned char i, j;
  ODE_Parameters GompertzParams;
  double N = NZero_par, t = t_i[0];
  double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
  unsigned char npar = n_par - 1;
  
  FILE* out;
  if ((out = fopen("solution.csv", "w")) == NULL) ExitError("could not open results.csv", 1);
  fprintf(out, "t;N\n");
  fprintf(out, "%.5f;%.5f\n", t, N);
  
  for(i=0; i < npar; i++){
    double tfin = t_i[i+1]; // Implementing treatment i
    GompertzParams.drift_i = 0.0;
    for(j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * (*(Cij++)>>4);
    
    while(t+h < tfin) {
      RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
      fprintf(out, "%.9f;%.0f\n", t, N);
      if(N > NMax_par) { fclose(out); return; }
    }
    
    do {
      h = tfin - t;
      RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
      fprintf(out, "%.5f;%.5f\n", t, N);
      if(N > NMax_par) { fclose(out); return; }
    } while (t < tfin);
    if(N < 1000) { fclose(out); return; }
  }
  fclose(out);
}


/* GENETIC ALGORITHM: FIND A CURATIVE OR PALIATIVE TREATMENT */
int GeneticAlgorithm(unsigned char objective) {
  
  int randomData = open("/dev/urandom", O_RDONLY);
  if (randomData < 0) ExitError("Could not open /dev/urandom", 1);
  
  individual * population = NULL;
  if ((population = (individual*) malloc(popsize*sizeof(individual))) == NULL) ExitError("could not allocate memory for the population", 1);
  
  /* Create a random population of candidate solutions: start by looking at possible curative treatment. */
  for (int i = 0; i < popsize; i++) {
    ssize_t result = read(randomData, population[i].Cij, (n_par-1)*d_par);
    if (result != (n_par-1)*d_par) ExitError("could not generate random initial population", 1);
    if (objective == 'c') population[i].fitness = Curative_Fitness(population[i].Cij);
    else if (objective == 'p') population[i].fitness = Paliative_Fitness(population[i].Cij);
    else ExitError("invalid objective", 1);
  }
  close(randomData);
  
  /* Find the fittest individual and store it in a copy. */
  individual* fittestInCurPop = find_fittest(population);
  individual fittestOverall;
  strcpy((char*) fittestOverall.Cij, (const char*) fittestInCurPop->Cij);
  fittestOverall.fitness = fittestInCurPop->fitness;
  
  /* MAIN LOOP */
  if (objective == 'c') printf("... Applying GA to find a curative treatment ...\n");
  if (objective == 'p') printf("... Applying GA to find a paliative treatment ...\n");
  
  int max_iter = max_update;
  int iter = 0;
  
  while (iter < max_iter) {
    individual * nextGen = NULL;
    if ((nextGen = (individual*) malloc(popsize*sizeof(individual))) == NULL) ExitError("could not allocate memory for the next generation", 1);
    /* make a new generation */
    for (int j = 0; j < popsize/2; j++) {
      individual * parent1 = TournamentSelection(population);
      individual * parent2 = TournamentSelection(population);
      OnePointCrossover(parent1, parent2, nextGen + 2*j, nextGen + 2*j + 1);
      Mutation(nextGen + 2*j); Mutation(nextGen + 2*j + 1);
      if (objective == 'c') { nextGen[2*j].fitness = Curative_Fitness(nextGen[2*j].Cij); nextGen[2*j+1].fitness = Curative_Fitness(nextGen[2*j+1].Cij);}
      else if (objective == 'p') { nextGen[2*j].fitness = Paliative_Fitness(nextGen[2*j].Cij); nextGen[2*j+1].fitness = Paliative_Fitness(nextGen[2*j+1].Cij);}
      else ExitError("invalid objective", 1);
    }
    /* replace previous population with new generation */
    individual * aux = population;
    population = nextGen;
    free(aux);
    /* find fittest inidividual in new population */
    fittestInCurPop = find_fittest(population);
    if (fittestInCurPop->fitness < fittestOverall.fitness) {
      strcpy((char*) fittestOverall.Cij, (const char*) fittestInCurPop->Cij);
      fittestOverall.fitness = fittestInCurPop->fitness;
      max_iter = iter + max_update + 1;
    }
    iter++ ;
  }
  free(population);
  
  if (objective == 'c') {
    if (fittestOverall.fitness < DBL_MAX) {
      printf("\nFound a solution for a curative treatment.\n\n");
      printf("Recommended dosages:\n");
      for (int i = 0; i < n_par-1; i++) { printf("Session %d: ", i+1); for (int j = 0; j < d_par; j++) printf("%2d ", fittestOverall.Cij[i*d_par + j]>>4); printf("\n"); }
      printf("\nIntegral N(t)dt = %.0f\n\n", fittestOverall.fitness);
      writeResult(fittestOverall.Cij);
      return 1;
    }
    else { printf("\nA feasible curative treatment could not be found.\n\n"); return 0; }
  }
  else if (objective == 'p') {
    if (fittestOverall.fitness < DBL_MAX) {
      printf("\nFound a solution for a paliative treatment.\n\n");
      printf("Recommended dosages:\n");
      for (int i = 0; i < n_par-1; i++) { printf("Session %d: ", i+1); for (int j = 0; j < d_par; j++) printf("%2d ", fittestOverall.Cij[i*d_par + j]>>4); printf("\n"); }
      printf("\nExpected lifespan: %.1f weeks.\n\n", 1/(fittestOverall.fitness));
      writeResult(fittestOverall.Cij);
      return 1;
    }
    else { printf("\nA feasible curative treatment could not be found.\n\n"); return 0; }
  }
  else ExitError("invalid objective", 1);
    
  return 0;
}
