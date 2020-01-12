#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "RKF78.c"
#include "auxiliary.c"
#include "GA_functions.c"


int main(int argc, char* argv[]) {
  
  int randomData = open("/dev/urandom", O_RDONLY);
  if (randomData < 0) return 1;
  
  individual * population = NULL;
  population = (individual*) malloc(popsize*sizeof(individual));
  individual fittest;
  
  /* Create a random population of candidate solution: start by looking at possible curative treatment. */
  ssize_t result = 0;
  for (int i = 0; i < popsize; i++) {
    result = read(randomData, population[i].Cij, sizeof(population[i].Cij));
    if (result < 0) return 1; // call ExitError function
    population[i].fitness = Curative_Fitness(population[i].Cij);
  }
  
  for (int i = 0; i < (n_par-1)*d_par; i++) printf("%u\n", population[1].Cij[i]>>4);
  for (int i = 0; i < popsize; i++) printf("%f\n", population[i].fitness);
  
  close(randomData);
  free(population);
    
  return 0;
}
