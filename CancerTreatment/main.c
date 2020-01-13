#include "RKF78.c"
#include "auxiliary.c"
#include "GA_functions.c"


int main(int argc, char* argv[]) {
  
  int cure = GeneticAlgorithm('c');
  
  if (!cure) GeneticAlgorithm('p');
  
  return 0;
}
