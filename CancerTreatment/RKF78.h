#include <stdio.h>
#include <stdlib.h>

double RKF78(double *, double *,
             double *, double, double,
             double,
             void *,
             void (*)(double, double, double *, void *));

void RKF782tfin(double *, double *, double, void *,
             void (*)(double, double, double *, void *));

void InitializeRKF78Sys(unsigned char);
double RKF78Sys(double *, double *,
             double *, double, double,
             double,
             void *,
             void (*)(double, double *, unsigned char, double *, void *));
void RKF78Sys2tfin(double *, double *, double, void *,
             void (*)(double, double *, unsigned char, double *, void *));
