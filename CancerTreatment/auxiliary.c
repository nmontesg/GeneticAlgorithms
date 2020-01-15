#include <float.h>

/* Tumour Parameters */
#define lambda_par 0.336
#define lambdalogTheta_par 9.284023094951992774680543277273261947643
#define NZero_par 20000
#define CCUMj_par 127
#define NMax_par 9.5e+11
#define m_par 4
#define d_par 10
#define n_par 11

static float k_j[d_par] = {0.12, 0.0502, 0.0637, 0.1347, 0.0902, 0.0546, 0.0767, 0.1121, 0.0971, 0.0403};
static float t_i[n_par] = {0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 33.0};
static float eta_kj_cur[m_par][d_par] = {
  {0.0036, 0.0098, 0.0061, 0.0009, 0.0003, 0.0108, 0.0045, 0.0021, 0.0096, 0.0125},
  {0.0063, 0.0082, 0.0062, 0.0062, 0.0083, 0.013, 0.0039, 0.0019, 0.0015, 0.005},
  {0.0129, 0.0018, 0.0116, 0.0021, 0.009, 0.0129, 0.0054, 0.0049, 0.0093, 0.0066},
  {0.0053, 0.0086, 0.0067, 0.0029, 0.0089, 0.0054, 0.0042, 0.0095, 0.0112, 0.0092}
};
static float eta_kj_pal[m_par][d_par] = {
  {0.00612, 0.01666, 0.01037, 0.00153, 0.00051, 0.01836, 0.00765, 0.00357, 0.01632, 0.02125},
  {0.01071, 0.01394, 0.01054, 0.01054, 0.01411, 0.0221, 0.00663, 0.00323, 0.00255, 0.0085},
  {0.02193, 0.00306, 0.01972, 0.00357, 0.0153, 0.02193, 0.00918, 0.00833, 0.01581, 0.01122},
  {0.00901, 0.01462, 0.01139, 0.00493, 0.01513, 0.00918, 0.00714, 0.01615, 0.01904, 0.01564}
};

typedef struct { double drift_i; } ODE_Parameters;
void Gompertz(double t, double N, double *der, void *Params) {
  *der = ((N < 1.e-16) ? 0.0 :
  N*(lambdalogTheta_par - lambda_par*log(N) - ((ODE_Parameters *) Params)->drift_i));
}


unsigned char TestIfConstraints2and3AreVerifiedCurative(unsigned char *Cij){
  register unsigned char i, j, k;
  unsigned char npar = n_par - 1, *Cijofi;
  for(j=0; j < d_par; j++) {
    unsigned int ccumj = 0U;
    for(i=0; i < npar; i++) ccumj += *(Cij + i*d_par + j)>>4;
    if(ccumj > CCUMj_par) return 0U;
  }
  
  for(i=0, Cijofi=Cij; i < npar; i++, Cijofi += d_par) {
    for(k=0; k < m_par; k++) {
      double Cseffk = 0.0;
      for(j=0; j < d_par; j++) Cseffk += eta_kj_cur[k][j] * (Cijofi[j]>>4);
      if(Cseffk > 1.0) return 0U;
    }
  }
  return 1U;
}

unsigned char TestIfConstraints2and3AreVerifiedPaliative(unsigned char *Cij){
  register unsigned char i, j, k;
  unsigned char npar = n_par - 1, *Cijofi;
  for(j=0; j < d_par; j++) {
    unsigned int ccumj = 0U;
    for(i=0; i < npar; i++) ccumj += *(Cij + i*d_par + j)>>4;
    if(ccumj > CCUMj_par) return 0U;
  }
  
  for(i=0, Cijofi=Cij; i < npar; i++, Cijofi += d_par) {
    for(k=0; k < m_par; k++) {
      double Cseffk = 0.0;
      for(j=0; j < d_par; j++) Cseffk += eta_kj_pal[k][j] * (Cijofi[j]>>4);
      if(Cseffk > 1.0) return 0U;
    }
  }
  return 1U;
}


double Curative_Fitness(unsigned char *Cij) {
  register unsigned char i, j;
  ODE_Parameters GompertzParams;
  double N = NZero_par, t = t_i[0];
  double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
  unsigned char curativecounter = 0U, npar = n_par - 1;
  double integral = 0.0, lastt = t, lastN = N;
  
  if(!TestIfConstraints2and3AreVerifiedCurative(Cij)) return DBL_MAX;
  
  for(i=0; i < npar; i++){
    double tfin = t_i[i+1]; // Implementing treatment i
    GompertzParams.drift_i = 0.0;
    for(j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * (Cij[i*d_par + j]>>4);
    
    while(t+h < tfin) {
      RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
      if(N > NMax_par) return DBL_MAX;
      integral += (lastN + N)*(t - lastt);
      lastt = t; lastN = N;
    }
    
    do {
      h = tfin - t;
      RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
      if(N > NMax_par) return DBL_MAX;
      integral += (lastN + N)*(t - lastt);
      lastt = t; lastN = N;
    } while (t < tfin);
    
    if(N < 1000) { curativecounter++; if(curativecounter > 2) return integral/2.0; }
    else curativecounter = 0U;
  }
  return DBL_MAX;
}


double Paliative_Fitness(unsigned char *Cij) {
  register unsigned char i, j;
  ODE_Parameters GompertzParams;
  double N = NZero_par, t = t_i[0];
  double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
  unsigned char npar = n_par - 1;
  
  if(!TestIfConstraints2and3AreVerifiedPaliative(Cij)) return DBL_MAX;
  
  for(i=0; i < npar; i++){
    double tfin = t_i[i+1]; // Implementing treatment i
    GompertzParams.drift_i = 0.0;
    for(j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * (Cij[i*d_par + j]>>4);
    while(t+h < tfin) {
      RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
      if(N > NMax_par) return 1.0/t;
    }
    
    do {
      h = tfin - t;
      RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
      if(N > NMax_par) return 1.0/t;
    } while (t < tfin);
  }
  return 1.0/t;
}

void setExcessDosagesToZero(unsigned char* Cij) {
  register unsigned char i, j;
  ODE_Parameters GompertzParams;
  double N = NZero_par, t = t_i[0];
  double hmin = 1.e-8, hmax = 1.0, h = 1.e-3, tol = 1.e-8;
  unsigned char curativecounter = 0U, npar = n_par - 1;
  
  if(!TestIfConstraints2and3AreVerifiedCurative(Cij)) return;
  
  for(i=0; i < npar; i++){
    double tfin = t_i[i+1]; // Implementing treatment i
    GompertzParams.drift_i = 0.0;
    for(j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * (Cij[i*d_par + j]>>4);
    
    while(t+h < tfin) {
      RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
      if(N > NMax_par) {
        for (int k = i+1; k < npar; k++) for (j = 0; j < d_par; j++) Cij[k*d_par+j] = 0x00;
        return;
      }
    }
    
    do {
      h = tfin - t;
      RKF78(&t, &N, &h, hmin, hmax, tol, &GompertzParams, Gompertz);
      if(N > NMax_par) {
        for (int k = i+1; k < npar; k++) for (j = 0; j < d_par; j++) Cij[k*d_par+j] = 0x00;
        return;
      }
    } while (t < tfin);
    
    if(N < 1000) {
      curativecounter++;
      if(curativecounter > 2) {
        for (int k = i+1; k < npar; k++) for (j = 0; j < d_par; j++) Cij[k*d_par+j] = 0x00;
        return;
      }
    }
    else curativecounter = 0U;
  } 
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
  fprintf(out, "%.7f;%.0f\n", t, N);
  
  for(i=0; i < npar; i++){
    double tfin = t_i[i+1]; // Implementing treatment i
    GompertzParams.drift_i = 0.0;
    for(j=0; j < d_par; j++) GompertzParams.drift_i += k_j[j] * (Cij[i*d_par + j]>>4);;
    
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
    if(N < 1000) { fclose(out); break; }
  }
  fclose(out);
}
