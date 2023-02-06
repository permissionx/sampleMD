#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.1415926

double randn2(double mu, double sigma) {
  double U1, U2, W, mult;
  double X1;
  
  do {
    U1 = ((double) rand () / RAND_MAX);
    U2 = ((double) rand () / RAND_MAX);
    W = U1 * U1 + U2 * U2;
  } while (W >= 1 || W == 0);
  X1 = sqrt(-2*log(U1)) * cos(2*PI*U2);
  return X1;
}

double randn(double mu, double sigma) {
  double U1, U2, W, mult;
  double X1, X2;
  
  do {
    U1 = -1 + ((double) rand () / RAND_MAX) * 2;
    U2 = -1 + ((double) rand () / RAND_MAX) * 2;
    W = U1 * U1 + U2 * U2;
  } while (W >= 1 || W == 0);
  
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
  
  
  return (mu + sigma *  X1);
}

int main()
{
    srand((unsigned)time(NULL)); 
    FILE *fp;
    fp = fopen("test.csv","w");
    fprintf(fp, "nums\n");
    int i;
    double mu = 0;
    double sigma = 1;
    double rand;
    for (i = 0; i<10000;i++)
    {
        rand = randn2(mu, sigma);
        fprintf(fp,"%f\n", rand);
    }
    fclose(fp);
}