#include <stdio.h>
#include <math.h>

void SteepestDescent_O2(double A[2][2], double b[2], double gCriteria, double x[2]);
void Direction_SD(double A[2][2], double b[2], double x[2], double d[2]);
double Lambda(double A[2][2], double g[2]);

void SteepestDescent_O2(double A[2][2], double b[2], double gCriteria, double x[2])
{
    double d[2];
    double gLength;
    double lambda;
    int nIter;

    nIter = 0;
    Direction_SD(A, b, x, d);
    printf("nIter x0 x1 gradient_length\n");
    do
    {
        printf("%d %f %f ", nIter, x[0], x[1]);
        lambda = Lambda(A, d);
        x[0] += lambda * d[0];
        x[1] += lambda * d[1];
        Direction_SD(A, b, x, d);
        gLength = sqrt(d[0] * d[0] + d[1] * d[1]);
        printf("%f\n", gLength);
        nIter += 1;
    } while (gLength > gCriteria);
}

void Direction_SD(double A[2][2], double b[2], double x[2], double d[2])
{
    d[0] = b[0] - (A[0][0] * x[0] + A[0][1] * x[1]);
    d[1] = b[1] - (A[1][0] * x[0] + A[1][1] * x[1]);
}

double Lambda(double A[2][2], double g[2])
{
    double lambda;
    lambda = g[0] * g[0] + g[1] * g[1];
    lambda /= g[0] * (A[0][0] * g[0] + A[0][1] * g[1]) + g[1] * (A[1][0] * g[0] + A[1][1] * g[1]);
    return lambda;
}

int main()
{
    double A[2][2];
    double b[2];
    double x[2];

    /* input */
    A[0][0] = 1;
    A[0][1] = 2;
    A[1][0] = 3;
    A[1][1] = 4;
    b[0] = 17;
    b[1] = 39;

    /* processing and output */
    x[0] = 2;
    x[1] = 3;
    printf("------Steepest descent demo------\n");
    printf("Equation:\n")
    printf("╭╰%f %fv╮╯")
    printf("╭%f %f╮")
    SteepestDescent_O2(A, b, 0.001, x);
}