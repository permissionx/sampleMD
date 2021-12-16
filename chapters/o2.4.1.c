#include <stdio.h>
#include <math.h>
#define MAX_ORDER 10

void SteepestDescent_O2(double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double gCriteria, double x[MAX_ORDER]);
void Direction_SD(double A[2][2], double b[2], double x[2], double d[2]);
double Lambda(double A[2][2], double g[2]);

void SteepestDescent_O2(double A[2][2], double b[2], double gCriteria, double x[2])
{
    double d[2];
    double gNorm;
    double lambda;
    int i;

    i = 0;
    Direction_SD(A, b, x, d);
    printf("i x1 x2 gradient_norm\n");
    do
    {
        printf("%d %f %f ", i, x[0], x[1]);
        lambda = Lambda(A, d);
        x[0] += lambda * d[0];
        x[1] += lambda * d[1];
        Direction_SD(A, b, x, d);
        gNorm = sqrt(d[0] * d[0] + d[1] * d[1]);
        printf("%f\n", gNorm);
        i += 1;
    } while (gNorm > gCriteria);
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

    /* parameters */
    double gCriteria = 0.01;
    A[0][0] = 3;
    A[0][1] = 2;
    A[1][0] = 2;
    A[1][1] = 6;
    b[0] = 2;
    b[1] = -8;

    /* processing and output */
    x[0] = 1.25;
    x[1] = -1.08073;
    printf("------Steepest descent demo------\n");
    printf("Equation:\n");
    printf("╭ %5.2f %5.2f ╮╭ x1 ╮ = ╭ %5.2f ╮\n", A[0][0], A[0][1], b[0]);
    printf("╰ %5.2f %5.2f ╯╰ x2 ╯   ╰ %5.2f ╯\n", A[1][0], A[1][1], b[1]);
    printf("\n");
    SteepestDescent_O2(A, b, gCriteria, x);
}