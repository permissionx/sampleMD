#include <stdio.h>
#include <math.h>
#define MO 10 // MAX ORDER

void SteepestDescent(int order, double A[MO][MO], double b[MO], double gCriteria, double x[MO]);
void Direction_SD(int order, double A[MO][MO], double b[MO], double x[MO], double d[MO]);
double Lambda(int order, double A[MO][MO], double g[MO]);
double Norm(int order, double vec[MO]);

void SteepestDescent(int order, double A[MO][MO], double b[MO], double gCriteria, double x[MO])
{
    double d[MO];
    double gNorm;
    double lambda;
    int i;
    int o;

    i = 0;
    Direction_SD(order, A, b, x, d);
    printf("i ");
    for (o = 0; o < order; o++)
    {
        printf("x%d ", o + 1);
    }
    printf("gradient_norm\n");
    do
    {
        printf("%d ", i);
        for (o = 0; o < order; o++)
        {
            printf("%f ", x[o]);
        }
        lambda = Lambda(order, A, d);
        for (o = 0; o < order; o++)
        {
            x[o] += lambda * d[o];
        }
        Direction_SD(order, A, b, x, d);
        gNorm = Norm(order, d);
        printf("%f\n", gNorm);
        i += 1;
    } while (gNorm > gCriteria);
}

void Direction_SD(int order, double A[MO][MO], double b[MO], double x[MO], double d[MO])
{
    int i, j;
    for (i = 0; i < order; i++)
    {
        d[i] = 0;
        for (j = 0; j < order; j++)
        {
            d[i] += A[i][j] * x[j];
        }
        d[i] = b[i] - d[i];
    }
}

double Lambda(int order, double A[MO][MO], double g[MO])
{
    double lambda;
    int i, j;
    double up, down, tmp;
    up = 0;
    down = 0;
    for (i = 0; i < order; i++)
    {
        tmp = 0;
        for (j = 0; j < order; j++)
        {
            tmp += A[i][j] * g[j];
        }
        down += g[i] * tmp;
        up += g[i] * g[i];
    }
    lambda = up / down;
    return lambda;
}

double Norm(int order, double vec[MO])
{
    double norm;
    int o;
    norm = 0;
    for (o = 0; o < order; o++)
    {
        norm += vec[o] * vec[o];
    }
    norm = sqrt(norm);
    return norm;
}

int main()
{
    double A[MO][MO];
    double b[MO];
    double x[MO];

    /* parameters */
    int order = 2;
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
    SteepestDescent(order, A, b, gCriteria, x);
}