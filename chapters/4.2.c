#include <stdio.h>
#include <math.h>
#define MAX_ORDER 10 // MAX ORDER

void Minimize_SD(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double tol, double x[MAX_ORDER]);
void Gradient(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER], double gradient[MAX_ORDER]);
double Lambda(int order, double A[MAX_ORDER][MAX_ORDER], double gradient[MAX_ORDER], double direction[MAX_ORDER]);
double f(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER]);

void Minimize_SD(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double tol, double x[MAX_ORDER])
{
    double gradient[MAX_ORDER];
    double lambda;
    double y, y1, delta;
    int i;

    y = f(order, A, b, x);
    do
    {
        Gradient(order, A, b, x, gradient);
        lambda = Lambda(order, A, gradient, gradient);
        for (i = 0; i < order; i++)
        {
            x[i] += lambda * gradient[i];
        }
    y1 = f(order, A, b ,x);
    delta = fabs(y1 - y);
    y = y1;
    } while (delta > tol);
}

void Gradient(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER], double gradient[MAX_ORDER])
{
    int i, j;
    for (i = 0; i < order; i++)
    {
        gradient[i] = 0;
        for (j = 0; j < order; j++)
        {
            gradient[i] += A[i][j] * b[j];
        }
    }
}

double Lambda(int order, double A[MAX_ORDER][MAX_ORDER], double gradient[MAX_ORDER], double direction[MAX_ORDER])
{
    int i,j;
    double ad_tmp, denominator_tmp, numerator_tmp;

    numerator_tmp = 0;
    denominator_tmp = 0;
    for (i = 0; i < order; i++)
    {
        ad_tmp = 0;
        for (j = 0; j < order; j++)
        {
            ad_tmp += A[i][j] * direction[j];
        }
        denominator_tmp += ad_tmp * gradient[i];
        numerator_tmp += gradient[i] * gradient[i];
    }
    return numerator_tmp / denominator_tmp;
}

double f(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER])
{
    int i,j;
    double result;
    double ax_tmp;
    result = 0;
    for (i = 0; i < order; i++)
    {
        ax_tmp = 0;
        for (j = 0; j < i; j++)
        {
            ax_tmp += A[i][j] * x[j];
        }
        result += x[i] * ax_tmp * 0.5 - b[i] * x[i];
    }
    return result;
}

int main()
{
    double A[MAX_ORDER][MAX_ORDER];
    double b[MAX_ORDER];
    double x[MAX_ORDER];
    double result;
    double order;

    x[0] = 2;
    x[1] = 2;
    A[0][0] = 3;
    A[0][1] = 2;
    A[1][0] = 2;
    A[1][1] = 6;
    b[0] = 2;
    b[1] = -8;
    order = 2;
    result = f(order, A, b, x);
    printf("%f\n", result);
}