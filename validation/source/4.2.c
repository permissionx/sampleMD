#include <stdio.h>
#include <math.h>
#define MAX_ORDER 10 // MAX ORDER

void Minimize_CG(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER]);
void Gradient(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER], double gradient[MAX_ORDER]);
double Lambda(int order, double A[MAX_ORDER][MAX_ORDER], double gradient[MAX_ORDER], double direction[MAX_ORDER]);
double f(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER]);
double Norm(int order, double vector[MAX_ORDER]);
void Direction_CG(int order, double gradient[MAX_ORDER], double lastGradient[MAX_ORDER], double direction[MAX_ORDER]);

void Minimize_CG(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER])
{
    double gradient[MAX_ORDER], lastGradient[MAX_ORDER], direction[MAX_ORDER];
    double gradientLength;
    double lambda;
    double y;
    int i;
    int nstep;

    printf("nstep ");
    for (i = 0; i < order; i++)
    {
        printf("x%d ", i + 1);
    }
    printf("gradient_length y\n");

    Gradient(order, A, b, x, gradient);
    gradientLength = Norm(order, gradient);
    for (i = 0; i < order; i++)
    {
        direction[i] = -gradient[i];
        lastGradient[i] = gradient[i];
    }
    printf("  0 ");
    for (i = 0; i < order; i++)
    {
        printf("%15.8f ",x[i]);
    }
    y = f(order, A, b, x);
    printf("%15.8f %15.8f\n", gradientLength, y);
    for (nstep = 1; nstep <= order; nstep++)
    {
        printf("%3d ", nstep);
        lambda = Lambda(order, A, gradient, direction);
        for (i = 0; i < order; i++)
        {
            x[i] += lambda * direction[i];
            printf("%15.8f ", x[i]);
        }
        Gradient(order, A, b, x, gradient);
        gradientLength = Norm(order, gradient);
        y = f(order, A, b, x);
        printf("%15.8f %15.8f\n", gradientLength, y);
        Direction_CG(order, gradient, lastGradient, direction);
        for (i = 0; i < order; i++)
        {
            lastGradient[i] = gradient[i];
        }
    }
}

void Gradient(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER], double gradient[MAX_ORDER])
{
    int i, j;
    for (i = 0; i < order; i++)
    {
        gradient[i] = 0;
        for (j = 0; j < order; j++)
        {
            gradient[i] += A[i][j] * x[j];
        }
        gradient[i] -= b[i];
    }
}

double Lambda(int order, double A[MAX_ORDER][MAX_ORDER], double gradient[MAX_ORDER], double direction[MAX_ORDER])
{
    int i, j;
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
        denominator_tmp += ad_tmp * direction[i];
        numerator_tmp += gradient[i] * gradient[i];
    }
    return numerator_tmp / denominator_tmp;
}

double f(int order, double A[MAX_ORDER][MAX_ORDER], double b[MAX_ORDER], double x[MAX_ORDER])
{
    int i, j;
    double result;
    double ax_tmp;
    result = 0;
    for (i = 0; i < order; i++)
    {
        ax_tmp = 0;
        for (j = 0; j < order; j++)
        {
            ax_tmp += A[i][j] * x[j];
        }
        result += x[i] * ax_tmp * 0.5 - b[i] * x[i];
    }
    return result;
}

double Norm(int order, double vector[MAX_ORDER])
{
    int i;
    double norm;
    norm = 0;
    for (i = 0; i < order; i++)
    {
        norm += vector[i] * vector[i];
    }
    return sqrt(norm);
}

void Direction_CG(int order, double gradient[MAX_ORDER], double lastGradient[MAX_ORDER], double direction[MAX_ORDER])
{
    int i;
    double numerator, denominator, preD;
    numerator = 0;
    denominator = 0;
    for (i = 0; i < order; i++)
    {
        numerator += gradient[i] * gradient[i];
        denominator += lastGradient[i] * lastGradient[i];
    }
    preD = numerator / denominator;
    for (i = 0; i < order; i++)
    {
        direction[i] = -gradient[i] + preD * direction[i];
    }
}

int main()
{
    double A[MAX_ORDER][MAX_ORDER];
    double b[MAX_ORDER];
    double x[MAX_ORDER];
    double result;
    int order;

    x[0] = 1.25;
    x[1] = -1.08073;
    A[0][0] = 3;
    A[0][1] = 2;
    A[1][0] = 2;
    A[1][1] = 6;
    b[0] = 2;
    b[1] = -8;
    order = 2;

    printf("------Conjugate gradient------\n");
    printf("Equation:\n");
    printf("╭ %5.2f %5.2f ╮╭ x1 ╮ = ╭ %5.2f ╮\n", A[0][0], A[0][1], b[0]);
    printf("╰ %5.2f %5.2f ╯╰ x2 ╯   ╰ %5.2f ╯\n", A[1][0], A[1][1], b[1]);
    printf("\n");

    Minimize_CG(order, A, b, x);
}
