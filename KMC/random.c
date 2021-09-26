#include <stdio.h>

void RandomNumber(int number, int a, int c, int m, int seed)
{
    int i, n;
    i = seed;
    for (n = 0; n < number; n++)
    {
        printf("%d\n",i);
        i = (i*a+c)%m;
    }

}

int main()
{
    int a = 107;
    int c = 31;
    int m = 10001;
    int seed = 5;
    int number = 1000;
    RandomNumber(number, a, c, m, seed);
}