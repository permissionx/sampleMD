#include <stdio.h>
    
void RandomNumber(int number, int a, int c, int m, int seed)
{
    int i, n;
    FILE *fp;
    fp = fopen("output/6.1_random-number.txt", "w");
    i = seed;
    for (n = 0; n < number; n++)
    {
        fprintf(fp, "%d\n", i);
        i = (i * a + c) % m;
    }
}

int main()
{
    int a = 1047;
    int c = 341;
    int m = 1000001;
    int seed = 7039;
    int number = 1000000;
    RandomNumber(number, a, c, m, seed);
}
