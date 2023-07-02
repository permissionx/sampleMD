/* function declarations */
void ZeroMomentum();

/* functions */
void ZeroMomentum()
{
    int i, d;
    double momentum[3] = {0, 0, 0};

    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            momentum[d] += atoms[i].velocity[d] * typeMasses[atoms[i].type];
        }
    }

    for (d = 0; d < 3; d++)
    {
        momentum[d] /= atomNumber;
    }

    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].velocity[d] -= momentum[d] / typeMasses[atoms[i].type];
        }
    }
}

/* main */
int main()
{
    /* parameters */
    typeMasses[1] = 183.85;
    InitMassUnit();
    unsigned int randomSeed;
    randomSeed = 1;
    srand(randomSeed);

    /* processing & output*/
    ConstructStdCrystal_BCC(3.14, 20);
    double temperature;
    temperature = 300;
    VelocityMaxwell(temperature);

    int i,d;
    double totalMomentum[3] = {0, 0, 0};
    for (i = 0; i < atomNumber; i++)
    {
        for(d = 0;d<3;d++)
        {
            totalMomentum[d] += atoms[i].velocity[d] * typeMasses[atoms[i].type];
        }
    }
    printf("Total totalMomentum before correction (x,y,z): %f %f %f\n", totalMomentum[0], totalMomentum[1], totalMomentum[2]);
    ZeroMomentum();
    totalMomentum[0] = 0; totalMomentum[1] = 0; totalMomentum[2] = 0;
    for (i = 0; i < atomNumber; i++)
    {
        for(d = 0;d<3;d++)
        {
            totalMomentum[d] += atoms[i].velocity[d] * typeMasses[atoms[i].type];
        }
    }
    printf("Total totalMomentum after correction (x,y,z): %f %f %f\n", totalMomentum[0], totalMomentum[1], totalMomentum[2]);

    return 0;
}
