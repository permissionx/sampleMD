void Dynamics(double stopTime, double timeStep)
{
    double time;
    int n, d;
    double temperature;
    double stress[6];
    char dumpName[30];
    time = 0;
    nStep = 0;
    printf("step time temperature stress_xx\n");
    temperature = ComputeTemperature();
    ComputeStress(stress);
    printf("%d %f %f %f\n",nStep, time, temperature, stress[0]);
    while (time <= stopTime)
    {
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
        if (nStep % 100 == 0)
        {
            temperature = ComputeTemperature();
            ComputeStress(stress);
            printf("%d %f %f %f\n",nStep, time, temperature, stress[0]);
        }
    }
}

/* main */
int main()
{
    /* parameters */
    double randomSeed;
    randomSeed = 1.0;
    srand(randomSeed);

    typeMasses[1] = 20.1797; // for Ne
    InitMassUnit();
    strcpy(potentialName, "LJ");
    potentialCutoff_LJ = 20;
    neighborCutoff = 20;
    neighborInterval = 100;
    strcpy(dynamicStyle, "VelocityVerlet");

    /* processing*/
    double temperature;
    for (temperature = 100; temperature <= 400; temperature += 50)
    {
        printf("----Case for initial temperature of %f K----\n", temperature);
        ConstructStdCrystal_FCC(4.23, 5);
        InitVelocity(temperature);
        Dynamics(1.0, 0.0005);
        printf("----------\n\n");
    }

    return 0;
}
