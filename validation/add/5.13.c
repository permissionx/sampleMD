void Dynamics(double stopTime, double timeStep)
{
    double time;
    int n, d;
    double temperature, pressure;
    double stress[6];
    time = 0;
    nStep = 0;
    printf("step time temperature pressure\n");
    temperature = ComputeTemperature();
    ComputeStress(stress);
    pressure = -(stress[0] + stress[1] + stress[2]) / 3;
    printf("%d %f %f %f\n", nStep, time, temperature, pressure);
    while (time <= stopTime)
    {
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
        if (nStep % 100 == 0)
        {
            temperature = ComputeTemperature();
            ComputeStress(stress);
            pressure = -(stress[0] + stress[1] + stress[2]) / 3;
            printf("%d %f %f %f\n", nStep, time, temperature, pressure);
        }
    }
}

/* main */
int main()
{
    /* parameters */
    unsigned int randomSeed;
    randomSeed = 1;
    srand(randomSeed);

    typeMasses[1] = 183.84; // for W
    InitMassUnit();
    strcpy(potentialName, "EAM");
    neighborCutoff = 6;
    neighborInterval = 100;
    strcpy(dynamicStyle, "VelocityVerlet");

    /* processing*/
    double temperature;
    for (temperature = 100; temperature < 401; temperature += 50)
    {
        printf("----Case for initial temperature of %f K----\n", temperature);
        ConstructStdCrystal_BCC(3.14, 10);
        InitVelocity(temperature);
        Dynamics(1.0, 0.0005);
        printf("----------\n\n");
    }

    return 0;
}
