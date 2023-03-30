void Dynamics(double stopTime, double timeStep)
{
    double time;
    int n, d;
    double temperature;
    double targetTemperature = 300;
    FILE *fp;
    char fileName[50] = "debug/thermostat/time-temperature.nh.csv";
    char dumpName[50] = "debug/thermostat/run.nh.dump";

    InitDynamic();

    time = 0;
    nStep = 0;
    fp = fopen(fileName, "w");
    fprintf(fp, "step time temperature\n");
    Dump_lammpstrj(dumpName, 1, nStep);
    while (time <= stopTime)
    {
        if (nStep % 100 == 0)
        {
            temperature = ComputeTemperature();
        }
        if (nStep % 100 == 0)
        {
            fprintf(fp, "%d %f %f\n", nStep, time, temperature);
            printf("%d %f %f\n", nStep, time, temperature);
            Dump_lammpstrj(dumpName, 0, nStep);
        }
        if (nStep >= 4000)
        {
            Thermostat(temperature, targetTemperature, 100, timeStep, "Berendsen");
        }
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
    }
    fclose(fp);
}

/* main */
int main()
{
    /* parameters */
    double randomSeed;
    randomSeed = 1.0;
    srand(randomSeed);

    typeMasses[1] = 183.84; // for W
    InitMassUnit();
    strcpy(potentialName, "EAM");
    neighborCutoff = 6;
    neighborInterval = 100;
    strcpy(dynamicStyle, "VelocityVerlet");

    /* processing*/
    ConstructStdCrystal_BCC(3.14, 10);
    InitVelocity(100);
    Dynamics(1, 0.0001);

    return 0;
}
