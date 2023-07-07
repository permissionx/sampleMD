/* function declarations*/
void Dynamics(double stopTime, double timeStep, double targetTemperature, double result[2])

/* functions */
void Dynamics(double stopTime, double timeStep, double targetTemperature, double result[2])
{
    double time;
    int n, d;
    int count;
    double temperature, pressure, ave_t, ave_p;
    double stress[6];

    InitDynamic();

    time = 0;
    nStep = 0;
    count = 0;
    ave_t = 0;
    ave_p = 0;

    while (time <= stopTime)
    {
        temperature = ComputeTemperature();
        if (time >= 8.0 && nStep%100==0)
        {
            ComputeStress(stress);
            pressure = -(stress[0] + stress[1] + stress[2])/3;
            ave_t += temperature;
            ave_p += pressure;
            count += 1;
        }
        Thermostat(temperature, targetTemperature, 1, timeStep, "Nose-Hoover");
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
    }
    result[0] = ave_t/count;
    result[1] = ave_p/count;
}

/* main */
int main()
{
    /* parameters */
    unsigned int randomSeed;
    randomSeed = 1;
    srand(randomSeed);
    typeMasses[1] = 20.1797;
    InitMassUnit();
    strcpy(potentialName, "LJ");
    potentialCutoff_LJ = 5;
    neighborCutoff = 6;
    neighborInterval = 100;
    strcpy(dynamicStyle, "VelocityVerlet");

    /* processing*/
    FILE *fp;
    char fileName[50];
    double targetTemperature, temperature;
    double scale, pressure, volume;
    double stress[6];
    double result[2]; // store average pressure and volume
    strcpy(fileName, "output/7.1_isotherm.txt");
    fp = fopen(fileName, "w");
    fprintf(fp, "targetTemperature volume pressure temperature\n");
    printf("targetTemperature volume pressure temperature\n");
    for (targetTemperature = 100; targetTemperature < 901; targetTemperature += 200)
    {
        for (scale = 7; scale < 14; scale += 1)
        {
            ConstructStdCrystal_FCC(4.23 * pow(scale, 1. / 3.), 5);
            InitVelocity(targetTemperature);
            Dynamics(10, 0.0001, targetTemperature, result);
            temperature = result[0];
            pressure = result[1];
            volume = ComputeBoxVolume();
            printf("%f %f %f %f\n", targetTemperature, volume, pressure, temperature);
            fprintf(fp, "%f %f %f %f\n", targetTemperature, volume, pressure, temperature);
        }          
    }
    fclose(fp);
    return 0;
}
