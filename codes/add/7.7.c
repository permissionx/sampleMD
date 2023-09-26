/* function declarations */
void Dynamics(double stopTime, double timeStep);

/* functions */
void Dynamics(double stopTime, double timeStep)
{
    double time;
    int n, d;
    double temperature;
    double targetTemperature = 1000;

    InitDynamic();

    time = 0;
    nStep = 0;
    while (time <= stopTime)
    {
        temperature = ComputeTemperature();
        if (nStep % 100 == 0)
        {
            Dump_lammpstrj("output/7.7_dynamic-recovery.lammpstrj", 0, nStep);
            printf("%d %f\n",nStep, temperature);
        }
        Thermostat(temperature, targetTemperature, 100, timeStep, "Nose-Hoover");
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
    }
}

/* main */
int main()
{
    int n;

    /* parameters */
    unsigned int randomSeed;
    randomSeed = 1;
    srand(randomSeed);
    typeMasses[1] = 183.84; // for W
    InitMassUnit();
    strcpy(potentialName, "EAM");
    strcpy(minimizeStyle, "CG");
    strcpy(dynamicStyle, "VelocityVerlet");
    neighborCutoff = 6;
    neighborInterval = 100;

    /* processing*/
    double siaPosition[3] = {5.25 * 3.14, 5.25 * 3.14, 5.25 * 3.14};
    ConstructStdCrystal_BCC(3.14, 10);
    InsertAtom(siaPosition, 1);
    Minimize();
    for (n = 0; n < atomNumber; n++)
    {
        if (atoms[n].id == 893)
        {
            DeleteAtomByIndex(n);
            atomNumber -= 1;
            break;
        }
    }
    Dump_lammpstrj("output/7.7_static-recovery.lammpstrj", 1, 0);
    Minimize();
    Dump_lammpstrj("output/7.7_static-recovery.lammpstrj", 0, 1);

    InitVelocity(2000);
    Dynamics(1, 0.0001);
    Minimize();
    Dump_lammpstrj("output/7.7_dynamic-recovery.lammpstrj", 0, nStep);

    return 0;
}
