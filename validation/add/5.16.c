void Dynamics(double stopTime, double timeStep)
{
    double time;
    int n, d;
    double temperature;
    double stress[6];
    double targetStress[3] = {-0.5, -0.5, -0.5};
    FILE *fp;
    char fileName[50] = "output/5.16_time-pressure.csv";
    char dumpName[50] = "output/5.16_run.dump";

    InitDynamic();

    time = 0;
    nStep = 0;
    fp = fopen(fileName, "w");
    fprintf(fp, "step time temperature strxx stryy strzz lx ly lz\n");
    printf("step time temperature strxx stryy strzz lx ly lz\n");
    Dump_lammpstrj(dumpName, 1, nStep);
    while (time <= stopTime)
    {
        if (nStep % 100 == 0)
        {
            temperature = ComputeTemperature();
            ComputeStress(stress);
        }
        if (nStep % 100 == 0)
        {
            fprintf(fp, "%d %f %f %f %f %f %f %f %f \n", nStep, time, temperature, stress[0], stress[1], stress[2], boxTranVecs[0][0], boxTranVecs[1][1], boxTranVecs[2][2]);
            printf("%d %f %f %f %f %f %f %f %f \n", nStep, time, temperature, stress[0], stress[1], stress[2], boxTranVecs[0][0], boxTranVecs[1][1], boxTranVecs[2][2]);
            Dump_lammpstrj(dumpName, 0, nStep);
        }
        Barostat(stress, targetStress, 100, timeStep, "Andersen");
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
    }
    fclose(fp);
}

/* main */
int main()

    /* parameters */
    unsigned int randomSeed;
    randomSeed = 1;
    srand(randomSeed);

    typeMasses[1] =   20.1797; // for Ne
    InitMassUnit();
    strcpy(potentialName, "LJ");
    potentialCutoff_LJ = 10;
    neighborCutoff = 15;
    neighborInterval = 100;
    strcpy(dynamicStyle, "VelocityVerlet");

    /* processing*/
    ConstructStdCrystal_BCC(3.15, 10);
    InitVelocity(300.0);
    Dynamics(10.0, 0.0001);

    return 0;
}
