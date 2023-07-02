/* function declarations*/
double ComputeTotalKineticEnergy();

/* functions */
double ComputeTotalKineticEnergy()
{
    int n, d;
    double e = 0;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            e += 0.5 * typeMasses[atoms[n].type] * atoms[n].velocity[d] * atoms[n].velocity[d];
        }
    }
    return e;
}

void Dynamics(double stopTime, double timestep)
{
    double time;
    int n, d;
    double totalEnergy;
    char dumpName[30];

    strcpy(dumpName, "output/5.9_bcc-run.lammpstrj");
    Dump_lammpstrj(dumpName, 1, 0);
    printf("step time totalEnergy\n");
    NeighborList(0);
    Potential(1, 0);
    totalEnergy = totalPotentialEnergy + ComputeTotalKineticEnergy();
    printf("%d %f %f\n",0, 0.0, totalEnergy);

    time = 0;
    nStep = 0;
    while (time <= stopTime)
    {
        IterRun(timestep);
        nStep += 1;
        time += timestep;
        if (nStep % 100 == 0)
        {
            Dump_lammpstrj(dumpName, 0, nStep);
            NeighborList(0);
            Potential(1, 0);
            totalEnergy = totalPotentialEnergy + ComputeTotalKineticEnergy();
            printf("%d %f %f\n",nStep, time, totalEnergy);
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

    typeMasses[1] = 183.85;
    InitMassUnit();
    strcpy(potentialName, "EAM");
    neighborCutoff = 6;
    neighborInterval = 100;
    strcpy(dynamicStyle, "Euler");

    /* processing*/
    ConstructStdCrystal_BCC(3.14, 10);
    InitVelocity(300);
    Dynamics(1, 0.001);
  
    return 0;
}

