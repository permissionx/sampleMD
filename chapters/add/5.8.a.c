/* functions */
void Dynamics(double stopTime, double timeStep)
{
    double time;
    int n, d;

    double dx, dv;
    FILE *fp;

    fp = fopen("dimer_Euler.txt", "w");
    fprintf(fp, "step time dx dv\n");
    dx = atoms[1].r[0] - atoms[0].r[0] - 3.076;
    dv = atoms[1].velocity[0] - atoms[0].velocity[0];
    fprintf(fp, "%d %f %f %f\n", nStep, time, dx, dv);

    time = 0;
    nStep = 0;
    while (time <= stopTime)
    {
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
        if (nStep % 1000 == 0)
        {
            dx = atoms[1].r[0] - atoms[0].r[0] - 3.076;
            dv = atoms[1].velocity[0] - atoms[0].velocity[0];
            fprintf(fp, "%d %f %f %f\n", nStep, time, dx, dv);
            printf("%d %f\n", nStep, time);
        }
    }
    fclose(fp);
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
    neighborCutoff = 5;
    neighborInterval = 50;
    strcpy(dynamicStyle, "Euler");

    /* processing*/
    ConstructStdCrystal_BCC(3, 10);
    atomNumber = 2;
    atoms[0].r[0] = 10;
    atoms[0].r[1] = 0;
    atoms[0].r[2] = 0;
    atoms[1].r[0] = 13.1;
    atoms[1].r[1] = 0;
    atoms[1].r[2] = 0;

    InitVelocity(0);
    Dynamics(30, 0.0005);

    return 0;
}
