/* classes */
struct Atom
{
    ...
    // Andersen Barostat
    double accelerationModify[3];
    double velocityModify[3];
};

/* function declarations */
void Barostat(double stress[6], double targetStress[3], int frequency, double timeStep, char algorithm[20]);
void Barostat_Andersen(double stress[6], double targetStress[3], int frequency, double timeStep);
void InitDynamic();
void ReviveVelocity();

/* functions */
void Barostat(double stress[3], double targetStress[3], int frequency, double timeStep, char algorithm[20])
{
    if (boxOrthogonal != 1 || !(boxStartPoint[0] == 0 && boxStartPoint[1] == 0 && boxStartPoint[2] == 0))
    {
        printf("Error: Barostat in wrong box, check function Barostat()\n");
        exit(1);
    }
    if (strcmp(algorithm, "Berendsen") == 0)
    {
        Barostat_Berendsen(stress, targetStress, frequency, timeStep);
    }
    else if (strcmp(algorithm, "Andersen") == 0)
    {
        Barostat_Andersen(stress, targetStress, frequency, timeStep);
    }
    else
    {
        printf("Barostat algorithm %s not found.", algorithm);
        exit(1);
    }
}

void Barostat_Andersen(double stress[6], double targetStress[3], int frequency, double timeStep)
{
    if (!(targetStress[0] == targetStress[1] && targetStress[0] == targetStress[2]))
    {
        printf("Error: Andersen barostat can only used for target sigma_xx == sigma_yy == sigma_zz\n");
        exit(1);
    }
    double M = 10; // piston mass, unit: unit: eV/(ps/A)^2

    static int count = 0;
    static double pistonVelocity[3] = {0, 0, 0};
    int n, i;
    double volume;
    double pistonAcceleration; // PistonPi/M
    double deltaTime;
    double AndersenVelocity, pistonVelocityPerLength; // PI/M/L

    deltaTime = frequency * timeStep;

    if (count == 0)
    {
        volume = ComputeBoxVolume();
        for (i = 0; i < 3; i++)
        {
            pistonAcceleration = (targetStress[i] - stress[i]) * volume / boxTranVecs[i][i] / M;
            pistonVelocity[i] += deltaTime * pistonAcceleration;
        }
    }
    count += 1;
    if (count == frequency)
    {
        count = 0;
    }

    for (i = 0; i < 3; i++)
    {
        boxTranVecs[i][i] += pistonVelocity[i] * timeStep;
        pistonVelocityPerLength = pistonVelocity[i] / boxTranVecs[i][i];
        for (n = 0; n < atomNumber; n++)
        {
            AndersenVelocity = atoms[n].r[i] * pistonVelocityPerLength;
            atoms[n].velocity[i] += AndersenVelocity;
            atoms[n].velocityModify[i] = AndersenVelocity;
            atoms[n].accelerationModify[i] = -atoms[n].velocity[i] * pistonVelocityPerLength;
        }
    }
}

void InitDynamic()
{
    int n, d;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].accelerationModify[d] = 0;
            atoms[n].velocityModify[d] = 0;
        }
    }
}

void ReviveVelocity()
{
    int n, d;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].velocity[d] -= atoms[n].velocityModify[d];
        }
    }
}

void IterRun(double timeStep)
{
    ...
    ReviveVelocity();
}

void ComputeAcceleration()
{
    ...
            atoms[n].acceleration[d] = atoms[n].force[d] / typeMasses[atoms[n].type] + atoms[n].accelerationModify[d];
    ...
}


