/* global variables */
char dynamicStyle[20];

/* function declarations */
void IterRun(double timeStep);
void Dynamics(double stopTime, double timeStep);
/* functions */
void Dynamics(double stopTime, double timeStep)
{
    double time;
    int n, d;

    time = 0;
    while (time <= stopTime)
    {
        PBC_r();
        NeighborList(0);
        Potential(0, 1);
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[n].acceleration[d] = atoms[n].force[d] / typeMasses[atoms[n].type];
            }
        }
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
    }
}

void IterRun(double timeStep)
{
    if (strcmp(dynamicStyle, "Euler") == 0)
        IterRun_Euler(timeStep);
    else if (strcmp(dynamicStyle, "Verlet") == 0)
        IterRun_Verlet(timeStep);
    else if (strcmp(dynamicStyle, "VelocityVerlet") == 0)
        IterRun_VelocityVerlet(timeStep);
}
