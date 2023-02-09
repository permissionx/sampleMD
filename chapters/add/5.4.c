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
        IterRun(timeStep);
        PBC_r();
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
