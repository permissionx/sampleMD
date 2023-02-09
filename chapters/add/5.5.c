/* function declarations */
void ComputeAcceleration();
void IterRun_Euler(double timeStep);

/* functions */
void ComputeAcceleration()
{
    int n,d;
    NeighborList(0);
    Potential(0, 1);
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].acceleration[d] = atoms[n].force[d] / typeMasses[atoms[n].type];
        }
    }
}

void IterRun_Euler(double timeStep)
{
    int n, d;
    ComputeAcceleration();
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].r[d] += atoms[n].velocity[d] * timeStep;
            atoms[n].velocity[d] += atoms[n].acceleration[d] * timeStep;
        }
    }
}