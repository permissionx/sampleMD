/* function declarations */
void IterRun_Euler(double timeStep);

/* functions */
void IterRun_Euler(double timeStep)
{
    int i, d;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].r[d] += atoms[i].velocity[d] * timeStep;
            atoms[i].velocity[d] += atoms[i].acceleration[d] * timeStep;
        }
    }
}
