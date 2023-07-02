/* classes */
struct Atom
{
    …
    // Velocity Verlet
    double lastA_vverlet[3];
};

/* function declarations */
void IterRun_VelocityVerlet(double timestep);

/* functions */
void IterRun_VelocityVerlet(double timestep)
{
    int n, d;
    if (nStep == 0)
    {
        ComputeAcceleration();
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[n].lastA_vverlet[d] = atoms[n].acceleration[d];
            }
        }
    }
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].r[d] += atoms[n].velocity[d] * timestep + 0.5 * atoms[n].acceleration[d] * timestep * timestep;
        }
    }
    PBC_r();
    ComputeAcceleration();
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].velocity[d] += 0.5 * (atoms[n].acceleration[d] + atoms[n].lastA_vverlet[d]) * timeStep;
            atoms[n].lastA_vverlet[d] = atoms[n].acceleration[d];
        }
    }
    
}

