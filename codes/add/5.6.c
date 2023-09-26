/* classes */
struct Atom
{
    …
    // Verlet
    double lastR_verlet[3];
};

/* function declarations */
void IterRun_Verlet(double timeStep);

/* functions */
void IterRun_Verlet(double timeStep)
{
    int n, d;
    double r_tmp;

    ComputeAcceleration();
    if (nStep == 0)
    {
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[n].lastR_verlet[d] = atoms[n].r[d];
                atoms[n].r[d] += atoms[n].velocity[d] * timeStep + 0.5 * atoms[n].acceleration[d] * timeStep * timeStep;
                atoms[n].velocity[d] = (atoms[n].r[d] - atoms[n].lastR_verlet[d]) / timeStep;
            }
        }
    }
    else
    {
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                r_tmp = atoms[n].r[d];
                atoms[n].r[d] = 2 * atoms[n].r[d] - atoms[n].lastR_verlet[d] + atoms[n].acceleration[d] * timeStep * timeStep;
                atoms[n].lastR_verlet[d] = r_tmp;
                atoms[n].velocity[d] = (atoms[n].r[d] - atoms[n].lastR_verlet[d]) / timeStep;
            }
        }
    }
    PBC_r();
}
