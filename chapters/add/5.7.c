/* classes */
struct Atom
{
    int id, type;
    double r[3];
    double reR[3];
    double boxReR[3];
    int neighborNumber;
    struct AtomNeighbor neighbors[MAX_NEIGHBOR_NUMBER];
    double force[3];
    double potentialEnergy;
    double rho_EAM;
    double af_EAM;
    double minDirection[3];
    double lastForce_CG[3];
    double startR_lineMin[3];
    double velocity[3];
    double acceleration[3];
    // Verlet
    double lastR_verlet[3];
    // Velocity Verlet
    double lastA_vverlet[3];
};

/* function declarations */
viod IterRun_VelocityVerlet(double timeStep);

/* functions */
void IterRun_VelocityVerlet(double timeStep)
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
            atoms[n].r[d] += atoms[n].velocity[d] * timeStep + 0.5 * atoms[n].acceleration[d] * timeStep * timeStep;
        }
    }
    PBC_r();
    ComputeAcceleration();
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].velocity[d] += 0.5 * (atoms[n].acceleration[d] + atoms[n].lastA_vverlet[d]) * timeStep;
        }
    }
    
}