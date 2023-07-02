/* function declarations */
void Barostat_Berendsen(double stress[6], double targetStress[3], int frequency, double timeStep);

/* functions */
void Barostat_Berendsen(double stress[6], double targetStress[3], int frequency, double timeStep)
{
    static int count = 0;
    double k_tau = 1; // parameter
    int n, d;
    double lambda;
    double deltaTime;
    deltaTime = frequency * timeStep;
    if (count == 0)
    {
        for (d = 0; d < 3; d++)
        {
            lambda = 1 + k_tau * deltaTime * (targetStress[d] - stress[d]);
            boxTranVecs[d][d] *= lambda;
            for (n = 0; n < atomNumber; n++)
            {
                atoms[n].r[d] *= lambda;
            }
        }
        PBC_r();
    }
    count += 1;
    if (count == frequency)
    {
        count = 0;
    }
}
