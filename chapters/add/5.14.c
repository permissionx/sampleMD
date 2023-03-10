/* function declarations */
void Barostat_Berendsen(double stress[6], double targetStress[3], double deltaTime);

/* functions */
void Barostat_Berendsen(double stress[6], double targetStress[3], double deltaTime)
{
    double k_tau = 1.0; // parameter
    int n, d;
    double lambda;
    for (d = 0; d < 3; d++)
    {
        lambda = 1 + k_tau * deltaTime * (targetStress[d] - stress[d]);
        boxTranVecs[d][d] *= lambda;
        for (n = 0; n < atomNumber; n++)
        {
            atoms[n].r[d] *= lambda;
        }
    }
}