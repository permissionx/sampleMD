/* function declarations */
double ComputeTemperature();
/* functions */
double ComputeTemperature()
{
    double Ek, T;
    Ek = ComputeTotalKineticEnergy();
    T = 2/3/K_B/atomNumber*Ek;
    return T;
}