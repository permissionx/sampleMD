/* function declarations */
void InitVelocity(double temperature);

/* functions */
void InitVelocity(double temperature)
{
    VelocityMaxwell(temperature);
    ZeroMomentum();
}