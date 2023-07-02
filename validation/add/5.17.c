/* function declarations*/
void Thermostat_Berendsen(double temperature, double targetTemperature, int frequency, double timeStep);

/* functions */
void Thermostat_Berendsen(double temperature, double targetTemperature, int frequency, double timeStep)
{
    static int count = 0;
    double lambda, deltaTime, deltaTime_tau; // deltaTime_tau: daltaTime/tau
    int n, d;
    double tau = 0.01;
    if (count == 0)
    {
        deltaTime_tau = frequency * timeStep / tau;
        lambda = sqrt(1 + deltaTime_tau * (targetTemperature / temperature - 1));
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[n].velocity[d] *= lambda;
            }
        }
    }
    count += 1;
    if (count == frequency)
    {
        count = 0;
    }
}
