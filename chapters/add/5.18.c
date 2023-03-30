/* global variables */
double xi_NoseHoover;

/* function declarations */
void Thermostat(double temperature, double targetTemperature, int frequency, double timeStep, char algorithm[20]);
void Thermostat_NoseHoover(double temperature, double targetTemperature, int frequency, double timeStep);

/* functions */
void Thermostat(double temperature, double targetTemperature, int frequency, double timeStep, char algorithm[20])
{
    if (strcmp(algorithm, "Berendsen") == 0)
    {
        Thermostat_Berendsen(temperature, targetTemperature, frequency, timeStep);
    }
    else if (strcmp(algorithm, "Nose-Hoover") == 0)
    {
        Thermostat_NoseHoover(temperature, targetTemperature, frequency, timeStep);
    }
    else
    {
        printf("Thermostat algorithm %s not found.", algorithm);
        exit(1);
    }
}

void Thermostat_NoseHoover(double temperature, double targetTemperature, int frequency, double timeStep)
{
    static int count = 0;
    double Q = 0.1; // parameter
    double xi_velocity;
    int n, d;
    if (count == 0)
    {
        xi_velocity = 3 * atomNumber * K_B / Q * (temperature - targetTemperature);
        xi_NoseHoover += xi_velocity * frequency * timeStep;
    }
    count += 1;
    if (count == frequency)
    {
        count = 0;
    }
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].accelerationModify[d] = -xi_NoseHoover * atoms[n].velocity[d];
        }
    }
}

void InitDynamic()
{
    ...
    xi_NoseHoover = 0;
}