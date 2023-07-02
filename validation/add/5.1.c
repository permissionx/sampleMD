/* constants */
// Math constants
#define PI 3.14159265358979323846
// Physical constants
#define K_B 0.0000861733          // unit: eV/K
#define UNIT_MASS 0.0001036426965 // unit: eV/(ps/A)^2

/* classes */
struct Atom
{
    ...
    double velocity[3];
};

/* global variables */
double typeMasses[5];

/* function declarations */
double GaussianRandom(double mu, double sigma);
void InitMassUnit();
void VelocityMaxwell(double temperature);

/* functions */
double GaussianRandom(double mu, double sigma)
{
    double u1, u2;
    double r;
    double z;
    do
    {
        u1 = -1 + ((double)rand() / RAND_MAX) * 2;
        u2 = -1 + ((double)rand() / RAND_MAX) * 2;
        r = u1 * u1 + u2 * u2;
    } while (r >= 1 || r == 0);
    z = sqrt(-2.0 * log(r) / r) * u1 * sigma + mu;
    return z;
}

void InitMassUnit()
{
    int i;
    for (i = 0; i < 5; i++)
    {
        typeMasses[i] *= UNIT_MASS;
    }
}

void VelocityMaxwell(double temperature)
{
    int i, d;
    double mass, sigma;
    for (i = 0; i < atomNumber; i++)
    {
        mass = typeMasses[atoms[i].type];
        for (d = 0; d < 3; d++)
        {
            sigma = sqrt(K_B * temperature / mass);
            atoms[i].velocity[d] = GaussianRandom(0.0, sigma);
        }
    }
}

/* main */
int main()
{
    /* parameters */
    typeMasses[1] = 183.85;
    InitMassUnit();
    unsigned int randomSeed;
    randomSeed = 1;
    srand(randomSeed);

    /* processing */
    ConstructStdCrystal_BCC(3.14, 20);
    double temperature;
    temperature = 300;
    VelocityMaxwell(temperature);

    /* Output */
    int i;
    double speed;
    FILE *fp;
    fp = fopen("output/5.1_velocity-maxwell.txt", "w");
    fprintf(fp, "velocity_x speed\n");
    for (i = 0; i < atomNumber; i++)
    {     
        speed = sqrt(atoms[i].velocity[0] * atoms[i].velocity[0] + atoms[i].velocity[1] * atoms[i].velocity[1] + atoms[i].velocity[2] * atoms[i].velocity[2]);
        fprintf(fp, "%f %f\n", atoms[i].velocity[0], speed);
    }
    fclose(fp);

    return 0;
}
