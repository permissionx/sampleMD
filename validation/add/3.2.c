/* constants */
// L-J parameters for Ne (energy: eV, length: A)
#define LJ_EPSILON 0.0031
#define LJ_SIGMA 2.74
#define LJ_2_EPSILON 0.0062                     // 4 * LJ_EPSILON
#define LJ_24_EPSILON_SIGMA_6 31.48301472289983 // 24 * LJ_EPSILON * pow(LJ_SIGMA, 6)
#define LJ_2_SIGMA_6 846.3176000779524          // 2 * pow(LJ_SIGMA, 6)

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
};

/* global variables */
double potentialCutoff_LJ;
double totalPotentialEnergy;

/* function declarations */
void Potential_LJ(int isEnergy, int isForce);

/* functions */
void Potential_LJ(int isEnergy, int isForce)
{
    int i, j, d;
    double distance;
    double tmpEnergy, tmpForce;

    if (isEnergy)
        totalPotentialEnergy = 0;
    for (i = 0; i < atomNumber; i++)
    {
        if (isEnergy)
            atoms[i].potentialEnergy = 0;
        if (isForce)
        {
            atoms[i].force[0] = 0;
            atoms[i].force[1] = 0;
            atoms[i].force[2] = 0;
        }
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            distance = atoms[i].neighbors[j].distance;
            if (distance > potentialCutoff_LJ)
            {
                continue;
            }
            if (isEnergy)
            {
                tmpEnergy = LJ_2_EPSILON * (pow((LJ_SIGMA / distance), 12.) - pow((LJ_SIGMA / distance), 6.));
                atoms[i].potentialEnergy += tmpEnergy;
                totalPotentialEnergy += tmpEnergy;
            }
            if (isForce)
            {
                tmpForce = LJ_2_SIGMA_6 / pow(distance, 14) - 1 / pow(distance, 8);
                for (d = 0; d < 3; d++)
                {
                    atoms[i].force[d] -= tmpForce * atoms[i].neighbors[j].dr[d];
                }
            }
        }
        if (isForce)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].force[d] *= LJ_24_EPSILON_SIGMA_6;
            }
        }
    }
}

/* main */
int main()
{
    /* parameters */
    neighborCutoff = 5.0;
    potentialCutoff_LJ = 5.0;

    atomNumber = 2;
    atoms[0].r[0] = 0;
    atoms[0].r[1] = 0;
    atoms[0].r[2] = 0;
    atoms[0].type = 1;
    atoms[0].id = 0;
    atoms[1].r[1] = 0;
    atoms[1].r[2] = 0;
    atoms[0].type = 1;
    atoms[1].id = 1;

    boxStartPoint[0] = 0;
    boxStartPoint[1] = 0;
    boxStartPoint[2] = 0;
    boxTranVecs[0][0] = 100;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = 1;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = 1;
    boxPerpendicular = 1;

    /* processing and output*/
    double r;

    printf("r(A) potentialEnergy(eV) force(eV/A)\n");
    for (r = 2.55; r < 5.0; r += 0.01)
    {
        atoms[1].r[0] = r;
        ConstructNeighborList();
        Potential_LJ(1, 1);
        printf("%f %f %f\n", r, totalPotentialEnergy, atoms[1].force[0]);
    }
    return 0;
}
