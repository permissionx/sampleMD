/* global variables */
int nStep;
int neighborInterval;

/* function declarations */
void NeighborList(int isForceConstruct);
void UpdateNeighborList();

/* functions */
void NeighborList(int isForceConstruct)
{
    if (isForceConstruct || nStep % neighborInterval == 0)
    {
        ConstructNeighborList();
    }
    else
    {
        UpdateNeighborList();
    }
}

void UpdateNeighborList()
{
    int i, j, d;
    int jAtomIndex;
    double dr[3];
    for (i = 0; i < atomNumber; i++)
    {
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            jAtomIndex = atoms[i].neighbors[j].index;
            PBC_dr(i, jAtomIndex, dr);
            atoms[i].neighbors[j].distance = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
            for (d = 0; d < 3; d++)
            {
                atoms[i].neighbors[j].dr[d] = dr[d];
            }
        }
    }
}

/* main */
int main()
{
    /* parameters */
    neighborInterval = 10000;
    neighborCutoff = 4.1 * 3.7;
    strcpy(potentialName, "LJ");

    /* processing and output*/
    double latticeConstant;

    printf("lattice_constant(A) pe(eV/atom)\n");
    nStep = 0;
    for (latticeConstant = 3.7; latticeConstant < 5.0; latticeConstant += 0.01)
    {
        ConstructStdCrystal_FCC(latticeConstant, 10);
        potentialCutoff_LJ = 4.1 * latticeConstant;
        NeighborList(0);
        Potential(1, 0);
        printf("%f %f\n", latticeConstant, totalPotentialEnergy / atomNumber);
        nStep += 1;
    }

    return 0;
}
