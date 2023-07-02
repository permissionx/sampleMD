/* constants */
#define MAX_LATTICE_NUMBER 20000 //maximum number of lattices

/* classes */
struct LatticePoint
{
    int reR[3];
};

/* global variables */
struct LatticePoint latticePoints[MAX_LATTICE_NUMBER];
int latticePointNumber;
int latticeSizes[3][2];

/* function declarations */
void ConstructReducedLattice();

/* functions */
void ConstructReducedLattice()
{
    int n, i, j, k;
    n = 0;
    for (i = latticeSizes[0][0]; i < latticeSizes[0][1]; i++)
    {
        for (j = latticeSizes[1][0]; j < latticeSizes[1][1]; j++)
        {
            for (k = latticeSizes[2][0]; k < latticeSizes[2][1]; k++)
            {
                latticePoints[n].reR[0] = i;
                latticePoints[n].reR[1] = j;
                latticePoints[n].reR[2] = k;
                n++;
            }
        }
    }
    latticePointNumber = n;
    if (latticePointNumber > MAX_LATTICE_NUMBER)
    {
        printf("Error: lattice point number exceeds the maximum number.\n");
        exit(1);
    }
}

/* main */
int main()
{
    /* parameters */
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 3;

    /* processing */
    ConstructReducedLattice();

    /* output */
    int n;
    printf("rho_x rho_y rho_z\n");
    for (n = 0; n < latticePointNumber; n++)
    {
        printf("%d %d %d\n", latticePoints[n].reR[0], latticePoints[n].reR[1], latticePoints[n].reR[2]);
    }

    return 0;
}
