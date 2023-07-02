/* classes */
struct LatticePoint
{
    int reR[3];
    double r[3];
};

/* global variables */
double priTranVecs[3][3]; // primitive translation vectors

/* function declarations */
void ConstructLattice();

/* functions */
void ConstructLattice()
{
    int n, d;
    for (n = 0; n < latticePointNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            latticePoints[n].r[d] = latticePoints[n].reR[0] * priTranVecs[0][d] + latticePoints[n].reR[1] * priTranVecs[1][d] + latticePoints[n].reR[2] * priTranVecs[2][d];
        }
    }
}

/* main */
int main()
{
    /* parameter */
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 3;

    priTranVecs[0][0] = 1;
    priTranVecs[1][0] = 0;
    priTranVecs[2][0] = 0;
    priTranVecs[0][1] = 0;
    priTranVecs[1][1] = 1;
    priTranVecs[2][1] = 0;
    priTranVecs[0][2] = 0;
    priTranVecs[1][2] = 0;
    priTranVecs[2][2] = 1;

    /* processings */
    ConstructReducedLattice();
    ConstructLattice();

    /* output */
    int n;
    printf("x y z\n");
    for (n = 0; n < latticePointNumber; n++)
    {
        printf("%f %f %f\n", latticePoints[n].r[0], latticePoints[n].r[1], latticePoints[n].r[2]);
    }


    return 0;
}
