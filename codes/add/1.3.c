/* main */
int main()
{
    /*parameters*/
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 3;

    priTranVecs[0][0] = 0.5;
    priTranVecs[0][1] = 0.5;
    priTranVecs[0][2] = -0.5;
    priTranVecs[1][0] = 0.5;
    priTranVecs[1][1] = -0.5;
    priTranVecs[1][2] = 0.5;
    priTranVecs[2][0] = -0.5;
    priTranVecs[2][1] = 0.5;
    priTranVecs[2][2] = 0.5;

    /*processing*/
    ConstructReducedLattice();
    ConstructLattice();

    /*output*/
    int n;
    printf("x y z\n");
    for (n = 0; n < latticePointNumber; n++)
    {
        printf("%f %f %f\n", latticePoints[n].r[0], latticePoints[n].r[1], latticePoints[n].r[2]);
    }

    return 0;
}
