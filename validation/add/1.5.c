/* main */
int main()
{
    /* parameters */
    double latticeConstant = 5.642;    //unit: angstrom
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 3;

    priTranVecs[0][0] = 0;
    priTranVecs[0][1] = 0.5 * latticeConstant;
    priTranVecs[0][2] = 0.5 * latticeConstant;
    priTranVecs[1][0] = 0.5 * latticeConstant;
    priTranVecs[1][1] = 0;
    priTranVecs[1][2] = 0.5 * latticeConstant;
    priTranVecs[2][0] = 0.5 * latticeConstant;
    priTranVecs[2][1] = 0.5 * latticeConstant;
    priTranVecs[2][2] = 0;

    cellAtomNumber = 2;
    cellAtomRs[0][0] = 0;
    cellAtomRs[0][1] = 0;
    cellAtomRs[0][2] = 0;
    cellAtomRs[1][0] = 0.5 * latticeConstant;
    cellAtomRs[1][1] = 0.5 * latticeConstant;
    cellAtomRs[1][2] = 0.5 * latticeConstant;
    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 2;

    /* processing */
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();

    /* output */
    int n;
    printf("id type x y z\n");
    for (n = 0; n < atomNumber; n++)
    {
        printf("%d %d %f %f %f\n", atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2]);
    }

    return 0;
}
