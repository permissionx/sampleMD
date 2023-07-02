/* classes */
struct Atom
{
    int id, type;
    double r[3];
    double reR[3];
};

/* function declarations */
void ComputeAtomReR();

/* functions */
void ComputeAtomReR()
{
    int n, d;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].reR[d] = recPriTranVecs[d][0] * atoms[n].r[0] + recPriTranVecs[d][1] * atoms[n].r[1] + recPriTranVecs[d][2] * atoms[n].r[2];
        }
    }
}

/* main */
int main()
{
    /* parameters */
    double latticeConstant = 3.14; //unit: angstrom
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 3;

    priTranVecs[0][0] = 0.5 * latticeConstant;
    priTranVecs[0][1] = 0.5 * latticeConstant;
    priTranVecs[0][2] = -0.5 * latticeConstant;
    priTranVecs[1][0] = 0.5 * latticeConstant;
    priTranVecs[1][1] = -0.5 * latticeConstant;
    priTranVecs[1][2] = 0.5 * latticeConstant;
    priTranVecs[2][0] = -0.5 * latticeConstant;
    priTranVecs[2][1] = 0.5 * latticeConstant;
    priTranVecs[2][2] = 0.5 * latticeConstant;

    cellAtomNumber = 1;
    cellAtomRs[0][0] = 0;
    cellAtomRs[0][1] = 0;
    cellAtomRs[0][2] = 0;
    cellAtomTypes[0] = 1;

    /* processing */
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    ComputeRecTranVecs(priTranVecs, recPriTranVecs);
    ComputeAtomReR();
    

    /* output */
    int n, d;
    printf("rho_x rho_y rho_z\n");
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            printf("%f ", atoms[n].reR[d]);
        }
        printf("\n");
    }

    return 0;
}
