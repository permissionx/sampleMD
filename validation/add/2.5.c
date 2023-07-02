/* function declarations */
void ConstructStdCrystal_BCC(double latticeConstant, int length);
void ConstructStdCrystal_FCC(double latticeConstant, int length);

/* functions */
void ConstructStdCrystal_BCC(double latticeConstant, int length)
{
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = length;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = length;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = length;

    priTranVecs[0][0] = latticeConstant;
    priTranVecs[0][1] = 0;
    priTranVecs[0][2] = 0;
    priTranVecs[1][0] = 0;
    priTranVecs[1][1] = latticeConstant;
    priTranVecs[1][2] = 0;
    priTranVecs[2][0] = 0;
    priTranVecs[2][1] = 0;
    priTranVecs[2][2] = latticeConstant;

    cellAtomNumber = 2;
    cellAtomRs[0][0] = 0;
    cellAtomRs[0][1] = 0;
    cellAtomRs[0][2] = 0;
    cellAtomRs[1][0] = 0.5 * latticeConstant;
    cellAtomRs[1][1] = 0.5 * latticeConstant;
    cellAtomRs[1][2] = 0.5 * latticeConstant;
    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 1;

    boxStartPoint[0] = 0;
    boxStartPoint[1] = 0;
    boxStartPoint[2] = 0;

    boxTranVecs[0][0] = latticeConstant * length;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * length;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * length;
    boxOrthogonal = 1;

    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
}

void ConstructStdCrystal_FCC(double latticeConstant, int length)
{
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = length;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = length;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = length;

    priTranVecs[0][0] = latticeConstant;
    priTranVecs[0][1] = 0;
    priTranVecs[0][2] = 0;
    priTranVecs[1][0] = 0;
    priTranVecs[1][1] = latticeConstant;
    priTranVecs[1][2] = 0;
    priTranVecs[2][0] = 0;
    priTranVecs[2][1] = 0;
    priTranVecs[2][2] = latticeConstant;

    cellAtomNumber = 4;
    cellAtomRs[0][0] = 0.0;
    cellAtomRs[0][1] = 0.0;
    cellAtomRs[0][2] = 0.0;
    cellAtomRs[1][0] = 0.5 * latticeConstant;
    cellAtomRs[1][1] = 0.5 * latticeConstant;
    cellAtomRs[1][2] = 0.0;
    cellAtomRs[2][0] = 0.5 * latticeConstant;
    cellAtomRs[2][1] = 0.0;
    cellAtomRs[2][2] = 0.5 * latticeConstant;
    cellAtomRs[3][0] = 0.0;
    cellAtomRs[3][1] = 0.5 * latticeConstant;
    cellAtomRs[3][2] = 0.5 * latticeConstant;
    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 1;
    cellAtomTypes[2] = 1;
    cellAtomTypes[3] = 1;

    boxStartPoint[0] = 0;
    boxStartPoint[1] = 0;
    boxStartPoint[2] = 0;

    boxTranVecs[0][0] = latticeConstant * length;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * length;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * length;
    boxOrthogonal = 1;

    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
}
