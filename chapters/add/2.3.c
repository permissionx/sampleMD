/* function declarations */
void PBC_dr_general(int i, int j, double dr[3]);

/* functions */
void PBC_dr_general(int i, int j, double dr[3])
{
    int d;
    double reDr[3];
    for (d = 0; d < 3; d++)
    {
        reDr[d] = atoms[j].boxReR[d] - atoms[i].boxReR[d];
        if (reDr[d] < -0.5)
        {
            reDr[d] += 1;
        }
        else if (reDr[d] > 0.5)
        {
            reDr[d] -= 1;
        }
    }
    for (d = 0; d < 3; d++)
    {
        dr[d] = boxTranVecs[0][d] * reDr[0] + boxTranVecs[1][d] * reDr[1] + boxTranVecs[2][d] * reDr[2];
    }
}

/* main */
int main()
{
    /* parameters */
    double latticeConstant = 5; 
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 4;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 4;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 4;

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

    boxTranVecs[0][0] = latticeConstant * 4;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * 4;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * 4;

    /* processing */
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    ComputeAtomBoxReR();
    printf("Coordinate of atom0: (%f,%f,%f)\n", atoms[0].r[0], atoms[0].r[1], atoms[0].r[2]);
    printf("Coordinate of atom6: (%f,%f,%f)\n", atoms[6].r[0], atoms[6].r[1], atoms[6].r[2]);
    double dr[3];
    PBC_dr_general(0, 6, dr);
    printf("displacement from atom0 to atom1: (%f,%f,%f)\n", dr[0], dr[1], dr[2]);
}