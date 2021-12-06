/* function declarations */
void EdgeDislocation_100(double latticeConstant);

/* functions */
void EdgeDislocation_100(double latticeConstant)
{
    int n;
    if (boxPerpendicular != 1)
    {
        printf("Error: EdgeDislocation_100() only works in cuboid.\n");
        exit(1);
    }
    double deleteBlock[3][2] = {
        {boxStartPoint[0] + boxTranVecs[0][0] - latticeConstant - 0.1, boxStartPoint[0] + boxTranVecs[0][0] + 0.1},
        {boxStartPoint[1] - 0.1, boxStartPoint[1] + boxTranVecs[1][1] + 0.1},
        {boxStartPoint[2] + boxTranVecs[2][2] / 2 - 0.1, boxStartPoint[2] + boxTranVecs[2][2] + 0.1}};

    DeleteAtomsByBlockRegion(deleteBlock);
    //shift atoms
    for (n = 0; n < atomNumber; n++)
    {
        if (atoms[n].r[2] > boxStartPoint[2] + boxTranVecs[2][2] / 2 - 0.1)
        {
            atoms[n].r[0] = boxStartPoint[0] + (atoms[n].r[0] - boxStartPoint[0]) * boxTranVecs[0][0] / (boxTranVecs[0][0] - latticeConstant);
        }
    }
}

/* main */
int main()
{
    /* parameters */
    double latticeConstant = 1; 
    latticeSizes[0][0] = -10;
    latticeSizes[0][1] = 10;
    latticeSizes[1][0] = -3;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = -5;
    latticeSizes[2][1] = 5;

    priTranVecs[0][0] = latticeConstant;
    priTranVecs[0][1] = 0;
    priTranVecs[0][2] = 0;
    priTranVecs[1][0] = 0;
    priTranVecs[1][1] = latticeConstant;
    priTranVecs[1][2] = 0;
    priTranVecs[2][0] = 0;
    priTranVecs[2][1] = 0;
    priTranVecs[2][2] = latticeConstant;

    cellAtomNumber = 1;
    cellAtomRs[0][0] = 0;
    cellAtomRs[0][1] = 0;
    cellAtomRs[0][2] = 0;
    cellAtomTypes[0] = 1;

    boxStartPoint[0] = -10 * latticeConstant;
    boxStartPoint[1] = -3 * latticeConstant;
    boxStartPoint[2] = -5 * latticeConstant;

    boxTranVecs[0][0] = latticeConstant * 20;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * 6;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * 10;
    boxPerpendicular = 1;

    /* processing*/
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();

    EdgeDislocation_100(latticeConstant);
    Dump_xyz("100_edge_dislocation.xyz");

    return 0;
}