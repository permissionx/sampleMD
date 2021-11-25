/* function declarations */
void DeleteAtomByIndex(int index);
void DeleteAtomsByShpereRegion(double center[3], double radius);
void DeleteAtomsByBlockRegion(double block[3][2]);

/* functions */
void DeleteAtomByIndex(int index)
{
    int i;
    for (i = index; i < atomNumber - 1; i++)
    {
        atoms[i] = atoms[i + 1];
    }
}

void DeleteAtomsByShpereRegion(double center[3], double radius)
{
    int n, d;
    double dr[3];
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            dr[d] = atoms[n].r[d] - center[d];
        }
        if (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] <= radius * radius)
        {
            DeleteAtomByIndex(n);
            n--;
            atomNumber--;
        }
    }
}

void DeleteAtomsByBlockRegion(double block[3][2])
{
    int n;
    for (n = 0; n < atomNumber; n++)
    {
        if (atoms[n].r[0] >= block[0][0] && atoms[n].r[0] <= block[0][1] &&
            atoms[n].r[1] >= block[1][0] && atoms[n].r[1] <= block[1][1] &&
            atoms[n].r[2] >= block[2][0] && atoms[n].r[2] <= block[2][1])
        {
            DeleteAtomByIndex(n);
            n--;
            atomNumber--;
        }
    }
}


/* main */
int main()
{
    /* parameters */
    double latticeConstant = 1.0;
    latticeSizes[0][0] = -3;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = -3;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = -3;
    latticeSizes[2][1] = 3;

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

    boxStartPoint[0] = -3;
    boxStartPoint[1] = -3;
    boxStartPoint[2] = -3;

    boxTranVecs[0][0] = latticeConstant * 6;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * 6;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * 6;
    boxVertical = 1;

    /* processing and output*/
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    
    Dump_xyz("origin.xyz");
    double center[3] = {0, 0, 0};
    DeleteAtomsByShpereRegion(center, 0.2);
    Dump_xyz("vacancy.xyz");

    return 0;
}