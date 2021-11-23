/* function declarations */
void Dump_xyz(char fileName[20]);

/* functions */
void Dump_xyz(char fileName[20])
{
    int n;
    FILE *fp;
    fp = fopen(fileName, "w");
    fprintf(fp, "%d\n", atomNumber);
    fprintf(fp, "id type x y z\n");
    for (n = 0; n < atomNumber; n++)
    {
        fprintf(fp, "%d %d %f %f %f\n", atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2]);
    }
    fclose(fp);
}

/*main*/
int main()
{
    /*parameters*/
    double latticeConstant = 5.642; //unit: angstrom
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = 0;
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

    cellAtomNumber = 8;
    cellAtomRs[0][0] = 0;
    cellAtomRs[0][1] = 0;
    cellAtomRs[0][2] = 0;
    cellAtomRs[1][0] = 0;
    cellAtomRs[1][1] = 0.5 * latticeConstant;
    cellAtomRs[1][2] = 0.5 * latticeConstant;
    cellAtomRs[2][0] = 0.5 * latticeConstant;
    cellAtomRs[2][1] = 0;
    cellAtomRs[2][2] = 0.5 * latticeConstant;
    cellAtomRs[3][0] = 0.5 * latticeConstant;
    cellAtomRs[3][1] = 0.5 * latticeConstant;
    cellAtomRs[3][2] = 0;
    cellAtomRs[4][0] = 0 + 0.5 * latticeConstant;
    cellAtomRs[4][1] = 0;
    cellAtomRs[4][2] = 0;
    cellAtomRs[5][0] = 0 + 0.5 * latticeConstant;
    cellAtomRs[5][1] = 0.5 * latticeConstant;
    cellAtomRs[5][2] = 0.5 * latticeConstant;
    cellAtomRs[6][0] = 0.5 * latticeConstant + 0.5 * latticeConstant;
    cellAtomRs[6][1] = 0;
    cellAtomRs[6][2] = 0.5 * latticeConstant;
    cellAtomRs[7][0] = 0.5 * latticeConstant + 0.5 * latticeConstant;
    cellAtomRs[7][1] = 0.5 * latticeConstant;
    cellAtomRs[7][2] = 0;

    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 1;
    cellAtomTypes[2] = 1;
    cellAtomTypes[3] = 1;
    cellAtomTypes[4] = 2;
    cellAtomTypes[5] = 2;
    cellAtomTypes[6] = 2;
    cellAtomTypes[7] = 2;

    /*processing*/
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();

    /*output*/
    Dump_xyz("NaCl.xyz");
}
