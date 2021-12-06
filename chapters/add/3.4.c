/* functions */
void Dump_lammpstrj(char fileName[20], int isNewFile, int nstep)
{
    int n;
    FILE *fp;
    if (boxPerpendicular != 1)
    {
        printf("Error: Dump_lammpstrj() only works in cuboid.\n");
        exit(1);
    }
    if (isNewFile)
    {
        fp = fopen(fileName, "w");
    }
    else
    {
        fp = fopen(fileName, "a");
    }
    fprintf(fp, "ITEM: TIMESTEP\n");
    fprintf(fp, "%d\n", nstep);
    fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fp, "%d\n", atomNumber);
    fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(fp, "%f %f\n", boxStartPoint[0], boxStartPoint[0] + boxTranVecs[0][0]);
    fprintf(fp, "%f %f\n", boxStartPoint[1], boxStartPoint[1] + boxTranVecs[1][1]);
    fprintf(fp, "%f %f\n", boxStartPoint[2], boxStartPoint[2] + boxTranVecs[2][2]);
    fprintf(fp, "ITEM: ATOMS id type x y z pe fx fy fz\n");
    for (n = 0; n < atomNumber; n++)
    {
        fprintf(fp, "%d %d %f %f %f %f %f %f %f\n", 
        atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2],
        atoms[n].potentialEnergy,
        atoms[n].force[0], atoms[n].force[1], atoms[n].force[2]);
    }
    fclose(fp);
}

/* main */
int main()
{
    /* parameters */
    double latticeConstant = 4.23;
    boxPerpendicular = 1;
    neighborCutoff = latticeConstant * 4.1;
    potentialCutoff_LJ = neighborCutoff;

    /* processing*/
    ConstructStdCrystal_FCC(4.23, 6);
    double vacancyPosition[3] = {3*latticeConstant, 3*latticeConstant, 3*latticeConstant};
    DeleteAtomsByShpereRegion(vacancyPosition, 0.1);
    ConstructNeighborList();
    Potential_LJ(1, 1);

    /* output */
    Dump_lammpstrj("vacancy_in_Ne_FCC.lammpstrj", 1, 1);


    return 0;
}
