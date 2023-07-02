/* function declarations */
void Dump_lammpstrj(char fileName[20], int isNewFile, int dumpStep);

/* functions */
void Dump_lammpstrj(char fileName[20], int isNewFile, int dumpStep)
{
    int n;
    FILE *fp;
    if (boxOrthogonal != 1)
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
    fprintf(fp, "%d\n", dumpStep);
    fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fp, "%d\n", atomNumber);
    fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(fp, "%f %f\n", boxStartPoint[0], boxStartPoint[0] + boxTranVecs[0][0]);
    fprintf(fp, "%f %f\n", boxStartPoint[1], boxStartPoint[1] + boxTranVecs[1][1]);
    fprintf(fp, "%f %f\n", boxStartPoint[2], boxStartPoint[2] + boxTranVecs[2][2]);
    fprintf(fp, "ITEM: ATOMS id type x y z\n");
    for (n = 0; n < atomNumber; n++)
    {
        fprintf(fp, "%d %d %f %f %f\n", atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2]);
    }
    fclose(fp);
}

/* main */
int main()
{
    /* processing and output*/
    ConstructStdCrystal_BCC(1.0, 10);
    Dump_lammpstrj("bcc.lammpstrj",1,1);
    ConstructStdCrystal_FCC(1.0, 10);
    Dump_lammpstrj("fcc.lammpstrj",1,1);

    return 0;
}