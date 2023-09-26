/* main */
int main()
{
    /* parameters */
    neighborCutoff = 6.0;
    strcpy(potentialName, "EAM");
    strcpy(minimizeStyle, "SD");
    neighborInterval = 100;
    energyTolerance_Minimize = 1E-9;

    char minimizeDumpFileName[20];
    strcpy(minimizeDumpFileName, "output/4.7_min-sia-SD.lammpstrj");

    /* processing & output*/
    ConstructStdCrystal_BCC(3.14, 10);
    double r[3] = {1, 1, 1};
    InsertAtom(r, 1);
    NeighborList(1);
    Potential(0, 1);
    Dump_lammpstrj(minimizeDumpFileName, 1, 1);
    Minimize();
    Dump_lammpstrj(minimizeDumpFileName, 0, 2);

    return 0;
}


