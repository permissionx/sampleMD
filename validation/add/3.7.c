/* main */
int main()
{
    /* parameters */
    neighborCutoff = 6.0;
    strcpy(potentialName, "EAM");

    /* processing*/
    ConstructStdCrystal_BCC(3.14, 10);
    double deleteBlock[3][2] = {{-0.1, 10.1 * 3.14}, {-0.1, 10.1 * 3.14}, {5.1 * 3.14, 10.1 * 3.14}};
    DeleteAtomsByBlockRegion(deleteBlock);
    ConstructNeighborList();
    Potential(1, 1);

    /* output */
    Dump_lammpstrj("output/3.7_surface-W-BCC.lammpstrj", 1, 1);
    
    return 0;
}
