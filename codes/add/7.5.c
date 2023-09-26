/* main */
int main()
{
    /* parameters */
    unsigned int randomSeed;
    randomSeed = 1;
    srand(randomSeed);
    typeMasses[1] = 183.84; // for W
    InitMassUnit();
    strcpy(potentialName, "EAM");
    strcpy(minimizeStyle, "CG");
    neighborCutoff = 6;
    neighborInterval = 100;

    /* processing*/
    double siaPosition[3] = {5.25*3.14,5.25*3.14,5.25*3.14};
    ConstructStdCrystal_BCC(3.14, 10);
    NeighborList(1);
    Potential(1,0);
    printf("Perfect system energy: %f\n", totalPotentialEnergy);
    InsertAtom(siaPosition, 1);
    Minimize();
    printf("System energy with a SIA: %f\n", totalPotentialEnergy);
    Dump_lammpstrj("output/7.5_sia.dump",1,0);
    return 0;
}
