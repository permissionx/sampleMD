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
    double vacancyCenter[3] = {5*3.14,5*3.14,5*3.14};
    ConstructStdCrystal_BCC(3.14, 10);
    NeighborList(1);
    Potential(1,0);
    printf("Perfect system energy: %f\n", totalPotentialEnergy);
    DeleteAtomsByShpereRegion(vacancyCenter, 0.1);
    Minimize();
    printf("System energy with a vacancy: %f\n", totalPotentialEnergy);

    return 0;
}
