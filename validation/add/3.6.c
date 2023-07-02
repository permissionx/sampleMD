/* global variables */
char potentialName[20];

/* function declarations */
void Potential(int isEnergy, int isForce);

/* functions */
void Potential(int isEnergy, int isForce)
{
    if (strcmp(potentialName, "LJ") == 0)

    {
        Potential_LJ(isEnergy, isForce);
    }
    else if (strcmp(potentialName, "EAM") == 0)
    {
        Potential_EAM(isEnergy, isForce);
    }
    else
    {
        printf("Error: Potential %s not found.\n", potentialName);
        exit(1);
    }
}

/* main */
int main()
{
    /* parameters */
    neighborCutoff = 6.0;
    strcpy(potentialName, "EAM");
    
    /* processing and output*/
    double latticeConstant;
    printf("lattice_constant pe(eV/atom)\n");
    for (latticeConstant = 2.9; latticeConstant < 3.4; latticeConstant += 0.01)
    {
        ConstructStdCrystal_BCC(latticeConstant, 7);
        ConstructNeighborList();
        Potential(1, 0);
        printf("%f %f\n", latticeConstant, totalPotentialEnergy / atomNumber);
    }

    return 0;
}
