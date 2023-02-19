/* main */
int main()
{
    /* parameters */
    double randomSeed;
    randomSeed = 1.0;
    srand(randomSeed);

    typeMasses[1] = 20.1797; // for Ne
    InitMassUnit();
    strcpy(potentialName, "LJ");
    potentialCutoff_LJ = 20;
    neighborCutoff = 20;
    neighborInterval = 100;

    /* processing*/
    double latticeConstant;
    double stress[6];
    int d;
    printf("lc str_xx str_yy str_zz str_xy str_xz str_yz potential\n");
    for (latticeConstant = 4.15; latticeConstant < 4.35; latticeConstant += 0.01)
    {
        ConstructStdCrystal_FCC(latticeConstant, 5);
        printf("%f ", latticeConstant);
        ComputeStress(stress);
        for (d = 0; d < 6; d++)
        {
            printf("%f ", stress[d]);
        }
        Potential(1, 0);
        printf("%f\n", totalPotentialEnergy/atomNumber);
        nStep += 1;
    }
    return 0;
}
