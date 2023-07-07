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
    neighborCutoff = 6;
    neighborInterval = 100;

    /* processing*/
    double latticeConstant;
    double stress[6];
    int d;
    printf("lc str_xx str_yy str_zz str_xy str_xz str_yz potential\n");
    for (latticeConstant =3.04; latticeConstant < 3.241; latticeConstant += 0.01)
    {
        ConstructStdCrystal_BCC(latticeConstant, 10);
        InitVelocity(0);
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
