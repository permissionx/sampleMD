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
    strcpy(dynamicStyle, "VelocityVerlet");

    /* processing*/
    
    double strain;
    int n;
    double stress[6];

    printf("strain_xx sigma_xx sigma_yy\n");
    for (strain=-0.001; strain<0.00101; strain+= 0.0001)
    {
        ConstructStdCrystal_BCC(3.14, 10);
        InitVelocity(0);
        for (n=0;n<atomNumber;n++)
        {
            atoms[n].r[0] *= 1+strain;
        }
        boxTranVecs[0][0] *= 1+strain;
        ComputeStress(stress);
        printf("%f%% %f %f\n", strain*100, stress[0], stress[1]);
    }

    return 0;
}
