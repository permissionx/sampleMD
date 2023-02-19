/* main */
int main()
{
    /* parameters */
    double randomSeed;
    randomSeed = 1.0;
    srand(randomSeed);

    typeMasses[1] = 183.85;
    InitMassUnit();
    strcpy(potentialName, "EAM");
    neighborCutoff = 6;
    neighborInterval = 100;
    strcpy(dynamicStyle, "VelocityVerlet");

    /* processing*/
    ConstructStdCrystal_BCC(3.14, 10);
    InitVelocity(300);
    Dynamics(1, 0.001);
  
    return 0;
}
