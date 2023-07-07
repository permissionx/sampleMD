/* main */
int main()
{
    int n;

    /* parameters */
    double randomSeed;
    randomSeed = 1.0;
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
    InsertAtom(siaPosition, 1);
    Minimize();
    for(n=0;n<atomNumber;n++)
    {
        if (atoms[n].id == 889)
        {
            DeleteAtomByIndex(n);
            atomNumber -= 1;
            break;
        }
    }
    Dump_lammpstrj("output/7.6_static-recovery.lammpstrj",1,0);
    Minimize();
    Dump_lammpstrj("output/7.6_static-recovery.lammpstrj",0,1);

    return 0;
}
