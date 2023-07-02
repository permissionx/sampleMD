/* global variables */
double energyTolerance_Minimize;
char minimizeStyle[20];

/* function declarations */
void Minimize();
void MinDirection();

/* functions */
void Minimize()
{
    double potentialEnergy_last;
    nStep = 0;
    printf("\n---Minimization start---\n");
    printf("iter pe dE\n");
    NeighborList(1);
    Potential(1, 1);
    printf("%d %20.10f\n", nStep, totalPotentialEnergy);
    do
    {
        potentialEnergy_last = totalPotentialEnergy;
        MinDirection();
        LineMinimize();
        Potential(0, 1);
        nStep += 1;
        printf("%d %30.20f  %30.20f\n", nStep, totalPotentialEnergy, totalPotentialEnergy - potentialEnergy_last);
    } while (fabs(potentialEnergy_last - totalPotentialEnergy) > energyTolerance_Minimize);
    printf("---Minimization end---\n");
}

void MinDirection()
{
    if (strcmp(minimizeStyle, "SD") == 0)
        MinDirection_SD();
    else if (strcmp(minimizeStyle, "CG") == 0)
        MinDirection_CG();
    else
        printf("\nMinimize style not found. Program terminated.\n"), exit(-1);
}
