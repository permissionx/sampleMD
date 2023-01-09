/* global variables */
double energyTolerance_Minimize;

/* function declarations */
void Minimize();

/* functions */
void Minimize()
{
    double potentialEnergy_begin;
    nStep = 0;
    printf("\n---Minimization start---\n");
    printf("iter pe dE\n");
    NeighborList();
    Potential(1, 0);
    do
    {
        potentialEnergy_begin = totalPotentialEnergy;
        MinDirection(nStep);
        LineMinimize();
        NeighborList();
        Potential(1, 0);
        nStep += 1;
        printf("%d %f %f %f\n", nStep, potentialEnergy_begin, totalPotentialEnergy, fabs(potentialEnergy_begin - totalPotentialEnergy));
    } while (fabs(potentialEnergy_begin - totalPotentialEnergy) > energyTolerance_Minimize);
    printf("---Minimization end---\n");
}