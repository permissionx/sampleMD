/* main */
int main()
{
    /* processing and output*/
    double latticeConstant;

    printf("lattice_constant(A) pe(eV)\n");
    for (latticeConstant = 3.7; latticeConstant < 5.0; latticeConstant += 0.01)
    {
        ConstructStdCrystal_FCC(latticeConstant, 7);
        neighborCutoff = 4.1 * latticeConstant;
        potentialCutoff_LJ = neighborCutoff;
        ConstructNeighborList();
        Potential_LJ(1, 0);
        printf("%f %f\n", latticeConstant, totalPotentialEnergy / atomNumber);
    }

    return 0;
}