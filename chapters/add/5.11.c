// stress 

/* function declarations */
void VirialPair(int n, double virial[6]);
void VirialPair_LJ(int n, double virial[6]);
double ComputeBoxVolume();
void ComputeStress(double stress[6]);

/* functions */
void VirialPair(int n, double virial[6])
{
    // result unit: eV
    NeighborList(0);
    if (strcmp(potentialName, "LJ") == 0)
    {
        VirialPair_LJ(n, virial);
    }
    else
    {
        printf("Error: Calculating Virial pair for potential %s is not supported.\n", potentialName);
        exit(1);
    }
}

void VirialPair_LJ(int n, double virial[6])
{
    int l, d;
    double sigma6;
    double distance, distance8, distance14;
    double x;
    double force[3];
    for (d = 0; d < 6; d++)
    {
        virial[d] = 0;
    }
    for (l = 0; l < atoms[n].neighborNumber; l++)
    {
        distance = atoms[n].neighbors[l].distance;
        distance8 = pow(distance, 8);
        distance14 = pow(distance, 14);
        for (d = 0; d < 3; d++)
        {
            x = atoms[n].neighbors[l].dr[d];
            force[d] = (LJ_2_SIGMA_6 * x / distance14 - x / distance8);
        }
        for (d = 0; d < 6; d++)
        {
            virial[0] += atoms[n].neighbors[l].dr[0] * force[0]; // xx
            virial[1] += atoms[n].neighbors[l].dr[1] * force[1]; // yy
            virial[2] += atoms[n].neighbors[l].dr[2] * force[2]; // zz
            virial[3] += atoms[n].neighbors[l].dr[0] * force[1]; // xy
            virial[4] += atoms[n].neighbors[l].dr[0] * force[2]; // xz
            virial[5] += atoms[n].neighbors[l].dr[1] * force[2]; // yz
        }
    }
    for (d = 0; d < 6; d++)
    {
        virial[d] *= -0.5 * LJ_24_EPSILON_SIGMA_6;
    }
}

double ComputeBoxVolume()
{
    double areaVector[3];
    VecCroMul(boxTranVecs[0], boxTranVecs[1], areaVector);
    return VecDotMul(areaVector, boxTranVecs[2]);
}

void ComputeStress(double stress[6])
{
    //result unit: GPa
    int n, l, d;
    double ke[6], virial[6];
    double volume;
    for (d = 0; d < 6; d++)
    {
        stress[d] = 0;
    }
    for (n = 0; n < atomNumber; n++)
    {
        ke[0] = -typeMasses[atoms[n].type] * atoms[n].velocity[0] * atoms[n].velocity[0]; // xx
        ke[1] = -typeMasses[atoms[n].type] * atoms[n].velocity[1] * atoms[n].velocity[1]; // yy
        ke[2] = -typeMasses[atoms[n].type] * atoms[n].velocity[2] * atoms[n].velocity[2]; // zz
        ke[3] = -typeMasses[atoms[n].type] * atoms[n].velocity[0] * atoms[n].velocity[1]; // xy
        ke[4] = -typeMasses[atoms[n].type] * atoms[n].velocity[0] * atoms[n].velocity[2]; // xz
        ke[5] = -typeMasses[atoms[n].type] * atoms[n].velocity[1] * atoms[n].velocity[2]; // yz
        VirialPair(n, virial);
        for (d = 0; d < 6; d++)
        {
            stress[d] += ke[d] + virial[d];
        }
    }
    volume = ComputeBoxVolume();
    for (d = 0; d < 6; d++)
    {
        stress[d] *= 160.21766208/volume;  // 1 eV/Angstrom3 = 160.21766208 GPa
    }
}