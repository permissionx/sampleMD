/* function declarations */
double ComputeBoxVolume();
void ComputeNonPBCForce(double nonPBCForce[3][3]);
void ComputeStress(double stress[6]);

/* functions */
double ComputeBoxVolume()
{
    double areaVector[3];
    VecCroMul(boxTranVecs[0], boxTranVecs[1], areaVector);
    return VecDotMul(areaVector, boxTranVecs[2]);
}


void ComputeNonPBCForce(double nonPBCForce[3][3])
{
    int n, i, j;
    for (i = 0; i < 3; i++)
    {
        boxTranVecs[i][i] *= 2;
        PBC_r();
        NeighborList(1);
        Potential(0, 1);
        for (j = 0; j < 3; j++)
        {
            nonPBCForce[i][j] = 0;
        }
        for (n = 0; n < atomNumber; n++)
        {
            if (atoms[n].r[i] >= boxTranVecs[i][i] / 4)
            {
                for (j = 0; j < 3; j++)
                {
                    nonPBCForce[i][j] += atoms[n].force[j];
                }
            }
        }
        boxTranVecs[i][i] /= 2;
    }
}

void ComputeStress(double stress[6])
{
    int n, d, i, j;
    double volume;
    double nonPBCForce[3][3], sumAtomForce[3][3];
    // only for orthogonal box with start point on (0,0,0)
    if (boxOrthogonal != 1)
    {
        printf("Error: computing stress in wrong box, check function ComputeStress()\n");
        exit(1);
    }

    // box virial
    ComputeNonPBCForce(nonPBCForce);
    PBC_r();
    NeighborList(1);
    Potential(0, 1);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            sumAtomForce[i][j] = 0;
        }
    }

    for (n = 0; n < atomNumber; n++)
    {
        for (i = 0; i < 3; i++)
        {
            if (atoms[n].r[i] >= boxTranVecs[i][i] / 2)
            {
                for (j = 0; j < 3; j++)
                {
                     sumAtomForce[i][j] += atoms[n].force[j];
                }
            }
        }
    }
    stress[0] = boxTranVecs[0][0] * (nonPBCForce[0][0] - sumAtomForce[0][0]);
    stress[1] = boxTranVecs[1][1] * (nonPBCForce[1][1] - sumAtomForce[1][1]);
    stress[2] = boxTranVecs[2][2] * (nonPBCForce[2][2] - sumAtomForce[2][2]);
    stress[3] = boxTranVecs[0][0] * (nonPBCForce[0][1] - sumAtomForce[0][1]);
    stress[4] = boxTranVecs[0][0] * (nonPBCForce[0][2] - sumAtomForce[0][2]);
    stress[5] = boxTranVecs[1][1] * (nonPBCForce[1][2] - sumAtomForce[1][2]);

    // add atom virial
    for (n = 0; n < atomNumber; n++)
    {
        stress[0] += atoms[n].r[0] * atoms[n].force[0] + atoms[n].velocity[0] * atoms[n].velocity[0] * typeMasses[atoms[n].type];
        stress[1] += atoms[n].r[1] * atoms[n].force[1] + atoms[n].velocity[1] * atoms[n].velocity[1] * typeMasses[atoms[n].type];
        stress[2] += atoms[n].r[2] * atoms[n].force[2] + atoms[n].velocity[2] * atoms[n].velocity[2] * typeMasses[atoms[n].type];
        stress[3] += atoms[n].r[0] * atoms[n].force[1] + atoms[n].velocity[0] * atoms[n].velocity[1] * typeMasses[atoms[n].type];
        stress[4] += atoms[n].r[0] * atoms[n].force[2] + atoms[n].velocity[0] * atoms[n].velocity[2] * typeMasses[atoms[n].type];
        stress[5] += atoms[n].r[1] * atoms[n].force[2] + atoms[n].velocity[1] * atoms[n].velocity[2] * typeMasses[atoms[n].type];
    }
    

    // stress
    volume = ComputeBoxVolume();
    for (d = 0; d < 6; d++)
    {
        stress[d] *= -160.21766208 / volume; // 1 eV/Angstrom3 = 160.21766208 GPa
    }
}
