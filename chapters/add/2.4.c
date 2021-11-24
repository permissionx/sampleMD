
/* global variables */
int boxVertical;

/* function declarations */
void PBC_r();
void PBC_dr();
void PBC_r_vertical();
void PBC_dr_vertical(int i, int j, double dr[3]);

/* functions */
void PBC_r_vertical()
{
    int n, d;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            if (atoms[n].r[d] < boxStartPoint[d])
            {
                atoms[n].r[d] += boxTranVecs[d][d];
            }
            else if (atoms[n].r[d] >= boxStartPoint[d] + boxTranVecs[d][d])
            {
                atoms[n].r[d] -= boxTranVecs[d][d];
            }
        }
    }
}

void PBC_dr_vertical(int i, int j, double dr[3])
{
    int d;
    for (d = 0; d < 3; d++)
    {
        dr[d] = atoms[j].r[d] - atoms[i].r[d];
        if (dr[d] < -boxTranVecs[d][d] / 2)
        {
            dr[d] += boxTranVecs[d][d];
        }
        else if (dr[d] > boxTranVecs[d][d] / 2)
        {
            dr[d] -= boxTranVecs[d][d];
        }
    }
}

void PBC_r()
{
    if (boxVertical == 1)
    {
        PBC_r_vertical();
    }
    else
    {
        ComputeAtomBoxReR();
        PBC_r_general();
    }
}

void PBC_dr(int i, int j, double dr[3])
{
    if (boxVertical == 1)
    {
        PBC_dr_vertical(i, j, dr);
    }
    else
    {
        ComputeAtomBoxReR();
        PBC_dr_general(i, j, dr);
    }
}

/* main */
int main()
{
    /* parameters */
    double latticeConstant = 5; 
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 4;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 4;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 4;

    priTranVecs[0][0] = latticeConstant;
    priTranVecs[0][1] = 0;
    priTranVecs[0][2] = 0;
    priTranVecs[1][0] = 0;
    priTranVecs[1][1] = latticeConstant;
    priTranVecs[1][2] = 0;
    priTranVecs[2][0] = 0;
    priTranVecs[2][1] = 0;
    priTranVecs[2][2] = latticeConstant;

    cellAtomNumber = 2;
    cellAtomRs[0][0] = 0;
    cellAtomRs[0][1] = 0;
    cellAtomRs[0][2] = 0;
    cellAtomRs[1][0] = 0.5 * latticeConstant;
    cellAtomRs[1][1] = 0.5 * latticeConstant;
    cellAtomRs[1][2] = 0.5 * latticeConstant;
    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 1;

    boxStartPoint[0] = 0;
    boxStartPoint[1] = 0;
    boxStartPoint[2] = 0;

    boxTranVecs[0][0] = latticeConstant * 4;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * 4;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * 4;
    boxVertical = 1;

    /* processing */
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    printf("Coordinate of atom0: (%f,%f,%f)\n", atoms[0].r[0], atoms[0].r[1], atoms[0].r[2]);
    printf("Coordinate of atom6: (%f,%f,%f)\n", atoms[6].r[0], atoms[6].r[1], atoms[6].r[2]);
    double dr[3];
    PBC_dr(0, 6, dr);
    printf("displacement from atom0 to atom1: (%f,%f,%f)\n", dr[0], dr[1], dr[2]);
}