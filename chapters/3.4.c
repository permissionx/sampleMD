/* include */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* constants */
#define MAX_LATTICE_NUMBER 2000  // maximum number of lattices
#define MAX_ATOM_NUMBER 20000    // maximum number of atoms
#define MAX_CELL_ATOM_NUMBER 10  // maximum number of atoms in a cell
#define MAX_NEIGHBOR_NUMBER 2000 // maximum number of neighbors
// L-J parameters for Ne
#define LJ_EPSILON 0.0031
#define LJ_SIGMA 2.74
#define LJ_2_EPSILON 0.0062                     // 4 * LJ_EPSILON
#define LJ_24_EPSILON_SIGMA_6 31.48301472289983 // 24 * LJ_EPSILON * pow(LJ_SIGMA, 6)
#define LJ_2_SIGMA_6 846.3176000779524          // 2 * pow(LJ_SIGMA, 6)

/* classes */
struct LatticePoint
{
    int reR[3];
    double r[3];
};

struct AtomNeighbor
{
    int index;
    double distance;
    double dr[3]; // from atom to neighbor
};

struct Atom
{
    int id, type;
    double r[3];
    double reR[3];
    double boxReR[3];
    int neighborNumber;
    struct AtomNeighbor neighbors[MAX_NEIGHBOR_NUMBER];
    double force[3];
    double potentialEnergy;
};

/* global variables */
struct LatticePoint latticePoints[MAX_LATTICE_NUMBER];
int latticePointNumber;
int latticeSizes[3][2];
double priTranVecs[3][3]; // primitive translation vectors
struct Atom atoms[MAX_ATOM_NUMBER];
int atomNumber;
int cellAtomNumber;
double cellAtomRs[MAX_CELL_ATOM_NUMBER][3];
int cellAtomTypes[MAX_CELL_ATOM_NUMBER];
double recPriTranVecs[3][3];

double boxStartPoint[3];
double boxTranVecs[3][3];    // box translation vectors
double boxRecTranVecs[3][3]; // box reciprocal translation vectors
int boxPerpendicular;

double neighborCutoff;
double potentialCutoff_LJ;
double totalPotentialEnergy;

/* function declarations */
void ConstructReducedLattice();
void ConstructLattice();
void ConstructCrystal();
void Dump_xyz(char fileName[20]);
double VecDotMul(double vec1[3], double vec2[3]);
void VecCroMul(double vec1[3], double vec2[3], double vecOut[3]);
void ComputeRecTranVecs(double tranVecs[3][3], double recTranVecs[3][3]);
void ComputeAtomReR();

void ComputeAtomBoxReR(int n);
void PBC_r();
void PBC_dr();
void PBC_r_general();
void PBC_dr_general(int i, int j, double dr[3]);
void PBC_r_orthogonal();
void PBC_dr_orthogonal(int i, int j, double dr[3]);
void ConstructStdCrystal_BCC(double latticeConstant, int length);
void ConstructStdCrystal_FCC(double latticeConstant, int length);
void Dump_lammpstrj(char fileName[20], int isNewFile, int dumpStep);

void DeleteAtomByIndex(int index);
void DeleteAtomsByShpereRegion(double center[3], double radius);
void DeleteAtomsByBlockRegion(double block[3][2]);
void InsertAtom(double r[3], int type);
void EdgeDislocation_100(double latticeConstant);

void ConstructNeighborList();
void Potential_LJ(int isEnergy, int isForce);

/* functions */
void ConstructReducedLattice()
{
    int n, i, j, k;
    n = 0;
    for (i = latticeSizes[0][0]; i < latticeSizes[0][1]; i++)
    {
        for (j = latticeSizes[1][0]; j < latticeSizes[1][1]; j++)
        {
            for (k = latticeSizes[2][0]; k < latticeSizes[2][1]; k++)
            {
                latticePoints[n].reR[0] = i;
                latticePoints[n].reR[1] = j;
                latticePoints[n].reR[2] = k;
                n++;
            }
        }
    }
    latticePointNumber = n;
    if (latticePointNumber > MAX_LATTICE_NUMBER)
    {
        printf("Error: lattice point number exceeds the maximum number.\n");
        exit(1);
    }
}

void ConstructLattice()
{
    int n, d;
    for (n = 0; n < latticePointNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            latticePoints[n].r[d] = priTranVecs[0][d] * latticePoints[n].reR[0] + priTranVecs[1][d] * latticePoints[n].reR[1] + priTranVecs[2][d] * latticePoints[n].reR[2];
        }
    }
}

void ConstructCrystal()
{
    int nLattice, nAtom, nCellAtom;
    nAtom = 0;

    for (nLattice = 0; nLattice < latticePointNumber; nLattice++)
    {
        for (nCellAtom = 0; nCellAtom < cellAtomNumber; nCellAtom++)
        {
            atoms[nAtom].r[0] = latticePoints[nLattice].r[0] + cellAtomRs[nCellAtom][0];
            atoms[nAtom].r[1] = latticePoints[nLattice].r[1] + cellAtomRs[nCellAtom][1];
            atoms[nAtom].r[2] = latticePoints[nLattice].r[2] + cellAtomRs[nCellAtom][2];

            atoms[nAtom].type = cellAtomTypes[nCellAtom];
            atoms[nAtom].id = nAtom + 1;
            nAtom++;
        }
    }
    atomNumber = nAtom;
    if (atomNumber > MAX_ATOM_NUMBER)
    {
        printf("Error: atom number exceeds the maximum number of atoms.\n");
        exit(1);
    }
}

void Dump_xyz(char fileName[20])
{
    int n;
    FILE *fp;
    fp = fopen(fileName, "w");
    fprintf(fp, "%d\n", atomNumber);
    fprintf(fp, "id type x y z\n");
    for (n = 0; n < atomNumber; n++)
    {
        fprintf(fp, "%d %d %f %f %f\n", atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2]);
    }
    fclose(fp);
}

double VecDotMul(double vec1[3], double vec2[3])
{
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

void VecCroMul(double vec1[3], double vec2[3], double vecOut[3])
{
    vecOut[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vecOut[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vecOut[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

void ComputeRecTranVecs(double tranVecs[3][3], double recTranVecs[3][3])
{
    int i, d;
    double cellVol;
    double tmpVec[3]; // temperary vectors

    VecCroMul(tranVecs[0], tranVecs[1], tmpVec);
    cellVol = VecDotMul(tmpVec, tranVecs[2]);
    for (i = 0; i < 3; i++)
    {
        VecCroMul(tranVecs[(i + 1) % 3], tranVecs[(i + 2) % 3], recTranVecs[i]);
        for (d = 0; d < 3; d++)
        {
            recTranVecs[i][d] /= cellVol; // 2pi factor ignored
        }
    }
}

void ComputeAtomReR()
{
    int n, d;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].reR[d] = recPriTranVecs[d][0] * atoms[n].r[0] + recPriTranVecs[d][1] * atoms[n].r[1] + recPriTranVecs[d][2] * atoms[n].r[2];
        }
    }
}

void ComputeAtomBoxReR(int n)
{
    int d;
    for (d = 0; d < 3; d++)
    {
        atoms[n].boxReR[d] = boxRecTranVecs[d][0] * (atoms[n].r[0] - boxStartPoint[0]) +
                             boxRecTranVecs[d][1] * (atoms[n].r[1] - boxStartPoint[1]) +
                             boxRecTranVecs[d][2] * (atoms[n].r[2] - boxStartPoint[2]);
        if (atoms[n].boxReR[d] < -1 || atoms[n].boxReR[d] >= 2)
        {
            printf("Lost atom %d!\n", atoms[n].id);
            exit(1);
        }
    }
}

void PBC_r_general()
{
    int n, d;
    int isDisplaced;
    for (n = 0; n < atomNumber; n++)
    {
        isDisplaced = 0;
        for (d = 0; d < 3; d++)
        {
            if (atoms[n].boxReR[d] < 0)
            {
                atoms[n].boxReR[d] += 1;
                isDisplaced = 1;
            }
            else if (atoms[n].boxReR[d] >= 1)
            {
                atoms[n].boxReR[d] -= 1;
                isDisplaced = 1;
            }
        }
        if (isDisplaced)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[n].r[d] = boxTranVecs[0][d] * atoms[n].boxReR[0] + boxTranVecs[1][d] * atoms[n].boxReR[1] + boxTranVecs[2][d] * atoms[n].boxReR[2] + boxStartPoint[d];
            }
        }
    }
}

void PBC_dr_general(int i, int j, double dr[3])
{
    int d;
    double reDr[3];
    for (d = 0; d < 3; d++)
    {
        reDr[d] = atoms[j].boxReR[d] - atoms[i].boxReR[d];
        if (reDr[d] < -0.5)
        {
            reDr[d] += 1;
        }
        else if (reDr[d] > 0.5)
        {
            reDr[d] -= 1;
        }
    }
    for (d = 0; d < 3; d++)
    {
        dr[d] = boxTranVecs[0][d] * reDr[0] + boxTranVecs[1][d] * reDr[1] + boxTranVecs[2][d] * reDr[2];
    }
}

void PBC_r_orthogonal()
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

void PBC_dr_orthogonal(int i, int j, double dr[3])
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
    int n;
    if (boxPerpendicular == 1)
    {
        PBC_r_orthogonal();
    }
    else
    {
        for (n = 0; n < atomNumber; n++)
        {
            ComputeAtomBoxReR(n);
        }
        PBC_r_general();
    }
}

void PBC_dr(int i, int j, double dr[3])
{
    if (boxPerpendicular == 1)
    {
        PBC_dr_orthogonal(i, j, dr);
    }
    else
    {
        ComputeAtomBoxReR(i);
        ComputeAtomBoxReR(j);
        PBC_dr_general(i, j, dr);
    }
}

void ConstructStdCrystal_BCC(double latticeConstant, int length)
{
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = length;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = length;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = length;

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

    boxTranVecs[0][0] = latticeConstant * length;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * length;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * length;
    boxPerpendicular = 1;

    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
}

void ConstructStdCrystal_FCC(double latticeConstant, int length)
{
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = length;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = length;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = length;

    priTranVecs[0][0] = latticeConstant;
    priTranVecs[0][1] = 0;
    priTranVecs[0][2] = 0;
    priTranVecs[1][0] = 0;
    priTranVecs[1][1] = latticeConstant;
    priTranVecs[1][2] = 0;
    priTranVecs[2][0] = 0;
    priTranVecs[2][1] = 0;
    priTranVecs[2][2] = latticeConstant;

    cellAtomNumber = 4;
    cellAtomRs[0][0] = 0.0;
    cellAtomRs[0][1] = 0.0;
    cellAtomRs[0][2] = 0.0;
    cellAtomRs[1][0] = 0.5 * latticeConstant;
    cellAtomRs[1][1] = 0.5 * latticeConstant;
    cellAtomRs[1][2] = 0.0;
    cellAtomRs[2][0] = 0.5 * latticeConstant;
    cellAtomRs[2][1] = 0.0;
    cellAtomRs[2][2] = 0.5 * latticeConstant;
    cellAtomRs[3][0] = 0.0;
    cellAtomRs[3][1] = 0.5 * latticeConstant;
    cellAtomRs[3][2] = 0.5 * latticeConstant;
    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 1;
    cellAtomTypes[2] = 1;
    cellAtomTypes[3] = 1;

    boxStartPoint[0] = 0;
    boxStartPoint[1] = 0;
    boxStartPoint[2] = 0;

    boxTranVecs[0][0] = latticeConstant * length;
    boxTranVecs[0][1] = 0;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * length;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * length;
    boxPerpendicular = 1;

    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
}

void Dump_lammpstrj(char fileName[20], int isNewFile, int dumpStep)
{
    int n;
    FILE *fp;
    if (boxPerpendicular != 1)
    {
        printf("Error: Dump_lammpstrj() only works in cuboid.\n");
        exit(1);
    }
    if (isNewFile)
    {
        fp = fopen(fileName, "w");
    }
    else
    {
        fp = fopen(fileName, "a");
    }
    fprintf(fp, "ITEM: TIMESTEP\n");
    fprintf(fp, "%d\n", dumpStep);
    fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fp, "%d\n", atomNumber);
    fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(fp, "%f %f\n", boxStartPoint[0], boxStartPoint[0] + boxTranVecs[0][0]);
    fprintf(fp, "%f %f\n", boxStartPoint[1], boxStartPoint[1] + boxTranVecs[1][1]);
    fprintf(fp, "%f %f\n", boxStartPoint[2], boxStartPoint[2] + boxTranVecs[2][2]);
    fprintf(fp, "ITEM: ATOMS id type x y z pe fx fy fz\n");
    for (n = 0; n < atomNumber; n++)
    {
        fprintf(fp, "%d %d %f %f %f %f %f %f %f\n", 
        atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2],
        atoms[n].potentialEnergy,
        atoms[n].force[0], atoms[n].force[1], atoms[n].force[2]);
    }
    fclose(fp);
}

void DeleteAtomByIndex(int index)
{
    int i;
    for (i = index; i < atomNumber - 1; i++)
    {
        atoms[i] = atoms[i + 1];
    }
}

void DeleteAtomsByShpereRegion(double center[3], double radius)
{
    int n, d;
    double dr[3];
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            dr[d] = atoms[n].r[d] - center[d];
        }
        if (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] <= radius * radius)
        {
            DeleteAtomByIndex(n);
            n--;
            atomNumber--;
        }
    }
}

void DeleteAtomsByBlockRegion(double block[3][2])
{
    int n;
    for (n = 0; n < atomNumber; n++)
    {
        if (atoms[n].r[0] >= block[0][0] && atoms[n].r[0] <= block[0][1] &&
            atoms[n].r[1] >= block[1][0] && atoms[n].r[1] <= block[1][1] &&
            atoms[n].r[2] >= block[2][0] && atoms[n].r[2] <= block[2][1])
        {
            DeleteAtomByIndex(n);
            n--;
            atomNumber--;
        }
    }
}

void InsertAtom(double r[3], int type)
{
    atoms[atomNumber].id = atoms[atomNumber - 1].id + 1;
    atoms[atomNumber].r[0] = r[0];
    atoms[atomNumber].r[1] = r[1];
    atoms[atomNumber].r[2] = r[2];
    atoms[atomNumber].type = type;
    atomNumber++;
}

void EdgeDislocation_100(double latticeConstant)
{
    int n;
    if (boxPerpendicular != 1)
    {
        printf("Error: EdgeDislocation_100() only works in cuboid.\n");
        exit(1);
    }
    double deleteBlock[3][2] = {
        {boxStartPoint[0] + boxTranVecs[0][0] - latticeConstant - 0.1, boxStartPoint[0] + boxTranVecs[0][0] + 0.1},
        {boxStartPoint[1] - 0.1, boxStartPoint[1] + boxTranVecs[1][1] + 0.1},
        {boxStartPoint[2] + boxTranVecs[2][2] / 2 - 0.1, boxStartPoint[2] + boxTranVecs[2][2] + 0.1}};

    DeleteAtomsByBlockRegion(deleteBlock);
    // shift atoms
    for (n = 0; n < atomNumber; n++)
    {
        if (atoms[n].r[2] > boxStartPoint[2] + boxTranVecs[2][2] / 2 - 0.1)
        {
            atoms[n].r[0] = boxStartPoint[0] + (atoms[n].r[0] - boxStartPoint[0]) * boxTranVecs[0][0] / (boxTranVecs[0][0] - latticeConstant);
        }
    }
}

void ConstructNeighborList()
{
    int i, j;
    int d;
    double distance;
    double dr[3];
    int neighborIndex_i, neighborIndex_j;
    int tmpNeighborNumber;

    for (i = 0; i < atomNumber; i++)
    {
        atoms[i].neighborNumber = 0;
    }
    for (i = 0; i < atomNumber; i++)
    {
        for (j = i + 1; j < atomNumber; j++)
        {
            PBC_dr(i, j, dr);
            distance = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
            if (distance < neighborCutoff)
            {
                neighborIndex_i = atoms[i].neighborNumber;
                neighborIndex_j = atoms[j].neighborNumber;
                atoms[i].neighbors[neighborIndex_i].index = j;
                atoms[j].neighbors[neighborIndex_j].index = i;
                atoms[i].neighbors[neighborIndex_i].distance = distance;
                atoms[j].neighbors[neighborIndex_j].distance = distance;
                for (d = 0; d < 3; d++)
                {
                    atoms[i].neighbors[neighborIndex_i].dr[d] = dr[d];
                    atoms[j].neighbors[neighborIndex_j].dr[d] = -dr[d];
                }
                tmpNeighborNumber = atoms[i].neighborNumber++;
                if (tmpNeighborNumber > MAX_NEIGHBOR_NUMBER)
                {
                    printf("Error: too many neighbors.\n");
                    exit(1);
                }
                tmpNeighborNumber = atoms[j].neighborNumber++;
                if (tmpNeighborNumber > MAX_NEIGHBOR_NUMBER)
                {
                    printf("Error: too many neighbors.\n");
                    exit(1);
                }
            }
        }
    }
}

void Potential_LJ(int isEnergy, int isForce)
{
    int i, j, d;
    double distance;
    double tmpEnergy, tmpForce;

    if (isEnergy)
        totalPotentialEnergy = 0;
    for (i = 0; i < atomNumber; i++)
    {
        if (isEnergy)
            atoms[i].potentialEnergy = 0;
        if (isForce)
        {
            atoms[i].force[0] = 0;
            atoms[i].force[1] = 0;
            atoms[i].force[2] = 0;
        }
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            distance = atoms[i].neighbors[j].distance;
            if (distance > potentialCutoff_LJ)
            {
                continue;
            }
            if (isEnergy)
            {
                tmpEnergy = LJ_2_EPSILON * (pow((LJ_SIGMA / distance), 12.) - pow((LJ_SIGMA / distance), 6.));
                atoms[i].potentialEnergy += tmpEnergy;
                totalPotentialEnergy += tmpEnergy;
            }
            if (isForce)
            {
                tmpForce = LJ_2_SIGMA_6 / pow(distance, 14) - 1 / pow(distance, 8);
                for (d = 0; d < 3; d++)
                {
                    atoms[i].force[d] -= tmpForce * atoms[i].neighbors[j].dr[d];
                }
            }
        }
        if (isForce)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].force[d] *= LJ_24_EPSILON_SIGMA_6;
            }
        }
    }
}

/* main */
int main()
{
    /* parameters */
    double latticeConstant = 4.23;
    boxPerpendicular = 1;
    neighborCutoff = latticeConstant * 4.1;
    potentialCutoff_LJ = neighborCutoff;

    /* processing*/
    ConstructStdCrystal_FCC(4.23, 6);
    double vacancyPosition[3] = {0, 0, 0};
    DeleteAtomsByShpereRegion(vacancyPosition, 0.1);
    ConstructNeighborList();
    Potential_LJ(1, 1);

    /* output */
    Dump_lammpstrj("vacancy_Ne_FCC.lammpstrj", 1, 1);


    return 0;
}
