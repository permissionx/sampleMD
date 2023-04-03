/* include */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* constants */
#define MAX_LATTICE_NUMBER 2000 // maximum number of lattices
#define MAX_ATOM_NUMBER 20000   // maximum number of atoms
#define MAX_CELL_ATOM_NUMBER 10 // maximum number of atoms in a cell

/* classes */
struct LatticePoint
{
    int reR[3];
    double r[3];
};

struct Atom
{
    int id, type;
    double r[3];
    double reR[3];
    double boxReR[3];
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

/* function declarations */
void ConstructReducedLattice();
void ConstructLattice();
void ConstructCrystal();
void Dump_xyz(char fileName[20]);
double VecDotMul(double vec1[3], double vec2[3]);
void VecCroMul(double vec1[3], double vec2[3], double vecOut[3]);
void ComputeRecTranVecs(double tranVecs[3][3], double recTranVecs[3][3]);
void ComputeAtomBoxReR(int n);

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
        if (atoms[n].boxReR[d] < -1 || atoms[n].boxReR[d] >= 2)
        {
            printf("Lost atom %d!\n", atoms[n].id);
            exit(1);
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

/* main */
int main()
{
    int n, d;
    /* parameters */
    double latticeConstant = 3.14; // unit: angstrom
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 4;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 4;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 4;

    priTranVecs[0][0] = 0.5 * latticeConstant;
    priTranVecs[0][1] = 0.5 * latticeConstant;
    priTranVecs[0][2] = -0.5 * latticeConstant;
    priTranVecs[1][0] = 0.5 * latticeConstant;
    priTranVecs[1][1] = -0.5 * latticeConstant;
    priTranVecs[1][2] = 0.5 * latticeConstant;
    priTranVecs[2][0] = -0.5 * latticeConstant;
    priTranVecs[2][1] = 0.5 * latticeConstant;
    priTranVecs[2][2] = 0.5 * latticeConstant;

    cellAtomNumber = 1;
    cellAtomRs[0][0] = 0;
    cellAtomRs[0][1] = 0;
    cellAtomRs[0][2] = 0;
    cellAtomTypes[0] = 1;

    boxStartPoint[0] = 0;
    boxStartPoint[1] = 0;
    boxStartPoint[2] = 0;

    boxTranVecs[0][0] = 0.5 * latticeConstant * 4;
    boxTranVecs[0][1] = 0.5 * latticeConstant * 4;
    boxTranVecs[0][2] = -0.5 * latticeConstant * 4;
    boxTranVecs[1][0] = 0.5 * latticeConstant * 4;
    boxTranVecs[1][1] = -0.5 * latticeConstant * 4;
    boxTranVecs[1][2] = 0.5 * latticeConstant * 4;
    boxTranVecs[2][0] = -0.5 * latticeConstant * 4;
    boxTranVecs[2][1] = 0.5 * latticeConstant * 4;
    boxTranVecs[2][2] = 0.5 * latticeConstant * 4;

    /* processing */
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    for (n = 0; n < atomNumber; n++)
    {
        ComputeAtomBoxReR(n);
    }

    /* output */
    printf("             rex              rey              rez\n");
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            printf("%16.6f ", atoms[n].boxReR[d]);
        }
        printf("\n");
    }
}

/* output
             rex              rey              rez
        0.000000         0.000000         0.000000
        0.000000         0.000000         0.250000
        0.000000         0.000000         0.500000
        0.000000         0.000000         0.750000
        0.000000         0.250000         0.000000
        0.000000         0.250000         0.250000
        0.000000         0.250000         0.500000
        0.000000         0.250000         0.750000
        0.000000         0.500000         0.000000
        0.000000         0.500000         0.250000
        0.000000         0.500000         0.500000
        0.000000         0.500000         0.750000
        0.000000         0.750000         0.000000
        0.000000         0.750000         0.250000
        0.000000         0.750000         0.500000
        0.000000         0.750000         0.750000
        0.250000         0.000000         0.000000
        0.250000         0.000000         0.250000
        0.250000         0.000000         0.500000
        0.250000         0.000000         0.750000
        0.250000         0.250000         0.000000
        0.250000         0.250000         0.250000
        0.250000         0.250000         0.500000
        0.250000         0.250000         0.750000
        0.250000         0.500000         0.000000
        0.250000         0.500000         0.250000
        0.250000         0.500000         0.500000
        0.250000         0.500000         0.750000
        0.250000         0.750000         0.000000
        0.250000         0.750000         0.250000
        0.250000         0.750000         0.500000
        0.250000         0.750000         0.750000
        0.500000         0.000000         0.000000
        0.500000         0.000000         0.250000
        0.500000         0.000000         0.500000
        0.500000         0.000000         0.750000
        0.500000         0.250000         0.000000
        0.500000         0.250000         0.250000
        0.500000         0.250000         0.500000
        0.500000         0.250000         0.750000
        0.500000         0.500000         0.000000
        0.500000         0.500000         0.250000
        0.500000         0.500000         0.500000
        0.500000         0.500000         0.750000
        0.500000         0.750000         0.000000
        0.500000         0.750000         0.250000
        0.500000         0.750000         0.500000
        0.500000         0.750000         0.750000
        0.750000         0.000000         0.000000
        0.750000         0.000000         0.250000
        0.750000         0.000000         0.500000
        0.750000         0.000000         0.750000
        0.750000         0.250000         0.000000
        0.750000         0.250000         0.250000
        0.750000         0.250000         0.500000
        0.750000         0.250000         0.750000
        0.750000         0.500000         0.000000
        0.750000         0.500000         0.250000
        0.750000         0.500000         0.500000
        0.750000         0.500000         0.750000
        0.750000         0.750000         0.000000
        0.750000         0.750000         0.250000
        0.750000         0.750000         0.500000
        0.750000         0.750000         0.750000
*/