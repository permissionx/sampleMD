/* include */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* constants */
#define MAX_LATTICE_NUMBER 20000 //maximum number of lattices
#define MAX_ATOM_NUMBER 20000  //maximum number of atoms
#define MAX_CELL_ATOM_NUMBER 10 //maximum number of atoms in a cell

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

/* function declarations */
void ConstructReducedLattice();
void ConstructLattice();
void ConstructCrystal();

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
            latticePoints[n].r[d] = latticePoints[n].reR[0] * priTranVecs[0][d] + latticePoints[n].reR[1] * priTranVecs[1][d] + latticePoints[n].reR[2] * priTranVecs[2][d];
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

/*main*/
int main()
{
    /*parameters*/
    double latticeConstant = 3.14; //unit: angstrom
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 3;

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

    /*processing*/
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();

    /*output*/
    int n;
    printf("    id   type                x                y                z\n");
    for (n = 0; n < atomNumber; n++)
    {
        printf("%6d %6d %16.6f %16.6f %16.6f\n", atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2]);
    }

    return 0;
}

/* output
    id   type                x                y                z
     1      1         0.000000         0.000000         0.000000
     2      1        -1.570000         1.570000         1.570000
     3      1        -3.140000         3.140000         3.140000
     4      1         1.570000        -1.570000         1.570000
     5      1         0.000000         0.000000         3.140000
     6      1        -1.570000         1.570000         4.710000
     7      1         3.140000        -3.140000         3.140000
     8      1         1.570000        -1.570000         4.710000
     9      1         0.000000         0.000000         6.280000
    10      1         1.570000         1.570000        -1.570000
    11      1         0.000000         3.140000         0.000000
    12      1        -1.570000         4.710000         1.570000
    13      1         3.140000         0.000000         0.000000
    14      1         1.570000         1.570000         1.570000
    15      1         0.000000         3.140000         3.140000
    16      1         4.710000        -1.570000         1.570000
    17      1         3.140000         0.000000         3.140000
    18      1         1.570000         1.570000         4.710000
    19      1         3.140000         3.140000        -3.140000
    20      1         1.570000         4.710000        -1.570000
    21      1         0.000000         6.280000         0.000000
    22      1         4.710000         1.570000        -1.570000
    23      1         3.140000         3.140000         0.000000
    24      1         1.570000         4.710000         1.570000
    25      1         6.280000         0.000000         0.000000
    26      1         4.710000         1.570000         1.570000
    27      1         3.140000         3.140000         3.140000
*/