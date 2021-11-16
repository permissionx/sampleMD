/* include */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* constants */
#define MAX_LATTICE_NUMBER 2000 //maximum number of lattices
#define MAX_ATOM_NUMBER 100000  //maximum number of atoms
#define MAX_CELL_ATOM_NUMBER 10 //maximum number of atoms in a cell

/* classes */
struct LatticePoint
{
    int id;
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
                latticePoints[n].id = n;
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
    double latticeConstant = 5.642; //unit: angstrom
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = 3;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = 3;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = 3;

    priTranVecs[0][0] = 0;
    priTranVecs[0][1] = 0.5 * latticeConstant;
    priTranVecs[0][2] = 0.5 * latticeConstant;
    priTranVecs[1][0] = 0.5 * latticeConstant;
    priTranVecs[1][1] = 0;
    priTranVecs[1][2] = 0.5 * latticeConstant;
    priTranVecs[2][0] = 0.5 * latticeConstant;
    priTranVecs[2][1] = 0.5 * latticeConstant;
    priTranVecs[2][2] = 0;

    cellAtomNumber = 2;
    cellAtomRs[0][0] = 0;
    cellAtomRs[0][1] = 0;
    cellAtomRs[0][2] = 0;
    cellAtomRs[1][0] = 0.5 * latticeConstant;
    cellAtomRs[1][1] = 0.5 * latticeConstant;
    cellAtomRs[1][2] = 0.5 * latticeConstant;
    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 2;

    /*processing*/
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();

    /*output*/
    int n;
    for (n = 0; n < atomNumber; n++)
    {
        printf("%d %d %f %f %f\n", atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2]);
    }

    return 0;
}

/* output
1 1 0.000000 0.000000 0.000000
2 2 2.821000 2.821000 2.821000
3 1 2.821000 2.821000 0.000000
4 2 5.642000 5.642000 2.821000
5 1 5.642000 5.642000 0.000000
6 2 8.463000 8.463000 2.821000
7 1 2.821000 0.000000 2.821000
8 2 5.642000 2.821000 5.642000
9 1 5.642000 2.821000 2.821000
10 2 8.463000 5.642000 5.642000
11 1 8.463000 5.642000 2.821000
12 2 11.284000 8.463000 5.642000
13 1 5.642000 0.000000 5.642000
14 2 8.463000 2.821000 8.463000
15 1 8.463000 2.821000 5.642000
16 2 11.284000 5.642000 8.463000
17 1 11.284000 5.642000 5.642000
18 2 14.105000 8.463000 8.463000
19 1 0.000000 2.821000 2.821000
20 2 2.821000 5.642000 5.642000
21 1 2.821000 5.642000 2.821000
22 2 5.642000 8.463000 5.642000
23 1 5.642000 8.463000 2.821000
24 2 8.463000 11.284000 5.642000
25 1 2.821000 2.821000 5.642000
26 2 5.642000 5.642000 8.463000
27 1 5.642000 5.642000 5.642000
28 2 8.463000 8.463000 8.463000
29 1 8.463000 8.463000 5.642000
30 2 11.284000 11.284000 8.463000
31 1 5.642000 2.821000 8.463000
32 2 8.463000 5.642000 11.284000
33 1 8.463000 5.642000 8.463000
34 2 11.284000 8.463000 11.284000
35 1 11.284000 8.463000 8.463000
36 2 14.105000 11.284000 11.284000
37 1 0.000000 5.642000 5.642000
38 2 2.821000 8.463000 8.463000
39 1 2.821000 8.463000 5.642000
40 2 5.642000 11.284000 8.463000
41 1 5.642000 11.284000 5.642000
42 2 8.463000 14.105000 8.463000
43 1 2.821000 5.642000 8.463000
44 2 5.642000 8.463000 11.284000
45 1 5.642000 8.463000 8.463000
46 2 8.463000 11.284000 11.284000
47 1 8.463000 11.284000 8.463000
48 2 11.284000 14.105000 11.284000
49 1 5.642000 5.642000 11.284000
50 2 8.463000 8.463000 14.105000
51 1 8.463000 8.463000 11.284000
52 2 11.284000 11.284000 14.105000
53 1 11.284000 11.284000 11.284000
54 2 14.105000 14.105000 14.105000
*/