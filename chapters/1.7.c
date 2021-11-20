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
void Dump(char fileName[20]);

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


void Dump(char fileName[20])
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

    /*processing*/
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();

    /*output*/
    Dump("W_BCC.xyz");
}

/* output in W_BCC.xyz 
54
id type x y z
1 1 0.000000 0.000000 0.000000
2 1 1.570000 1.570000 1.570000
3 1 0.000000 0.000000 3.140000
4 1 1.570000 1.570000 4.710000
5 1 0.000000 0.000000 6.280000
6 1 1.570000 1.570000 7.850000
7 1 0.000000 3.140000 0.000000
8 1 1.570000 4.710000 1.570000
9 1 0.000000 3.140000 3.140000
10 1 1.570000 4.710000 4.710000
11 1 0.000000 3.140000 6.280000
12 1 1.570000 4.710000 7.850000
13 1 0.000000 6.280000 0.000000
14 1 1.570000 7.850000 1.570000
15 1 0.000000 6.280000 3.140000
16 1 1.570000 7.850000 4.710000
17 1 0.000000 6.280000 6.280000
18 1 1.570000 7.850000 7.850000
19 1 3.140000 0.000000 0.000000
20 1 4.710000 1.570000 1.570000
21 1 3.140000 0.000000 3.140000
22 1 4.710000 1.570000 4.710000
23 1 3.140000 0.000000 6.280000
24 1 4.710000 1.570000 7.850000
25 1 3.140000 3.140000 0.000000
26 1 4.710000 4.710000 1.570000
27 1 3.140000 3.140000 3.140000
28 1 4.710000 4.710000 4.710000
29 1 3.140000 3.140000 6.280000
30 1 4.710000 4.710000 7.850000
31 1 3.140000 6.280000 0.000000
32 1 4.710000 7.850000 1.570000
33 1 3.140000 6.280000 3.140000
34 1 4.710000 7.850000 4.710000
35 1 3.140000 6.280000 6.280000
36 1 4.710000 7.850000 7.850000
37 1 6.280000 0.000000 0.000000
38 1 7.850000 1.570000 1.570000
39 1 6.280000 0.000000 3.140000
40 1 7.850000 1.570000 4.710000
41 1 6.280000 0.000000 6.280000
42 1 7.850000 1.570000 7.850000
43 1 6.280000 3.140000 0.000000
44 1 7.850000 4.710000 1.570000
45 1 6.280000 3.140000 3.140000
46 1 7.850000 4.710000 4.710000
47 1 6.280000 3.140000 6.280000
48 1 7.850000 4.710000 7.850000
49 1 6.280000 6.280000 0.000000
50 1 7.850000 7.850000 1.570000
51 1 6.280000 6.280000 3.140000
52 1 7.850000 7.850000 4.710000
53 1 6.280000 6.280000 6.280000
54 1 7.850000 7.850000 7.850000
*/