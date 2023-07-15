/* constants */
#define MAX_ATOM_NUMBER 20000  //maximum number of atoms
#define MAX_CELL_ATOM_NUMBER 10 //maximum number of atoms in a cell

/* classes */
struct Atom
{
    int id, type;
    double r[3];
};

/* global variables */
struct Atom atoms[MAX_ATOM_NUMBER];
int atomNumber;
int cellAtomNumber;
double cellAtomRs[MAX_CELL_ATOM_NUMBER][3];
int cellAtomTypes[MAX_CELL_ATOM_NUMBER];

/* function declarations */
void ConstructCrystal();

/* functions */
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
    for (n = 0; n < atomNumber; n++)
    {
        printf("%d %d %f %f %f\n", atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2]);
    }

    return 0;
}

/* output
1 1 0.000000 0.000000 0.000000
2 1 -1.570000 1.570000 1.570000
3 1 -3.140000 3.140000 3.140000
4 1 1.570000 -1.570000 1.570000
5 1 0.000000 0.000000 3.140000
6 1 -1.570000 1.570000 4.710000
7 1 3.140000 -3.140000 3.140000
8 1 1.570000 -1.570000 4.710000
9 1 0.000000 0.000000 6.280000
10 1 1.570000 1.570000 -1.570000
11 1 0.000000 3.140000 0.000000
12 1 -1.570000 4.710000 1.570000
13 1 3.140000 0.000000 0.000000
14 1 1.570000 1.570000 1.570000
15 1 0.000000 3.140000 3.140000
16 1 4.710000 -1.570000 1.570000
17 1 3.140000 0.000000 3.140000
18 1 1.570000 1.570000 4.710000
19 1 3.140000 3.140000 -3.140000
20 1 1.570000 4.710000 -1.570000
21 1 0.000000 6.280000 0.000000
22 1 4.710000 1.570000 -1.570000
23 1 3.140000 3.140000 0.000000
24 1 1.570000 4.710000 1.570000
25 1 6.280000 0.000000 0.000000
26 1 4.710000 1.570000 1.570000
27 1 3.140000 3.140000 3.140000
*/