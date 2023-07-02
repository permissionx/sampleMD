/*constant*/
#define MAX_LATTICE_NUMBER 10000
#define MAX_ATOM_NUMBER 100000
#define MAX_CELL_ATOM_NUMBER 10


/*class definition*/
struct LatticePoint
{
    int id;
    double reR[3], r[3];
};

struct Atom
{
    double r[3]; 
    int id, type;
};


/*global variable declaration*/
struct LatticePoint latticePoints[MAX_LATTICE_NUMBER];
int latticePointNumber;
int latticeSizes[3][2];
double priTranVecs[3][3]; // primitive translation vectors
struct Atom atoms[MAX_ATOM_NUMBER];
int atomNumber;
int cellAtomNumber; 
double cellAtomRs[MAX_CELL_ATOM_NUMBER][3];
int cellAtomTypes[MAX_CELL_ATOM_NUMBER];


/*function declaration*/
void ConstructReducedLattice();
void Mul_3_1(double p[3], double a[3][3], double b[3]); // matrix multiplication
void ConstructLattice();
void ConstructCrystal();


/*function*/
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


void Mul_3_1(double p[3], double a[3][3], double b[3])
{
    int i;
    for (i = 0; i < 3; i++)
        p[i] = a[0][i] * b[0] + a[1][i] * b[1] + a[2][i] * b[2];
}


void ConstructLattice()
{
    int n;
    for (n = 0; n < latticePointNumber; n++)
        Mul_3_1(latticePoints[n].r, priTranVecs, latticePoints[n].reR);
}

void ConstructCrystal()
{
    int nLattice, nAtom, nCellAtom;
    nAtom = 0;

    for (nLattice = 0;nLattice < latticePointNumber; nLattice++)
    {
        for ( nCellAtom = 0; nCellAtom < cellAtomNumber; nCellAtom++)
        {
            atoms[nAtom].r[0] = latticePoints[nLattice].r[0] + cellAtomRs[nCellAtom][0];
            atoms[nAtom].r[1] = latticePoints[nLattice].r[1] + cellAtomRs[nCellAtom][1];
            atoms[nAtom].r[2] = latticePoints[nLattice].r[2] + cellAtomRs[nCellAtom][2];

            atoms[nAtom].type = cellAtomTypes[nCellAtom];
            atoms[nAtom].id = nCellAtom+1;
            nAtom++;
        }
    }
    atomNumber = nAtom;
}

/*main*/
int main() //
{
    double latticeConstant = 3.165; //unit: Angstrom
    /*parameter*/
    latticeSizes[0][0] = 0;  latticeSizes[0][1] = 10;
    latticeSizes[1][0] = 0;  latticeSizes[1][1] = 10;
    latticeSizes[2][0] = 0;  latticeSizes[2][1] = 10;
    
    priTranVecs[0][0] = 0.5*latticeConstant; priTranVecs[1][0] = 0.5*latticeConstant; priTranVecs[2][0] = -0.5*latticeConstant;
    priTranVecs[0][1] = 0.5*latticeConstant; priTranVecs[1][1] = -0.5*latticeConstant; priTranVecs[2][1] = 0.5*latticeConstant;
    priTranVecs[0][2] = -0.5*latticeConstant; priTranVecs[1][2] = 0.5*latticeConstant; priTranVecs[2][2] = 0.5*latticeConstant;
    
    cellAtomNumber = 1; 
    cellAtomRs[0][0] = 0; cellAtomRs[0][1] = 0; cellAtomRs[0][2] = 0;
    cellAtomTypes[0] = 1;


    /*process*/
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    PBC_r();
}
