/*include*/
#include <stdio.h>

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
    double reR_box[3];
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
double boxStartPoint[3];
double boxTranVecs[3][3]; // box translation vectors
double boxRecTranVecs[3][3]; //box reciprocal translation vectors
double boxRecTranVecs_inv[3][3];


/*function declaration*/
void ConstructReducedLattice();
void Mul_3_1(double a[3][3], double b[3], double p[3]); // matrix multiplication 3x1
void ConstructLattice();
void ConstructCrystal();
void DumpSingle(char fileName[20]);
double VecDotMul(double vec0[3], double vec1[3]); // vector dot multiplication
void VecCroMul(double vec0[3], double vec1[3], double vecOut[3]); // vector cross multiplication
void ComputeRecTranVecs(tranVecs[3][3], recTranVecs[3][3]);
void ComputeReR(recTranVecs[3][3]);
void MatInv(double matIn[3][3], double matOut[3][3]);    //Matrix Inversion
void ComputeReR(recTranVecs_inv[3][3], boxStartPoint[3]);
void ComputeRecTranVecs();
void ComputeAtomBoxReR();


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


void Mul_3_1(double a[3][3], double b[3], double p[3])
{
    int i;
    for (i = 0; i < 3; i++)
        p[i] = a[0][i] * b[0] + a[1][i] * b[1] + a[2][i] * b[2];
}


void ConstructLattice()
{
    int n;
    for (n = 0; n < latticePointNumber; n++)
        Mul_3_1(priTranVecs, latticePoints[n].reR, latticePoints[n].r,);
}


void ConstructCrystal()
{
    int nLattice, nAtom, nCellAtom;
    nAtom = 0;

    for (nLattice = 0; nLattice < latticePointNumber; nLattice++)
    {
        for ( nCellAtom = 0; nCellAtom < cellAtomNumber; nCellAtom++)
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


void DumpSingle(char fileName[20])
{
    int i;
    FILE *file;
    file = fopen(fileName, "w");
    fprintf(file, "%d\n", atomNumber);
    fprintf(file, "id type x y z\n");
    for (i = 0; i < atomNumber; i++)
    {
        fprintf(file, "%d %d %f %f %f\n", atoms[i].id, atoms[i].type, atoms[i].r[0], atoms[i].r[1], atoms[i].r[2]);
    }
    fclose(file);
}

double VecDotMul(double vec0[3], double vec1[3])
{
    result = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
    return result;
}

void VecCroMul(double vec0[3], double vec1[3], double vecOut[3])
{
    vecOut[0] = vec0[1] * vec1[2] - vec0[2] * vec1[1];
    vecOut[1] = vec0[2] * vec1[0] - vec0[0] * vec1[2];
    vecOut[2] = vec0[0] * vec1[1] - vec0[1] * vec1[0];
}

void ComputeRecTranVecs(tranVecs[3][3], recTranVecs[3][3])
{
    int i, d; // d: direction
    double cellVol;
    double tmpVec[3]; // temperary vectors

    VecCroMul(tmpVec, tranVecs[0], tranVecs[1]);
    cellVol = VecDotMul(tmpVec, tranVecs[2]);
    for (i = 0; i < 3; i++)
    {
        VecCroMul(recTranVecs[i], tranVecs[(i + 1) % 3], tranVecs[(i + 2) % 3]);
        for (d = 0; d < 3; d++)
        {
            recTranVecs[i][d] /= cellVol;  // 2pi factor ignored
        }
    }
}

void MatInv(double matIn[3][3], double matOut[3][3])
{
    int i, j;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            matOut[i][j] = matIn[j][i];
        }
    }
}

void ComputeRecTranVecs()
{
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    MatInv(boxRecTranVecs, boxRecTranVecs_inv);
}

void ComputeReR(r[3], recTranVecs_inv[3][3], startPoint[3], reR[3])
{
    double relativeR[3];
    int i;
    for (i = 0; i < 3; i++)
    {
        relativeR[i] = r[i] - startPoint[i];
    }
    Mul_3_1(recTranVecs_inv, relativeR, reR);
}

void ComputeAtomBoxReR()
{
    int n;
    for (n=0;n<atomNumber;n++)
    {
        ComputeReR(atoms[n].r, boxRecTranVecs_inv, boxStartPoint, atoms[n].reR);
    }
}



/*main*/
int main() //
{
    double latticeConstant = 3.165; //unit: Angstrom
    /*parameter*/
    latticeSizes[0][0] = 0;  latticeSizes[0][1] = 10;
    latticeSizes[1][0] = 0;  latticeSizes[1][1] = 10;
    latticeSizes[2][0] = 0;  latticeSizes[2][1] = 10;

    priTranVecs[0][0] = latticeConstant; priTranVecs[1][0] = 0; priTranVecs[2][0] = 0;
    priTranVecs[0][1] = 0; priTranVecs[1][1] = latticeConstant; priTranVecs[2][1] = 0;
    priTranVecs[0][2] = 0; priTranVecs[1][2] = 0; priTranVecs[2][2] = latticeConstant;

    cellAtomNumber = 2;
    cellAtomRs[0][0] = 0; cellAtomRs[0][1] = 0; cellAtomRs[0][2] = 0;
    cellAtomRs[1][0] = 0.5 * latticeConstant; cellAtomRs[1][1] = 0.5 * latticeConstant; cellAtomRs[1][2] = 0.5 * latticeConstant;
    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 1;

    startPoint[0] = 0; startPoint[1] = 0; startPoint[2] = 0;
    boxTranVecs[0][0] = latticeConstant*10; boxTranVecs[1][0] = 0; boxTranVecs[2][0] = 0;
    boxTranVecs[0][1] = 0; boxTranVecs[1][1] = latticeConstant*10; boxTranVecs[2][1] = 0;
    boxTranVecs[0][2] = 0; boxTranVecs[1][2] = 0; boxTranVecs[2][2] = latticeConstant*10;




    /*process*/
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    DumpSingle("W_BCC.xyz");


}
