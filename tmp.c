/*class definition*/
struct LatticePoint
{
    int id;
    double reR[3], r[3];
};

/*global variable declaration*/
double boxStartPoint[3];
double boxTranVecs[3][3]; // box translation vectors
double boxRecTranVecs[3][3]; //box reciprocal translation vectors
double boxRecTranVecs_inv[3][3];

/*function declaration*/
void ComputeRecTranVecs();
void ComputeAtomBoxReR();

/*function*/
void ComputeRecTranVecs()
{
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    MatInv(boxRecTranVecs, boxRecTranVecs_inv);
}

void ComputeAtomBoxReR()
{
    int n;
    for (n=0;n<atomNumber;n++)
    {
        ComputeReR(atoms[n].r, boxRecTranVecs_inv, boxStartPoint, atoms[n].reR_box);
    }
}

void PBC_r(atom)
{
    
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
    ComputeRecTranVecs();
    ComputeAtomBoxReR();
}