/*constant*/
#define MAX_LATTICE_NUMBER 10000


/*class definition*/
struct LatticePoint
{
    int id;
    double reR[3], r[3];
};


/*global variable declaration*/
struct LatticePoint latticePoints[MAX_LATTICE_NUMBER];
int latticePointNumber;
int latticeSizes[3][2];
double priTranVecs[3][3]; // primitive translation vectors



/*function declaration*/
void ConstructReducedLattice();
void Mul_3_1(double p[3], double a[3][3], double b[3]); // matrix multiplication
void ConstructLattice();


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



/*main*/
//1.3
int main()
{
    /*parameter*/
    latticeSizes[0][0] = 0;  latticeSizes[0][1] = 10;
    latticeSizes[1][0] = 0;  latticeSizes[1][1] = 10;
    latticeSizes[2][0] = 0;  latticeSizes[2][1] = 10;
    priTranVecs[0][0] = 0.5; priTranVecs[1][0] = 0.5; priTranVecs[2][0] = -0.5;
    priTranVecs[0][1] = 0.5; priTranVecs[1][1] = -0.5; priTranVecs[2][1] = 0.5;
    priTranVecs[0][2] = -0.5; priTranVecs[1][2] = 0.5; priTranVecs[2][2] = 0.5;
    /*process*/
    ConstructReducedLattice();
    ConstructLattice();
}
