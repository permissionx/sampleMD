/*include*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* constants */
// chapter 1
#define MAX_LATTICE_NUMBER 20000 //maximum number of lattices
#define MAX_ATOM_NUMBER 200000  //maximum number of atoms
#define MAX_CELL_ATOM_NUMBER 10 //maximum number of atoms in a cell

// chapter 2
#define MAX_NEIGHBOR_NUMBER 2000 // maximum number of neighbors

// chapter 3 
// L-J parameters for Ne (energy: eV, length: A)
#define LJ_EPSILON 0.0031
#define LJ_SIGMA 2.74
#define LJ_2_EPSILON 0.0062                     // 4 * LJ_EPSILON
#define LJ_24_EPSILON_SIGMA_6 31.48301472289983 // 24 * LJ_EPSILON * pow(LJ_SIGMA, 6)
#define LJ_2_SIGMA_6 846.3176000779524          // 2 * pow(LJ_SIGMA, 6)

// EAM parameters for W
// data source: M. C. Marinica, et al., J. Phys. Condens. Matter 25, (2013).
const int n_phi_EAM = 15;
const int n_rho_EAM = 4;
const double rc_EAM = 2.002970124727;
const double pho_rc_EAM = 1.193547157792;
double a_phi_EAM[15] =
    {
        0.960851701343041e2,
        -0.184410923895214e3,
        0.935784079613550e2,
        -0.798358265041677e1,
        0.747034092936229e1,
        -0.152756043708453e1,
        0.125205932634393e1,
        0.163082162159425e1,
        -0.141854775352260e1,
        -0.819936046256149e0,
        0.198013514305908e1,
        -0.696430179520267e0,
        0.304546909722160e-1,
        -0.163131143161660e1,
        0.138409896486177e1};
double a_rho_EAM[4] =
    {
        -0.420429107805055e1,
        0.518217702261442e0,
        0.562720834534370e-1,
        0.344164178842340e-1};

const double a1_f_EAM = -5.946454472402710;
const double a2_f_EAM = -0.049477376935239;

double delta_phi_EAM[15] =
    {
        2.5648975000,
        2.6297950000,
        2.6946925000,
        2.8663175000,
        2.9730450000,
        3.0797725000,
        3.5164725000,
        3.8464450000,
        4.1764175000,
        4.7008450000,
        4.8953000000,
        5.0897550000,
        5.3429525000,
        5.4016950000,
        5.4604375000};
double delta_rho_EAM[4] =
    {
        2.500000000000000,
        3.100000000000000,
        3.500000000000000,
        4.900000000000000};

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
    double dr[3];  // from atom to neighbor
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
    double rho_EAM;
    double af_EAM;
};




/* global variables */
// chapter 1
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

// chapter 2
double boxStartPoint[3];
double boxTranVecs[3][3]; // box translation vectors
double boxRecTranVecs[3][3]; //box reciprocal translation vectors
int boxPerpendicular; 

// chapter 3
double neighborCutoff;
double potentialCutoff_LJ;
double totalPotentialEnergy;
char potentialName[20];



/* function declarations */
// chapter 1 
void ConstructReducedLattice();
void ConstructLattice();
void ConstructCrystal();
void Dump_xyz(char fileName[20]);

double VecDotMul(double vec1[3], double vec2[3]);
void VecCroMul(double vec1[3], double vec2[3], double vecOut[3]);
void ComputeRecTranVecs(double tranVecs[3][3], double recTranVecs[3][3]);

void ComputeAtomReR();

// chapter 2 
void ComputeAtomBoxReR();
void PBC_r_general();
void PBC_dr_general(int i, int j, double dr[3]);
void PBC_r();
void PBC_dr();
void PBC_r_vertical();
void PBC_dr_vertical(int i, int j, double dr[3]);

void ConstructStdCrystal_BCC(double latticeConstant, int length);
void ConstructStdCrystal_FCC(double latticeConstant, int length);

void Dump_lammpstrj(char fileName[20], int isNewFile, int dumpStep);

void DeleteAtomByIndex(int index);
void DeleteAtomsByShpereRegion(double center[3], double radius);
void DeleteAtomsByBlockRegion(double block[3][2]);

void InsertAtom(double r[3], int type);

void EdgeDislocation_100(double latticeConstant);

// chatper 3
void ConstructNeighborList();
void Potential_LJ(int isEnergy, int isForce);
void Potential_EAM(int isEnergy, int isForce);
void Potential(int isEnergy, int isForce);

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
            recTranVecs[i][d] /= cellVol;  // 2pi factor ignored
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


void ComputeAtomBoxReR()
{
    int n, d;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].boxReR[d] = boxRecTranVecs[d][0] * (atoms[n].r[0] - boxStartPoint[0]) + boxRecTranVecs[d][1] * (atoms[n].r[1] - boxStartPoint[1])+ boxRecTranVecs[d][2] * (atoms[n].r[2] - boxStartPoint[2]);
            if (atoms[n].boxReR[d] < -1 || atoms[n].boxReR[d] >= 2)
            {
                printf("Lost atom %d!\n", atoms[n].id);
                exit(1);
            }
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

void PBC_r_vertical()
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

void PBC_dr_vertical(int i, int j, double dr[3])
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
    if (boxPerpendicular == 1)
    {
        PBC_r_vertical();
    }
    else
    {
        ComputeAtomBoxReR();
        PBC_r_general();
    }
}

void PBC_dr(int i, int j, double dr[3])
{
    if (boxPerpendicular == 1)
    {
        PBC_dr_vertical(i, j, dr);
    }
    else
    {
        ComputeAtomBoxReR();
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
    //shift atoms
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
                atoms[i].neighborNumber++;
                atoms[j].neighborNumber++;
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

void Potential_EAM(int isEnergy, int isForce)
{
    int i, j;
    double distance;
    double energyPhi, energyRho;
    double forcePhi, forceRho, tmpForce;
    double tmpRho, sqrtTmpRho;
    double tmpVariable; // very short lifetime: two lines
    int n;
    int d;
    int jAtomIndex;
    if (isEnergy)
    {
        totalPotentialEnergy = 0.0;
    }
    // compute atom rho and rho self-related expressions
    for (i = 0; i < atomNumber; i++)
    {
        tmpRho = 0.0;
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            distance = atoms[i].neighbors[j].distance;
            if (distance > rc_EAM)
            {

                for (n = 0; n < n_rho_EAM; n++)
                {
                    if (distance < delta_rho_EAM[n])
                    {
                        tmpVariable = delta_rho_EAM[n] - distance;
                        tmpRho += a_rho_EAM[n] * tmpVariable * tmpVariable * tmpVariable;
                    }
                }
            }
            else
            {
                tmpRho += pho_rc_EAM;
            }
        }
        atoms[i].rho_EAM = tmpRho;
        sqrtTmpRho = sqrt(tmpRho);
        if (isEnergy)
        {
            // rho component for energy
            energyRho = a1_f_EAM * sqrtTmpRho + a2_f_EAM * tmpRho * tmpRho;
            totalPotentialEnergy += energyRho;
            atoms[i].potentialEnergy = energyRho;
        }
        if (isForce)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].force[d] = 0;
            }
            atoms[i].af_EAM = 0.5 * a1_f_EAM / sqrtTmpRho + 2 * a2_f_EAM * tmpRho;
        }
    }

    for (i = 0; i < atomNumber; i++)
    {
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            distance = atoms[i].neighbors[j].distance;
            if (isEnergy)
            // phi component for energy
            {
                energyPhi = 0.0;
                for (n = 0; n < n_phi_EAM; n++)
                {
                    if (distance < delta_phi_EAM[n])
                    {
                        tmpVariable = delta_phi_EAM[n] - distance;
                        energyPhi += a_phi_EAM[n] * tmpVariable * tmpVariable * tmpVariable;
                    }
                }
                energyPhi *= 0.5;
                atoms[i].potentialEnergy += energyPhi;
                totalPotentialEnergy += energyPhi;
            }
            if (isForce)
            {
                // phi component for force
                forcePhi = 0.0;
                for (n = 0; n < n_phi_EAM; n++)
                {
                    if (distance < delta_phi_EAM[n])
                    {
                        tmpVariable = delta_phi_EAM[n] - distance;
                        forcePhi += -3 * a_phi_EAM[n] * tmpVariable * tmpVariable;
                    }
                }

                // rho component for force
                forceRho = 0.0;
                for (n = 0; n < n_rho_EAM; n++)
                {
                    if (distance < delta_rho_EAM[n])
                    {
                        tmpVariable = delta_rho_EAM[n] - distance;
                        forceRho += -3 * a_rho_EAM[n] * tmpVariable * tmpVariable;
                    }
                }

                // sum force
                jAtomIndex = atoms[i].neighbors[j].index;
                tmpForce = (((atoms[i].af_EAM + atoms[jAtomIndex].af_EAM) * forceRho) + forcePhi) / distance;
                for (d = 0; d < 3; d++)
                {
                    atoms[i].force[d] += tmpForce * atoms[i].neighbors[j].dr[d];
                }
            }
        }
    }
}

void Potential(int isEnergy, int isForce)
{
    if (strcmp(potentialName, "LJ") == 0)

    {
        Potential_LJ(isEnergy, isForce);
    }
    else if (strcmp(potentialName, "EAM") == 0)
    {
        Potential_EAM(isEnergy, isForce);
    }
    else
    {
        printf("Error: Potential %s not found.\n", potentialName);
        exit(1);
    }
}

/* main */
int main()
{
    /* parameters */
    neighborCutoff = 6.0;
    strcpy(potentialName, "EAM");

    /* processing*/
    ConstructStdCrystal_BCC(3.14, 10);
    double deleteBlock[3][2] = {{-0.1, 10.1 * 3.14}, {-0.1, 10.1 * 3.14}, {5.1 * 3.14, 10.1 * 3.14}};
    DeleteAtomsByBlockRegion(deleteBlock);
    ConstructNeighborList();
    Potential(1, 1);

    /* output */
    Dump_lammpstrj("output/3.7_surface-W-BCC.lammpstrj", 1, 1);

    return 0;
}
