/* include */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* constants */
#define MAX_LATTICE_NUMBER 20000 // maximum number of lattices
#define MAX_ATOM_NUMBER 200000   // maximum number of atoms
#define MAX_CELL_ATOM_NUMBER 10  // maximum number of atoms in a cell
#define MAX_NEIGHBOR_NUMBER 2000 // maximum number of neighbors
// L-J parameters for Ne
#define LJ_EPSILON 0.0031
#define LJ_SIGMA 2.74
#define LJ_2_EPSILON 0.0062                     // 4 * LJ_EPSILON
#define LJ_24_EPSILON_SIGMA_6 31.48301472289983 // 24 * LJ_EPSILON * pow(LJ_SIGMA, 6)
#define LJ_2_SIGMA_6 846.3176000779524          // 2 * pow(LJ_SIGMA, 6)
// Math constants
#define PI 3.14159265358979323846
// Physical constants
#define K_B 0.0000861733          // unit: eV/K
#define UNIT_MASS 0.0001036426965 // unit: eV/(ps/A)^2
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
    double dr[3]; // from atom to neighbor
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
    double minDirection[3];
    double lastForce_CG[3];
    double startR_lineMin[3];
    double velocity[3];
    double acceleration[3];
    // Verlet
    double lastR_verlet[3];
    // Velocity Verlet
    double lastA_vverlet[3];
    double stress_Volume[6];

    // Andersen Barostat
    double accelerationModify[3];
    double velocityModify[3];
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
int boxOrthogonal;

double neighborCutoff;
double potentialCutoff_LJ;
double totalPotentialEnergy;
char potentialName[20];
int nStep;
int neighborInterval;

double energyTolerance_Minimize;
char minimizeStyle[20];

double delta_lineMinimize;

double typeMasses[5];
char dynamicStyle[20];

double xi_NoseHoover;

/* function declarations */
void ConstructReducedLattice();
void ConstructLattice();
void ConstructCrystal();
void Dump_xyz(char fileName[20]);
double VecDotMul(double vec1[3], double vec2[3]);
void VecCroMul(double vec1[3], double vec2[3], double vecOut[3]);
void ComputeRecTranVecs(double tranVecs[3][3], double recTranVecs[3][3]);
void ComputeAtomReR();

void ComputeAtomBoxReR(int n);
void PBC_r();
void PBC_dr();
void PBC_r_general();
void PBC_dr_general(int i, int j, double dr[3]);
void PBC_r_orthogonal();
void PBC_dr_orthogonal(int i, int j, double dr[3]);
void ConstructStdCrystal_BCC(double latticeConstant, int length);
void ConstructStdCrystal_FCC(double latticeConstant, int length);
void Dump_lammpstrj(char fileName[20], int isNewFile, int dumpStep);

void DeleteAtomByIndex(int index);
void DeleteAtomsByShpereRegion(double center[3], double radius);
void DeleteAtomsByBlockRegion(double block[3][2]);
void InsertAtom(double r[3], int type);
void EdgeDislocation_100(double latticeConstant);

void ConstructNeighborList();
void Potential_LJ(int isEnergy, int isForce);
void Potential_EAM(int isEnergy, int isForce);
void Potential(int isEnergy, int isForce);
void NeighborList(int isForceConstruct);
void UpdateNeighborList();

void Minimize();
void MinDirection();
void MinDirection_SD();
void MinDirection_CG();
void LineMinimize();

// chpater 5
double GaussianRandom(double mu, double sigma);
void InitMassUnit();
void VelocityMaxwell(double temperature);
void ZeroMomentum();
void InitVelocity(double temperature);
void Dynamics(double stopTime, double timeStep);
void IterRun();
void ComputeAcceleration();
void IterRun_Euler(double timeStep);
void IterRun_Verlet(double timeStep);
void IterRun_VelocityVerlet(double timeStep);
double ComputeTotalKineticEnergy();
double ComputeTemperature();
double ComputeBoxVolume();
void ComputeNonPBCForce(double nonPBCForce[3][3]);
void ComputeStress(double stress[6]);
void Barostat(double stress[6], double targetStress[3], int frequency, double timeStep, char algorithm[20]);
void Barostat_Berendsen(double stress[6], double targetStress[3], int frequency, double timeStep);
void Barostat_Andersen(double stress[6], double targetStress[3], int frequency, double timeStep);
void InitDynamic();
void ReviveVelocity();

void Thermostat(double temperature, double targetTemperature, int frequency, double timeStep, char algorithm[20]);
void Thermostat_Berendsen(double temperature, double targetTemperature, int frequency, double timeStep);
void Thermostat_NoseHoover(double temperature, double targetTemperature, int frequency, double timeStep);

void ConstructStdCrystal_BCC_Shear(double latticeConstant, int length, double xy);

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
    if (latticePointNumber > MAX_LATTICE_NUMBER)
    {
        printf("Error: lattice point number exceeds the maximum number.\n");
        exit(1);
    }
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
    if (atomNumber > MAX_ATOM_NUMBER)
    {
        printf("Error: atom number exceeds the maximum number of atoms.\n");
        exit(1);
    }
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

void PBC_r_orthogonal()
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

void PBC_dr_orthogonal(int i, int j, double dr[3])
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
    int n;
    if (boxOrthogonal == 1)
    {
        PBC_r_orthogonal();
    }
    else
    {
        for (n = 0; n < atomNumber; n++)
        {
            ComputeAtomBoxReR(n);
        }
        PBC_r_general();
    }
}

void PBC_dr(int i, int j, double dr[3])
{
    if (boxOrthogonal == 1)
    {
        PBC_dr_orthogonal(i, j, dr);
    }
    else
    {
        ComputeAtomBoxReR(i);
        ComputeAtomBoxReR(j);
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
    boxOrthogonal = 1;

    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    PBC_r();
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
    boxOrthogonal = 1;

    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    PBC_r();
}

void Dump_lammpstrj(char fileName[20], int isNewFile, int dumpStep)
{
    int n;
    FILE *fp;
    if (0) // boxOrthogonal != 1)
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
    fprintf(fp, "ITEM: ATOMS id type x y z pe fx fy fz vx vy vz\n");
    for (n = 0; n < atomNumber; n++)
    {
        fprintf(fp, "%d %d %30.20f %30.20f %30.20f %30.20f %30.20f %30.20f %30.20f  %30.20f  %30.20f  %30.20f\n",
                atoms[n].id, atoms[n].type, atoms[n].r[0], atoms[n].r[1], atoms[n].r[2],
                atoms[n].potentialEnergy,
                atoms[n].force[0], atoms[n].force[1], atoms[n].force[2],
                atoms[n].velocity[0], atoms[n].velocity[1], atoms[n].velocity[2]);
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
    if (boxOrthogonal != 1)
    {
        printf("Error: EdgeDislocation_100() only works in cuboid.\n");
        exit(1);
    }
    double deleteBlock[3][2] = {
        {boxStartPoint[0] + boxTranVecs[0][0] - latticeConstant - 0.1, boxStartPoint[0] + boxTranVecs[0][0] + 0.1},
        {boxStartPoint[1] - 0.1, boxStartPoint[1] + boxTranVecs[1][1] + 0.1},
        {boxStartPoint[2] + boxTranVecs[2][2] / 2 - 0.1, boxStartPoint[2] + boxTranVecs[2][2] + 0.1}};

    DeleteAtomsByBlockRegion(deleteBlock);
    // shift atoms
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
    int tmpNeighborNumber;

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
                tmpNeighborNumber = atoms[i].neighborNumber++;
                if (tmpNeighborNumber > MAX_NEIGHBOR_NUMBER)
                {
                    printf("Error: too many neighbors.\n");
                    exit(1);
                }
                tmpNeighborNumber = atoms[j].neighborNumber++;
                if (tmpNeighborNumber > MAX_NEIGHBOR_NUMBER)
                {
                    printf("Error: too many neighbors.\n");
                    exit(1);
                }
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
        double all_phi = 0;
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

                all_phi += energyPhi;
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

void NeighborList(int isForceConstruct)
{
    if (isForceConstruct)
    {
        ConstructNeighborList();
    }
    else if (nStep % neighborInterval == 0)
    {
        ConstructNeighborList();
    }
    else
    {
        UpdateNeighborList();
    }
}

void UpdateNeighborList()
{
    int i, j, d;
    int jAtomIndex;
    double dr[3];
    for (i = 0; i < atomNumber; i++)
    {
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            jAtomIndex = atoms[i].neighbors[j].index;
            PBC_dr(i, jAtomIndex, dr);
            atoms[i].neighbors[j].distance = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
            for (d = 0; d < 3; d++)
            {
                atoms[i].neighbors[j].dr[d] = dr[d];
            }
        }
    }
}

void Minimize()
{
    double potentialEnergy_last;
    nStep = 0;
    printf("\n---Minimization start---\n");
    printf("iter pe dE\n");
    NeighborList(1);
    Potential(1, 1);
    printf("%d %20.10f\n", nStep, totalPotentialEnergy);
    do
    {
        potentialEnergy_last = totalPotentialEnergy;
        MinDirection();
        LineMinimize();
        Potential(0, 1);
        nStep += 1;
        printf("%d %30.20f  %30.20f\n", nStep, totalPotentialEnergy, totalPotentialEnergy - potentialEnergy_last);
    } while (fabs(potentialEnergy_last - totalPotentialEnergy) > energyTolerance_Minimize);
    printf("---Minimization end---\n");
}

void MinDirection()
{
    if (strcmp(minimizeStyle, "SD") == 0)
        MinDirection_SD();
    else if (strcmp(minimizeStyle, "CG") == 0)
        MinDirection_CG();
    else
        printf("\nMinimize style not found. Program terminated.\n"), exit(-1);
}

void MinDirection_SD()
{
    int i, d;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].minDirection[d] = atoms[i].force[d];
        }
    }
}

void MinDirection_CG()
{
    int i, d;
    double beta_up, beta_down, beta;
    double directionCheck;
    if (nStep == 0)
    {
        for (i = 0; i < atomNumber; i++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].minDirection[d] = atoms[i].force[d];
                atoms[i].lastForce_CG[d] = atoms[i].force[d];
            }
        }
    }
    else
    {
        beta_up = 0;
        beta_down = 0;
        for (i = 0; i < atomNumber; i++)
        {
            beta_up += VecDotMul(atoms[i].force, atoms[i].force);
            beta_down += VecDotMul(atoms[i].lastForce_CG, atoms[i].lastForce_CG);
        }
        beta = beta_up / beta_down;
        directionCheck = 0;
        for (i = 0; i < atomNumber; i++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].minDirection[d] = atoms[i].force[d] + beta * atoms[i].minDirection[d];
                atoms[i].lastForce_CG[d] = atoms[i].force[d];
                directionCheck += atoms[i].minDirection[d] * atoms[i].force[d];
            }
        }
        if (directionCheck < 0)
        {
            for (i = 0; i < atomNumber; i++)
            {
                for (d = 0; d < 3; d++)
                {
                    atoms[i].minDirection[d] = atoms[i].force[d];
                }
            }
        }
    }
}

void LineMinimize()
{
    double startPotentialEnergy;
    double k, c;
    double slop;
    double lambda_t;
    int i, d;

    k = 0.5;
    c = 0.5;
    slop = 0;
    lambda_t = 100;

    startPotentialEnergy = totalPotentialEnergy;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            slop += atoms[i].minDirection[d] * atoms[i].minDirection[d];
            atoms[i].startR_lineMin[d] = atoms[i].r[d];
        }
    }
    slop = sqrt(slop);

    lambda_t /= k;
    do
    {
        lambda_t *= k;
        for (i = 0; i < atomNumber; i++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].r[d] = atoms[i].startR_lineMin[d] + lambda_t * atoms[i].minDirection[d];
            }
        }
        PBC_r();
        NeighborList(1);
        Potential(1, 0);
    } while (totalPotentialEnergy > startPotentialEnergy - c * slop * lambda_t);
}

double GaussianRandom(double mu, double sigma)
{
    double u1, u2;
    double r;
    double z;
    do
    {
        u1 = -1 + ((double)rand() / RAND_MAX) * 2;
        u2 = -1 + ((double)rand() / RAND_MAX) * 2;
        r = u1 * u1 + u2 * u2;
    } while (r >= 1 || r == 0);
    z = sqrt(-2.0 * log(r) / r) * u1 * sigma + mu;
    return z;
}

void InitMassUnit()
{
    int i;
    for (i = 0; i < 5; i++)
    {
        typeMasses[i] *= UNIT_MASS;
    }
}

void VelocityMaxWell(double temperature)
{
    int i, d;
    double mass, sigma;
    for (i = 0; i < atomNumber; i++)
    {
        mass = typeMasses[atoms[i].type];
        for (d = 0; d < 3; d++)
        {
            sigma = sqrt(K_B * temperature / mass);
            atoms[i].velocity[d] = GaussianRandom(0.0, sigma);
        }
    }
}

void ZeroMomentum()
{
    int i, d;
    double momentum[3] = {0, 0, 0};

    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            momentum[d] += atoms[i].velocity[d] * typeMasses[atoms[i].type];
        }
    }

    for (d = 0; d < 3; d++)
    {
        momentum[d] /= atomNumber;
    }

    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].velocity[d] -= momentum[d] / typeMasses[atoms[i].type];
        }
    }
}

void InitVelocity(double temperature)
{
    VelocityMaxWell(temperature);
    ZeroMomentum();
}

void IterRun(double timeStep)
{
    if (strcmp(dynamicStyle, "Euler") == 0)
        IterRun_Euler(timeStep);
    else if (strcmp(dynamicStyle, "Verlet") == 0)
        IterRun_Verlet(timeStep);
    else if (strcmp(dynamicStyle, "VelocityVerlet") == 0)
        IterRun_VelocityVerlet(timeStep);
    else
    {
        printf("Eorror: Dynamic style %s not found.", dynamicStyle);
        exit(1);
    }
    ReviveVelocity();
}

void ComputeAcceleration()
{
    int n, d;
    NeighborList(0);
    Potential(0, 1);
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].acceleration[d] = atoms[n].force[d] / typeMasses[atoms[n].type] + atoms[n].accelerationModify[d];
        }
    }
}

void IterRun_Euler(double timeStep)
{
    int n, d;
    ComputeAcceleration();
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].r[d] += atoms[n].velocity[d] * timeStep;
            atoms[n].velocity[d] += atoms[n].acceleration[d] * timeStep;
        }
    }
    PBC_r();
}

void IterRun_Verlet(double timeStep)
{
    int n, d;
    double r_tmp;

    ComputeAcceleration();
    if (nStep == 0)
    {
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[n].lastR_verlet[d] = atoms[n].r[d];
                atoms[n].r[d] += atoms[n].velocity[d] * timeStep + 0.5 * atoms[n].acceleration[d] * timeStep * timeStep;
                atoms[n].velocity[d] = (atoms[n].r[d] - atoms[n].lastR_verlet[d]) / timeStep;
            }
        }
    }
    else
    {
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                r_tmp = atoms[n].r[d];
                atoms[n].r[d] = 2 * atoms[n].r[d] - atoms[n].lastR_verlet[d] + atoms[n].acceleration[d] * timeStep * timeStep;
                atoms[n].lastR_verlet[d] = r_tmp;
                atoms[n].velocity[d] = (atoms[n].r[d] - atoms[n].lastR_verlet[d]) / timeStep;
            }
        }
    }
    PBC_r();
}

void IterRun_VelocityVerlet(double timeStep)
{
    int n, d;
    if (nStep == 0)
    {
        ComputeAcceleration();
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[n].lastA_vverlet[d] = atoms[n].acceleration[d];
            }
        }
    }

    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].r[d] += atoms[n].velocity[d] * timeStep + 0.5 * atoms[n].acceleration[d] * timeStep * timeStep;
        }
    }
    PBC_r();
    ComputeAcceleration();
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].velocity[d] += 0.5 * (atoms[n].acceleration[d] + atoms[n].lastA_vverlet[d]) * timeStep;
            atoms[n].lastA_vverlet[d] = atoms[n].acceleration[d];
        }
    }
}

double ComputeTotalKineticEnergy()
{
    int n, d;
    double e = 0;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            e += 0.5 * typeMasses[atoms[n].type] * atoms[n].velocity[d] * atoms[n].velocity[d];
        }
    }
    return e;
}

double ComputeTemperature()
{
    double Ek, T;
    Ek = ComputeTotalKineticEnergy();
    T = 2. / 3. / K_B / atomNumber * Ek;
    return T;
}

double ComputeBoxVolume()
{
    double areaVector[3];
    VecCroMul(boxTranVecs[0], boxTranVecs[1], areaVector);
    return VecDotMul(areaVector, boxTranVecs[2]);
}

void ComputeNonPBCForce(double nonPBCForce[3][3])
{
    int n, i, j;
    for (i = 0; i < 3; i++)
    {
        boxTranVecs[i][i] *= 2;
        NeighborList(1);
        Potential(0, 1);
        for (j = 0; j < 3; j++)
        {
            nonPBCForce[i][j] = 0;
        }
        for (n = 0; n < atomNumber; n++)
        {
            if (atoms[n].r[i] >= boxTranVecs[i][i] / 4)
            {
                for (j = 0; j < 3; j++)
                {
                    nonPBCForce[i][j] += atoms[n].force[j];
                }
            }
        }
        boxTranVecs[i][i] /= 2;
    }
}

void ComputeStress(double stress[6])
{
    int n, d, i, j;
    double volume;
    double nonPBCForce[3][3], sumAtomForce[3][3];
    // only for orthogonal box with start point on (0,0,0)
    if (boxOrthogonal != 1)
    {
        printf("Error: computing stress in wrong box, check function ComputeStress()\n");
        exit(1);
    }

    // box virial
    ComputeNonPBCForce(nonPBCForce);
    NeighborList(1);
    Potential(0, 1);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            sumAtomForce[i][j] = 0;
        }
    }

    for (n = 0; n < atomNumber; n++)
    {
        for (i = 0; i < 3; i++)
        {
            if (atoms[n].r[i] >= boxTranVecs[i][i] / 2)
            {
                for (j = 0; j < 3; j++)
                {
                    sumAtomForce[i][j] += atoms[n].force[j];
                }
            }
        }
    }
    stress[0] = boxTranVecs[0][0] * (nonPBCForce[0][0] - sumAtomForce[0][0]);
    stress[1] = boxTranVecs[1][1] * (nonPBCForce[1][1] - sumAtomForce[1][1]);
    stress[2] = boxTranVecs[2][2] * (nonPBCForce[2][2] - sumAtomForce[2][2]);
    stress[3] = boxTranVecs[0][0] * (nonPBCForce[0][1] - sumAtomForce[0][1]);
    stress[4] = boxTranVecs[0][0] * (nonPBCForce[0][2] - sumAtomForce[0][2]);
    stress[5] = boxTranVecs[1][1] * (nonPBCForce[1][2] - sumAtomForce[1][2]);

    // add atom virial and ke
    for (n = 0; n < atomNumber; n++)
    {
        stress[0] += atoms[n].r[0] * atoms[n].force[0] + atoms[n].velocity[0] * atoms[n].velocity[0] * typeMasses[atoms[n].type];
        stress[1] += atoms[n].r[1] * atoms[n].force[1] + atoms[n].velocity[1] * atoms[n].velocity[1] * typeMasses[atoms[n].type];
        stress[2] += atoms[n].r[2] * atoms[n].force[2] + atoms[n].velocity[2] * atoms[n].velocity[2] * typeMasses[atoms[n].type];
        stress[3] += atoms[n].r[0] * atoms[n].force[1] + atoms[n].velocity[0] * atoms[n].velocity[1] * typeMasses[atoms[n].type];
        stress[4] += atoms[n].r[0] * atoms[n].force[2] + atoms[n].velocity[0] * atoms[n].velocity[2] * typeMasses[atoms[n].type];
        stress[5] += atoms[n].r[1] * atoms[n].force[2] + atoms[n].velocity[1] * atoms[n].velocity[2] * typeMasses[atoms[n].type];
    };

    // stress
    volume = ComputeBoxVolume();
    for (d = 0; d < 6; d++)
    {
        stress[d] *= -160.21766208 / volume; // 1 eV/Angstrom3 = 160.21766208 GPa
    }
}

void Barostat(double stress[3], double targetStress[3], int frequency, double timeStep, char algorithm[20])
{
    if (boxOrthogonal != 1 || !(boxStartPoint[0] == 0 && boxStartPoint[1] == 0 && boxStartPoint[2] == 0))
    {
        printf("Error: Barostat in wrong box, check function Barostat()\n");
        exit(1);
    }
    if (strcmp(algorithm, "Berendsen") == 0)
    {
        Barostat_Berendsen(stress, targetStress, frequency, timeStep);
    }
    else if (strcmp(algorithm, "Andersen") == 0)
    {
        Barostat_Andersen(stress, targetStress, frequency, timeStep);
    }
    else
    {
        printf("Barostat algorithm %s not found.", algorithm);
        exit(1);
    }
}

void Barostat_Berendsen(double stress[6], double targetStress[3], int frequency, double timeStep)
{
    static int count = 0;
    double k_tau = 1; // parameter
    int n, d;
    double lambda;
    double deltaTime;
    deltaTime = frequency * timeStep;
    if (count == 0)
    {
        for (d = 0; d < 3; d++)
        {
            lambda = 1 + k_tau * deltaTime * (targetStress[d] - stress[d]);
            boxTranVecs[d][d] *= lambda;
            for (n = 0; n < atomNumber; n++)
            {
                atoms[n].r[d] *= lambda;
            }
        }
        PBC_r();
    }
    count += 1;
    if (count == frequency)
    {
        count = 0;
    }
}

void Barostat_Andersen(double stress[6], double targetStress[3], int frequency, double timeStep)
{
    if (!(targetStress[0] == targetStress[1] && targetStress[0] == targetStress[2]))
    {
        printf("Error: Andersen barostat can only used for target sigma_xx == sigma_yy == sigma_zz\n");
        exit(1);
    }
    double M = 10; // piston mass, unit: unit: eV/(ps/A)^2

    static int count = 0;
    static double pistonVelocity[3] = {0, 0, 0};
    int n, i;
    double volume;
    double pistonAcceleration; // PistonPi/M
    double deltaTime;
    double AndersenVelocity, pistonVelocityPerLength; // PI/M/L

    deltaTime = frequency * timeStep;

    if (count == 0)
    {
        volume = ComputeBoxVolume();
        for (i = 0; i < 3; i++)
        {
            pistonAcceleration = (targetStress[i] - stress[i]) * volume / boxTranVecs[i][i] / M;
            pistonVelocity[i] += deltaTime * pistonAcceleration;
        }
    }
    count += 1;
    if (count == frequency)
    {
        count = 0;
    }

    for (i = 0; i < 3; i++)
    {
        boxTranVecs[i][i] += pistonVelocity[i] * timeStep;
        pistonVelocityPerLength = pistonVelocity[i] / boxTranVecs[i][i];
        for (n = 0; n < atomNumber; n++)
        {
            AndersenVelocity = atoms[n].r[i] * pistonVelocityPerLength;
            atoms[n].velocity[i] += AndersenVelocity;
            atoms[n].velocityModify[i] = AndersenVelocity;
            atoms[n].accelerationModify[i] = -atoms[n].velocity[i] * pistonVelocityPerLength;
        }
    }
}

void InitDynamic()
{
    int n, d;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].accelerationModify[d] = 0;
            atoms[n].velocityModify[d] = 0;
        }
    }
    xi_NoseHoover = 0;
}

void ReviveVelocity()
{
    int n, d;
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].velocity[d] -= atoms[n].velocityModify[d];
        }
    }
}

void Thermostat(double temperature, double targetTemperature, int frequency, double timeStep, char algorithm[20])
{
    if (strcmp(algorithm, "Berendsen") == 0)
    {
        Thermostat_Berendsen(temperature, targetTemperature, frequency, timeStep);
    }
    else if (strcmp(algorithm, "Nose-Hoover") == 0)
    {
        Thermostat_NoseHoover(temperature, targetTemperature, frequency, timeStep);
    }
    else
    {
        printf("Thermostat algorithm %s not found.", algorithm);
        exit(1);
    }
}

void Thermostat_Berendsen(double temperature, double targetTemperature, int frequency, double timeStep)
{
    static int count = 0;
    double lambda, deltaTime, deltaTime_tau; // deltaTime_tau: daltaTime/tau
    int n, d;
    double tau = 0.01;
    if (count == 0)
    {
        deltaTime_tau = frequency * timeStep / tau;
        lambda = sqrt(1 + deltaTime_tau * (targetTemperature / temperature - 1));
        for (n = 0; n < atomNumber; n++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[n].velocity[d] *= lambda;
            }
        }
    }
    count += 1;
    if (count == frequency)
    {
        count = 0;
    }
}

void Thermostat_NoseHoover(double temperature, double targetTemperature, int frequency, double timeStep)
{
    static int count = 0;
    double Q = 0.1; // parameter
    double xi_velocity;
    int n, d;
    if (count == 0)
    {
        xi_velocity = 3 * atomNumber * K_B / Q * (temperature - targetTemperature);
        xi_NoseHoover += xi_velocity * frequency * timeStep;
    }
    count += 1;
    if (count == frequency)
    {
        count = 0;
    }
    for (n = 0; n < atomNumber; n++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[n].accelerationModify[d] = -xi_NoseHoover * atoms[n].velocity[d];
        }
    }
}

void Dynamics(double stopTime, double timeStep)
{
    double time;
    int n, d;
    double temperature;
    double targetTemperature = 300;
    FILE *fp;
    char fileName[50] = "debug/thermostat/time-temperature.nh.csv";
    char dumpName[50] = "debug/thermostat/run.nh.dump";

    InitDynamic();

    time = 0;
    nStep = 0;
    fp = fopen(fileName, "w");
    fprintf(fp, "step time temperature\n");
    Dump_lammpstrj(dumpName, 1, nStep);
    while (time <= stopTime)
    {
        if (nStep % 100 == 0)
        {
            temperature = ComputeTemperature();
        }
        if (nStep % 100 == 0)
        {
            fprintf(fp, "%d %f %f\n", nStep, time, temperature);
            printf("%d %f %f\n", nStep, time, temperature);
            Dump_lammpstrj(dumpName, 0, nStep);
        }
        if (nStep >= 4000)
        {
            Thermostat(temperature, targetTemperature, 100, timeStep, "Nose-Hoover");
        }
        IterRun(timeStep);
        nStep += 1;
        time += timeStep;
    }
    fclose(fp);
}

void ConstructStdCrystal_BCC_Shear(double latticeConstant, int length, double xy)
{
    latticeSizes[0][0] = 0;
    latticeSizes[0][1] = length;
    latticeSizes[1][0] = 0;
    latticeSizes[1][1] = length;
    latticeSizes[2][0] = 0;
    latticeSizes[2][1] = length;

    priTranVecs[0][0] = latticeConstant;
    priTranVecs[0][1] = latticeConstant * xy;
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
    cellAtomRs[1][1] = 0.5 * latticeConstant + 0.5 * latticeConstant * xy;
    cellAtomRs[1][2] = 0.5 * latticeConstant;
    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 1;

    boxStartPoint[0] = 0;
    boxStartPoint[1] = 0;
    boxStartPoint[2] = 0;

    boxTranVecs[0][0] = latticeConstant * length;
    boxTranVecs[0][1] = latticeConstant * length * xy;
    boxTranVecs[0][2] = 0;
    boxTranVecs[1][0] = 0;
    boxTranVecs[1][1] = latticeConstant * length;
    boxTranVecs[1][2] = 0;
    boxTranVecs[2][0] = 0;
    boxTranVecs[2][1] = 0;
    boxTranVecs[2][2] = latticeConstant * length;
    boxOrthogonal = 0;
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();
    PBC_r();
}

/* main */
int main()
{
    /* parameters */
    unsigned int randomSeed;
    randomSeed = 1;
    srand(randomSeed);

    typeMasses[1] = 183.84; // for W
    InitMassUnit();
    strcpy(potentialName, "EAM");
    neighborCutoff = 6;
    neighborInterval = 100;
    strcpy(dynamicStyle, "VelocityVerlet");

    /* processing*/

    double strain;
    int n;
    double stress[6];

    printf("strain_xy energy\n");
    for (strain = -0.001; strain < 0.00101; strain += 0.0001)
    {
        ConstructStdCrystal_BCC_Shear(3.14, 10, strain);
        InitVelocity(0);
        NeighborList(1);
        Potential(1, 0);
        printf("%f%% %f\n", strain * 100, totalPotentialEnergy);
    }

    return 0;
}
