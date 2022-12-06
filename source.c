/*include*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*constant*/
#define K_B 8.61733326E-5
#define PI 3.1415926
#define EBASE 2.71828182846

#define MAX_LATTICE_NUMBER 2000
#define MAX_ATOM_NUMBER 20000
#define MAX_CELL_ATOM_NUMBER 10
#define MAX_TYPES 10
#define MAX_NEIGHBOR_NUMBER 2000
#define UNITMASS 0.0001040459571284739

//LJ parameters
const double LJ_epsilon = 0.0031;
const double LJ_sigma = 2.74;


//EAM parameters
const  int n_phi_EAM = 15;
const  int n_rho_EAM = 4;
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
    0.138409896486177e1
};
double a_rho_EAM[4] =
{
    -0.420429107805055e1,
    0.518217702261442e0,
    0.562720834534370e-1,
    0.344164178842340e-1
};

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
    5.4604375000
};
double delta_rho_EAM[4] =
{
    2.500000000000000,
    3.100000000000000,
    3.500000000000000,
    4.900000000000000
};

/*class definition*/
struct LatticePoint
{
    int id;
    double reR[3], r[3];
};

struct Neighbor
{
    int index;
    double distance;
    double dr[3];
};

struct Atom
{
    double r[3];
    double reR_box[3];
    int id, type;
    double potentialEnergy;
    double rho_EAM;
    double force[3];
    int neighborNumber;
    struct Neighbor neighborList[MAX_NEIGHBOR_NUMBER]; // error if exceed
    double lastForce_CG[3];
    double minDirection[3];
    double lineMinStartR[3];
    double v[3], a[3];
    double m;
    double r_last_Verlet[3];
};


/*global variable declaration*/
struct LatticePoint latticePoints[MAX_LATTICE_NUMBER];
int latticePointNumber;
int latticeSizes[3][2];
double priTranVecs[3][3]; // primitive translation vectors
struct Atom atoms[MAX_ATOM_NUMBER];
int atomNumber;
int maxID;
int cellAtomNumber;
double cellAtomRs[MAX_CELL_ATOM_NUMBER][3];
int cellAtomTypes[MAX_CELL_ATOM_NUMBER];
double typeMasses[MAX_TYPES];
double boxStartPoint[3];
double boxTranVecs[3][3]; // box translation vectors
double boxRecTranVecs[3][3]; //box reciprocal translation vectors
double boxRecTranVecs_inv[3][3];
int boxOrthogonalFlag;

char potentialName[20];
double neighborCutoff;

double LJ_24_epsilon_sigma6;
double LJ_2_sigma6;
double LJ_4_epsilon;
double totalPotentialEnergy;

char minimizeStyle[20];
double energyCriterion_min;
double lambda_min;
double c_min;
double rho_min;

char dynamicStyle[20];
double totalTime;
double timeStep;
double time_2_Verlet;
double time_m_2_Verlet;

char dumpName[20];

double totalKineticEnergy;


double tau_NoseHoover;
double tau2_r_NoseHoover; // 1/tau^2
double xi_NoseHoover;
double tau_Berendsen;
char thermostatStyle[20];

/*function declaration*/
void ConstructReducedLattice();
void Mul_3_1(double a[3][3], double b[3], double p[3]); // matrix multiplication 3x1
void ConstructLattice();
void ConstructCrystal();
void Dump(int step);
double VecDotMul(double vec1[3], double vec2[3]); // vector dot multiplication
void VecCroMul(double vec1[3], double vec2[3], double vecOut[3]); // vector cross multiplication
void ComputeRecTranVecs(double tranVecs[3][3], double recTranVecs[3][3]);
void MatInv(double matIn[3][3], double matOut[3][3]);    //Matrix Inversion
void ComputeReR(double r[3], double recTranVecs_inv[3][3], double startPoint[3], double reR[3]);

void ComputeBoxRecTranVecs();
void ComputeAtomBoxReR();
void PBC_r();
void PBC_dr(int i, int j, double dr[3]);
void PBC_r_nonorthogonal();
void PBC_dr_nonorthogonal(int i, int j, double dr[3]);
void PBC_r_orthogonal();
void PBC_dr_orthogonal(int i, int j, double dr[3]);

void DeleteAtomByIndex(int index);
void DeleteAtomsByRegion(double block[3][2]);

void FindNeighbors();
void InitLJ();
void Potential_LJ(int energyFlag, int forceFlag);
void Potential_EAM(int energyFlag, int forceFlag);

void Minimize();
void MinDirection(int step);
void MinDirection_SD();
void MinDirection_CG(int step);
double LineMinimize();


double GenerateSpeed(double randomNumber, double temperature, double prefactor1, double prefactor2);
double MaxwellProbabilityDensity(double v, double prefactor1, double prefactor2);
void DistributeVelocity(double temperature);
void RandomSpherePoint(double radius, double vector[3]);
void ZeroMomentum();

void Dynamics();
void IterRun();
void LaunchRun();
void IterRun_Euler();
void IterRun_Verlet();
void LaunchRun_Verlet();
void DynamicsProcessing(int step, double time);

void ComputeTotalKineticEnergy();
double ComputeTemperature();
double ComputeStress(double stress[6]);
void VirialPairs(int i, double virial[6]);
void VirialPairs_LJ(int i, double virial[6]);
void VirialPairs_EAM(int i, double virial[6]);

void Thermostat_NoseHoover(int step, double targetTemperature);
void Friction_NoseHoover();
void Thermostat_Berendsen(double targetTemperature);

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
        Mul_3_1(priTranVecs, latticePoints[n].reR, latticePoints[n].r);
}


void ConstructCrystal()
{
    int nLattice, nAtom, nCellAtom;
    nAtom = 0;

    ConstructReducedLattice();
    ConstructLattice();
    for (nLattice = 0; nLattice < latticePointNumber; nLattice++)
    {
        for (nCellAtom = 0; nCellAtom < cellAtomNumber; nCellAtom++)
        {
            atoms[nAtom].r[0] = latticePoints[nLattice].r[0] + cellAtomRs[nCellAtom][0];
            atoms[nAtom].r[1] = latticePoints[nLattice].r[1] + cellAtomRs[nCellAtom][1];
            atoms[nAtom].r[2] = latticePoints[nLattice].r[2] + cellAtomRs[nCellAtom][2];

            atoms[nAtom].type = cellAtomTypes[nCellAtom];
            atoms[nAtom].id = nAtom;
            atoms[nAtom].m = typeMasses[atoms[nAtom].type] * UNITMASS;
            nAtom++;
        }
    }
    maxID = nAtom - 1;
    atomNumber = nAtom;
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
    int i, d; // d: direction
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

void solveInverseMatrix(double mat[3][3], double invMat[3][3])
{
    double det;
    int i, j;
    det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
          mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
          mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    if (fabs(det) < 1e-10)
    {
        printf("Singular matrix!\n");
        exit(1);
    }
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            invMat[i][j] = (mat[(i + 1) % 3][(j + 1) % 3] * mat[(i + 2) % 3][(j + 2) % 3] -
                            mat[(i + 1) % 3][(j + 2) % 3] * mat[(i + 2) % 3][(j + 1) % 3]) / det;
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

void ComputeReR(double r[3], double recTranVecs_inv[3][3], double startPoint[3], double reR[3])
{
    double relativeR[3];
    int i;
    for (i = 0; i < 3; i++)
    {
        relativeR[i] = r[i] - startPoint[i];
    }
    Mul_3_1(recTranVecs_inv, relativeR, reR);
}

void ComputeBoxRecTranVecs()
{
    ComputeRecTranVecs(boxTranVecs, boxRecTranVecs);
    MatInv(boxRecTranVecs, boxRecTranVecs_inv);
}

void ComputeAtomBoxReR()
{
    int n;
    for (n = 0; n < atomNumber; n++)
    {
        ComputeReR(atoms[n].r, boxRecTranVecs_inv, boxStartPoint, atoms[n].reR_box);
    }
}

void PBC_r()
{
    if (boxOrthogonalFlag)
    {
        PBC_r_orthogonal();
    }
    else
    {
        ComputeAtomBoxReR();
        PBC_r_nonorthogonal();
    }
}

void PBC_dr(int i, int j, double dr[3])
{
    if (boxOrthogonalFlag)
    {
        PBC_dr_orthogonal(i, j, dr);
    }
    else
    {
        ComputeAtomBoxReR();
        PBC_dr_nonorthogonal(i, j, dr);
    }
}


void PBC_r_nonorthogonal()
{
    int n, d;
    int outFlag;

    for (n = 0; n < atomNumber; n++)
    {
        outFlag = 0;
        for (d = 0; d < 3; d++)
        {
            if (atoms[n].reR_box[d] < 0)
            {
                atoms[n].reR_box[d] += 1;
                outFlag = 1;
            }
            else if (atoms[n].reR_box[d] >= 1)
            {
                atoms[n].reR_box[d] -= 1;
                outFlag = 1;
            }
        }
        if (outFlag == 1)
        {
            Mul_3_1(boxTranVecs, atoms[n].reR_box, atoms[n].r);
        }
    }
}



void PBC_dr_nonorthogonal(int i, int j, double dr[3]) //vector: 1->2
{
    int d;
    double reDr[3];
    for (d = 0; d < 3; d++)
    {
        reDr[d] = atoms[j].reR_box[d] - atoms[i].reR_box[d];
        if (reDr[d] < -0.5)
        {
            reDr[d] += 1;
        }
        else if (reDr[d] > 0.5)
        {
            reDr[d] -= 1;
        }
    }
    Mul_3_1(boxTranVecs, reDr, dr);
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
                atoms[n].r_last_Verlet[d] += boxTranVecs[d][d];
            }
            else if (atoms[n].r[d] >= boxStartPoint[d] + boxTranVecs[d][d])
            {
                atoms[n].r[d] -= boxTranVecs[d][d];
                atoms[n].r_last_Verlet[d] -= boxTranVecs[d][d];
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

void DeleteAtomByIndex(int index)
{
    int i;
    for (i = index; i < atomNumber; i++)
    {
        atoms[i] = atoms[i + 1];
    }
    atomNumber--;
}

void DeleteAtomsByRegion(double block[3][2])
{
    int i, d;
    int dFlag[3];
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            if (atoms[i].r[d] > block[d][0] && atoms[i].r[d] < block[d][1])
            {
                dFlag[d] = 1;
            }
            else
            {
                dFlag[d] = 0;
            }
        }
        if (dFlag[0] && dFlag[1] && dFlag[2])
        {
            DeleteAtomByIndex(i);
            i--;
        }
    }
}


void CreateAtom(double r[3], int type)
{
    int d;
    for (d = 0; d < 3; d++)
    {
        atoms[atomNumber].r[d] = r[d];
    }
    atoms[atomNumber].type = type;
    maxID++;
    atoms[atomNumber].id = maxID;
    atomNumber++;
}


void FindNeighbors()
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
            distance = sqrt(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2));
            if (distance < neighborCutoff)
            {
                neighborIndex_i = atoms[i].neighborNumber;
                neighborIndex_j = atoms[j].neighborNumber;
                atoms[i].neighborList[neighborIndex_i].index = j;
                atoms[j].neighborList[neighborIndex_j].index = i;
                atoms[i].neighborList[neighborIndex_i].distance = distance;
                atoms[j].neighborList[neighborIndex_j].distance = distance;
                for (d = 0; d < 3; d++)
                {
                    atoms[i].neighborList[neighborIndex_i].dr[d] = dr[d];
                    atoms[j].neighborList[neighborIndex_j].dr[d] = -dr[d];
                }
                atoms[i].neighborNumber++;
                atoms[j].neighborNumber++;
            }
        }
    }
}


void InitLJ()
{
    LJ_4_epsilon = 4 * LJ_epsilon;
    LJ_24_epsilon_sigma6 = 24 * LJ_epsilon * pow(LJ_sigma, 6);
    LJ_2_sigma6 = 2 * pow(LJ_sigma, 6);
}


void Potential_LJ(int energyFlag, int forceFlag)
{
    int i, j, d;
    double distance;
    double energy;
    double force_inner;

    if (energyFlag) totalPotentialEnergy = 0;
    for (i = 0; i < atomNumber; i++)
    {
        if (energyFlag) atoms[i].potentialEnergy = 0;
        if (forceFlag)
        {
            atoms[i].force[0] = 0;
            atoms[i].force[1] = 0;
            atoms[i].force[2] = 0;
        }
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            distance = atoms[i].neighborList[j].distance;
            if (energyFlag)
            {
                energy = LJ_4_epsilon * (pow((LJ_sigma / distance), 12.) - pow((LJ_sigma / distance), 6.));
                atoms[i].potentialEnergy += energy / 2;
                totalPotentialEnergy += energy / 2;
            }
            if (forceFlag)
            {
                force_inner = LJ_2_sigma6 / pow(distance, 14) - 1 / pow(distance, 8);
                for (d = 0; d < 3; d++)
                {
                    atoms[i].force[d] -= force_inner * atoms[i].neighborList[j].dr[d];
                }
            }
        }
        if (forceFlag)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].force[d] *= LJ_24_epsilon_sigma6;
            }
        }
    }
}


void Potential_EAM(int energyFlag, int forceFlag)
{
    int i, j;
    double distance;
    double energyPhi, energyRho;
    double forcePhi, forceRho, forceInner, forceRho_ai;
    int n;
    int d;
    int jAtomIndex;
    // compute atom rho
    for (i = 0; i < atomNumber; i++)
    {
        atoms[i].rho_EAM = 0.0;
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            distance = atoms[i].neighborList[j].distance;
            if (distance > 2.002970124727)
            {

                for (n = 0; n < n_rho_EAM; n++)
                {
                    if (distance < delta_rho_EAM[n])
                    {
                        atoms[i].rho_EAM += a_rho_EAM[n] * pow((delta_rho_EAM[n] - distance), 3);
                    }
                }
            }
            else
            {
                atoms[i].rho_EAM += 1.193547157792;
            }
        }
    }

    if (energyFlag)
    {
        totalPotentialEnergy = 0.0;
    }
    for (i = 0; i < atomNumber; i++)
    {
        if (energyFlag)
        {
            //rho component for energy
            energyRho = a1_f_EAM * pow(atoms[i].rho_EAM, 0.5) + a2_f_EAM * pow(atoms[i].rho_EAM, 2);
            totalPotentialEnergy += energyRho;
            atoms[i].potentialEnergy = energyRho;
            //phi component for energy
        }
        if (forceFlag)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].force[d] = 0;
            }
            forceRho_ai = 0.5 * a1_f_EAM * pow(atoms[i].rho_EAM, -0.5) + 2 * a2_f_EAM * atoms[i].rho_EAM;
        }
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            distance = atoms[i].neighborList[j].distance;
            if (energyFlag)
            {
                energyPhi = 0.0;
                for (n = 0; n < n_phi_EAM; n++)
                {
                    if (distance < delta_phi_EAM[n])
                    {
                        energyPhi += a_phi_EAM[n] * pow((delta_phi_EAM[n] - distance), 3);
                    }
                }
                energyPhi *= 0.5;
                atoms[i].potentialEnergy += energyPhi;
                totalPotentialEnergy += energyPhi;
            }
            if (forceFlag)
            {
                //phi component for force
                forcePhi = 0.0;
                for (n = 0; n < n_phi_EAM; n++)
                {
                    if (distance < delta_phi_EAM[n])
                    {
                        forcePhi += -3 * a_phi_EAM[n] * pow(delta_phi_EAM[n] - distance, 2);
                    }
                }

                //rho component for force
                forceRho = 0.0;
                for (n = 0; n < n_rho_EAM; n++)
                {
                    if (distance < delta_rho_EAM[n])
                    {
                        forceRho += -3 * a_rho_EAM[n] * pow(delta_rho_EAM[n] - distance, 2);
                    }
                }

                // sum force
                jAtomIndex = atoms[i].neighborList[j].index;
                forceInner = ((((0.5 * a1_f_EAM * pow(atoms[jAtomIndex].rho_EAM, -0.5) + 2 * a2_f_EAM * atoms[jAtomIndex].rho_EAM)
                                + forceRho_ai)
                               * forceRho) + forcePhi) / distance;
                for (d = 0; d < 3; d++)
                {
                    atoms[i].force[d] += forceInner * atoms[i].neighborList[j].dr[d];
                }
            }
        }
    }
}


void Potential(int energyFlag, int forceFlag)
{
    FindNeighbors();
    if (strcmp(potentialName, "LJ") == 0) Potential_LJ(energyFlag, forceFlag);
    else if (strcmp(potentialName, "EAM") == 0) Potential_EAM(energyFlag, forceFlag);
    else printf("\nPotential not found. Program terminated.\n"), exit(-1);
}


void Minimize()
{
    double potentialEnergy_begin;
    int step;
    step = 0;
    printf("\n---Minimization start---\n");
    printf("iter pe dE\n");
    Potential(1, 0);
    do
    {
        potentialEnergy_begin = totalPotentialEnergy;
        MinDirection(step);
        LineMinimize();
        Potential(1, 0);
        step += 1;
        printf("%d %f %f %f\n", step, potentialEnergy_begin, totalPotentialEnergy, fabs(potentialEnergy_begin - totalPotentialEnergy));
    } while (fabs(potentialEnergy_begin - totalPotentialEnergy) > energyCriterion_min);
    printf("---Minimization end---\n");
}


void MinDirection(int step)
{
    if (strcmp(minimizeStyle, "SD") == 0) MinDirection_SD();
    else if (strcmp(minimizeStyle, "CG") == 0) MinDirection_CG(step);
    else printf("\nMinimize style not found. Program terminated.\n"), exit(-1);
}


void MinDirection_SD()
{
    int i, d;
    Potential(0, 1);
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {   atoms[i].minDirection[d] = atoms[i].force[d];
        }
    }
}   

void MinDirection_CG(int step)
{
    int i, d;
    double beta_up, beta_down, beta;
    double directionCheck;
    if (step == 0)
    {    Potential(0, 1);

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
                    atoms[i].minDirection[d] *= -1;
                }
            }
        }
    }
}

double delta_lineMin = 0.02;

double LineMinimize()
{
    int i, d;
    double alpha_new, alpha_1, alpha_2; //alpha = delta_lineMin*alpha_1/(alpha_2-alpha_1)
    alpha_1 = 0;
    alpha_2 = 0;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            alpha_1 += atoms[i].force[d] * atoms[i].minDirection[d];
        }
    }
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].r[d] += delta_lineMin * atoms[i].minDirection[d];
        }
    }
    PBC_r();
    Potential(0, 1);
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            alpha_2 += atoms[i].force[d] * atoms[i].minDirection[d];
        }
    }
    alpha_new = -delta_lineMin * alpha_1 / (alpha_2 - alpha_1) - delta_lineMin;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].r[d] += alpha_new * atoms[i].minDirection[d];
        }
    }
    PBC_r();
}

/*
double LineMinimize_Backtrace()
{
    int i, d;
    double startEnergy, acceptableEnergy;
    double lambda;
    int debug = 0;

    Potential(1, 0);
    startEnergy = totalPotentialEnergy;
    lambda = lambda_min;

    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].lineMinStartR[d] = atoms[i].r[d];
        }
    }
    do
    {
        acceptableEnergy = startEnergy;
        for (i = 0; i < atomNumber; i++)
        {
            for (d = 0; d < 3; d++)
            {
                atoms[i].r[d] = atoms[i].lineMinStartR[d] + atoms[i].minDirection[d] * lambda;
                acceptableEnergy -= atoms[i].minDirection[d] * atoms[i].minDirection[d] * lambda * c_min;
            }
        }
        Potential(1, 1);

        strcpy(dumpName, "debug.dump");
        Dump(debug);
        debug ++;

        lambda *= rho_min;
    } while (acceptableEnergy <= totalPotentialEnergy);
    return startEnergy - totalPotentialEnergy;
}
*/

void Dynamics()
{
    int step;
    double time;
    step = 0;
    time = 0;
    LaunchRun();
    printf("\n---Dynamics start---\n");
    while (time < totalTime)
    {
        Potential(0, 1);
        DynamicsProcessing(step, time);
        IterRun();
        time += timeStep;
        step++;
    }
    printf("---Dynamics end---\n");
}


void IterRun()
{
    if (strcmp(dynamicStyle, "Euler") == 0) IterRun_Euler();
    else if (strcmp(dynamicStyle, "Verlet") == 0) IterRun_Verlet();
    else printf("\nDynamic style not found. Program terminated.\n"), exit(-1);
    PBC_r();
}


void LaunchRun()
{
    if (strcmp(dynamicStyle, "Euler") == 0) 1;
    else if (strcmp(dynamicStyle, "Verlet") == 0) LaunchRun_Verlet();
    else printf("\nDynamic style not found. Program terminated.\n"), exit(-1); //fix: exit programe
    PBC_r();
}



void IterRun_Euler()
{
    int i, d;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].a[d] = atoms[i].force[d] / atoms[i].m;
            atoms[i].r[d] += atoms[i].v[d] * timeStep;
            atoms[i].v[d] += atoms[i].a[d] * timeStep;
        }
    }
}


void IterRun_Verlet()
{
    int i, d;
    double rTmp;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].a[d] = atoms[i].force[d] / atoms[i].m;
            rTmp = atoms[i].r[d];
            atoms[i].r[d] = 2 * atoms[i].r[d] - atoms[i].r_last_Verlet[d] + atoms[i].a[d] * time_2_Verlet;
            atoms[i].v[d] = (atoms[i].r[d] - atoms[i].r_last_Verlet[d]) / time_m_2_Verlet;
            atoms[i].r_last_Verlet[d] = rTmp;
        }
    }
}


void LaunchRun_Verlet()
{
    int i, d;
    time_2_Verlet = timeStep * timeStep;
    time_m_2_Verlet = timeStep * 2.0;
    Potential(0, 1);
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].r_last_Verlet[d] = atoms[i].r[d];
            atoms[i].a[d] = atoms[i].force[d] / atoms[i].m;
            atoms[i].r[d] += atoms[i].v[d] * timeStep + 0.5 * time_2_Verlet * atoms[i].a[d];
        }
    }
}


void Dump(int step) //only suitable for orthogonal box
{
    int i, d;
    FILE *file;
    if (step != 0)
    {
        file = fopen(dumpName, "a");
    }
    else
    {
        file = fopen(dumpName, "w");
    }
    fprintf(file, "ITEM: TIMESTEP\n");
    fprintf(file, "%d\n", step);
    fprintf(file, "ITEM: NUMBER OF ATOMS\n");
    fprintf(file, "%d\n", atomNumber);
    fprintf(file, "ITEM: BOX BOUNDS pp pp pp\n");
    for (d = 0; d < 3; d++)
    {
        fprintf(file, "%f %f\n", boxStartPoint[d], boxStartPoint[d] + boxTranVecs[d][d]);
    }
    fprintf(file, "ITEM: ATOMS id type x y z fx fy fz pe mx my mz\n");
    {
        for (i = 0; i < atomNumber; i++)
        {
            fprintf(file, "%d %d %f %f %f %f %f %f %f %f %f %f\n", atoms[i].id, atoms[i].type, atoms[i].r[0], atoms[i].r[1], atoms[i].r[2],
                    atoms[i].force[0], atoms[i].force[1], atoms[i].force[2], atoms[i].potentialEnergy,
                    atoms[i].minDirection[0], atoms[i].minDirection[1], atoms[i].minDirection[2]);
        }
    }
    fclose(file);
}


void DistributeVelocity(double temperature) // need to check correctness
{
    double randomNumber;
    int i;
    double speed;
    double prefactor1Base;
    double prefactor2Base;
    double prefactor1;
    double prefactor2;

    prefactor1Base = 4.0 * PI * pow(1.0 / 2.0 / PI / K_B / temperature, 3.0 / 2.0);
    prefactor2Base = -1.0 / 2.0 / K_B / temperature;


    for (i = 0; i < atomNumber; i++)
    {
        randomNumber = rand() / (RAND_MAX + 1.0);
        prefactor1 = prefactor1Base * pow(atoms[i].m, 3.0 / 2.0);
        prefactor2 = prefactor2Base * atoms[i].m;
        speed = GenerateSpeed(randomNumber, temperature, prefactor1, prefactor2);
        RandomSpherePoint(speed, atoms[i].v);
    }
    ZeroMomentum();
    //ZeroRotate();
}


double GenerateSpeed(double randomNumber, double temperature, double prefactor1, double prefactor2)
{
    double sigma;
    double x;

    x = 0.0;
    sigma = 0.0;
    while (sigma <= randomNumber)
    {
        sigma += MaxwellProbabilityDensity(x, prefactor1, prefactor2) * 0.01;
        x += 0.01;
    }
    return x;
}


double MaxwellProbabilityDensity(double v, double prefactor1, double prefactor2)
{
    return prefactor1 * pow(EBASE, pow(v, 2) * prefactor2) * pow(v, 2);
}


void RandomSpherePoint(double radius, double vector[3])
{
    double u, phi;
    u = (rand() / (RAND_MAX + 1.0) - 0.5) * 2;
    phi = (rand() / (RAND_MAX + 1.0)) * 2 * PI;

    vector[0] = sqrt(1 - u * u) * cos(phi) * radius;
    vector[1] = sqrt(1 - u * u) * sin(phi) * radius;
    vector[2] = u * radius;
}


void ZeroMomentum()
{
    int i, d;
    double totalMomentum[3];
    double deltaMomentumPerAtom[3];

    totalMomentum[0] = 0, totalMomentum[1] = 0, totalMomentum[2] = 0;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            totalMomentum[d] += atoms[i].v[d] * atoms[i].m;
        }
    }

    for (d = 0; d < 3; d++) deltaMomentumPerAtom[d] = -totalMomentum[d] / atomNumber;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].v[d] += deltaMomentumPerAtom[d] / atoms[i].m;
        }
    }
}


double ComputeTemperature()
{
    ComputeTotalKineticEnergy();
    return totalKineticEnergy * 2.0 / 3.0 / atomNumber / K_B;
}


void ComputeTotalKineticEnergy()
{
    int i, d;
    totalKineticEnergy = 0.0;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            totalKineticEnergy += 0.5 * atoms[i].m * atoms[i].v[d] * atoms[i].v[d];
        }
    }
}


double ComputeVolume()
{
    double s[3];
    VecCroMul(boxTranVecs[0], boxTranVecs[1], s);
    return VecDotMul(s, boxTranVecs[2]);
}


double ComputeStress(double stress[6]) //need to check correctness. need to fix the unit.
{
    int i, d;
    double ke[6];
    double virial[6];
    double volume;
    for (d = 0; d < 6; d++)
    {
        stress[d] = 0;
    }
    for (i = 0; i < atomNumber; i++)   // need to be divided by dimension
    {
        ke[0] = atoms[i].m * atoms[i].v[0] * atoms[i].v[0];
        ke[1] = atoms[i].m * atoms[i].v[1] * atoms[i].v[1];
        ke[2] = atoms[i].m * atoms[i].v[2] * atoms[i].v[2];
        ke[3] = atoms[i].m * atoms[i].v[0] * atoms[i].v[1];
        ke[4] = atoms[i].m * atoms[i].v[0] * atoms[i].v[2];
        ke[5] = atoms[i].m * atoms[i].v[1] * atoms[i].v[2];
        VirialPairs(i, virial);
        for (d = 0; d < 6; d++)
        {
            stress[d] -= ke[d] + virial[d];
        }
    }
    volume = ComputeVolume();
    for (d = 0; d < 6; d++)
    {
        stress[d] /= volume;
    }
    return  -(stress[0] + stress[1] + stress[2]) / 3;
}


void VirialPairs(int i, double virial[6])
{
    if (strcmp(potentialName, "LJ") == 0) VirialPairs_LJ(i, virial);
    else if (strcmp(potentialName, "EAM") == 0) VirialPairs_EAM(i, virial);
    else printf("\nPotential not found. Program terminated.\n"), exit(-1);
}


void VirialPairs_EAM(int i, double virial[6])
{
    int j, d, n;
    double force[3];
    double distance;
    double forcePhi, forceRho, forceInner, forceRho_ai;
    int jAtomIndex;
    for (d = 0; d < 6; d++)
    {
        virial[d] = 0;
    }
    forceRho_ai = 0.5 * a1_f_EAM * pow(atoms[i].rho_EAM, -0.5) + 2 * a2_f_EAM * atoms[i].rho_EAM;
    for (j = 0; j < atoms[i].neighborNumber; j++)
    {
        distance = atoms[i].neighborList[j].distance;
        //phi component for force
        forcePhi = 0;
        for (n = 0; n < n_phi_EAM; n++)
        {
            if (distance < delta_phi_EAM[n])
            {
                forcePhi += -3 * a_phi_EAM[n] * pow(delta_phi_EAM[n] - distance, 2);
            }
        }

        //rho component for force
        forceRho = 0;
        for (n = 0; n < n_rho_EAM; n++)
        {
            if (distance < delta_rho_EAM[n])
            {
                forceRho += -3 * a_rho_EAM[n] * pow(delta_rho_EAM[n] - distance, 2);
            }
        }

        // sum force
        jAtomIndex = atoms[i].neighborList[j].index;
        forceInner = ((((0.5 * a1_f_EAM * pow(atoms[jAtomIndex].rho_EAM, -0.5) + 2 * a2_f_EAM * atoms[jAtomIndex].rho_EAM)
                        + forceRho_ai)
                       * forceRho) + forcePhi) / distance;
        for (d = 0; d < 3; d++)
        {
            force[d] = forceInner * atoms[i].neighborList[j].dr[d];
        }
        virial[0] -= atoms[i].neighborList[j].dr[0] * force[0];
        virial[1] -= atoms[i].neighborList[j].dr[1] * force[1];
        virial[2] -= atoms[i].neighborList[j].dr[2] * force[2];
        virial[3] -= atoms[i].neighborList[j].dr[0] * force[1];
        virial[4] -= atoms[i].neighborList[j].dr[0] * force[2];
        virial[5] -= atoms[i].neighborList[j].dr[1] * force[2];
    }
}


void VirialPairs_LJ(int i, double virial[6])
{
    int j;
    int d;
    double forceInner;
    double force[3];
    double distance;
    for (d = 0; d < 6; d++)
    {
        virial[d] = 0;
    }
    for (j = 0; j < atoms[i].neighborNumber; j++)
    {
        distance = atoms[i].neighborList[j].distance;
        forceInner = LJ_2_sigma6 / pow(distance, 14) - 1 / pow(distance, 8);
        for (d = 0; d < 3; d++)
        {
            force[d] = forceInner * atoms[i].neighborList[j].dr[d];
        }
        virial[0] += atoms[i].neighborList[j].dr[0] * force[0];
        virial[1] += atoms[i].neighborList[j].dr[1] * force[1];
        virial[2] += atoms[i].neighborList[j].dr[2] * force[2];
        virial[3] += atoms[i].neighborList[j].dr[0] * force[1];
        virial[4] += atoms[i].neighborList[j].dr[0] * force[2];
        virial[5] += atoms[i].neighborList[j].dr[1] * force[2];
    }
    for (d = 0; d < 6; d++)
    {
        virial[d] *= LJ_24_epsilon_sigma6;
    }
}

void Thermostat_NoseHoover(int step, double targetTemperature)
{
    double temperature;
    if (step == 0)
    {
        xi_NoseHoover = 0;
        tau2_r_NoseHoover = 1.0 / tau_NoseHoover / tau_NoseHoover;
    }
    else
    {
        temperature = ComputeTemperature();
        xi_NoseHoover += tau2_r_NoseHoover * (temperature / targetTemperature - 1) * timeStep;
    }
    Friction_NoseHoover();
}

void Friction_NoseHoover()
{
    int i, d;
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].force[d] -= atoms[i].m * xi_NoseHoover * atoms[i].v[d];
        }
    }
}

void Thermostat_Berendsen(double targetTemperature)
{
    double lambda;
    double temperature;
    int i, d;
    temperature = ComputeTemperature();
    lambda = pow(1 + timeStep / tau_Berendsen * (targetTemperature / temperature - 1), 0.5);
    for (i = 0; i < atomNumber; i++)
    {
        for (d = 0; d < 3; d++)
        {
            atoms[i].v[d] *= lambda;
        }
    }
}

void Thermostat(int step, double targetTemperature)
{
    if (strcmp(thermostatStyle, "NoseHoover") == 0) Thermostat_NoseHoover(step, targetTemperature);
    else if (strcmp(thermostatStyle, "Berendsen") == 0) Thermostat_Berendsen(targetTemperature);
    else printf("\nThermo style not found. Program terminated.\n"), exit(-1);
}

void DynamicsProcessing(int step, double time)
{
    double temperature;
    double stress[6];
    if (step == 0)
    {
        printf("%d atoms are in the dynamic simulation.\n", atomNumber);
        printf("step time temperature stress[0] stress[1] stress[2] stress[3] stress[4] stress[5] energy\n");
    }
    if (step % 100 == 0)
    {
        temperature = ComputeTemperature();
        ComputeStress(stress);
        printf("%d %f %f %f %f %f %f %f %f %f \n",
               step, time, temperature, stress[0], stress[1], stress[2], stress[3], stress[4], stress[5], totalPotentialEnergy + totalKineticEnergy);
        Potential(1, 0);
        Dump(step);
    }
    Thermostat(step, 300);
}

void DisturbAtoms(double maxDistance)
{
    int i, d;
    double distance;
    double dr[3];
    for (i = 0; i < atomNumber; i++)
    {
        distance = (rand() / (RAND_MAX + 1.0)) * maxDistance;
        RandomSpherePoint(distance, dr);
        for (d = 0; d < 3; d++)
        {
            atoms[i].r[d] += dr[d];
        }
    }
}

/*
void MC_Metropolis(int totalSteps)
{
    int step;
    struct Atom storedAtoms[MAX_ATOM_NUMBER];
    double storedEnergy;
    int acceptFlag;
    storedEnergy = Potential(1,0);
    for (step = 0; step < totalSteps; step++)
    {
        newEnergy = Attempt_MC(storedAtoms);
        MetropolisCriterion(storedEnergy, newEnergy, temperature)
    }
}


void MetropolisCriterion(double storedEnergy, double newEnergy, double temperature)
{
    int randomNumber;
    randomNumber = rand() / (RAND_MAX + 1);
    if randomNumber <= pow(EBASE, (storedEnergy - newEnergy) / K_B / temperature)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void Recover(struct Atom storedAtoms[MAX_ATOM_NUMBER], int storedAtomNumber)
{
    int n;
    for (n = 0; n < storedAtomNumber; n++)
    {
        atoms[n] = storedAtoms[n];
    }
    atomNumber = storedAtomNumber;
}

void Attempt_MC(storedAtoms)
// todo list: 压强控制 径向分布函数  弛豫曲线 等概率球形模型
// 通过noose hoover 介绍系综原理
// 介绍数值积分方法
*/


/*main*/
int main() //
{
    double latticeConstant; //unit: Angstrom
    latticeConstant = 3.14;
    //box and atoms
    latticeSizes[0][0] = 0;  latticeSizes[0][1] = 5;
    latticeSizes[1][0] = 0;  latticeSizes[1][1] = 5;
    latticeSizes[2][0] = 0;  latticeSizes[2][1] = 5;

    priTranVecs[0][0] = latticeConstant; priTranVecs[1][0] = 0; priTranVecs[2][0] = 0;
    priTranVecs[0][1] = 0; priTranVecs[1][1] = latticeConstant; priTranVecs[2][1] = 0;
    priTranVecs[0][2] = 0; priTranVecs[1][2] = 0; priTranVecs[2][2] = latticeConstant;
    boxOrthogonalFlag = 1;
    if (!boxOrthogonalFlag)
    {
        ComputeBoxRecTranVecs();
    }

    cellAtomNumber = 2;
    cellAtomRs[0][0] = 0; cellAtomRs[0][1] = 0; cellAtomRs[0][2] = 0;
    cellAtomRs[1][0] = 0.5 * latticeConstant; cellAtomRs[1][1] = 0.5 * latticeConstant; cellAtomRs[1][2] = 0.5 * latticeConstant;

    cellAtomTypes[0] = 0; //start from 0
    cellAtomTypes[1] = 0;

    typeMasses[0] = 183;  // relative atomic mass


    boxStartPoint[0] = 0; boxStartPoint[1] = 0; boxStartPoint[2] = 0;

    boxTranVecs[0][0] = latticeConstant * 5; boxTranVecs[1][0] = 0;                   boxTranVecs[2][0] = 0;
    boxTranVecs[0][1] = 0;                   boxTranVecs[1][1] = latticeConstant * 5; boxTranVecs[2][1] = 0;
    boxTranVecs[0][2] = 0;                   boxTranVecs[1][2] = 0;                   boxTranVecs[2][2] = latticeConstant * 5;


    ConstructCrystal();
    double block[3][2] = { { 7.8, 7.9}, { 7.8, 7.9}, { 7.8, 7.9} };
    //DeleteAtomsByRegion(block);
    DisturbAtoms(1.0);

    //process
    //--------------minimize--------------
    neighborCutoff = 5.5;

    strcpy(potentialName, "EAM");
    //InitLJ();
    strcpy(minimizeStyle, "CG");
    energyCriterion_min = 1e-9;
    lambda_min = 1;
    c_min = 0.5;
    rho_min = 0.5;
       
    strcpy(dumpName, "min.dump");
    Dump(0);

    //-------------dynamic----------------
    srand(1);
    strcpy(dynamicStyle, "Euler");
    strcpy(thermostatStyle, "Berendsen");
    tau_NoseHoover = 0.01;
    tau_Berendsen = 0.01;
    totalTime = 0;
    timeStep = 0.0001;
    DistributeVelocity(2000);
    strcpy(dumpName, "force.dump");
    Dynamics();

    //FindNeighbors();
    //Potential(1, 1);
    //Dump(0);
    //printf("Pe: %f\n", totalPotentialEnergy);
    /*
    double r0[3] = {0,0,0};
    double r1[3] = {0,0,3};
    CreateAtom(r0,1);
    CreateAtom(r1,1);

    FindNeighbors();
    neighborCutoff = 4 * latticeConstant + 0.1;
    InitLJ();
    Potential_LJ(1, 1);
    */
}

