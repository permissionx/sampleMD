/*include*/
#include <stdio.h>
#include <math.h>

/*constant*/
#define MAX_LATTICE_NUMBER 2000
#define MAX_ATOM_NUMBER 20000
#define MAX_CELL_ATOM_NUMBER 10
#define MAX_NEIGHBOR_NUMBER 2000

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
    double distance;
    double dr[3];
};

struct Atom
{
    double r[3];
    double reR_box[3];
    int id, type;
    double potentialEnergy;
    double force[3];
    int neighborNumber;
    struct Neighbor neighborList[MAX_NEIGHBOR_NUMBER]; // error if exceed
    double minDirection[3];
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
double boxStartPoint[3];
double boxTranVecs[3][3]; // box translation vectors
double boxRecTranVecs[3][3]; //box reciprocal translation vectors
double boxRecTranVecs_inv[3][3];
int boxOrthogonalFlag;

char interactionName[20];

double cutoff;
double LJ_sigma;
double LJ_epsilon;
double LJ_24_epsilon_sigma6;
double LJ_2_sigma6;
double LJ_4_epsilon;
double totalPotentialEnergy;

char minimizeStyle[20];


/*function declaration*/
void ConstructReducedLattice();
void Mul_3_1(double a[3][3], double b[3], double p[3]); // matrix multiplication 3x1
void ConstructLattice();
void ConstructCrystal();
void DumpSingle(char fileName[20]);
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
void Interaction_LJ(int energyFlag, int forceFlag);
void Interaction_EAM(int energyFlag, int forceFlag);

void Minimize(double energyCriterion);
void MinimizeSD(double energyCriterion);
void LineMinimize();

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

    for (nLattice = 0; nLattice < latticePointNumber; nLattice++)
    {
        for (nCellAtom = 0; nCellAtom < cellAtomNumber; nCellAtom++)
        {
            atoms[nAtom].r[0] = latticePoints[nLattice].r[0] + cellAtomRs[nCellAtom][0];
            atoms[nAtom].r[1] = latticePoints[nLattice].r[1] + cellAtomRs[nCellAtom][1];
            atoms[nAtom].r[2] = latticePoints[nLattice].r[2] + cellAtomRs[nCellAtom][2];

            atoms[nAtom].type = cellAtomTypes[nCellAtom];
            atoms[nAtom].id = nAtom;
            nAtom++;
        }
    }
    maxID = nAtom - 1;
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
        fprintf(file, "%d %d %f %f %f %f %f %f\n", atoms[i].id, atoms[i].type, atoms[i].r[0], atoms[i].r[1], atoms[i].r[2], atoms[i].force[0], atoms[i].force[1], atoms[i].force[2]);
    }
    fclose(file);
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
    int i, d;
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
            if (distance < cutoff)
            {
                neighborIndex_i = atoms[i].neighborNumber;
                neighborIndex_j = atoms[j].neighborNumber;
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

void Interaction_LJ(int energyFlag, int forceFlag)
{
    int i, j, d;
    double distance, distance8, distance14, dr[3];
    double energy;
    double force_temp;

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
                distance14 = pow(distance, 14);
                distance8 = pow(distance, 8);
                for (d = 0; d < 3; d++)
                {
                    dr[d] = atoms[i].neighborList[j].dr[d];
                    force_temp = LJ_24_epsilon_sigma6 * (LJ_2_sigma6 * dr[d] / distance14 - dr[d] / distance8);
                    atoms[i].force[d] += force_temp;
                }
            }
        }
    }
}

void Interaction_EAM(int energyFlag, int forceFlag)
{
    int i, j;
    double distance;
    double rho;
    double atomicEnergy;
    double dr_distance; // dr[d]/distance
    double forcePhiTmp[3], forceRhoTmp[3], forceRhoTmp_inner[3], forcePhiTmp_inner[3];
    int n;
    int d;
    if (energyFlag) totalPotentialEnergy = 0;
    for (i = 0; i < atomNumber; i++)
    {
        rho = 0.0;
        if (energyFlag) atomicEnergy = 0;
        if (forceFlag)
        {
            for (d = 0; d < 3; d++)
            {
                forcePhiTmp[d] = 0;
                forceRhoTmp[d] = 0;
                atoms[i].force[d] = 0;
            }
        }
        for (j = 0; j < atoms[i].neighborNumber; j++)
        {
            distance = atoms[i].neighborList[j].distance;
            for (n = 0; n < n_rho_EAM; n++)
            {
                if (distance < delta_rho_EAM[n])
                {
                    rho += a_rho_EAM[n] * pow((delta_rho_EAM[n] - distance), 3);
                }
            }
            //phi component for energy
            if (energyFlag)
            {
                for (n = 0; n < n_phi_EAM; n++)
                {
                    if (distance < delta_phi_EAM[n])
                    {
                        atomicEnergy += 0.5 * a_phi_EAM[n] * pow((delta_phi_EAM[n] - distance), 3);
                    }
                }
            }

            if (forceFlag)
            {
                for (d = 0; d < 3; d++)
                {
                    forcePhiTmp_inner[d] = 0;
                    forceRhoTmp_inner[d] = 0;
                    //phi component for force
                    for (n = 0; n < n_phi_EAM; n++)
                    {
                        if (distance < delta_phi_EAM[n])
                        {
                            forcePhiTmp_inner[d] += -3 * a_phi_EAM[n] * pow(delta_phi_EAM[n] - distance, 2);
                        }
                    }

                    //rho component for force
                    for (n = 0; n < n_rho_EAM; n++)
                    {
                        if (distance < delta_rho_EAM[n])
                        {
                            forceRhoTmp_inner[d] += -3 * a_rho_EAM[n] * pow(delta_rho_EAM[n] - distance, 2);
                        }
                    }
                    dr_distance = atoms[i].neighborList[j].dr[d];
                    forcePhiTmp[d] += dr_distance * forcePhiTmp_inner[d];
                    forceRhoTmp[d] += dr_distance * forceRhoTmp_inner[d];
                }
            }
        }

        //rho component for energy and total energy
        if (energyFlag)
        {
            atomicEnergy += a1_f_EAM * pow(rho, 0.5) + a2_f_EAM * pow(rho, 2);
            atoms[i].potentialEnergy = atomicEnergy;
            totalPotentialEnergy += atomicEnergy;
        }

        // total force
        if (forceFlag)
        {
            for (d = 0; d < 3; d++)
            {
                forceRhoTmp[d] *= a1_f_EAM * pow(rho, -0.5) + 4 * a2_f_EAM * rho;
                atoms[i].force[d] = -(forcePhiTmp[d] + forceRhoTmp[d]);
            }
        }
    }
}

void Interaction(int energyFlag, int forceFlag)
{
    if (interactionName == "LJ") Interaction_LJ(energyFlag, forceFlag);
    else if (interactionName == "EAM") Interaction_EAM(energyFlag, forceFlag);
    else printf("Potential not found.\n");
}

void minimizeSD(double energyCriterion)
{
    int i, d;
    double deltaEnergy, lastEnergy;

    FindNeighbors();
    Interaction(0, 1);
    do
    {
        for (i = 0; i < atomNumber; i++)
        {
            for (d = 0; d < 3; d++) atoms[i].minDirection[d] = atoms[i].force[d];
        }
        deltaEnergy = LineMinimize();
    } while (deltaEnergy < energyCriterion);
}

void LineMinimize()
{
    int i,d;
    double startEnergy, endEnergy;
    double lambda;
    lambda = 100;
    Interaction(1, 0);
    startEnergy = totalPotentialEnergy;
    do
    {
        for(i=0;i<atomNumber;i++)
        {
            for (d=0;d<3;d++)
            {
                atoms[i].r[d] += atoms[i].minDirection[d]*lambda;
            }
        }
    }
}

/*main*/
int main() //
{
    double latticeConstant; //unit: Angstrom
    latticeConstant = 3.14;
    //parameter
    latticeSizes[0][0] = 0;  latticeSizes[0][1] = 0;
    latticeSizes[1][0] = 0;  latticeSizes[1][1] = 0;
    latticeSizes[2][0] = 0;  latticeSizes[2][1] = 0;

    priTranVecs[0][0] = latticeConstant; priTranVecs[1][0] = 0; priTranVecs[2][0] = 0;
    priTranVecs[0][1] = 0; priTranVecs[1][1] = latticeConstant; priTranVecs[2][1] = 0;
    priTranVecs[0][2] = 0; priTranVecs[1][2] = 0; priTranVecs[2][2] = latticeConstant;
    boxOrthogonalFlag = 1;

    cellAtomNumber = 2;
    cellAtomRs[0][0] = 0; cellAtomRs[0][1] = 0; cellAtomRs[0][2] = 0;
    cellAtomRs[1][0] = 0.5 * latticeConstant; cellAtomRs[1][1] = 0.5 * latticeConstant; cellAtomRs[1][2] = 0.5 * latticeConstant;
    cellAtomRs[2][0] = 0.5 * latticeConstant; cellAtomRs[2][1] = 0; cellAtomRs[2][2] = 0.5 * latticeConstant;
    cellAtomRs[3][0] = 0.5 * latticeConstant; cellAtomRs[3][1] = 0.5 * latticeConstant; cellAtomRs[3][2] = 0;

    cellAtomTypes[0] = 1;
    cellAtomTypes[1] = 1;
    cellAtomTypes[2] = 1;
    cellAtomTypes[3] = 1;

    boxStartPoint[0] = 0; boxStartPoint[1] = 0; boxStartPoint[2] = 0;

    boxTranVecs[0][0] = latticeConstant * 10; boxTranVecs[1][0] = 0;                   boxTranVecs[2][0] = 0;
    boxTranVecs[0][1] = 0;                   boxTranVecs[1][1] = latticeConstant * 10; boxTranVecs[2][1] = 0;
    boxTranVecs[0][2] = 0;                   boxTranVecs[1][2] = 0;                   boxTranVecs[2][2] = latticeConstant * 10;


    cutoff = 4 * latticeConstant + 0.1;


    //process
    ConstructReducedLattice();
    ConstructLattice();
    ConstructCrystal();

    if (!boxOrthogonalFlag)
    {
        ComputeBoxRecTranVecs();
        ComputeAtomBoxReR();
    }

    //double block[3][2] = {{ -0.1, 0.1}, { -0.1, 0.1}, { -0.1, 0.1}};
    //DeleteAtomsByRegion(block);
    double r1[3] = { 0,0,0 };
    double r2[3] = { 1,1,1 };
    CreateAtom(r1, 1);
    CreateAtom(r2, 1);
    FindNeighbors();
    Interaction_EAM(1, 1);
    /*
    LJ_epsilon = 0.0031;
    LJ_sigma = 2.74;
    cutoff = 4 * latticeConstant + 0.1;
    InitLJ();
    Interaction_LJ(1, 1);
    */
    DumpSingle("force.xyz");
    printf("%f\n", totalPotentialEnergy);

}
