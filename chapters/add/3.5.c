/* constants */
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

/* function declarations */
void Potential_EAM(int isEnergy, int isForce);

/* functions */
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

/* main */
int main()
{
    /* parameters */
    neighborCutoff = 6.0;
    
    /* processing and output*/
    double latticeConstant;
    printf("lattice_constant pe(eV/atom)\n");
    for (latticeConstant = 2.9; latticeConstant < 3.4; latticeConstant += 0.01)
    {
        ConstructStdCrystal_BCC(latticeConstant, 7);
        ConstructNeighborList();
        Potential_EAM(1, 0);
        printf("%f %f\n", latticeConstant, totalPotentialEnergy / atomNumber);
    }

    return 0;
}