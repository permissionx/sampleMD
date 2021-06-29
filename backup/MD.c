#pragma warning(disable:4996)  
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MAX_ATOM_NUMBER 50000
#define MAX_LATTICE_NUMBER 10000
#define MAX_CELL_ATOM_NUMBER 20
#define MAX_TYPE_NUMBER 3
#define PI 3.1415926
#define k_B 8.6173324e-5
#define MAX_TEMP_NODES 20




struct LatPntT  // lattice point type，晶格格点类型
{
	double r[3], reR[3];  //格点的空间坐标和约化坐标
	int id;  //格点id
};

struct AtomT
{
	double r[3], reR[3], supCellReR[3], force[3], lastForce[3], minDirection[3], lastMinDirection[3], velocity[3], ft[3], dr[3];
	double backtrackingStartR[3];
	double energy;
	double mass;
	int id;
	int ttype;
};

struct AtomT atom[MAX_ATOM_NUMBER];
int nAtom;
int cellAtomN;                                         //number of atoms in a cell
double cellAtomR[MAX_CELL_ATOM_NUMBER][3];             //each atom relative position in a cell
int cellAtomT[MAX_CELL_ATOM_NUMBER];                   //each atom type in a cell
double typeMass[MAX_TYPE_NUMBER];

struct LatPntT latPnt[MAX_LATTICE_NUMBER];             //定义格点数组lattice point,最多存储10000个格点
double tranVec[3][3];                                  //primitive translation vector,三个初基平移矢量的三个分量
int nLatPnt;                                           //总格点数
double latSize[3];                                     //晶格模型在三个方向上的格点数
double lattice_constant = 0.0;						   //晶格常数
double rcpVec[3][3];                                   //rcpiprocal lattice vectors

double supCellVec[3][3];
double rcpSupCellVec[3][3];
//L-J
char* pairStyle;
double potentialEnergy;
double epsilons[MAX_TYPE_NUMBER][MAX_TYPE_NUMBER];
double sigmas[MAX_TYPE_NUMBER][MAX_TYPE_NUMBER];
double cutoff[MAX_TYPE_NUMBER][MAX_TYPE_NUMBER];

int iterTime;
double pressure; //into local 
double pressureAndvolume;
double volume = 1.0;
const double R_const = 8.31;
const double NA = 6.02214076e23;

int ReLat();
int MatMul(double a[3], double b[3][3], double c[3]);
int Lat();
int Crystal();
int Output(char filename[20], int timeStep);
void InitOut(char filename[20]);
int CRcpVec();
int CroPro(double vecOut[3], double vecIn1[3], double vecIn2[3]);
int DotPro(double valueOut[1], double vecIn1[3], double vecIn2[3]);
int CReR();
int MatInv(double out[3][3], double in[3][3]);
int CVec_supCell();
int CRcpVec_supCell();
int CR_supCell();
int PBC_r();
int PBC_dr(double dr[3], int index1, int index2);
int CreateAtom(double block[3], int type); // type to ttype
int DeleteAtoms(double block[3][2]);
int EnergyAndForce_2Body(int ifEnergy, int ifForce);
double Energy_LJ(double distance, double epsilon, double sigma);
int Force_LJ(double force[3], double dr[3], double distance, double epsilon, double sigma);
int Minimize_SD(double energyCrit, double alpha0, double rho, double c, char filename[20], int outSteps);
int Minimize_SD_EAM(double energyCrit, double alpha0, double rho, double c, char filename[20], int outSteps);
int Minimize_CG(double energyCrit, double alpha0, double rho, double c, char filename[20], int outSteps);
double LineSearch_Backtracking(double alpha0, double rho, double c, char filename[20], int outSteps);
double LineSearch_Backtracking_EAM(double alpha0, double rho, double c, char filename[20], int outSteps);
double ComputeBeta_CG();
int ComputeMinDirection_CG(double beta);
double ComputeTemperature();
void ComputeVolume();
void ControlTempeture(double temperature, double targetTempeture);
double GetTargetTempeture(double tempratureNodes[MAX_TEMP_NODES], double temperatureTimeNodes[MAX_TEMP_NODES], double time, double total_time);
int VelocityVerlet(double dt);
double EnergyAndForce_EAM_Cu(double cutoff);

//EAM 参数
const  int EAM_n_phi = 15;
const  int EAM_n_rho = 4;
double EAM_parameter_a_phi[15] =
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
double EAM_parameter_a_rho[4] =
{
-0.420429107805055e1,
0.518217702261442e0,
0.562720834534370e-1,
0.344164178842340e-1
};
double EAM_parameter_a_F[2] =
{
-5.946454472402710,
-0.049477376935239
};
double EAM_parameter_delta_phi[15] =
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
double EAM_parameter_delta_rho[4] =
{
 2.500000000000000,
 3.100000000000000,
 3.500000000000000,
 4.900000000000000
};
/*
double EAM_parameter_a_phi[15] =
{
 0.954071477542914e2,
-0.181161004118916e3,
 0.930215233132627e2,
-0.112027468539573e2,
 0.112027468539573e2 ,//* 3.95425,
-0.312459176640227e1,//*7.908515,
 0.123028140617302e1,
 0.154767467307454e1,
-0.128861387780439e1,
-0.843327493551467e0,
 0.214009882965042e1,
-0.102898314532388e1,
 0.138163259601912e1,
-0.360872433001551e1,
 0.217655968740690e1
};
double EAM_parameter_a_rho[4] =
{
-0.420429107805055e1,
0.518217702261442e0,
0.562720834534370e-1,
0.344164178842340e-1
};
double EAM_parameter_a_F[2] =
{
-5.553986589859130,
-0.045691157657292
};
double EAM_parameter_delta_phi[15] =
{
2.564897500000000,
2.629795000000000,
2.694692500000000,
2.866317500000000,
2.973045000000000,
3.079772500000000,
3.516472500000000,
3.846445000000000,
4.176417500000000,
4.700845000000000,
4.895300000000000,
5.089755000000000,
5.342952500000000,
5.401695000000000,
5.460437000000000
};
double EAM_parameter_delta_rho[4] =
{
 2.500000000000000,
 3.100000000000000,
 3.500000000000000,
 4.900000000000000
};
*/
int ReLat()                                            //reduced lattice
{
	int n, i, j, k;
	n = 0;
	for (i = 0; i < latSize[0]; i++)
	{
		for (j = 0; j < latSize[1]; j++)
		{
			for (k = 0; k < latSize[2]; k++)
			{
				latPnt[n].reR[0] = i;
				latPnt[n].reR[1] = j;
				latPnt[n].reR[2] = k;
				latPnt[n].id = n + 1;
				n++;
			}
		}
	}
	nLatPnt = n;
	return 0;
}
int MatMul(double p[3], double a[3][3], double b[3])
{
	int i;
	for (i = 0; i < 3; i++)
		p[i] = a[0][i] * b[0] + a[1][i] * b[1] + a[2][i] * b[2];
	return 0;
}
int Lat()
{
	int n;
	ReLat();
	for (n = 0; n < nLatPnt; n++)
		MatMul(latPnt[n].r, tranVec, latPnt[n].reR);
	return 0;
}
int Crystal()
{
	int	n, na, i, j;
	na = 0;

	Lat();
	for (n = 0; n < nLatPnt; n++)
		for (i = 0; i < cellAtomN; i++)
		{
			for (j = 0; j < 3; j++)
			{
				atom[na].r[j] = latPnt[n].r[j] + cellAtomR[i][j];
			}
			atom[na].ttype = cellAtomT[i];
			atom[na].id = na + 1;
			atom[na].mass = typeMass[atom[na].ttype];
			na++;
		}
	nAtom = na;
	return 0;
}

void InitOut(char filename[20])
{
	FILE *p;
	p = fopen(filename, "w");
	fclose(p);
}

int Output(char filename[20], int timeStep)
{
	int n;
	FILE *p;
	p = fopen(filename, "a");
	fprintf(p, "ITEM: TIMESTEP\n");
	fprintf(p, "%d\n", timeStep);
	fprintf(p, "ITEM: NUMBER OF ATOMS\n");
	fprintf(p, "%d\n", nAtom);
	fprintf(p, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(p, "0\t%f\n", latSize[0] * lattice_constant); // where lc comes from
	fprintf(p, "0\t%f\n", latSize[1] * lattice_constant);
	fprintf(p, "0\t%f\n", latSize[2] * lattice_constant);
	fprintf(p, "ITEM: ATOMS id type x y z force_x force_y force_z\n");
	for (n = 0; n < nAtom; n++)
	{
		fprintf(p, "%d %d %f %f %f %f %f %f %f\n", atom[n].id, atom[n].ttype, atom[n].r[0], atom[n].r[1], atom[n].r[2], atom[n].force[0], atom[n].force[1], atom[n].force[2], atom[n].energy);

	}
	fclose(p);
	return 0;
}

int CRcpVec()											//compute reciprocal lattice vectors
{
	int i, j, k, d;
	double cellVol[1];
	double tempVec[3];									//temporary vector

	CroPro(tempVec, tranVec[0], tranVec[1]);
	DotPro(cellVol, tempVec, tranVec[2]);

	for (i = 0; i < 3; i++)
	{
		j = i + 1;
		k = j + 1;
		if (j > 2)
		{
			j -= 3;
		}
		if (k > 2)
		{
			k -= 3;
		}
		CroPro(rcpVec[i], tranVec[j], tranVec[k]);
		for (d = 0; d < 3; d++)
		{
			rcpVec[i][d] *= 2 * PI / cellVol[0];
		}
	}
	return 0;
}
int CroPro(double vecOut[3], double vecIn1[3], double vecIn2[3])	//cross product
{
	vecOut[0] = vecIn1[1] * vecIn2[2] - vecIn1[2] * vecIn2[1];
	vecOut[1] = vecIn1[2] * vecIn2[0] - vecIn1[0] * vecIn2[2];
	vecOut[2] = vecIn1[0] * vecIn2[1] - vecIn1[1] * vecIn2[0];
	return 0;
}
int DotPro(double valueOut[1], double vecIn1[3], double vecIn2[3])		//dot product
{
	valueOut[0] = vecIn1[0] * vecIn2[0] + vecIn1[1] * vecIn2[1] + vecIn1[2] * vecIn2[2];
	return 0;
}
int CReR()													//compute reR
{
	int n, d;
	double rcpVec_inv[3][3];				//rcperted lattice vector matrix of inversion		

	CRcpVec();
	MatInv(rcpVec_inv, rcpVec);
	for (n = 0; n < nAtom; n++)
	{
		MatMul(atom[n].reR, rcpVec_inv, atom[n].r);
		for (d = 0; d < 3; d++)
		{
			atom[n].reR[d] /= 2 * PI;
		}
	}
	return 0;
}
int MatInv(double out[3][3], double in[3][3])							//Matrix Inversion
{
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			out[i][j] = in[j][i];
		}
	}
	return 0;
}

int CVec_supCell()
{
	int i, d;
	for (i = 0; i < 3; i++)
	{
		for (d = 0; d < 3; d++)
		{
			supCellVec[i][d] = latSize[i] * tranVec[i][d];
		}
	}
	return 0;
}
int CRcpVec_supCell()						//compute rcpiprical lattice vectors
{
	int i, d;
	double cellVol[1];
	double tempVec[3];						//temporary vector

	CVec_supCell();
	CroPro(tempVec, supCellVec[0], supCellVec[1]);
	DotPro(cellVol, tempVec, supCellVec[2]);

	for (i = 0; i < 3; i++)
	{

		CroPro(rcpSupCellVec[i], supCellVec[(i + 1) % 3], supCellVec[(i + 2) % 3]);
		for (d = 0; d < 3; d++)
		{
			rcpSupCellVec[i][d] *= 2 * PI / cellVol[0];
		}
	}
	return 0;
}
int CReR_supCell()						//compute reR
{
	int n, d;
	double rcpSupCellVec_inv[3][3];			//inverse rcpiprocal lattive vector matrix 

	CRcpVec_supCell();
	MatInv(rcpSupCellVec_inv, rcpSupCellVec);
	for (n = 0; n < nAtom; n++)
	{
		MatMul(atom[n].supCellReR, rcpSupCellVec_inv, atom[n].r);
		for (d = 0; d < 3; d++)
		{
			atom[n].supCellReR[d] /= 2 * PI;
		}
	}
	return 0;
}
int CR_supCell()								// compute R of atoms by supercell vector
{
	int n;

	for (n = 0; n < nAtom; n++)
	{
		MatMul(atom[n].r, supCellVec, atom[n].supCellReR);
	}
	return 0;
}
int PBC_r()						//Periodic boundary condition Positions 
{
	int n, d;
	int CReR_supCell();
	int CR_supCell();

	CReR_supCell();
	for (n = 0; n < nAtom; n++)
	{
		for (d = 0; d < 3; d++)
		{
			if (atom[n].supCellReR[d] < 0)
			{
				atom[n].supCellReR[d] += 1;
			}
			else if (atom[n].supCellReR[d] >= 1)
			{
				atom[n].supCellReR[d] -= 1;
			}
		}
	}
	CR_supCell();
	return 0;
}
int PBC_dr(double dr[3], int index1, int index2)			//距离矢量2->1
{
	int d;
	double reDr[3];

	for (d = 0; d < 3; d++)
	{
		if (atom[index1].supCellReR[d] - atom[index2].supCellReR[d] < -0.5)
		{
			reDr[d] = atom[index1].supCellReR[d] - (atom[index2].supCellReR[d] - 1);
		}
		else if (atom[index1].supCellReR[d] - atom[index2].supCellReR[d] > 0.5)
		{
			reDr[d] = atom[index1].supCellReR[d] - (atom[index2].supCellReR[d] + 1);
		}
		else
		{
			reDr[d] = atom[index1].supCellReR[d] - atom[index2].supCellReR[d];
		}
	}
	MatMul(dr, supCellVec, reDr);
	return 0;
}

int CreateAtom(double block[3], int type)
{
	atom[nAtom].r[0] = block[0];
	atom[nAtom].r[1] = block[1];
	atom[nAtom].r[2] = block[2];
	atom[nAtom].mass = typeMass[type];
	atom[nAtom].ttype = cellAtomT[type];
	atom[nAtom].id = nAtom + 1;
	nAtom++;
	return 0;
}
int DeleteAtoms(double block[3][2])
{
	int n, d, i;
	int flag[3];

	for (n = 0; n < nAtom; n++)
	{
		for (d = 0; d < 3; d++)
		{
			flag[d] = 0;
			if (atom[n].r[d] >= block[d][0] && atom[n].r[d] <= block[d][1])
			{
				flag[d] = 1;
			}
		}
		//改了索引范围
		if (flag[0] && flag[1] && flag[2])
		{
			for (i = n; i < nAtom - 1; i++)
			{
				atom[i] = atom[i + 1];
			}
			n--;
			nAtom -= 1;
		}
	}
	return 0;
}
int RandomDisplace()
{
	int n, j;

	for (n = 0; n < nAtom; n++)
	{
		for (j = 0; j < 3; j++)
		{
			atom[n].r[j] += 0.5 * (rand() / (RAND_MAX + 1.0));
		}
	}
	return 0;
}

int EnergyAndForce_2Body(int ifEnergy, int ifForce)
{
	int i;
	int j;
	double distance, energy;
	double dr[3], force[3];

	if (ifEnergy == 1)
	{
		potentialEnergy = 0.0;
	}
	for (i = 0; i < nAtom; i++)
	{
		if (ifEnergy == 1)
		{
			atom[i].energy = 0;
		}
		if (ifForce == 1)
		{
			atom[i].force[0] = 0;
			atom[i].force[1] = 0;
			atom[i].force[2] = 0;
		}
	}
	PBC_r();
	for (i = 0; i < nAtom; i++)
	{
		for (j = i + 1; j < nAtom; j++)
		{
			PBC_dr(dr, i, j);
			distance = sqrt(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2));
			if (distance <= cutoff[atom[i].ttype][atom[j].ttype])
			{

				if (ifEnergy == 1)
				{
					energy = Energy_LJ(distance, epsilons[atom[i].ttype][atom[j].ttype], sigmas[atom[i].ttype][atom[j].ttype]);
				}
				if (ifForce == 1)
				{
					Force_LJ(force, dr, distance, epsilons[atom[i].ttype][atom[j].ttype], sigmas[atom[i].ttype][atom[j].ttype]);
				}

				if (ifEnergy == 1)
				{
					potentialEnergy += energy;
					atom[i].energy += energy / 2;
					atom[j].energy += energy / 2;
				}
				if (ifForce == 1)
				{
					atom[i].force[0] += force[0];
					atom[i].force[1] += force[1];
					atom[i].force[2] += force[2];
					atom[j].force[0] += -force[0];
					atom[j].force[1] += -force[1];
					atom[j].force[2] += -force[2];
				}
			}
		}
	}
	return 0;
}
double Energy_LJ(double distance, double epsilon, double sigma)
{
	double  energy;

	energy = 4.*epsilon*(pow((sigma / distance), 12.) - pow((sigma / distance), 6.));
	return energy;
}
int Force_LJ(double force[3], double dr[3], double distance, double epsilon, double sigma)
{
	int d;
	double  distance_p8, distance_p14;


	distance_p14 = pow(distance, 14.);
	distance_p8 = pow(distance, 8.);

	for (d = 0; d != 3; d++)
	{
		force[d] = 24 * epsilon * pow(sigma, 6.) * (2 * pow(sigma, 6.) * dr[d] / distance_p14 - dr[d] / distance_p8);
	}
	return 0;
}
double EnergyAndForce_EAM()
{
	double energy = 0;
	double dr[3];
	PBC_r();
	for (int i = 0; i != nAtom; i++)
	{
		double rho = 0;
		for (int j = 0; j != nAtom; j++)
		{
			if (i != j)
			{
				PBC_dr(dr, i, j);
				double r = sqrt(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2));
				for (int t = 0; t != EAM_n_phi; t++)
					if (EAM_parameter_delta_phi[t] >= r)
						energy += EAM_parameter_a_phi[t] * pow(EAM_parameter_delta_phi[t] - r, 3) / 2.;
				for (int t = 0; t != EAM_n_rho; t++)
				{
					if (r < 2.002970124727)
					{
						rho += EAM_parameter_a_rho[t] * pow(EAM_parameter_delta_rho[t] - 2.002970124727, 3);
					}
					else if (EAM_parameter_delta_rho[t] >= r)
					{
						rho += EAM_parameter_a_rho[t] * pow(EAM_parameter_delta_rho[t] - r, 3);
					}
				}
			}
		}
		energy += (EAM_parameter_a_F[0] * sqrt(rho) + EAM_parameter_a_F[1] * rho * rho);
		for (int d = 0; d != 3; d++)
		{
			atom[i].force[d] = 0;
			for (int j = 0; j != nAtom; j++)
			{
				if (i != j)
				{
					PBC_dr(dr, i, j);
					double r = sqrt(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2));
					for (int t = 0; t != EAM_n_phi; t++)
						if (EAM_parameter_delta_phi[t] >= r)
							atom[i].force[d] += 3 * EAM_parameter_a_phi[t] * pow(EAM_parameter_delta_phi[t] - r, 2) * (dr[d] / r);
					for (int t = 0; t != EAM_n_rho; t++)
						if (EAM_parameter_delta_rho[t] >= r)
							atom[i].force[d] += (EAM_parameter_a_F[0] * (1 / (2 * sqrt(rho))) + EAM_parameter_a_F[1] * 2 * rho)
							* EAM_parameter_a_rho[t] * 3 * pow(EAM_parameter_delta_rho[t] - r, 2) * (dr[d] / r);
				}
			}
		}
	}
	potentialEnergy = energy;
	return energy;
}
double EnergyAndForce_EAM_Cu(double cutoff)
{
	static const double F0 = 3.54 - 1.30, F1 = 1.0241, alpha = 0.3902, beta = 6.0641, ra = 2.3051, Xb = 3.00, re = 2.3051;
	double A1, A2, B1;
	double energy = 0.0, rou = 0.0, rou_e = 6.05;
	//rou_e = 0.015;
	double dr[3];
	PBC_r();
	//two-body potential
	for (int i = 0; i != nAtom; i++)
	{
		for (int j = i + 1; j != nAtom; j++)
		{
			PBC_dr(dr, i, j);
			double r = sqrt(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2));
			if (r < cutoff)
			{
				energy += -alpha * (1 + beta * ((r / ra) - 1)) * exp(-beta * ((r / ra) - 1));
			}
		}
	}
	//embedding energy
	for (int i = 0; i != nAtom; i++)
	{
		rou = 0;
		for (int j = 0; j != nAtom; j++)
		{
			if (i != j)
			{
				PBC_dr(dr, i, j);
				double r = sqrt(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2));
				if (r < cutoff)
				{
					rou += exp(Xb * re - Xb * r);
				}
			}
		}
		energy += -F0 * (1 - log(pow(rou / rou_e, 0.5))) * pow((rou / rou_e), 0.5) + F1 * (rou / rou_e);
	}
	potentialEnergy = energy;
	A1 = F0 * Xb * 0.5 *exp(Xb * ra) / rou;
	A2 = F1 * exp(Xb * ra) / rou;
	B1 = alpha * beta * exp(beta) / rou;
	for (int i = 0; i != nAtom; i++)
	{

		for (int j = 0; j != nAtom; j++)
		{
			if (i != j)
			{
				PBC_dr(dr, i, j);
				double r = sqrt(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2));
				if (r < cutoff)
				{
					for (int k = 0; k != 3; k++)
					{
						atom[i].force[k] += (A1 * Xb * dr[k] - A1 - A2 * dr[k]) * exp(-Xb * dr[k]) + (B1 - B1 * beta * dr[k] / ra) * exp(-beta * dr[k] / ra);
					}
				}
			}
		}
	}
	return energy;
}
int Minimize_SD(double energyCrit, double alpha0, double rho, double c, char filename[20], int outSteps)
{
	double dE;
	double thisEnergy, lastEnergy;
	char filenameStep[20];
	int n, d;

	EnergyAndForce_2Body(1, 0);
	lastEnergy = potentialEnergy;
	iterTime = 0;
	do
	{
		EnergyAndForce_2Body(0, 1);
		for (n = 0; n < nAtom; n++)
		{
			for (d = 0; d < 3; d++)
			{
				atom[n].minDirection[d] = atom[n].force[d];
			}
		}
		thisEnergy = LineSearch_Backtracking(alpha0, rho, c, filename, outSteps);
		dE = lastEnergy - thisEnergy;
		lastEnergy = thisEnergy;
	} while (dE > energyCrit);
	printf("%d %30.25f\n", iterTime, potentialEnergy);
	sprintf(filenameStep, "min/%d.%s", iterTime, filename);
	Output(filename, 1);

	return 0;
}
int Minimize_SD_EAM(double energyCrit, double alpha0, double rho, double c, char filename[20], int outSteps)
{
	double dE;
	double thisEnergy, lastEnergy;
	char filenameStep[20];
	int n, d;

	potentialEnergy = EnergyAndForce_EAM();
	lastEnergy = potentialEnergy;
	iterTime = 0;
	do
	{
		potentialEnergy = EnergyAndForce_EAM();
		for (n = 0; n < nAtom; n++)
		{
			for (d = 0; d < 3; d++)
			{
				atom[n].minDirection[d] = atom[n].force[d];
			}
		}
		thisEnergy = LineSearch_Backtracking_EAM(alpha0, rho, c, filename, outSteps);
		dE = lastEnergy - thisEnergy;
		lastEnergy = thisEnergy;
	} while (dE > energyCrit);
	printf("%d %30.25f\n", iterTime, potentialEnergy);
	sprintf(filenameStep, "min/%d.%s", iterTime, filename);
	Output(filename, 1);
	return 0;
}
int Minimize_CG(double energyCrit, double alpha0, double rho, double c, char filename[20], int outSteps)
{
	double dE;
	double thisEnergy, lastEnergy;
	double beta;
	int n, d;
	char filenameStep[20];

	EnergyAndForce_2Body(1, 0);
	lastEnergy = potentialEnergy;
	iterTime = 0;

	EnergyAndForce_2Body(0, 1);
	for (n = 0; n < nAtom; n++)
	{
		for (d = 0; d < 3; d++)
		{
			atom[n].minDirection[d] = atom[n].force[d];
		}
	}
	do
	{
		thisEnergy = LineSearch_Backtracking(alpha0, rho, c, filename, outSteps);
		dE = lastEnergy - thisEnergy;
		lastEnergy = thisEnergy;


		// Compute CG direction
		for (n = 0; n < nAtom; n++)
		{
			for (d = 0; d < 3; d++)
			{
				atom[n].lastForce[d] = atom[n].force[d];
			}
		}
		EnergyAndForce_2Body(0, 1);
		beta = ComputeBeta_CG();
		for (n = 0; n < nAtom; n++)
		{
			for (d = 0; d < 3; d++)
			{
				atom[n].lastMinDirection[d] = atom[n].minDirection[d];
			}
		}
		ComputeMinDirection_CG(beta);
	} while (dE > energyCrit);

	printf("%d %30.25f\n", iterTime, potentialEnergy);
	sprintf(filenameStep, "min/%d.%s", iterTime, filename);
	Output(filenameStep, 1);

	return 0;
}
int Minimize_CG_EAM(double energyCrit, double alpha0, double rho, double c, char filename[20], int outSteps)
{
	double dE;
	double thisEnergy, lastEnergy;
	double beta;
	int n, d;
	char filenameStep[20];

	EnergyAndForce_EAM();
	lastEnergy = potentialEnergy;
	iterTime = 0;

	EnergyAndForce_EAM();
	for (n = 0; n < nAtom; n++)
	{
		for (d = 0; d < 3; d++)
		{
			atom[n].minDirection[d] = atom[n].force[d];
		}
	}
	do
	{
		thisEnergy = LineSearch_Backtracking_EAM(alpha0, rho, c, filename, outSteps);
		dE = lastEnergy - thisEnergy;
		lastEnergy = thisEnergy;

		// Compute CG direction
		for (n = 0; n < nAtom; n++)
		{
			for (d = 0; d < 3; d++)
			{
				atom[n].lastForce[d] = atom[n].force[d];
			}
		}
		EnergyAndForce_EAM();
		beta = ComputeBeta_CG();
		for (n = 0; n < nAtom; n++)
		{
			for (d = 0; d < 3; d++)
			{
				atom[n].lastMinDirection[d] = atom[n].minDirection[d];
			}
		}
		ComputeMinDirection_CG(beta);
	} while (dE > energyCrit);

	printf("%d %30.25f\n", iterTime, potentialEnergy);
	sprintf(filenameStep, "min/%d.%s", iterTime, filename);
	Output(filename, 1);

	return 0;
}
double ComputeBeta_CG()
{
	int n, d;
	double numerator = 0;
	double denominator = 0;
	double beta;
	for (n = 0; n < nAtom; n++)
	{
		for (d = 0; d < 3; d++)
		{
			numerator += (atom[n].force[d] - atom[n].lastForce[d]) * atom[n].force[d];
			denominator += atom[n].lastForce[d] * atom[n].lastForce[d];
		}
	}
	beta = numerator / denominator;
	if (beta > 0)
	{
		return beta;
	}
	else
	{
		return 0;
	}
}
int ComputeMinDirection_CG(double beta)
{
	int n, d;
	for (n = 0; n < nAtom; n++)
	{
		for (d = 0; d < 3; d++)
		{
			atom[n].minDirection[d] = atom[n].force[d] + beta * atom[n].lastMinDirection[d];
		}
	}
	return 0;
}
double LineSearch_Backtracking(double alpha0, double rho, double c, char filename[20], int outSteps)
{
	double startEnergy, endEnergy;
	double alpha;
	int i, d;
	char filenameStep[20];

	alpha = alpha0;

	//EnergyAndForce_2Body(1, 0);
	startEnergy = potentialEnergy;

	for (i = 0; i < nAtom; i++)
	{
		for (d = 0; d < 3; d++)
		{
			atom[i].backtrackingStartR[d] = atom[i].r[d];
		}
	}

	do
	{
		endEnergy = startEnergy;
		for (i = 0; i < nAtom; i++)
		{
			for (d = 0; d < 3; d++)
			{
				atom[i].r[d] = atom[i].backtrackingStartR[d] + atom[i].minDirection[d] * alpha0;
				endEnergy -= c * atom[i].minDirection[d] * atom[i].minDirection[d] * alpha0;
			}

		}
		//printf("%f,end:%f\n", startEnergy, endEnergy);
		EnergyAndForce_2Body(1, 0);
		alpha0 *= rho;
		if (iterTime % outSteps == 0)
		{
			printf("%d %30.25f\n", iterTime, potentialEnergy);
			sprintf(filenameStep, "min/%d.%s", iterTime, filename);
			Output(filename, 1);
		}
		iterTime++;

		if (iterTime > 500)
		{
			//printf("%30.25f, %30.25f, %30.25f, %30.25f\n", alpha0, startEnergy, endEnergy, potentialEnergy-endEnergy);
		}

		//printf("%40.35f, %40.35f, %d\n", potentialEnergy, endEnergy, potentialEnergy > endEnergy);
	} while (potentialEnergy > endEnergy);

	return potentialEnergy;
}

double LineSearch_Backtracking_EAM(double alpha0, double rho, double c, char filename[20], int outSteps)
{
	double startEnergy, endEnergy;
	double alpha;
	int i, d;
	char filenameStep[20];

	alpha = alpha0;

	startEnergy = potentialEnergy;

	for (i = 0; i < nAtom; i++)
	{
		for (d = 0; d < 3; d++)
		{
			atom[i].backtrackingStartR[d] = atom[i].r[d];
		}
	}

	do
	{
		endEnergy = startEnergy;
		for (i = 0; i < nAtom; i++)
		{
			for (d = 0; d < 3; d++)
			{
				atom[i].r[d] = atom[i].backtrackingStartR[d] + atom[i].minDirection[d] * alpha0;
				endEnergy -= c * atom[i].minDirection[d] * atom[i].minDirection[d] * alpha0;
			}

		}
		//printf("%f,end:%f\n", startEnergy, endEnergy);
		potentialEnergy = EnergyAndForce_EAM();;
		alpha0 *= rho;
		if (iterTime % outSteps == 0)
		{
			printf("%d %30.25f\n", iterTime, potentialEnergy);
			sprintf(filenameStep, "min/%d.%s", iterTime, filename);
			Output(filename, iterTime);
		}
		iterTime++;

		if (iterTime > 500)
		{
			//printf("%30.25f, %30.25f, %30.25f, %30.25f\n", alpha0, startEnergy, endEnergy, potentialEnergy-endEnergy);
		}

		//printf("%40.35f, %40.35f, %d\n", potentialEnergy, endEnergy, potentialEnergy > endEnergy);
	} while (potentialEnergy > endEnergy);

	return potentialEnergy;
}


double Gaussrand(double V, double E) // 方差V,期望E
{
	static double V1, V2, S;
	static int phase = 0;
	double X;


	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X * sqrt(V) + E;
}
void InitVelocity(double T)
{
	int i, a, b, temp;
	int M[MAX_ATOM_NUMBER];
	for (i = 0; i != nAtom; i++)
		M[i] = i;
	for (i = 0; i != nAtom; i++)
	{
		a = rand() % nAtom;
		b = rand() % nAtom;
		temp = M[a];
		M[a] = M[b];
		M[b] = temp;
	}
	for (i = 0; i < nAtom;)
	{
		if (i == nAtom - 1) {
			atom[i].velocity[0] = atom[i].velocity[1] = atom[i].velocity[2] = 0;
			break;
		}
		atom[i].velocity[0] = Gaussrand(k_B * T / atom[0].mass, 0);
		atom[i + 1].velocity[0] = -atom[i].velocity[0];
		atom[i].velocity[1] = Gaussrand(k_B * T / atom[0].mass, 0);
		atom[i + 1].velocity[1] = -atom[i].velocity[1];
		atom[i].velocity[2] = Gaussrand(k_B * T / atom[0].mass, 0);
		atom[i + 1].velocity[2] = -atom[i].velocity[2];
		i += 1;
	}
}
void Dynamic(double total_time, double dt, int nplot, int nProb, int nControl, char filename[20], double tempratureNodes[MAX_TEMP_NODES], double temperatureTimeNodes[MAX_TEMP_NODES])
{
	//后续更新：检测nControl是否整除nProb
	int timeStep = 0;
	double time = 0;
	double temperature, targetTempeture;
	while (time < total_time)
	{
		VelocityVerlet(dt);
		time += dt;
		timeStep += 1;

		if (timeStep % nProb == 0)
		{
			temperature = ComputeTemperature();
		}

		if (timeStep % nControl == 0)
		{
			targetTempeture = GetTargetTempeture(tempratureNodes, temperatureTimeNodes, time, total_time); //global variable
			ControlTempeture(temperature, targetTempeture);
		}

		if (timeStep % nplot == 0)
		{
			printf("%d\t%f\t%f\n", timeStep, time, temperature);
			Output(filename, timeStep);
		}
	}

}
double ComputeTemperature()
{
	int i, j;
	double kineticEnergy = 0.0;
	double temperature;
	for (i = 0; i != nAtom; i++)
	{
		for (j = 0; j != 3; j++)
		{
			kineticEnergy += 0.5 * atom[i].mass * atom[i].velocity[j] * atom[i].velocity[j];
		}
	}
	temperature = kineticEnergy * 2.0 / 3.0 / nAtom / k_B;
	return temperature;
}
void ComputeVolume()
{
	volume = lattice_constant * lattice_constant * lattice_constant * latSize[0] * latSize[1] * latSize[2];
}
void ControlTempeture(double temperature, double targetTempeture)
{
	double scaleFactor;
	int i, d;
	//printf("T:%f  %f \n", targetTempeture, temperature);
	scaleFactor = sqrt((targetTempeture / temperature));
	for (i = 0; i < nAtom; i++)
	{
		for (d = 0; d < 3; d++)
		{
			atom[i].velocity[d] = atom[i].velocity[d] * scaleFactor;

		}
		//printf("@@@@@@@%f\n", atom[i].velocity[0]);
	}
}
double GetTargetTempeture(double tempratureNodes[MAX_TEMP_NODES], double temperatureTimeNodes[MAX_TEMP_NODES], double time, double total_time)
{
	double targetTempeture, temp_o, temp_i, timeInter, timeGone;
	int i = 1;

	while (1)
	{
		if (temperatureTimeNodes[i] >= time)
		{
			temp_o = tempratureNodes[i];
			temp_i = tempratureNodes[i - 1];
			timeInter = temperatureTimeNodes[i] - temperatureTimeNodes[i - 1];
			timeGone = time - temperatureTimeNodes[i - 1];
			break;
		}
		i++;
	}

	targetTempeture = temp_i + (temp_o - temp_i) / timeInter * timeGone;
	return targetTempeture;
}
int VelocityVerlet(double dt)
{
	int i, d;
	EnergyAndForce_EAM();
	for (i = 0; i < nAtom; i++)
	{
		for (d = 0; d < 3; d++)
		{
			atom[i].ft[d] = atom[i].force[d];
			atom[i].dr[d] = dt * atom[i].velocity[d] + atom[i].force[d] / (2 * atom[i].mass) * (dt * dt);
			atom[i].r[d] += atom[i].dr[d];
		}
	}
	PBC_r();
	EnergyAndForce_EAM();
	for (i = 0; i < nAtom; i++)
	{
		for (d = 0; d < 3; d++)
		{
			atom[i].velocity[d] += (atom[i].force[d] + atom[i].ft[d]) * dt / (2 * atom[i].mass);
		}
	}
	return 0;
}

/*int main()
{
	double distance;
	int i = 0;
	for (distance =3.7; distance <= 5.0; distance += 0.001)
	{
		potentialEnergy = 0;
		//---------------INPUT-----------------
		tranVec[0][0] = distance; tranVec[1][0] = 0; tranVec[2][0] = 0;
		tranVec[0][1] = 0; tranVec[1][1] = distance; tranVec[2][1] = 0;
		tranVec[0][2] = 0; tranVec[1][2] = 0; tranVec[2][2] = distance;
		latSize[0] =5; latSize[1] =5; latSize[2] = 5;


		cellAtomN = 4;
		cellAtomR[0][0] = 0.0; cellAtomR[0][1] = 0.0;  cellAtomR[0][2] = 0.0;
		cellAtomR[1][0] = distance / 2; cellAtomR[1][1] = distance / 2;  cellAtomR[1][2] = 0.0;
		cellAtomR[2][0] = distance / 2; cellAtomR[2][1] = 0.0;  cellAtomR[2][2] = distance / 2;
		cellAtomR[3][0] = 0.0; cellAtomR[3][1] = distance / 2;  cellAtomR[3][2] = distance / 2;

		cellAtomT[0] = 0;
		cellAtomT[1] = 0;
		cellAtomT[2] = 0;
		cellAtomT[3] = 0;

		epsilons[0][0] = 0.0031;
		sigmas[0][0] = 2.74;
		cutoff[0][0] = 4.01*distance;
		//------------RUN---------------------
		Crystal();
		PBC_r();
		EnergyAndForce_2Body(1, 1);
		printf("%f  %lf  %lf\n", distance, potentialEnergy, potentialEnergy / nAtom);
	}
	return 0;
}*/
/*
int main()
{
	double distance;
	int i = 0;
	for (distance = 2.5; distance <= 4.0; distance += 0.01)
	{
		potentialEnergy = 0;
		//---------------INPUT-----------------
		tranVec[0][0] = 100; tranVec[1][0] = 0; tranVec[2][0] = 0;
		tranVec[0][1] = 0; tranVec[1][1] = 100; tranVec[2][1] = 0;
		tranVec[0][2] = 0; tranVec[1][2] = 0; tranVec[2][2] = 100;
		latSize[0] = 1; latSize[1] = 1; latSize[2] = 1;


		cellAtomN = 2;
		cellAtomR[0][0] = 0.0; cellAtomR[0][1] = 0.0;  cellAtomR[0][2] = 0.0;
		cellAtomR[1][0] = distance; cellAtomR[1][1] = 0.0;  cellAtomR[1][2] = 0.0;

		cellAtomT[0] = 0;
		cellAtomT[1] = 0;

		epsilons[0][0] = 0.0031;
		sigmas[0][0] = 2.74;
		cutoff[0][0] = 4.01 * distance;
		//------------RUN---------------------
		Crystal();
		PBC_r();
		EnergyAndForce_2Body(1, 1);
		printf("%f  %lf  %lf\n", distance, potentialEnergy, atom[0].force[0]);
	}
	return 0;
}
*/
/*
int main()
{
	double distance;
	int i = 0;
	for (distance = 3.0; distance <= 3.5; distance += 0.01)
	{
		potentialEnergy = 0;
		//---------------INPUT-----------------
		tranVec[0][0] = distance; tranVec[1][0] = 0; tranVec[2][0] = 0;
		tranVec[0][1] = 0; tranVec[1][1] = distance; tranVec[2][1] = 0;
		tranVec[0][2] = 0; tranVec[1][2] = 0; tranVec[2][2] = distance;
		latSize[0] = 10; latSize[1] = 10; latSize[2] = 10;


		cellAtomN = 2;
		cellAtomR[0][0] = 0.0; cellAtomR[0][1] = 0.0;  cellAtomR[0][2] = 0.0;
		cellAtomR[1][0] = distance / 2; cellAtomR[1][1] = distance / 2;  cellAtomR[1][2] = distance / 2;
		cellAtomR[2][0] = distance / 2; cellAtomR[2][1] = 0.0;  cellAtomR[2][2] = distance / 2;
		cellAtomR[3][0] = 0.0; cellAtomR[3][1] = distance / 2;  cellAtomR[3][2] = distance / 2;

		cellAtomT[0] = 0;
		cellAtomT[1] = 0;
		cellAtomT[2] = 0;
		cellAtomT[3] = 0;

		epsilons[0][0] = 0.0031;
		sigmas[0][0] = 2.74;
		cutoff[0][0] = 1.65 * distance;
		//------------RUN---------------------
		Crystal();
		PBC_r();
		EnergyAndForce_2Body(1, 1);
		printf("%f  %lf  %lf\n", distance, potentialEnergy, atom[1].force[0]);
	}
	return 0;
}*/
/*
int main()
{
	double  Vacancy [3][2];
	double initialEnergy, endEnergy;
	lattice_constant = 3.14;
	tranVec[0][0] = lattice_constant; tranVec[1][0] = 0; tranVec[2][0] = 0;
	tranVec[0][1] = 0; tranVec[1][1] = lattice_constant; tranVec[2][1] = 0;
	tranVec[0][2] = 0; tranVec[1][2] = 0; tranVec[2][2] = lattice_constant;
	latSize[0] =5; latSize[1] = 5; latSize[2] =5;

	cellAtomN = 2;
	cellAtomR[0][0] = 0; cellAtomR[0][1] = 0.0;  cellAtomR[0][2] = 0.0;
	cellAtomR[1][0] = 0.5 * lattice_constant; cellAtomR[1][1] = 0.5 * lattice_constant;  cellAtomR[1][2] = 0.5 * lattice_constant;
	cellAtomT[0] = 0;
	cellAtomT[1] = 0;

	Vacancy[0][0] = 1.9 * lattice_constant; Vacancy[0][1] = 2.1 * lattice_constant;
	Vacancy[1][0] = 1.9 * lattice_constant; Vacancy[1][1] = 2.1 * lattice_constant;
	Vacancy[2][0] = 1.9 * lattice_constant; Vacancy[2][1] = 2.1 * lattice_constant;
	cutoff[0][0] = 10;
	//Minimize
	double Backtracking_alpha0 =1;
	double Backtracking_rho = 0.5;
	double Backtracking_c = 0.001;
	double energyCrit = 1E-5;

	char filename[20] = "vacancy.xyz";
	int outSteps = 1;
	InitOut(filename);
	double cutoff_cu = 2.01 * lattice_constant;

	Crystal();
	PBC_r();
	//初始计算能量
	initialEnergy = EnergyAndForce_EAM();
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	Output(filename, 1);
	//删除原子后的能量
	DeleteAtoms(Vacancy);
	potentialEnergy = EnergyAndForce_EAM();
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	Output(filename, 1);
	//驰豫后的能量
	Minimize_SD_EAM(energyCrit, Backtracking_alpha0, Backtracking_rho, Backtracking_c, filename, outSteps);
	endEnergy = EnergyAndForce_EAM();
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	printf("%lf\n", endEnergy - (nAtom ) * initialEnergy / (1+nAtom));
	return 0;
}*/

int main()
{
	double  Vacancy[3][2];
	double initialEnergy, endEnergy;
	lattice_constant = 3.14;
	tranVec[0][0] = lattice_constant; tranVec[1][0] = 0; tranVec[2][0] = 0;
	tranVec[0][1] = 0; tranVec[1][1] = lattice_constant; tranVec[2][1] = 0;
	tranVec[0][2] = 0; tranVec[1][2] = 0; tranVec[2][2] = lattice_constant;
	latSize[0] = 5; latSize[1] = 5; latSize[2] = 5;


	cellAtomN = 2;
	cellAtomR[0][0] = 0.0; cellAtomR[0][1] = 0.0;  cellAtomR[0][2] = 0.0;
	cellAtomR[1][0] = lattice_constant / 2; cellAtomR[1][1] = lattice_constant / 2;  cellAtomR[1][2] = lattice_constant / 2;
	cellAtomR[2][0] = lattice_constant / 2; cellAtomR[2][1] = 0.0;  cellAtomR[2][2] = lattice_constant / 2;
	cellAtomR[3][0] = 0.0; cellAtomR[3][1] = lattice_constant / 2;  cellAtomR[3][2] = lattice_constant / 2;

	cellAtomT[0] = 0;
	cellAtomT[1] = 0;
	cellAtomT[2] = 0;
	cellAtomT[3] = 0;

	Vacancy[0][0] = -0.1 * lattice_constant; Vacancy[0][1] = 0.1 * lattice_constant;
	Vacancy[1][0] =-0.1 * lattice_constant; Vacancy[1][1] = 0.1 * lattice_constant;
	Vacancy[2][0] =-0.1 * lattice_constant; Vacancy[2][1] = 0.1 * lattice_constant;
	//cutoff[0][0] = 10;
	//Minimize 
	double Backtracking_alpha0 = 0.01;
	double Backtracking_rho = 0.5;
	double Backtracking_c = 0.001;
	double energyCrit = 1E-5;

	char filename[20] = "V.xyz";
	int outSteps = 1;
	InitOut(filename);
	double cutoff_cu = 1.65 * lattice_constant;
	Crystal();
	PBC_r();
	//初始计算能量
	initialEnergy = EnergyAndForce_EAM();
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	printf("%d", nAtom);
	//删除原子后的能量
	DeleteAtoms(Vacancy);
	potentialEnergy = EnergyAndForce_EAM();
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	printf("%d", nAtom);
	//filename = "V.xyz"
	Output(filename, 1);

	/*
	//驰豫后的能量
	Minimize_SD_EAM(energyCrit, Backtracking_alpha0, Backtracking_rho, Backtracking_c, filename, outSteps);
	endEnergy = EnergyAndForce_EAM_Cu(cutoff_cu);
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	printf("%lf\n", endEnergy - (nAtom)*initialEnergy / (1 + nAtom));
	*/
	return 0;
}

/*int main()
{
	double  Create[3];
	double initialEnergy, endEnergy;
	lattice_constant = 3.14;
	tranVec[0][0] = lattice_constant; tranVec[1][0] = 0; tranVec[2][0] = 0;
	tranVec[0][1] = 0; tranVec[1][1] = lattice_constant; tranVec[2][1] = 0;
	tranVec[0][2] = 0; tranVec[1][2] = 0; tranVec[2][2] = lattice_constant;
	latSize[0] = 15; latSize[1] = 15; latSize[2] = 15;

	cellAtomN = 2;
	cellAtomR[0][0] = 0; cellAtomR[0][1] = 0.0;  cellAtomR[0][2] = 0.0;
	cellAtomR[1][0] = 0.5 * lattice_constant; cellAtomR[1][1] = 0.5 * lattice_constant;  cellAtomR[1][2] = 0.5 * lattice_constant;
	cellAtomT[0] = 0;
	cellAtomT[1] = 0;

	cutoff[0][0] = 7;

	Create[0] = 1.8* lattice_constant;
	Create[1] = 1.8 * lattice_constant;
	Create[2] = 1.8 * lattice_constant;

	//Minimize
	double Backtracking_alpha0 = 0.01;
	double Backtracking_rho = 0.5;
	double Backtracking_c = 0.001;
	double energyCrit = 1E-5;

	char filename[20] = "sia.xyz";
	int outSteps = 1;
	InitOut(filename);


	pairStyle = "LJ";
	epsilons[0][0] = 0.0031;
	sigmas[0][0] = 2.74;
	cutoff[0][0] = 12;
	Crystal();
	PBC_r();
	//初始计算能量
	initialEnergy = EnergyAndForce_EAM();
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	Output(filename, 1);
	//增加SIA后的能量
	atom[124].r[0] += 0.628; atom[124].r[1] += 0.628; atom[124].r[2] += 0.628;
	CreateAtom(Create, 0);
	PBC_r();
	Output(filename, 1);
	potentialEnergy = EnergyAndForce_EAM();
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	Output(filename, 1);
	//驰豫后的能量
	Minimize_SD_EAM(energyCrit, Backtracking_alpha0, Backtracking_rho, Backtracking_c, filename, outSteps);
	endEnergy = EnergyAndForce_EAM();
	printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
	printf("%lf\n", endEnergy - nAtom*initialEnergy/(nAtom-1));
	return 0;
}*/
/*
int main()
{
	double  Vacancy[3][2], Create[3];
	double initialEnergy, endEnergy;
	int flag_V_SIA = 2, flag_static_or_dynamic = 1;
	lattice_constant = 3.14;
	tranVec[0][0] = lattice_constant; tranVec[1][0] = 0; tranVec[2][0] = 0;
	tranVec[0][1] = 0; tranVec[1][1] = lattice_constant; tranVec[2][1] = 0;
	tranVec[0][2] = 0; tranVec[1][2] = 0; tranVec[2][2] = lattice_constant;
	latSize[0] = 5; latSize[1] = 5; latSize[2] = 5;

	cellAtomN = 2;
	cellAtomR[0][0] = 0; cellAtomR[0][1] = 0.0;  cellAtomR[0][2] = 0.0;
	cellAtomR[1][0] = 0.5 * lattice_constant; cellAtomR[1][1] = 0.5 * lattice_constant;  cellAtomR[1][2] = 0.5 * lattice_constant;
	cellAtomT[0] = 0;
	cellAtomT[1] = 0;

	if (flag_V_SIA == 1)
	{
		//case 1:V and SIA annihilate (Static relaxation)
		Vacancy[0][0] = 2.9 * lattice_constant; Vacancy[0][1] = 3.1 * lattice_constant;
		Vacancy[1][0] = 2.9 * lattice_constant; Vacancy[1][1] = 3.1 * lattice_constant;
		Vacancy[2][0] = 2.9 * lattice_constant; Vacancy[2][1] = 3.1 * lattice_constant;
	}
	else if (flag_V_SIA == 2)
	{
		//case 2:V and SIA (Static relaxation) and case 3:V and SIA (dynamic relaxation)
		Vacancy[0][0] = 3.9 * lattice_constant; Vacancy[0][1] = 4.1 * lattice_constant;
		Vacancy[1][0] = 0.9 * lattice_constant; Vacancy[1][1] = 1.1 * lattice_constant;
		Vacancy[2][0] = 3.9 * lattice_constant; Vacancy[2][1] = 4.1 * lattice_constant;
	}

	Create[0] = 1.8 * lattice_constant;
	Create[1] = 1.8 * lattice_constant;
	Create[2] = 1.8 * lattice_constant;
	typeMass[0] = 74.0;
	cutoff[0][0] = 5.01*lattice_constant;
	//Minimize
	double Backtracking_alpha0 = 0.001;
	double Backtracking_rho = 0.5;
	double Backtracking_c = 0.001;
	double energyCrit = 1E-5;

	char filename[20] = "static2.xyz";
	int outSteps = 1;
	InitOut(filename);

	Crystal();
	PBC_r();
	if (flag_static_or_dynamic == 1)
	{
		//初始计算能量
		initialEnergy = EnergyAndForce_EAM();
		printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
		Output(filename, 1);
		//制造空位和SIA后的能量
		DeleteAtoms(Vacancy);
		atom[124].r[0] += 0.628; atom[124].r[1] += 0.628; atom[124].r[2] += 0.628;
		CreateAtom(Create, 0);
		potentialEnergy = EnergyAndForce_EAM();
		printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
		Output(filename, 1);
		//驰豫后的能量
		Minimize_SD_EAM(energyCrit, Backtracking_alpha0, Backtracking_rho, Backtracking_c, filename, outSteps);
		endEnergy = EnergyAndForce_EAM();
		printf("%f  %lf  %lf \n", lattice_constant, potentialEnergy, potentialEnergy / nAtom);
		printf("%lf\n", endEnergy - (nAtom)*initialEnergy / (1 + nAtom));
	}
	else if (flag_static_or_dynamic == 2)
	{
		//初始化必要参数
		double total_time = 300;
		double dt = 0.1;
		int nplot = 10;
		int nProb = 1;
		int nControl = 10;
		double initTemperature = 800;

		double tempratureNodes[MAX_TEMP_NODES];
		double temperatureTimeNodes[MAX_TEMP_NODES];
		tempratureNodes[0] = 800; tempratureNodes[1] = 800;
		temperatureTimeNodes[0] = 0; temperatureTimeNodes[1] = total_time;
		//构建V和SIA,输出构型
		DeleteAtoms(Vacancy);
		atom[124].r[0] += 0.628; atom[124].r[1] += 0.628; atom[124].r[2] += 0.628;
		CreateAtom(Create, 0);
		Output(filename, 1);
		printf("%d\n", nAtom);

		//Minmize
		Minimize_SD_EAM(energyCrit, Backtracking_alpha0, Backtracking_rho, Backtracking_c, filename, outSteps);
		//设定初始温度
		InitVelocity(initTemperature);
		//Dynamic relax
		Dynamic(total_time, dt, nplot, nProb, nControl, filename, tempratureNodes, temperatureTimeNodes);
		Minimize_SD_EAM(energyCrit, Backtracking_alpha0, Backtracking_rho, Backtracking_c, filename, outSteps);

	}
	return 0;
}
*/