#ifndef PML_H
#define PML_H
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string>
using namespace std;
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148
class PML  {
public:
	// result files
	ofstream  file1, file2;
	// partition number for time
	long long Nt;
	// parameters for physical domain
	int Npx, Npy;
	double apx, bpx, apy, bpy, at, bt;
	double dt;
	double *dx, *dy;
	double rv;	// R.V.
	int reg;	// The region where is current computing
	int flag;	// always denote the array which is solving. take value 0 or 1.
	// parameters for PML
	int nPML;
	double dd, err, cv;
	// common parameters
	double we, wm, ge, gm, ax, bx, ay, by;
	int Nx, Ny;
	double *px, *py;
	double pt(long long);
	// meta
	double xM0, xM1, yM0, yM1;
	// partition coefficeint. 0<p_co. partition size is p_co*dx in meta region.
	// p_co < 1 means finer.
	double p_co; 
	// source
	double xf, yf1, yf2;
	int Nyf0, Nyf1, Nxf;
	double sw_f0, sw_w0, sw_m, sw_k;
	// solutions
	double ***Hx, ***Hy, ***Kx, ***Ky, ***E, ***J, **Em;	// order: y,x,t
	// parameter for write solution
	int lapse;	// every [lapse] time step, write a solution
	// utility functions
	double Hx_fcn(double, double, double);
	double Hy_fcn(double, double, double);
	double Kx_fcn(double, double, double);
	double Ky_fcn(double, double, double);
	double E_fcn(double, double, double);
	double J_fcn(double, double, double);
	double g1_fcn(double, double, double);
	double g2_fcn(double, double, double);
	double g3_fcn(double, double, double);
	double g4_fcn(double, double, double);
	double g5_fcn(double, double, double);
	double g6_fcn(double, double, double);
	double e0(double, double);
	double mu0(double, double);
	// member functions
	PML(int, int, long long, double, char*, char*); 	// 4th parameter is R.V.
	void Allo();		// initialize arrays
	void Initialize();
	void Partition();
	void GetSource(int,long long);
	void Match();
	void InitValue();	// write initial value, use for test
	void BondValue(long long);	// write boundary value, use for test
	void Iterate();		// iterate the solution
	void CloseFile();
	void WriteFile();
	void WritePara();
	void Solve();
	void test();
	int WhichReg(double, double);
	int* xy2N(double, double);
	double sgmx(double);
	double sgmy(double);
	double source(double, double, double);
	double pJE(double, double);
	double pJ(double, double);
	double pKxHx(double, double);
	double pKx(double, double);
	double pKyHy(double, double);
	double pKy(double, double);

};
#endif