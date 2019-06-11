#include "PML.h"
PML::PML(int nx, int ny, long long nt, double rv0, char* fn1, char* fn2){
	Npx = nx;
	Npy = ny;
	Nt = nt;
	rv = rv0;
	reg = 1;
	flag = 1;
	file1.open(fn1);
	file2.open(fn2);
	apx = 0.0; bpx = 0.07; apy = 0.0; bpy = 0.064;
	at = 0.0; bt = 1e-9;
	err = 1e-7;
	cv = 3.0e8;
	nPML = 12;
	we = 2 * PI*sqrt(2)*30e9;
	wm = we;
	ge = 1e8;
	gm = ge;
	xM0 = 0.024; xM1 = 0.054;
	yM0 = 0.002; yM1 = 0.062;
	xf = 0.004;
	yf1 = 0.027; yf2 = 0.037;
	sw_f0 = 30e9;
	sw_w0 = 2.0*PI*sw_f0; sw_m = 2.0; sw_k = 100.0;
	lapse = 10;
	p_co = 0.5;
}
void PML::Initialize(){
	dt = (bt - at) / (Nt - 1);
	Nx = Npx + 2 * nPML; Ny = Npy + 2 * nPML;
}
void PML::Allo(){
	// Hx
	Hx = new double**[Ny];
	for (int j = 0; j < Ny; j++)
	{
		Hx[j] = new double*[Nx - 1];
	}
	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx - 1; i++)
		{
			Hx[j][i] = new double[2]();
		}
	}
	// Hy
	Hy = new double**[Ny - 1];
	for (int j = 0; j < Ny - 1; j++)
	{
		Hy[j] = new double*[Nx];
	}
	for (int j = 0; j < Ny - 1; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			Hy[j][i] = new double[2]();
		}
	}
	// Kx
	Kx = new double**[Ny];
	for (int j = 0; j < Ny; j++)
	{
		Kx[j] = new double*[Nx - 1];
	}
	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx - 1; i++)
		{
			Kx[j][i] = new double[2]();
		}
	}
	// Ky
	Ky = new double**[Ny - 1];
	for (int j = 0; j < Ny - 1; j++)
	{
		Ky[j] = new double*[Nx];
	}
	for (int j = 0; j < Ny - 1; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			Ky[j][i] = new double[2]();
		}
	}
	// E
	E = new double**[Ny - 1];
	for (int j = 0; j < Ny - 1; j++)
	{
		E[j] = new double*[Nx - 1];
	}
	for (int j = 0; j < Ny - 1; j++)
	{
		for (int i = 0; i < Nx - 1; i++)
		{
			E[j][i] = new double[2]();
		}
	}
	// J
	J = new double**[Ny - 1];
	for (int j = 0; j < Ny - 1; j++)
	{
		J[j] = new double*[Nx - 1];
	}
	for (int j = 0; j < Ny - 1; j++)
	{
		for (int i = 0; i < Nx - 1; i++)
		{
			J[j][i] = new double[2]();
		}
	}
	// Em
	Em = new double*[Npy - 1];
	for (int j = 0; j < Npy - 1; j++)
	{
		Em[j] = new double[Npx - 1]();
	}
}
void PML::Partition(){
	px = new double[Nx]();
	py = new double[Ny]();
	dx = new double[Nx-1]();
	dy = new double[Ny-1]();
	double L1, L2, L3;	// length of left, meta and right respectively
	int N1, N2, N3;		// interval number of left, meta and right respectively
	double dx1, dx2, dx3, dy1, dy2, dy3;
	// x direction
	L1 = xM0 - apx;
	L2 = xM1 - xM0;
	L3 = bpx - xM1;
	N1 = (int)ceil((Npx-1)/(L1+L2/p_co+L3)*L1);
	N3 = (int)ceil((Npx-1)/(L1+L2/p_co+L3)*L3);
	N2 = Npx - 1 - N1 - N3;
	dx1 = L1/N1;
	dx2 = L2/N2;
	dx3 = L3/N3;

	dd = nPML*dx1;
	ax = apx - dd; bx = bpx + dd;
	ay = apy - dd; by = bpy + dd;

	px[0] = ax;
	// iterate by interval
	for (int i = 1; i <= N1 + nPML; i++)
	{
		px[i] = ax + dx1*i;
		dx[i-1] = dx1;
	}
	for (int i = 1; i <= N2; i++)
	{
		px[i + N1 + nPML] = xM0 + dx2*i;
		dx[i + N1 + nPML - 1] = dx2;
	}
	for (int i = 1; i <= N3 + nPML; i++)
	{
		px[i + N1 + N2 + nPML] = xM1 + dx3*i;
		dx[i + N1 + N2 + nPML - 1] = dx3;
	}

	// y direction
	L1 = yM0 - apy;
	L2 = yM1 - yM0;
	L3 = bpy - yM1;
	N1 = (int)ceil((Npy-1)/(L1+L2/p_co+L3)*L1);
	N3 = (int)ceil((Npy-1)/(L1+L2/p_co+L3)*L3);
	N2 = Npy - 1 - N1 - N3;
	dy1 = L1/N1;
	dy2 = L2/N2;
	dy3 = L3/N3;

	py[0] = ay;
	// iterate by interval
	for (int i = 1; i <= N1 + nPML; i++)
	{
		py[i] = ay + dy1*i;
		dy[i-1] = dy1;
	}
	for (int i = 1; i <= N2; i++)
	{
		py[i + N1 + nPML] = yM0 + dy2*i;
		dy[i + N1 + nPML - 1] = dy2;
	}
	for (int i = 1; i <= N3 + nPML; i++)
	{
		py[i + N1 + N2 + nPML] = yM1 + dy3*i;
		dy[i + N1 + N2 + nPML - 1] = dy3;
	}

	// set source
	int *nxy;
	nxy = xy2N(xf, yf1);
	Nxf = nxy[0]; Nyf0 = nxy[1];
	nxy = xy2N(xf, yf2);
	Nyf1 = nxy[1];
	delete[] nxy;
	WritePara();
}
void PML::GetSource(int n,long long nt){	// put in E[][][n] and time step is nt
	for (int j = Nyf0; j < Nyf1; j++)
	{
		E[j][Nxf][n] = source(px[Nxf], py[j], pt(nt));
	}
}
void PML::Solve(){
	Initialize();
	Allo();
	Partition();
	//InitValue();	// for test
	GetSource(0,0);	// for simulation
	Match();
	//test();		// for test
	CloseFile();
	cout << "done!" << endl;
}
double PML::pt(long long nt){
	return nt * dt;
}
void PML::CloseFile(){
	file1.close();
	file2.close();
}
int* PML::xy2N(double x, double y){
	int *Nxy = new int[2]();
	int temp = 0;
	for (int i = 0; i < Nx; i++)
	{
		if(px[i]<=x & px[i+1]>x)
		{
			Nxy[0] = i;
			break;
		}
	}
	for (int i = 0; i < Ny; i++)
	{
		if(py[i]<=y & py[i+1]>y)
		{
			Nxy[1] = i;
			break;
		}
	}
	return Nxy;
	// need delete[] Nxy after use
}
void PML::WriteFile(){
	// write Em
	for (int i = 0; i < Npx - 1; i++)
	{
		for (int j = 0; j < Npy - 1; j++)
		{
			file2 << Em[j][i];
			if (j != Npy - 2)
				file2 << " ";	
		}
		file2 << endl;
	}
	file2 << endl;
}
void PML::WritePara(){
	// [apx,bpx,apy,bpy,at,bt,Npx,Npy,Nt,lapse]=vars
	// partition info
	file1 << apx << endl;
	file1 << bpx << endl;
	file1 << apy << endl;
	file1 << bpy << endl;
	file1 << at << endl;
	file1 << bt << endl;
	file1 << Npx << endl;
	file1 << Npy << endl; 
	file1 << Nt << endl;
	// lapse
	file1 << lapse << endl;
	// partition
	for (int i = 0; i < Nx; i++)
	{
		file1 << px[i] << " ";
	}
	file1 << endl;
	for (int i = 0; i < Ny; i++)
	{
		file1 << py[i] << " ";
	}
	file1 << endl;
	file1 << nPML << endl;
}