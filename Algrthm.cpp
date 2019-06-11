#include "PML.h"
void PML::Match(){
	int i, j;
	long long it;
	double c0, c1, c2;	// temporary parameters
	double ct, cx, cy;
	for (it = 1; it < Nt; it++)
	{
		GetSource(flag, it);	// for simulation
		//BondValue(it); // use for test
		for (i = 0; i < Nx - 1; i++)
		{
			for (j = 0; j < Ny - 1; j++)
			{
				// loop for E
				// (1) in paper, point of eqn : (i + 1 / 2, j + 1 / 2, n + 1 / 2)
				if ((i != Nxf)|(j >= Nyf1)|(j <= Nyf0))		// do not update source
				{
					ct = pt(it) - 0.5*dt;
					cx = 0.5*(px[i] + px[i + 1]);
					cy = 0.5*(py[j] + py[j + 1]);
					reg = WhichReg(cx, cy);
					c0 = 1 / (e0(cx, cy) / dt + 0.5*e0(cx, cy)*(sgmx(cx) + sgmy(cy)));
					c1 = e0(cx, cy) / dt - 0.5*e0(cx, cy)*(sgmx(cx) + sgmy(cy));
					E[j][i][flag] = c0*((Hy[j][i + 1][!flag] - Hy[j][i][!flag]) / (0.5*dx[i]+0.5*dx[i+1]) - (Hx[j + 1][i][!flag] - Hx[j][i][!flag]) / (0.5*dy[j]+0.5*dy[j+1]) \
						- J[j][i][!flag] + c1*E[j][i][!flag] + g3_fcn(cx, cy, ct));
				}
				// loop for J
				// (2) in paper, point of eqn : (i + 1 / 2, j + 1 / 2, n + 1)
				ct = pt(it);
				cx = 0.5*(px[i] + px[i + 1]);
				cy = 0.5*(py[j] + py[j + 1]);
				reg = WhichReg(cx, cy);
				c0 = 1 / (1 / dt + 0.5*pJ(cx, cy));
				c1 = e0(cx, cy)*pJE(cx, cy);
				c2 = 1 / dt - 0.5*pJ(cx, cy);
				J[j][i][flag] = c0*(c1*E[j][i][flag] + c2*J[j][i][!flag] + g4_fcn(cx, cy, ct));
				if (j != 0)
				{
					// loop for Kx
					// (3) in paper, point of eqn : (i + 1 / 2, j, n + 1 / 2)
					ct = pt(it) - 0.5*dt;
					cx = 0.5*(px[i] + px[i + 1]);
					cy = py[j];
					reg = WhichReg(cx, cy);
					c0 = 1 / (1 / dt + 0.5*pKx(cx, cy));
					c1 = mu0(cx, cy)*pKxHx(cx, cy);
					c2 = 1 / dt - 0.5*pKx(cx, cy);
					Kx[j][i][flag] = c0*(c1*Hx[j][i][!flag] + c2*Kx[j][i][!flag] + g5_fcn(cx, cy, ct));
				}
				if (i != 0)
				{
					// loop for Ky
					// (4) in paper, point of eqn : (i, j + 1 / 2, n + 1 / 2)
					ct = pt(it) - 0.5*dt;
					cx = px[i];
					cy = 0.5*(py[j] + py[j + 1]);
					reg = WhichReg(cx, cy);
					c0 = 1 / (1 / dt + 0.5*pKy(cx, cy));
					c1 = mu0(cx, cy)*pKyHy(cx, cy);
					c2 = 1 / dt - 0.5*pKy(cx, cy);
					Ky[j][i][flag] = c0*(c1*Hy[j][i][!flag] + c2*Ky[j][i][!flag] + g6_fcn(cx, cy, ct));
				}
			}
		}
		
		for (i = 0; i < Nx - 1; i++)
		{
			for (j = 0; j < Ny - 1; j++)
			{
				if (j != 0)
				{	
					// loop for Hx
					// (5) in paper, point of eqn : (i + 1 / 2, j, n + 1)
					ct = pt(it);
					cx = 0.5*(px[i] + px[i + 1]);
					cy = py[j];
					reg = WhichReg(cx, cy);
					c0 = 1 / (mu0(cx, cy) / dt - 0.5*mu0(cx, cy)*(sgmx(cx) - sgmy(cy)));
					c1 = mu0(cx, cy) / dt + 0.5*mu0(cx, cy)*(sgmx(cx) - sgmy(cy));
					Hx[j][i][flag] = c0*(-(E[j][i][flag] - E[j - 1][i][flag]) / dy[j - 1] - Kx[j][i][flag] + c1*Hx[j][i][!flag] + g1_fcn(cx, cy, ct));
				}

				if (i != 0)
				{
					// loop for Hy
					// (6) in paper, point of eqn : (i, j + 1 / 2, n + 1)
					ct = pt(it);
					cx = px[i];
					cy = 0.5*(py[j] + py[j + 1]);
					reg = WhichReg(cx, cy);
					c0 = 1 / (mu0(cx, cy) / dt + 0.5*mu0(cx, cy)*(sgmx(cx) - sgmy(cy)));
					c1 = mu0(cx, cy) / dt - 0.5*mu0(cx, cy)*(sgmx(cx) - sgmy(cy));
					Hy[j][i][flag] = c0*((E[j][i][flag] - E[j][i - 1][flag]) / dx[i - 1] - Ky[j][i][flag] + c1*Hy[j][i][!flag] + g2_fcn(cx, cy, ct));
				}
			}
		}
		
		for (i = nPML; i < Npx + nPML - 1; i++)	// copy the solution E to Em
			for (j = nPML; j < Npy + nPML - 1; j++)
				Em[j - nPML][i - nPML] = E[j][i][flag];

		if (it % lapse == 0)
		{
			WriteFile();
		}
		flag = !flag;	// iterate
		//Iterate();
		/*test();
		int aaaaaaaaaaaaaa = 1;*/
		if (it % 10 == 0)
			cout << "time  " << it << endl;
	}

}
void PML::Iterate(){
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
		{
			if (i != Nx - 1) {
				Hx[j][i][0] = Hx[j][i][1];
				Hx[j][i][1] = 0.0;
				Kx[j][i][0] = Kx[j][i][1];
				Kx[j][i][1] = 0.0;
				if (j != Ny - 1){
					E[j][i][0] = E[j][i][1];
					E[j][i][1] = 0.0;
					J[j][i][0] = J[j][i][1];
					J[j][i][1] = 0.0;
				}
			}
			if (j != Ny - 1)
			{
				Hy[j][i][0] = Hy[j][i][1];
				Hy[j][i][1] = 0.0;
				Ky[j][i][0] = Ky[j][i][1];
				Ky[j][i][1] = 0.0;
			}
		}
}
void PML::InitValue()	// write initial value, use for test
{
	double cx, cy;
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
		{
			if (i != Nx - 1) {
				cx = 0.5*(px[i] + px[i + 1]);
				cy = py[j];
				Hx[j][i][0] = Hx_fcn(cx, cy, at + 0.5*dt);
				Kx[j][i][0] = Kx_fcn(cx,cy,at);
				if (j != Ny - 1){
					cx = 0.5*(px[i] + px[i + 1]);
					cy = 0.5*(py[j] + py[j + 1]);
					E[j][i][0] = E_fcn(cx,cy,at);
					J[j][i][0] = J_fcn(cx, cy, at + 0.5*dt);
				}
			}
			if (j != Ny - 1)
			{
				cx = px[i];
				cy = 0.5*(py[j] + py[j + 1]);
				Hy[j][i][0] = Hy_fcn(cx, cy, at + 0.5*dt);
				Ky[j][i][0] = Ky_fcn(cx, cy, at);
			}
		}
}
void PML::BondValue(long long nt)	// write boundary value, use for test
{
	// always for 2nd layer variables
	double cx, cy;
	double ct = pt(nt);
	for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny; j++)
	{
		if (i != Nx - 1){
			cx = 0.5*(px[i] + px[i + 1]);
			Hx[0][i][1] = Hx_fcn(cx, py[0], ct + 0.5*dt);
			Hx[Ny - 1][i][1] = Hx_fcn(cx, py[Ny - 1], ct + 0.5*dt);
			Kx[0][i][1] = Kx_fcn(cx, py[0], ct);
			Kx[Ny - 1][i][1] = Kx_fcn(cx, py[Ny - 1], ct);
		}
		if (j != Ny - 1){
			cy = 0.5*(py[j] + py[j + 1]);
			Hy[j][0][1] = Hy_fcn(px[0], cy, ct + 0.5*dt);
			Hy[j][Nx - 1][1] = Hy_fcn(px[Nx - 1], cy, ct + 0.5*dt);
			Ky[j][0][1] = Ky_fcn(px[0], cy, ct);
			Ky[j][Nx - 1][1] = Ky_fcn(px[Nx - 1], cy, ct);
		}
	}
}
void PML::test(){
	// test soluition
	double t = -1.0;
	double er = 0.0;
	double cx, cy, ct;
	ct = bt;
	// test Hx
	for (int i = 0; i < Nx - 1; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			cx = 0.5*(px[i] + px[i + 1]);
			cy = py[j];
			er = fabs(Hx[j][i][0] - Hx_fcn(cx, cy, ct + 0.5*dt));
			if (t < er)
				t = er;
			cout << Hx[j][i][0] << " ";
		}
		cout << endl;
	}
	//// test Kx
	/*for (int i = 0; i < Nx - 1; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			cx = 0.5*(px[i] + px[i + 1]);
			cy = py[j];
			er = fabs(Kx[j][i][0] - Kx_fcn(cx, cy, ct));
			if (t < er)
				t = er;
			cout << Kx[j][i][0] << " ";
		}
		cout << endl;
	}*/
	//// test E
	//for (int i = 0; i < Nx - 1; i++)
	//{
	//	for (int j = 0; j < Ny - 1; j++)
	//	{
	//		cx = 0.5*(px[i] + px[i + 1]);
	//		cy = 0.5*(py[j] + py[j + 1]);
	//		er = fabs(E[j][i][0] - E_fcn(cx, cy, ct));
	//		if (t < er)
	//			t = er;
	//		cout << E[j][i][0] << " ";
	//	}
	//	cout << endl;
	//}
	//// test J
	//for (int i = 0; i < Nx - 1; i++)
	//{
	//	for (int j = 0; j < Ny - 1; j++)
	//	{
	//		cx = 0.5*(px[i] + px[i + 1]);
	//		cy = 0.5*(py[j] + py[j + 1]);
	//		er = fabs(J[j][i][0] - J_fcn(cx, cy, ct + 0.5*dt));
	//		if (t < er)
	//			t = er;
	//		cout << J[j][i][0] << " ";
	//	}
	//	cout << endl;
	//}
	cout << "error" << "  " << t << endl;
}