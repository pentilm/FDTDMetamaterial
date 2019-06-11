#include "PML.h"
//double PML::Hx_fcn(double x, double y, double t){
//	return exp(t)*sin(x + y);
//}
//double PML::Hy_fcn(double x, double y, double t){
//	return exp(t)*sin(x + y);
//}
//double PML::Kx_fcn(double x, double y, double t){
//	return exp(t)*sin(x + y);
//}
//double PML::Ky_fcn(double x, double y, double t){
//	return exp(t)*sin(x + y);
//}
//double PML::E_fcn(double x, double y, double t){
//	return exp(t)*sin(x + y);
//}
//double PML::J_fcn(double x, double y, double t){
//	return exp(t)*sin(x + y);
//}
//double PML::g1_fcn(double x, double y, double t){
//	return -exp(t)*((-1 + (sgmx(x) - sgmy(y) - 1)*mu0(x, y))*sin(x + y) - cos(x + y));
//}
//double PML::g2_fcn(double x, double y, double t){
//	return exp(t)*((1 + (sgmx(x) - sgmy(y) + 1)*mu0(x, y))*sin(x + y) - cos(x + y));
//}
//double PML::g3_fcn(double x, double y, double t){
//	return exp(t)*sin(x + y)*(1 + (sgmx(x) + sgmy(y) + 1)*e0(x, y));
//}
//double PML::g4_fcn(double x, double y, double t){
//	return -exp(t)*sin(x + y)*(e0(x, y)*pJE(x, y) - pJ(x, y) - 1);
//}
//double PML::g5_fcn(double x, double y, double t){
//	return -exp(t)*sin(x + y)*(mu0(x, y)*pKxHx(x, y) - pKx(x, y) - 1);
//}
//double PML::g6_fcn(double x, double y, double t){
//	return -exp(t)*sin(x + y)*(mu0(x, y)*pKyHy(x, y) - pKy(x, y) - 1);
//}
double PML::Hx_fcn(double x, double y, double t){
	return 0.0;
}
double PML::Hy_fcn(double x, double y, double t){
	return 0.0;
}
double PML::Kx_fcn(double x, double y, double t){
	return 0.0;
}
double PML::Ky_fcn(double x, double y, double t){
	return 0.0;
}
double PML::E_fcn(double x, double y, double t){
	return 0.0;
}
double PML::J_fcn(double x, double y, double t){
	return 0.0;
}
double PML::g1_fcn(double x, double y, double t){
	return 0.0;
}
double PML::g2_fcn(double x, double y, double t){
	return 0.0;
}
double PML::g3_fcn(double x, double y, double t){
	return 0.0;
}
double PML::g4_fcn(double x, double y, double t){
	return 0.0;
}
double PML::g5_fcn(double x, double y, double t){
	return 0.0;
}
double PML::g6_fcn(double x, double y, double t){
	return 0.0;
}
double PML::e0(double x, double y){
	//int reg = WhichReg(x, y);
	//if ((reg == 1) | (reg == 2))
		return 8.85*1e-12;
	//else
	//	return 1.0;
}
double PML::mu0(double x, double y){
	//int reg = WhichReg(x, y);
	//if ((reg == 1) | (reg == 2))
		return 4*PI*1e-7;
	//else
	//	return 1.0;
}
int PML::WhichReg(double x, double y){
	if ((xM0 < x) & (x < xM1) & (yM0 < y) & (y < yM1))
	{
		return 3;
	}
	else if ((apx < x) & (x < bpx) & (apy < y) & (y < bpy))
	{
		return 2;
	}
	else
	{
		return 1;
	}
}
double PML::sgmx(double x){
	//double sm = -log(err)*5.0*bpx*cv/(2.0*dd);
	double sm = -log(err)*5.0*cv / (2.0*dd);
	if (x >= bpx)
	{
		return sm*pow((x - bpx) / dd, 4);
	}
	else if (x <= apx)
	{
		return sm*pow(fabs(x / dd), 4);
	}
	else
		return 0.0;

}
double PML::sgmy(double y){
	//double sm = -log(err)*5.0*bpy*cv / (2.0*dd);
	double sm = -log(err)*5.0*cv / (2.0*dd);
	if (y >= bpy)
	{
		return sm*pow((y - bpy) / dd, 4);
	}
	else if (y <= apy)
	{
		return sm*pow(fabs(y / dd), 4);
	}
	else
		return 0.0;
}
double PML::source(double x, double y, double t){
	double tm, sp;
	double Tp = 1 / sw_f0;
	double x1 = t / sw_m / Tp;
	double x2 = (t - (sw_m + sw_k)*Tp) / sw_m / Tp;
	double g1 = 10 * x1*x1*x1 - 15 * x1*x1*x1*x1 + 6 * x1*x1*x1*x1*x1;
	double g2 = 1 - (10 * x2*x2*x2 - 15 * x2*x2*x2*x2 + 6 * x2*x2*x2*x2*x2*x2);
	if ((0 < t)&(t < sw_m*Tp))
	{
		tm = g1*sin(sw_w0*t);
	}
	else if ((sw_m*Tp < t) & (t < (sw_m + sw_k)*Tp))
	{
		tm = sin(sw_w0*t);
	}
	else if (((sw_m + sw_k)*Tp < t) & (t < (2 * sw_m + sw_k)*Tp))
	{
		tm = g2*sin(sw_w0*t);
	}
	else
	{
		tm = 0.0;
	}
	sp = exp(-(y - 0.032)*(y - 0.032) / 250.0 / dx[0] / dx[0]);
	return tm*sp;
}
double PML::pJE(double x, double y){
	if ((reg == 1) | (reg == 2))
		return sgmx(x)*sgmy(y);
	else if (reg==3)
	{
		return we*we;
	}
	return 0.0;
}
double PML::pJ(double x, double y){
	if (reg == 3)
	{
		return ge;
	}
	return 0.0;
}
double PML::pKxHx(double x, double y){
	if ((reg == 1) | (reg == 2))
		return sgmx(x)*(sgmx(x)-sgmy(y));
	else if (reg == 3)
	{
		return wm*wm;
	}
	return 0.0;
}
double PML::pKx(double x, double y){
	if ((reg == 1) | (reg == 2))
		return sgmx(x);
	else if (reg == 3)
	{
		return gm;
	}
	return 0.0;
}
double PML::pKyHy(double x, double y){
	if ((reg == 1) | (reg == 2))
		return -sgmy(y)*(sgmx(x) - sgmy(y));
	else if (reg == 3)
	{
		return wm*wm;
	}
	return 0.0;
}
double PML::pKy(double x, double y){
	if ((reg == 1) | (reg == 2))
		return sgmy(y);
	else if (reg == 3)
	{
		return gm;
	}
	return 0.0;
}
